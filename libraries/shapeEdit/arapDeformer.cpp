/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "shapeEdit" library , Copyright (C) 2018 USC                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Hongyi Xu, Koki Nagano, Yijing Li, Jernej Barbic        *
 * http://www.jernejbarbic.com/vega                                      *
 *                                                                       *
 * Research: Jernej Barbic, Hongyi Xu, Yijing Li,                        *
 *           Danyong Zhao, Bohan Wang,                                   *
 *           Fun Shing Sin, Daniel Schroeder,                            *
 *           Doug L. James, Jovan Popovic                                *
 *                                                                       *
 * Funding: National Science Foundation, Link Foundation,                *
 *          Singapore-MIT GAMBIT Game Lab,                               *
 *          Zumberge Research and Innovation Fund at USC,                *
 *          Sloan Foundation, Okawa Foundation,                          *
 *          USC Annenberg Foundation                                     *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#include "arapDeformer.h"
#include <float.h>
#include "vec3d.h"
//#include "matrix.h"
#include "polarDecomposition.h"
#include "objMeshOrientable.h"
#include "triKey.h"
#include "tetKey.h"
#include "rectKey.h"
#include "geometryQuery.h"
#include "generateLaplacian.h"
#include "triple.h"
#include "basicAlgorithms.h"
#include <algorithm>
#include <iostream>
#ifdef USE_TBB
  #include <tbb/tbb.h>
#else
  #include "range.h"
#endif
//#include "matrixBLAS.h"
using namespace std;

#define USE_DIRICHLET

// ARAP Energy: E = \sum_i \sum_j\in N(i) wij | (pi'-pj') - Ri (pi-pj) |^2,
//                  i: vtx index, j: nbr index in i's neighborhood N(i), pi: rest position of vtx i, pi': deformed position of i
//                  Ri: the current rotation around vtx i
//                  wij = wji = 0.5 * \sum_k cot(angle_ikj), the weight on edge (i,j) is the 0.5 * sum of each angle /_ikj which
//                  opposites the edge (i,j) on a triangle (i,j,k)
// partial E/ partial pi' = 4 \sum_j\in N(i) wij ( (pi'-pj') - 0.5 * (Ri+Rj)(pi-pj)
//                        = 4 \sum_j wij [ (pi'-pj') - \sum_j 0.5 * wij (Ri+Rj)(pi-pj) ]
// To set pE = 0, we have: L p = b, where
//               L_ij = wij is a laplacian matrix, b_i = \sum_j 0.5 * wij (Ri+Rj)(pi-pj)

ARAPModel::ARAPModel(const ObjMesh * objMesh)
{
  n = objMesh->getNumVertices();
  for(int i = 0; i < n; i++)
    restVertices.push_back(objMesh->getPosition(i));

  L = generateLaplacian(objMesh, true);
  assert(L->HasInfOrNaN() == false);
  W = nullptr;
  initialize();
}

ARAPModel::ARAPModel(const TetMesh * tetMesh, const double * vtxWeights)
{
  n = tetMesh->getNumVertices();
  for(int i = 0; i < n; i++)
    restVertices.push_back(tetMesh->getVertex(i));

  SparseMatrixOutline LOutline(n), WOutline(n);
  const double min_angle = acos(0.9999);
  const double max_angle = acos(-0.9999);
  vector<UTriKey> visitedFaces;
  for(int ele = 0; ele < tetMesh->getNumElements(); ele++)
  {
    OTetKey tet(tetMesh->getVertexIndices(ele));
    for(size_t j = 0; j < 4; j++)
    {
      UTriKey uface = tet.uFaceKey(j);
      visitedFaces.push_back(uface);
    }
  }
  sortAndDeduplicate(visitedFaces);

#ifdef USE_TBB
  tbb::enumerable_thread_specific<vector<triple<int,int,double>>> threadLocalWeightBuffer;
  tbb::parallel_for(tbb::blocked_range<int>(0, visitedFaces.size()), [&](const tbb::blocked_range<int> & rng)
  {
    auto & weightBuffer = threadLocalWeightBuffer.local();
#else
    vector<triple<int,int,double>> weightBuffer;
    Range<int> rng(0, visitedFaces.size());
#endif
    for(int face = rng.begin(); face != rng.end(); ++face)
    {
      UTriKey oface = visitedFaces[face];
      assert(oface[0] != oface[1] && oface[1] != oface[2] && oface[2] != oface[0]);
      for(int k = 0; k < 3; k++)
      {
        int v0 = oface[k];
        int v1 = oface[(k+1)%3];
        int v2 = oface[(k+2)%3];
        Vec3d p0 = restVertices[v0];
        Vec3d p1 = restVertices[v1];
        Vec3d p2 = restVertices[v2];

        //compute 0.5 * cot (angle(p1p0p2))
        double w = 0.5 / tan(clamp(getTriangleAngleRobust(p0, p1, p2), min_angle, max_angle));
        // cout << "computing vtx " << v0 << " " << v1 << " " << v2 << " " << w << endl;
        //clamp
        if (w < 0) w = 0;

        // formula from the paper: wij = 0.5 (cot alpha_ij + cot beta_ij)
        // alpha_ij and beta_ij are the angles opposite of the mesh edge (i, j)

        weightBuffer.emplace_back(v1, v2, w);
      }
    }
#ifdef USE_TBB
  }, tbb::static_partitioner());

  for(const auto & weightBuffer : threadLocalWeightBuffer)
#endif
  {
    for(const auto & t : weightBuffer)
    {
      int v1 = t.first, v2 = t.second;
      double w = t.third;
      WOutline.AddEntry(v1, v2, w);
      WOutline.AddEntry(v2, v1, w);
      if (vtxWeights) // default vtxWeight is 1.0 for each vtx.
        w *= (vtxWeights[v1] + vtxWeights[v2]) / 2.0;
      LOutline.AddEntry(v1, v2, -w);
      LOutline.AddEntry(v2, v1, -w);
      LOutline.AddEntry(v1, v1, w);
      LOutline.AddEntry(v2, v2, w);
    }
  }

  //we have to do this to ensure every row has at least one element
  //otherwise the solver will complain
  for(int row = 0; row < n; row++)
    LOutline.AddEntry(row, row, 0);

  L = new SparseMatrix(&LOutline);
  W = new SparseMatrix(&WOutline);

//  W->Save("W.sparse",0);
  initialize();
}

ARAPModel::ARAPModel(const CubicMesh * cubicMesh, const double * vtxWeights)
{
  n = cubicMesh->getNumVertices();
  for(int i = 0; i < n; i++)
    restVertices.push_back(cubicMesh->getVertex(i));

  SparseMatrixOutline LOutline(n), WOutline(n);

  // face vtx indices
  int faceVtx[6][4] = {{0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {1,2,6,5}, {0,4,7,3}};

  set<URectKey> visitedFaces;
  for(int ele = 0; ele < cubicMesh->getNumElements(); ele++)
  {
    for(int f = 0; f < 6; f++)
    {
      int v[4];
      v[0] = cubicMesh->getVertexIndex(ele, faceVtx[f][0]);
      v[1] = cubicMesh->getVertexIndex(ele, faceVtx[f][1]);
      v[2] = cubicMesh->getVertexIndex(ele, faceVtx[f][2]);
      v[3] = cubicMesh->getVertexIndex(ele, faceVtx[f][3]);
      URectKey key(v[0],v[1],v[2],v[3]);
      if (visitedFaces.find(key) != visitedFaces.end())
        continue;
      visitedFaces.insert(key);

      for(int i = 0; i < 4; i++)
      {
        int j = (i+1)%4;
        double w = 0.5;
        WOutline.AddEntry(v[i], v[j], w);
        WOutline.AddEntry(v[j], v[i], w);

        if (vtxWeights) // default vtxWeight is 1.0 for each vtx.
          w *= (vtxWeights[v[i]] + vtxWeights[v[j]]) / 2.0;

        LOutline.AddEntry(v[i], v[j], -w);
        LOutline.AddEntry(v[j], v[i], -w);
        LOutline.AddEntry(v[i], v[i], w);
        LOutline.AddEntry(v[j], v[j], w);
      }
    }
  }
  //we have to do this to ensure every row has at least one element
  //otherwise the solver will complain
  for(int row = 0; row < n; row++)
    LOutline.AddEntry(row, row, 0);

  L = new SparseMatrix(&LOutline);
  W = new SparseMatrix(&WOutline);

//  W->Save("W.sparse",0);
  initialize();
}

void ARAPModel::initialize()
{
  R.resize(n, Mat3d(1.0));
  vector<double> p0(n); //Lp = L * p0
  //LP stores the multiplication of L and p. By using Lp, we can rewrite the original formula in terms of u not p
  for(int i = 0; i < 3; i++)
  {
    Lp[i].resize(n);
    for(int j = 0; j < n; j++)
      p0[j] = restVertices[j][i];
    L->MultiplyVector(&p0[0],&Lp[i][0]);
  }
}

ARAPModel::~ARAPModel()
{
  delete L;
  delete W;
}

void ARAPModel::updateRotations(const double * disp)
{
  double Si[9], U[9], Ri[9];

  for(int vtxIdx = 0; vtxIdx < n; vtxIdx++)
  {
    Vec3d pos_i_def, pos_j_def, pos_i_rest, pos_j_rest;
    pos_i_rest = restVertices[vtxIdx];
    pos_i_def = pos_i_rest + (disp + vtxIdx*3);

    int rowLength = L->GetRowLength(vtxIdx);

    memset(Si, 0, sizeof(double) * 9);

    for(int col = 0; col < rowLength; col++)
    {
      int neighbor = L->GetColumnIndex(vtxIdx, col);
      if (neighbor == vtxIdx) continue;
      pos_j_rest = restVertices[neighbor];
      pos_j_def = pos_j_rest + (disp + neighbor*3);

      double w = -1 * (L->GetEntry(vtxIdx, col));
      Vec3d e_rest = pos_i_rest - pos_j_rest;
      Vec3d e_def = pos_i_def - pos_j_def;

      //w_ij * eij * (eij_deformed)T
      Mat3d cov = w * tensorProduct(e_rest, e_def);
      cov.addToArray(Si);
    }

    double tolerance = 1E-6;
    int forceRotation = 1;

    PolarDecomposition::Compute(Si, Ri, U, tolerance, forceRotation);

    Mat3d RiM(Ri);
    R[vtxIdx] = trans(RiM);
  }
}


void ARAPModel::buildPosRHS(std::vector<double> rhs[3])
{
  for(int i = 0; i < 3; i++)
    rhs[i].assign(n, 0.0);

  for(int vtxIdx = 0; vtxIdx < n; vtxIdx++)
  {
    int rowLength = L->GetRowLength(vtxIdx);
    Vec3d pos_i = restVertices[vtxIdx];
    // iterate all the neighbors of current vertex
    for(int entryIndex = 0; entryIndex < rowLength; entryIndex++)
    {
      // w_ij
      double w = -1 * (L->GetEntry(vtxIdx, entryIndex)); // L stores negative wij
      int nbrIdx = L->GetColumnIndex(vtxIdx, entryIndex);
      if (nbrIdx == vtxIdx)
        continue;

      // R_i + R_j
      Mat3d Rsum = R[vtxIdx] + R[nbrIdx];
      Vec3d pos_j = restVertices[nbrIdx];
      // p_i - p_j
      Vec3d pos_diff = pos_i - pos_j;
      // (w_ij / 2) * (R_i + R_j) * (p_i - p_j)
      Vec3d res = (w / 2) * Rsum * pos_diff;

      //seperate the dofs
      for(int d = 0; d < 3; d++)
        rhs[d][vtxIdx] += res[d];
    }
  }

//  vector<double> nrhs[3];
//  for(int d = 0; d < 3; d++)
//    nrhs[d].resize(n);
//  for(int i = 0; i < n; i++)
//    buildOneVertexRHS(i, nrhs);
//
//  for(int i = 0; i < n; i++)
//    nrhs[1][i] -= rhs[1][i];
}

void ARAPModel::buildDispRHS(std::vector<double> rhs[3])
{
  buildPosRHS(rhs);
  for(int d = 0; d < 3; d++)
    for(int i = 0; i < n; i++)
      rhs[d][i] -= Lp[d][i];
}

void ARAPModel::buildPosRHSOnOneRotation(int vtxIdx, std::vector<double> rhs[3])
{
  int rowLength = L->GetRowLength(vtxIdx);
  Vec3d pos_i = restVertices[vtxIdx];
  // iterate all the neighbors of current vertex
  for(int entryIndex = 0; entryIndex < rowLength; entryIndex++)
  {
    // w_ij
    double w = -1 * (L->GetEntry(vtxIdx, entryIndex)); // L stores negative wij
    int nbrIdx = L->GetColumnIndex(vtxIdx, entryIndex);
    if (nbrIdx == vtxIdx)
      continue;

    Vec3d pos_j = restVertices[nbrIdx];

    // p_i - p_j
    Vec3d pos_diff = pos_i - pos_j;

    // (w_ij / 2) * R_i * (p_i - p_j)
    Vec3d res = (w / 2) * R[vtxIdx] * pos_diff;

    //seperate the dofs
    for(int d = 0; d < 3; d++)
    {
      rhs[d][vtxIdx] += res[d];
      rhs[d][nbrIdx] -= res[d];
    }
  }
}

void ARAPModel::buildEnergyHessian(SparseMatrix ** sparseMatrix) const
{
  SparseMatrixOutline outline(3*n);

  outline.AddBlockMatrix(0,0, L, 1., 3);
  *sparseMatrix = new SparseMatrix(&outline);
}

// E = \sum_i \sum_j\in N(i) wij | (pi'-pj') - Ri (pi-pj) |^2,
// partial E/ partial pi' = 4 \sum_j\in N(i) wij ( (pi'-pj') - 0.5 * (Ri+Rj)(pi-pj) )
void ARAPModel::getEnergyAndGradient(const double * disp, double * energy, double * gradient)
{
  if (gradient)
    memset(gradient, 0, n * 3);
  updateRotations(disp);

  if (energy)
    *energy = 0.;

  for(int vtxIdx = 0; vtxIdx < n; vtxIdx++)
  {
    int rowLength = L->GetRowLength(vtxIdx);
    Vec3d resPi = restVertices[vtxIdx];
    Vec3d defPi = resPi + Vec3d(disp + 3*vtxIdx);
    // iterate all the neighbors of current vertex
    for(int j = 0; j < rowLength; j++)
    {
      // w_ij
      double w = -1 * (L->GetEntry(vtxIdx, j)); // L stores negative wij
      int nbrIdx = L->GetColumnIndex(vtxIdx, j);
      if (nbrIdx == vtxIdx)
        continue;
      Vec3d resPj = restVertices[nbrIdx];
      Vec3d defPj = resPj + Vec3d(disp + 3*nbrIdx);

      // p_i - p_j
      Vec3d resDif = resPi - resPj;
      Vec3d defDif = defPi - defPj;

      // (R_i) * (p_i - p_j)
      Vec3d rotResDif = R[vtxIdx] * resDif;

      if (energy)
        *energy += w * len2(defDif - rotResDif);

      if (gradient)
      {
        Vec3d g = defDif - 0.5 * rotResDif;
        g *= 4 * w;
        g.addToArray(gradient + 3 * vtxIdx);
      }
    }
  }
}

void ARAPModel::setRotations(const Mat3d * rotations)
{
  for(int i = 0; i < n; i++)
    R[i] = rotations[i];
//  memcpy(R.data(), rotations, sizeof(Mat3d) * n);
}

void ARAPModel::assembleDimension(const std::vector<double> input[3], double * output)
{
  int n = input[0].size();
  for(int i = 0; i < n; i++)
  {
    output[3 * i + 0] = input[0][i];
    output[3 * i + 1] = input[1][i];
    output[3 * i + 2] = input[2][i];
  }
}

ARAPDeformer::ARAPDeformer(const ObjMesh * objMesh, int numFixedVertices_, const int * fixedVertices_, int numThreads_) :
    arapModel(objMesh), numThreads(numThreads_)
{
  //assert(objMesh->isTriangularMesh());
  n = objMesh->getNumVertices();
  n3 = 3*n;

  fixedVertices.resize(numFixedVertices_);
  memcpy(fixedVertices.data(), fixedVertices_, sizeof(int) * numFixedVertices_);
  fixedVerticesSet.insert(fixedVertices.begin(), fixedVertices.end());
}

ARAPDeformer::ARAPDeformer(const TetMesh * tetMesh, int numFixedVertices_, const int * fixedVertices_, int numThreads_) :
    arapModel(tetMesh), numThreads(numThreads_)
{
  n = tetMesh->getNumVertices();
  n3 = 3*n;

  fixedVertices.resize(numFixedVertices_);
  memcpy(fixedVertices.data(), fixedVertices_, sizeof(int) * numFixedVertices_);
  fixedVerticesSet.insert(fixedVertices.begin(), fixedVertices.end());
}

ARAPDeformer::ARAPDeformer(const CubicMesh * cubicMesh, int numFixedVertices_, const int * fixedVertices_, int numThreads_) :
    arapModel(cubicMesh), numThreads(numThreads_)
{
  n = cubicMesh->getNumVertices();
  n3 = 3*n;

  fixedVertices.resize(numFixedVertices_);
  memcpy(fixedVertices.data(), fixedVertices_, sizeof(int) * numFixedVertices_);
  fixedVerticesSet.insert(fixedVertices.begin(), fixedVertices.end());
}

ARAPDeformer::~ARAPDeformer()
{
  delete solver[0];
  delete solver[1];
  delete solver[2];
}

bool ARAPDeformer::checkHandles(const std::map<int, Vec3d> & newHandles)
{
  bool rebuild = false;
  if (newHandles.size() != handles.size())
    rebuild = true;
  else
  {
    // checking handle vertex indices and updating handle displacements
    map<int, Vec3d>::iterator oldIt = handles.begin();
    map<int, Vec3d>::const_iterator newIt = newHandles.begin();
    for(; newIt != newHandles.end(); oldIt++, newIt++)
      if (oldIt->first != newIt->first)
      {
        rebuild = true;
        break;
      }
      else
      {
        oldIt->second = newIt->second;
        //cout << "move handle " << oldIt->first  << " to " << oldIt->second << endl;
      }
  }
  return rebuild;
}

void ARAPDeformer::updateHandles(const std::map<int, Vec3d> & newHandles)
{
  bool rebuild = checkHandles(newHandles);
  if (rebuild)
    rebuildSolver(newHandles);
}

void ARAPDeformer::rebuildSolver(const std::map<int, Vec3d> & newHandles)
{
  for(int i = 0; i < 3; i++)
    delete (solver[i]);
  handles = newHandles;

  // no constraints, so return
  if (newHandles.size() == 0)
  {
    for(int i = 0; i < 3; i++)
      solver[i] = nullptr;
    return;
  }

  cout << "Build ARAPDeformer solver with " << handles.size() << " handles" << endl;

  // form constraints
#ifdef USE_DIRICHLET
  set<int> constrainedVertices(fixedVertices.begin(), fixedVertices.end());
  for(map<int, Vec3d>::const_iterator it = newHandles.begin(); it != newHandles.end(); it++)
  {
    constrainedVertices.insert(it->first);
    assert(it->first < n);
  }
  vector<int> constrainedVerticesVec(constrainedVertices.begin(), constrainedVertices.end());
  for(size_t i = 0; i < constrainedVerticesVec.size(); i++)
    assert(constrainedVerticesVec[i] >= 0 && constrainedVerticesVec[i] < n);

  for(int i = 0; i < 3; i++)
  {
    //reconstruct the solver and do Cholesky decomposition
    int updatable = 0, addDirichlet = 1;
    solver[i] = new LagrangeMultiplierSolver(arapModel.getL(), nullptr, nullptr, constrainedVerticesVec.size(), constrainedVerticesVec.data(),
        numThreads, updatable, addDirichlet);
  }
#else
  SparseMatrixOutline JOutline(newHandles.size());
  int row = 0;
  for(map<int, Vec3d>::const_iterator it = newHandles.begin(); it != newHandles.end(); it++)
  {
    JOutline.AddEntry(row, it->first, constrainedScale);
    row++;
  }
  JOutline.AddEntry(0, n-1, 0);
  SparseMatrix J(&JOutline);

  for(int i = 0; i < 3; i++)
  {
    //reconstruct the solver and do Cholesky decomposition
    solver[i] = new LagrangeMultiplierSolver(L, &J, nullptr, fixedVertices.size(), fixedVertices.data(), numThreads);
  }
#endif
}

void ARAPDeformer::deform(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration)
{
  if (handles.size() == 0)
  {
    memset(disp, 0, sizeof(double) * 3 * n);
    return;
  }

  double ep2 = epsilon * epsilon;


  vector<double> buffer(n3);
  memcpy(buffer.data(), dispLast, sizeof(double) * n3);

  double error2 = 0;
  unsigned int iter = 0;
  double dispLastLen2 = 0;
  do
  {
    dispLastLen2 = 0;
    for(int i = 0; i < n3; i++)
      dispLastLen2 += buffer[i] * buffer[i];

    deformOneIter(&buffer[0], disp);
    error2 = 0;
    for(int i = 0; i < n3; i++)
    {
      error2 += (buffer[i] - disp[i]) * (buffer[i] - disp[i]);
    }
    cout << sqrt(error2 / dispLastLen2) << " ";
    iter++;
    if (iter >= maxIteration)
      break;
    memcpy(buffer.data(), disp, sizeof(double) * n3);
  } while(error2 / dispLastLen2 > ep2);
  cout << iter << " " << endl;
}

void ARAPDeformer::deformOneIter(const double * dispLast, double * disp)
{
  if (handles.size() == 0) //if there is no constraints, set to rest configuration
  {
    memset(disp, 0, sizeof(double) * 3 * n);
    return;
  }

  arapModel.updateRotations(dispLast);
  optimizePositions(disp);
}

void ARAPDeformer::deformOneIter(double * disp)
{
  deformOneIter(disp, disp);
}

void ARAPDeformer::optimizePositions(double * disp)
{
  vector<double> rhs[3], x[3];
  for(int i = 0; i < 3; i++)
  {
    rhs[i].resize(n);
    x[i].resize(n);
  }

  arapModel.buildDispRHS(rhs);

  for(std::map<int, Vec3d>::iterator it = handles.begin(); it != handles.end(); it++)
  {
    Vec3d disp = it->second;

    int vtx = it->first;
    for(int d = 0; d < 3; d++)
      x[d][vtx] = disp[d];
  }

  // 3 solves
  for(int i = 0; i < 3; i++)
    solver[i]->SolveLinearSystem(&x[i][0], &rhs[i][0]);

  arapModel.assembleDimension(x, disp);
}

void ARAPDeformer::setNumThreads(int threads)
{
  if (solver[0])
    for(int i = 0; i < 3; i++)
      solver[i]->SetNumThreads(threads);
}

//void ARAPDeformer::setConstraintMatrix(const SparseMatrix * C[3], const SparseMatrix * B[3])
//{
//  vector<int> oldFixedDOFs(solver[0]->GetFixedDOFs(), solver[0]->GetFixedDOFs() + solver[0]->GetNumFixedDOFs());
//
//  if (C[0] == nullptr)
//    consRhs.clear();
//  else
//    consRhs.assign(C[0]->GetNumRows(), 0.0);
//  for(int i = 0; i < 3; i++)
//  {
//    assert((C[i] == nullptr && consRhs.size() == 0) || (C[i] && C[i]->GetNumRows() == (int)consRhs.size()));
//    //reconstruct the solver and do Cholesky decomposition
//    int updatable = 0, addDirichlet = 1;
//    solver[i] = new LagrangeMultiplierSolver(L, C[i], B[i], oldFixedDOFs.size(), oldFixedDOFs.data(),
//        numThreads, updatable, addDirichlet);
//  }
//}
//
//void ARAPDeformer::setConstraintRhs(const double * constraintRhs)
//{
//  memcpy(consRhs.data(), constraintRhs, consRhs.size() * sizeof(double));
//}
//
//void ARAPDeformer::clearConstraints()
//{
//  vector<const SparseMatrix *> C(3, nullptr);
//  setConstraintMatrix(&C[0]);
//}

