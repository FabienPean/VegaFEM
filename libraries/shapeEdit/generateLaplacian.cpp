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

#include <iostream>
#include <vector>
#include <set>
#include <cassert>
#include <limits>
#include "generateLaplacian.h"
#include "tetMeshManifold.h"
#include "matrixBLAS.h"
#include "triMeshGeo.h"
#include "geometryQuery.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
using namespace std;

// the direction that is perpendicular to the line (start, dir) and go from the line to an external point p
static Vec3d getNormalOfLine(Vec3d start, Vec3d dir, Vec3d p)
{
  dir.normalize();
  double projOnLine = dot(p - start, dir);
  Vec3d closestPoint = start + (projOnLine / len2(dir)) * dir;
  Vec3d disp = p - closestPoint;
  disp.normalize();
  return disp;
}

static void addVolumetricLaplacianTerm(const TetMesh * tetMesh, int vi, int vj, int op0, int op1,
    SparseMatrixOutline & outline, int lineID)
{
  set<int> as;
  as.insert(vi);
  as.insert(vj);
  as.insert(op0);
  as.insert(op1);
  assert(as.size() == 4);

  Vec3d pi = tetMesh->getVertex(vi);
  Vec3d pj = tetMesh->getVertex(vj);
  Vec3d p0 = tetMesh->getVertex(op0);
  Vec3d p1 = tetMesh->getVertex(op1);
  Vec3d dir = p1 - p0;
  double lij = len(dir);
  Vec3d ni = getNormalOfLine(p0, dir, pi);
  Vec3d nj = getNormalOfLine(p0, dir, pj);
  double cosij = dot(ni, nj);
  //  cout << vi << " " << vj << " angle: " << acos(cosij) * 180 / M_PI << endl;
  double sinij = sqrt(1 - cosij * cosij);
  double cotij = cosij / sinij;
  assert(Vec3d::isNaN(cotij) == false);
  double w = lij * cotij / 6.;

  //if (lineID == vj && w < 0)
  //  w = 0;

  outline.AddEntry(lineID, vi, -w);
  outline.AddEntry(lineID, vj, w);
}

SparseMatrix * generateLaplacian(const TetMesh * tetMesh, bool LinearlyPrecise)
{
  int anotherIndex[4][3] = { {1,2,3}, {0,2,3}, {0,1,3}, {0,1,2} };
  int oppositeIndex[4][3][2] = { { {2,3}, {1,3}, {1,2} }, { {2,3}, {0,3}, {0,2} }, { {1,3}, {0,3}, {0,1} }, { {1,2}, {0,2}, {0,1} } };
  int n = tetMesh->getNumVertices();
  int t = tetMesh->getNumElements();
  SparseMatrixOutline outline(n);
  for(int e = 0; e < t; e++)
    for(int i = 0; i < 4; i++)
    {
      int vi = tetMesh->getVertexIndex(e, i);
      for(int j = 0; j < 3; j++)
      {
        int vj = tetMesh->getVertexIndex(e, anotherIndex[i][j]);
        assert(vi != vj);
        int op0 = tetMesh->getVertexIndex(e, oppositeIndex[i][j][0]);
        int op1 = tetMesh->getVertexIndex(e, oppositeIndex[i][j][1]);
        addVolumetricLaplacianTerm(tetMesh, vi, vj, op0, op1, outline, vi);
      }
    }

  if (LinearlyPrecise)
  {
    TetMeshManifold m;
    for(int i = 0; i < t; i++)
      m.add(tetMesh->getVertexIndices(i), i);
    const map<OTriKey, TetMeshManifold::Tetrahedron*> & surface = m.getSurfaceMap();
    const TetMeshManifold::TriMap & triMap = m.getTriMap();
    int triangleOtherIndex[3][2] = { {1,2}, {0,2}, {0,1} };
    for(map<OTriKey, TetMeshManifold::Tetrahedron *>::const_iterator it = surface.begin(); it != surface.end(); it++)
    {
      const OTriKey & tri = it->first;
      //cout << "tri: " << tri[0] << " " << tri[1] << " " << tri[2] << endl;
      UTriKey utri(tri.indices());
      TetMeshManifold::TriCIter triIter = triMap.find(utri);
      const TetMeshManifold::Tetrahedron * tet = triIter->second->getTet(0);
      assert(triIter->second->getTet(1) == NULL);
      int vf = tet->getOppositeVtx(triIter->second);
      assert(vf >= 0);
      for(int i = 0; i < 3; i++)
      {
        int vi = tri[i];
        assert(vi != vf);
        for(int j = 0; j < 3; j++)
        {
          int vj = tri[j];
          assert(vj != vf);
          int op0 = tri[triangleOtherIndex[j][0]];
          int op1 = tri[triangleOtherIndex[j][1]];
          addVolumetricLaplacianTerm(tetMesh, vf, vj, op0, op1, outline, vi);
        }
      }
    }
  }

  return new SparseMatrix(&outline);
}

// compute cotangent for angle aob
static double computeCotangent(const ObjMesh * objMesh, int va, int vo, int vb)
{
  // clamp cos values to remove the effect of nearly degenerate triangles
  const double min_angle = 0.015;
  const double max_angle = M_PI - min_angle;

  Vec3d pa = objMesh->getPosition(va);
  Vec3d pb = objMesh->getPosition(vb);
  Vec3d po = objMesh->getPosition(vo);

  double angle = getTriangleAngleRobust(po, pa, pb);
  angle = clamp(angle, min_angle, max_angle);
  double tanAngle = tan(angle);
  assert(tanAngle != 0.0);
  assert(tanAngle != numeric_limits<double>::infinity());
  assert(tanAngle != - numeric_limits<double>::infinity());
  return 1.0 / tanAngle;
}

SparseMatrix * generateLaplacian(const ObjMesh * objMesh, bool clampCotangent)
{
  bool triangulated = false;
  if (!objMesh->isTriangularMesh())
  {
    ObjMesh * triMesh = new ObjMesh(*objMesh);
    triMesh->triangulate();
    objMesh = triMesh;
    triangulated = true;
  }

  //assert(objMesh->isTriangularMesh());
  int n = objMesh->getNumVertices();

  //add entries to the
  SparseMatrixOutline outline(n);

  //ObjMesh * objMesh = sceneObjectDeformable->GetMesh();

  int pi0[] = { 2, 0, 1 };
  int pi1[] = { 1, 2, 0 };
  int pi2[] = { 0, 1, 2 };
  //w_ij = w_ji = 1/2 [ cot(alpha_ij) + cot(beta_ij) ],   alpha_ij, beta_ij are angles opposite of the edge (i,j)
  //typedef pair<int, int> Edge;
  for(size_t i = 0; i < objMesh->getNumGroups(); i++)
  {
    const ObjMesh::Group* g = objMesh->getGroupHandle(i);
    for(size_t j = 0; j < g->getNumFaces(); j++)
    {
      const ObjMesh::Face * face = g->getFaceHandle(j);
      //TODO: add a function to compute face angle
      assert(face->getVertexHandle(0) != face->getVertexHandle(1) && face->getVertexHandle(1) != face->getVertexHandle(2)
          && face->getVertexHandle(2) != face->getVertexHandle(0));
      for(int k = 0; k < 3; k++)
      {
        unsigned int v0 = face->getVertexHandle(pi0[k])->getPositionIndex();
        unsigned int v1 = face->getVertexHandle(pi1[k])->getPositionIndex();
        unsigned int v2 = face->getVertexHandle(pi2[k])->getPositionIndex();

        double w = 0.5 * computeCotangent(objMesh, v0, v2, v1);

        if (clampCotangent && w < 0) //clamp
          w = 0.0;
        assert(v0 != v1);
        outline.AddEntry(v0, v1, -w);
        outline.AddEntry(v1, v0, -w);

        outline.AddEntry(v0, v0, w);
        outline.AddEntry(v1, v1, w);
      }
    }
  }

  //we have to do this to ensure every row has at least one element
  //otherwise the solver will complain
  for(int row = 0; row < n; row++)
    outline.AddEntry(row, row, 0);


  // NOTE: the following commented algorithm only works on 2D. So it's disabled

//  if (linearlyPrecise)
//  {
//    typedef pair<int, int> edge;
//#define SORTED_EDGE(a, b) ((a) > (b) ? edge((b),(a)) : edge((a),(b)))
//    map<edge, std::vector<const ObjMesh::Face*> > faceMap;
//
//    for(size_t i = 0; i < objMesh->getNumGroups(); i++)
//    {
//      const ObjMesh::Group* g = objMesh->getGroupHandle(i);
//      for(size_t j = 0; j < g->getNumFaces(); j++)
//      {
//        const ObjMesh::Face * face = g->getFaceHandle(j);
//        for(int k = 0; k < 3; k++)
//        {
//          unsigned int v0 = face->getVertexHandle(pi0[k])->getPositionIndex();
//          unsigned int v1 = face->getVertexHandle(pi1[k])->getPositionIndex();
//          edge e = SORTED_EDGE(v0, v1);
//          faceMap[e].push_back(face);
//          assert(faceMap[e].size() <= 2);
//        }
//      }
//    }
//    // find all boundary edges
//    int numBoundaryEdges = 0;
//    for(map<edge, std::vector<const ObjMesh::Face*> >::iterator it = faceMap.begin(); it != faceMap.end(); it++)
//    {
//      edge e = it->first;
//      if (it->second.size() == 2) // skip non-boundary edges
//        continue;
//      numBoundaryEdges++;
//      const ObjMesh::Face* face = it->second[0];
//
//      for(int i = 0; i < 2; i++)
//      {
//        int vi = (i == 0 ? e.first : e.second);
//        int vj = (i == 1 ? e.first : e.second);
//        int vk = vi;
//        for(int k = 0; k < 3; k++)
//        {
//          vk = face->getVertexHandle(k)->getPositionIndex();
//          if (vk != vi && vk != vj)
//            break;
//        }
//        assert(vk != vi);
//
//        double w = 0.5 * computeCotangent(objMesh, vi, vk, vj);
//        if (clampCotangent && w < 0) //clamp
//          w = 0.0;
//
//        outline.AddEntry(vi, vk, -w);
//        outline.AddEntry(vi, vi, w);
//
//        outline.AddEntry(vj, vk, -w);
//        outline.AddEntry(vj, vj, w);
//      }
//    }
//    cout << "numBoundaryEdges: " << numBoundaryEdges << endl;
//  }

  if (triangulated)
    delete objMesh;

  return new SparseMatrix(&outline);
}

// return ||A 1_n||_2
double testConstantPrecision(const SparseMatrix * A)
{
  int n = A->GetNumRows();
  vector<double> ones(n, 1.0);
  vector<double> result(n, 0.0);
  A->MultiplyVector(&ones[0], &result[0]);
  return EuclideanNorm(n, &result[0]);
}

// return ||A \bar{V}||_2, \bar{V} is a nx3 matrix storing rest positions of the mesh
double testLinearPrecision(const SparseMatrix * A, const VolumetricMesh * mesh)
{
  int n = A->GetNumRows();
  assert(n == mesh->getNumVertices());
  vector<double> V(3*n, 0.0);
  vector<double> result(3*n, 0.0);
  for(int i = 0; i < n; i++)
  {
    const Vec3d & v = mesh->getVertex(i);
    for(int j = 0; j < 3; j++)
      V[i + j * n] = v[j];
  }

  for(int i = 0; i < 3; i++)
  {
    A->MultiplyVector(&V[i*n], &result[i*n]);
  }
  return EuclideanNorm(3*n, &result[0]);
}

