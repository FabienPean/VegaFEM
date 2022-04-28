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

#ifndef ARAPDEFORMER_H
#define ARAPDEFORMER_H

#include "sparseMatrix.h"
#include "objMesh.h"
#include "LagrangeMultiplierSolver.h"
#include <set>
#include <map>
#include <vector>
#include "tetMesh.h"
#include "cubicMesh.h"

/*
  As-rigid-as-possible model for obj meshes, tet meshes and cubic meshes.

  1. For obj meshes, this class implements the model for fast geometric shape deformation described in the following paper:

    Olga Sorkine and Marc Alexa: As-Rigid-As-Possible Surface Modeling
    In Proc. of Eurographics Symposium on Geometry Processing (2007), Vol. 4, p. 30

  2. For tet and cubic meshes, the class implements the extension for volumetric meshes described in the same paper.
*/

class ARAPModel
{
public:
  ARAPModel(const ObjMesh * objMesh);
  ARAPModel(const TetMesh * tetMesh, const double * vtxWeights = nullptr);
  ARAPModel(const CubicMesh * tetMesh, const double * vtxWeights = nullptr);
  virtual ~ARAPModel();

  void updateRotations(const double * disp);

  void buildDispRHS(std::vector<double> rhs[3]);
  void buildPosRHS(std::vector<double> rhs[3]);

  const SparseMatrix * getL() const { return L; }
  const SparseMatrix * getW() const { return W; }
  const Vec3d & getRestVertex(int i) const { return restVertices[i]; }
  const std::vector<Vec3d> & getRestVertices() const { return restVertices; }

  const Mat3d & getRotation(int i) const { return R[i]; }
  const std::vector<Mat3d> & getRotations() const { return R; }

  void setRotations(const Mat3d * rotations);

  // we omit the derivative of rotation Ri, so energy hessian is always the same: a weighted Laplacian matrix
  // the hessian is of size 3*#vtx
  void buildEnergyHessian(SparseMatrix ** sparseMatrix) const;

  // compute energy and its gradient; the derivative of rotation is omitted in calculation
  // either energy or gradient can be nullptr, dim is of size 3*#vtx
  void getEnergyAndGradient(const double * disp, double * energy, double * gradient = nullptr);

  void buildPosRHSOnOneRotation(int vtxIdx, std::vector<double> rhs[3]);

  static void assembleDimension(const std::vector<double> input[3], double * output);

protected:
  void initialize();

  int n;
  SparseMatrix * L;
  SparseMatrix * W; // W is the wij matrix
  std::vector<Vec3d> restVertices;
  std::vector<double> Lp[3];
  std::vector<Mat3d> R;
};


class ARAPDeformer
{
public:
  ARAPDeformer(const ObjMesh * objMesh, int numFixedVertices, const int * fixedVertices, int numThreads);
  ARAPDeformer(const TetMesh * tetMesh, int numFixedVertices, const int * fixedVertices, int numThreads);
  ARAPDeformer(const CubicMesh * tetMesh, int numFixedVertices, const int * fixedVertices, int numThreads);
  virtual ~ARAPDeformer();

  void setNumThreads(int threads);

  virtual void updateHandles(const std::map<int, Vec3d> & newHandles);

  // do one iteration of rotation and position optimization, return a rest shape (disp = 0) if no constraints
  void deformOneIter(double * disp); //disp serves as input and output
  virtual void deformOneIter(const double * dispLast, double * disp);

  // do rotation and position optimization until relative error < epsilon or maxIteration is reached
  void deform(const double * dispLast, double * disp, double epsilon, unsigned int maxIteration);

  const SparseMatrix * getL() const { return arapModel.getL(); }

protected:
  // update handle displacements. return true if handle vertex indices don't match
  bool checkHandles(const std::map<int, Vec3d> & newHandles);

  void rebuildSolver(const std::map<int, Vec3d> & newHandles);

  //Solve 3 linear systems (n * n)
  //call Lagrange multiplier solver SolveLinearSystem()
  //store the position to sceneObjectDeformable (it uses displacement. instead of position)
  void optimizePositions(double * disp);

  ARAPModel arapModel;

  std::map<int, Vec3d> handles;

  std::vector<int> fixedVertices;
  std::set<int> fixedVerticesSet;

  int numThreads;

  //per-edge Laplacian matrix: n * n

  LagrangeMultiplierSolver * solver[3] = {nullptr, nullptr, nullptr};
  int n, n3;
};

#endif
