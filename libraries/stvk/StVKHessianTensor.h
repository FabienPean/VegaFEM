/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "StVK" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
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

#ifndef VEGAFEM_STVKHESSIANTENSOR_H
#define VEGAFEM_STVKHESSIANTENSOR_H

/*
  The second derivative (Hessian tensor) of internal elastic forces.
  See also StVKInternalForces.h .
*/

#include "triple.h"
#include "sparseMatrix.h"
#include "volumetricMesh.h"
#include "StVKElementABCD.h"
#include "StVKStiffnessMatrix.h"

namespace vegafem
{

class StVKHessianTensor
{
public:

  StVKHessianTensor(StVKStiffnessMatrix * stVKStiffnessMatrix);
  virtual ~StVKHessianTensor();

  inline int numVertices() { return numVertices_; }
  inline int numElements() { return numElements_; }

  // uses the Hessian tensor to approximate the change in the tangent stiffness matrix, if configuration is changed from u to u + du (du is assumed small)
  // any required Hessian tensor quantities (in configuration u) are computed automatically
  virtual void ComputeStiffnessMatrixCorrection(double * u, double * du, SparseMatrix * dK);

  // compute the vector result=(H:u)v, where H is the Hessian in the zero deformation configuration
  // high-memory version; must call ComputeHessianAtZero before calling EvaluateHessianQuadraticForm
  int ComputeHessianAtZero(int verbose=1); 
  int SaveHessianAtZeroToFile(const char * filename);
  void EvaluateHessianQuadraticForm(const double * u, const double * v, double * result); 

  // compute the vector result=(H:u)v, where H is the Hessian in the zero deformation configuration
  // low-memory version; no need to call ComputeHessianAtZero, but longer computation times for smaller meshes
  void EvaluateHessianQuadraticFormDirect(double * u, double * v, double * result);

  // computes (H:u)v, for all pairs of columns of the matrix Ulin; k is the number of columns of Ulin
  // the first numRigidModes columns will not be used for this computation
  // result is the output matrix; must be pre-allocated with (k-numRigidModes) * (k-numRigidModes+1) / 2 columns
  // low-memory version; no need to call ComputeHessianAtZero
  void EvaluateHessianQuadraticFormDirectAll(double * Ulin, int k, double * result, int numRigidModes=0, int verbose=1);
  
  // low-level routines (advanced use)
  void AddQuadraticTermsContribution(double * u, double * du, SparseMatrix * dK, int elementLow=-1, int elementHigh=-1);
  void AddCubicTermsContribution(double * u, double * du, SparseMatrix * dK, int elementLow=-1, int elementHigh=-1);

protected:
  int numVertices_;
  int numElements_;

  StVKStiffnessMatrix * stVKStiffnessMatrix;
  VolumetricMesh * volumetricMesh;
  StVKElementABCD * precomputedIntegrals;

  int numElementVertices;
 
  typedef triple<int,int,int> triIndex;
  typedef std::map<triIndex,double*> hessianType;
  hessianType hessian;
  void AddTensor3x3x3Block(int v1, int v2, int v3, Vec3d & vec, int type);

  double * lambdaLame;
  double * muLame;
};


}//namespace vegafem

#endif


