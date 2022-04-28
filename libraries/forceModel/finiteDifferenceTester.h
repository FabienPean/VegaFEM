/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "forceModel" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC     *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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

/*
   Tests a Vega FEM force model, using two-point, or five-point finite differences.
*/

#ifndef FINITEDIFFERENCETESTER_H
#define FINITEDIFFERENCETESTER_H

#include "forceModel.h"
#include "sparseMatrix.h"
#include <vector>

class FiniteDifferenceTester
{
public:
  enum Mode
  {
    TWO_POINT, // f'(x) <- (f(x+h) - f(x)) / h
    FIVE_POINT // f'(x) <- (-f(x+2h) + 8 f(x+h) - 8 f(x-h) + f(x-2h)) / (12h)
  };

  // Warning: parallelism enabled by numThreads > 1 only works if forceModel->GetElasticEnergy() and
  //   forceModel->GetInternalForce() are parallel-safe
  FiniteDifferenceTester(ForceModel * forceModel, double timestep, Mode mode, int numThreads);
  virtual ~FiniteDifferenceTester();

  void setTimestep(double timestep) { h = timestep; }

  // use ForceModel::GetInternalForce() to compute analtical force f_a
  // use ForceModel::GetElasticEnergy() to compute finite difference force f_d
  // return relative error: ||f_a - f_d|| / ||f_a||, if ||f_a|| == 0, return ||f_a - f_d||
  double testInternalForce(const double * u);

  // use ForceModel::GetInternalForceAndMatrix() to compute analtical tangent stiffness matrix K_a
  // use ForceModel::GetInternalForce() to compute finite difference matrix K_d
  // return relative error: ||K_a - K_d||_F / ||K_a||_F,  || ||_F is frobNorm, if ||K_a|| == 0, return ||K_a - K_d||_F
  // to save time, we assume the topology (non-zero entry locations) are correct in the stiffness matrix. So when
  //   computing ||K_a - K_d||_F, we only compare those non-zero entries, ignoring other entries in K_d which is computed densely
  // relativeErrorOnUnsymmetry returns ||K_a - K_a^T||_F / ||K_a||_F, if ||K_a||_F == 0, return ||K_a - K_a^T||_F
  double testStiffnessMatrix(const double * u, double * relativeErrorOnUnsymmetry = nullptr);

protected:
  ForceModel * forceModel = nullptr;
  double h = 0.0; // timestep
  Mode mode = TWO_POINT;
  int numThreads = 0;
  int r = 0;

  std::vector<double> analyticForce, finiteForce;
  SparseMatrix * stiffnessMatrix0 = nullptr, * transposedStiffness0 = nullptr;
  SparseMatrix * unsymmetricStiffness = nullptr;
};

#endif /* FINITEDIFFERENCETESTER_H_ */
