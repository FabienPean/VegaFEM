/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "sparseSolver" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC   *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li, Hongyi Xu                     *
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

#ifndef _LAGRANGEMULTIPLIERSOLVER_H_
#define _LAGRANGEMULTIPLIERSOLVER_H_

/*

  Solves

        [ A  J^T ] [y]      = [rhs0]
        [ J  B   ] [lambda]   [rhs1]

  where A is sparse, symmetric and square (m x m), J is sparse and rectangular (c x m) (often a constraint gradient), and B is square and symmetric (c x c) (often zero)

  This system arises when forward-simulating a dynamical system with constraints (enforced via Lagrange multipliers).

  The Pardiso solver is used to solve the system.

*/

#include "PardisoSolver.h"

class LagrangeMultiplierSolver
{
public:

  // the constructor will also factor the matrix
  // if J is passed as NULL, empty matrix is assumed for it (no constraints)
  // if B is passed as NULL, zero matrix is assumed for B
  // if updatable=1, one can use function UpdateAJB to update matrices A, J, B (note: unlike in the derived Updatable class, the entire matrix will be re-factored, even if only updating J or B)
  // if addDirichletBoundaryCondition=1, fixedDOFs can be updated by assigning correct fixed values in x when calling SolveLinearSystem
  LagrangeMultiplierSolver(const SparseMatrix * A, const SparseMatrix * J, const SparseMatrix * B = NULL,
      int numFixedDOFs = 0, const int * fixedDOFs = NULL, int numThreads = 0, int updatable = 0, int addDirichletBoundaryCondition = 0);
  LagrangeMultiplierSolver(const LagrangeMultiplierSolver &);
  virtual ~LagrangeMultiplierSolver();

  // updates the A, J, B matrix, assuming equal topology as previous A, J, B. 
  // Pass NULL if you do not wish to update the corresponding matrix. 
  // If you pass non-NULL for J (or B), then J (or B) MUST have been set in the constructor.
  // updatable=1 must have been set in the constructor to use this function.
  void UpdateAJB(const SparseMatrix * A, const SparseMatrix * J = NULL, const SparseMatrix * B = NULL); // updates the A, J, B matrix, assuming equal topology as previous A, J, B. Pass NULL if you do not wish to update the corresponding matrix. If you pass non-NULL for B, then B MUST have been set in the constructor.

  // solve the linear system, using Pardiso
  // the routine does not modify rhs
  // if addDirichletBoundaryCondition is 0 when initializing,
  // x is not an initial guess (it is ignored as input, and used only as output parameter)
  // otherwise, the fixedDOFs in x are used for Dirichlet boundary conditions
  // returns Pardiso exit code 
  MKL_INT SolveLinearSystem(double * x, const double * rhs);

  double CheckLinearSystemSolution(const double * x, const double * rhs, int verbose=1);

  const SparseMatrix * GetSystemMatrix() const { return systemMatrix; }

  int GetNumFixedDOFs() const { return numFixedDOFs; }
  const int * GetFixedDOFs() const { return fixedDOFs.data(); }

  int GetNumARows() const { return m; }
  int GetNumJRows() const { return c; }

  void SetNumThreads(int numThreads); // call at any time to change the number of threads

protected:
  void AddDirichletOnRhs(const double * x); // modify rhsConstrained to satisfy Dirichlet boundary condition
  bool decompositionDone;
  int m, c;
  int numFixedDOFs; 
  std::vector<int> fixedDOFs;
  int numThreads;
  SparseMatrix * systemMatrix;
  PardisoSolver * pardisoSolver;
  std::vector<double> xConstrained, rhsConstrained;
  SparseMatrix * Acons;
  SparseMatrix * Jcons;
  SparseMatrix * JTcons;
  SparseMatrix * ADirichlet, * JDirichlet; // used for updating Dirichlet conditions
  std::vector<double> DirichletBuffer, DirichletFreeBuffer;
  int mFree;
};

#endif


