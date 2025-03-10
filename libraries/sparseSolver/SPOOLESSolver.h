/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "sparseSolver" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC   *
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

#ifndef VEGAFEM_SPOOLESSOLVER_H
#define VEGAFEM_SPOOLESSOLVER_H

/*
  Solves A * x = rhs, where A is sparse, usually large, and symmetric.
  
  The solution is obtained using the SPOOLES library (which is free software).
  The solution method is direct (not iterative). As such, convergence
  is often very robust, and there is no need to tune convergence parameters,
  unlike, say, in the Conjugate gradient method.
  Memory requirements are minimized by re-ordering the matrix before applying
  Cholesky decomposition.
  However, for very large systems (e.g. 200,000 x 200,000 matrices on a 
  2Gb machine), the Cholesky decomposition might run out of memory.
*/

#include <cstdio>
#include <cstdlib>
#include "linearSolver.h"
#include "sparseMatrix.h"

namespace vegafem
{

class SPOOLESSolver : public LinearSolver
{
public:

  // this constructor will re-order A (in an internal copy), and then perform complete Cholesky factorization (via SPOOLES)
  // A is not modified
  SPOOLESSolver(const SparseMatrix * A, int verbose=0);
  virtual ~SPOOLESSolver();

  // solve: A * x = rhs, using SPOOLES
  // uses the Cholesky factors obtained in the constructor
  // rhs is not modified
  virtual int SolveLinearSystem(double * x, const double * rhs);

protected:
  int n;
  void * bridgePointer;
  void * mtx_xPointer;
  void * mtx_rhsPointer;
  void * APointer;
  FILE * msgFile;
  int verbose;

  static void DisabledSolverError();
};


}//namespace vegafem

#endif


