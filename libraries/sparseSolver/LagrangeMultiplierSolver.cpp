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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "LagrangeMultiplierSolver.h"
#include "constrainedDOFs.h"
#include "sparseMatrix.h"
#include "performanceCounter.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <set>
#include <algorithm>
using namespace std;

LagrangeMultiplierSolver::LagrangeMultiplierSolver(const SparseMatrix * A, const SparseMatrix * J, const SparseMatrix * B, int numFixedDOFs_, const int * fixedDOFs_, int numThreads_, int updatable, int addDirichlet) : numFixedDOFs(numFixedDOFs_), numThreads(numThreads_)
{
  decompositionDone = false;
  fixedDOFs.assign(fixedDOFs_, fixedDOFs_ + numFixedDOFs);
  // create system matrix
  m = A->GetNumRows();
  mFree = m - numFixedDOFs;
  assert(mFree >= 0);
  if (J != NULL)
    c = J->GetNumRows();
  else 
    c = 0;

  if ((B != NULL) && (B->GetNumRows() != c)) 
  {
    printf("Error: matrix size mismatch in B.\n");
    throw 1;
  }

  if (B != NULL) 
  {
    SparseMatrixOutline outline(m+c);
    for(int i=0; i<m; i++)
      for(int j=0; j<A->GetRowLength(i); j++)
        outline.AddEntry(i, A->GetColumnIndex(i,j), A->GetEntry(i,j));

    for(int i=0; i<c; i++)
      for(int j=0; j<J->GetRowLength(i); j++)
      {
        outline.AddEntry(m+i, J->GetColumnIndex(i,j), J->GetEntry(i,j));
        outline.AddEntry(J->GetColumnIndex(i,j), m+i, J->GetEntry(i,j));
      }

    for(int i=0; i<c; i++)
      for(int j=0; j<B->GetRowLength(i); j++)
        outline.AddEntry(m+i, m+B->GetColumnIndex(i,j), B->GetEntry(i,j));

    // must add this because Pardiso requires diagonal elements to be set, even if zero
    for(int i=0; i<c; i++)
      outline.AddEntry(m+i, m+i, 0.0);

    systemMatrix = new SparseMatrix(&outline);
  }
  else 
  {
    systemMatrix = new SparseMatrix(*A);
    if (c > 0)
      systemMatrix->AppendRowsColumns(J);
  }
  systemMatrix->RemoveRowsColumns(numFixedDOFs, &fixedDOFs[0]);

  ADirichlet = JDirichlet = NULL;
  if (addDirichlet)
  {
    // [ S11 S12 ] [x1]   [b1]
    // [ S21 S22 ] [x2] = [b2], x2 is the fixedDOF
    // then Dirichlet boundary condition is to give fixed values on x2
    // To solve x1, we have:
    // S11 x1 = b1 - S12 x2
    // S12 is formed by
    // S12 = [ A12 ],  A12 is ADirichlet
    //       [ J12 ],  J12 is JDirichlet
    DirichletBuffer.resize(numFixedDOFs);
    DirichletFreeBuffer.resize(max(mFree, c));
    ADirichlet = new SparseMatrix(*A);
    ADirichlet->RemoveRows(numFixedDOFs, &fixedDOFs[0]);

    vector<int> remainedDOFs(mFree);
    ConstrainedDOFs::FindFreeDOFs(m, &remainedDOFs[0], numFixedDOFs, &fixedDOFs[0]);

    ADirichlet->RemoveColumns(remainedDOFs.size(), remainedDOFs.data());

    assert(ADirichlet->GetNumColumns() <= numFixedDOFs);
    if (updatable)
      ADirichlet->BuildSuperMatrixIndices(numFixedDOFs, &fixedDOFs[0], remainedDOFs.size(), remainedDOFs.data(), A);

    if (c > 0)
    {
      JDirichlet = new SparseMatrix(*J);
      JDirichlet->RemoveColumns(remainedDOFs.size(), remainedDOFs.data());
      if (updatable)
        JDirichlet->BuildSuperMatrixIndices(0, NULL, remainedDOFs.size(), remainedDOFs.data(), J);
    }
  }

  pardisoSolver = new PardisoSolver(systemMatrix, numThreads, PardisoSolver::REAL_SYM_INDEFINITE);
  //pardisoSolver->ComputeCholeskyDecomposition(systemMatrix);

  xConstrained.resize(m+c - numFixedDOFs);
  rhsConstrained = xConstrained;

  Acons = NULL;
  Jcons = NULL;
  JTcons = NULL;
  if (updatable)
  {
    Acons = new SparseMatrix(*A);
    Acons->RemoveRowsColumns(numFixedDOFs, &fixedDOFs[0]);
    Acons->BuildSuperMatrixIndices(numFixedDOFs, &fixedDOFs[0], A);

    if (c > 0)
    {
      Jcons = new SparseMatrix(*J);
      Jcons->RemoveColumns(numFixedDOFs, &fixedDOFs[0]);
      Jcons->BuildSuperMatrixIndices(0, NULL, numFixedDOFs, &fixedDOFs[0], J);
      Jcons->BuildTranspositionIndices();
      JTcons = Jcons->Transpose(m-numFixedDOFs);
      systemMatrix->BuildSubMatrixIndices(*JTcons, 0, 0, m-numFixedDOFs);
    }

    if (B != NULL)
      systemMatrix->BuildSubMatrixIndices(*B, 1, m-numFixedDOFs, m-numFixedDOFs);
  }
//  cout << "LagMulSolver: " << endl << profiler.toString() << endl;
}

LagrangeMultiplierSolver::LagrangeMultiplierSolver(const LagrangeMultiplierSolver & other)
{
  decompositionDone = false;
  m = other.m;
  c = other.c;
  numFixedDOFs = other.numFixedDOFs;
  fixedDOFs = other.fixedDOFs;
  numThreads = other.numThreads;
  systemMatrix = new SparseMatrix(*other.systemMatrix);
  xConstrained = other.xConstrained;
  rhsConstrained = other.rhsConstrained;
  Acons = (other.Acons ? new SparseMatrix(*other.Acons) : NULL);
  Jcons = (other.Jcons ? new SparseMatrix(*other.Jcons) : NULL);
  JTcons = (other.JTcons ? new SparseMatrix(*other.JTcons) : NULL);
  ADirichlet = (other.ADirichlet ? new SparseMatrix(*other.ADirichlet) : NULL);
  JDirichlet = (other.JDirichlet ? new SparseMatrix(*other.JDirichlet) : NULL);
  DirichletBuffer = other.DirichletBuffer;
  DirichletFreeBuffer = other.DirichletFreeBuffer;
  mFree = other.mFree;
  pardisoSolver = new PardisoSolver(systemMatrix, numThreads, PardisoSolver::REAL_SYM_INDEFINITE);
}

LagrangeMultiplierSolver::~LagrangeMultiplierSolver()
{
  delete systemMatrix;
  delete pardisoSolver;
  delete Acons;
  delete Jcons;
  delete JTcons;
  delete ADirichlet;
  delete JDirichlet;
}

void LagrangeMultiplierSolver::AddDirichletOnRhs(const double * x)
{
  if (DirichletBuffer.size() > 0)
  {
    for(int i = 0; i < numFixedDOFs; i++)
      DirichletBuffer[i] = x[fixedDOFs[i]];
    assert(ADirichlet);
    ADirichlet->MultiplyVector(&DirichletBuffer[0], &DirichletFreeBuffer[0]);
    for(int i = 0; i < mFree; i++)
      rhsConstrained[i] -= DirichletFreeBuffer[i];

    if (JDirichlet)
    {
      JDirichlet->MultiplyVector(&DirichletBuffer[0], &DirichletFreeBuffer[0]);
      for(int i = 0; i < c; i++)
        rhsConstrained[mFree + i] -= DirichletFreeBuffer[i];
    }
  }
}

MKL_INT LagrangeMultiplierSolver::SolveLinearSystem(double * x, const double * rhs)
{
  if (decompositionDone == false)
  {
    pardisoSolver->FactorMatrix(systemMatrix);
    decompositionDone = true;
  }

  ConstrainedDOFs::RemoveDOFs(m+c, &rhsConstrained[0], rhs, numFixedDOFs, &fixedDOFs[0]);

  AddDirichletOnRhs(x);

  MKL_INT code =  pardisoSolver->SolveLinearSystem(&xConstrained[0], &rhsConstrained[0]);
  ConstrainedDOFs::InsertDOFs(m+c, &xConstrained[0], x, numFixedDOFs, &fixedDOFs[0]);

  // add back Dirichlet boundary conditions on x since those values have been set to zero in ConstrainedDOFs::InsertDOFs
  if (DirichletBuffer.size() > 0)
  {
    for(int i = 0; i < numFixedDOFs; i++)
      x[fixedDOFs[i]] = DirichletBuffer[i];
  }
  return code;
}

double LagrangeMultiplierSolver::CheckLinearSystemSolution(const double * x, const double * rhs, int verbose)
{
  ConstrainedDOFs::RemoveDOFs(m+c, &rhsConstrained[0], rhs, numFixedDOFs, &fixedDOFs[0]);
  ConstrainedDOFs::RemoveDOFs(m+c, &xConstrained[0], x, numFixedDOFs, &fixedDOFs[0]);
  AddDirichletOnRhs(x);
  return systemMatrix->CheckLinearSystemSolution(&xConstrained[0], &rhsConstrained[0], verbose);
}

void LagrangeMultiplierSolver::UpdateAJB(const SparseMatrix * A, const SparseMatrix * J, const SparseMatrix * B)
{
  if (A != NULL)
  {
    assert(Acons);
    Acons->AssignSuperMatrix(*A);
    for(int i = 0; i < Acons->Getn(); i++) 
    {
      int rowLen = Acons->GetRowLength(i);
      for(int j = 0; j < rowLen; j++)
        systemMatrix->SetEntry(i,j,Acons->GetEntry(i,j));
    }
    //systemMatrix->AssignSubMatrix(*Acons, 0);

    if (ADirichlet)
      ADirichlet->AssignSuperMatrix(*A);
  }

  if (J != NULL && J->GetNumRows() > 0)
  {
    assert(Jcons);
    assert(JTcons);
    Jcons->AssignSuperMatrix(*J);
    for(int i = 0; i < Jcons->Getn(); i++) 
    {
      int rowLen = Jcons->GetRowLength(i);
      int row = i + Acons->Getn();
      for(int j = 0; j < rowLen; j++)
        systemMatrix->SetEntry(row,j,Jcons->GetEntry(i,j));
    }
    JTcons->AssignTransposedMatrix(*Jcons);
    systemMatrix->AssignSubMatrix(*JTcons,0);

    if (JDirichlet)
      JDirichlet->AssignSuperMatrix(*J);
  }

  if (B != NULL)
    systemMatrix->AssignSubMatrix(*B, 1);

  decompositionDone = false;
}

void LagrangeMultiplierSolver::SetNumThreads(int numThreads_)
{ 
  numThreads = numThreads_; 
  if (pardisoSolver != NULL)
    pardisoSolver->SetNumThreads(numThreads);
}

