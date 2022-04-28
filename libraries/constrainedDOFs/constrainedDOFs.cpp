/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "constrainedDOFs" library , Copyright (C) 2007 CMU, 2009 MIT          *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic, Yijing Li                                 *
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

#include "constrainedDOFs.h"
#include "stdio.h"
#include "stdlib.h"
#include <cassert>

void ConstrainedDOFs::InsertDOFs(int mFull, const double * xConstrained, double * x, int numFixedRows, const int * fixedRows, int oneIndexed)
{
  int destRow = 0; // current row in the dest matrix
  int sourceRow = 0; // in source

  for(int i=0; i<numFixedRows; i++)
  {
    int index = fixedRows[i] + 1 - oneIndexed;
    assert((index <= mFull) && (index >= 1));
    //if ((index > mFull) || (index < 1))
    index--;

    while (destRow < index)
    {
      // while row index smaller than index, keep on copying from source
      x[destRow] = xConstrained[sourceRow];
      destRow++;
      sourceRow++;
    }

    // insert zero row
    x[destRow] = 0.0;
    destRow++;
  }
  
  while (destRow < mFull)
  { 
    x[destRow] =  xConstrained[sourceRow];
    destRow++;
    sourceRow++;
  }
}

void ConstrainedDOFs::RemoveDOFs(int mFull, double * xConstrained, const double * x, int numFixedRows, const int * fixedRows, int oneIndexed)
{
  int numrows = 0;
  int row = 0;

  for(int i=0; i<numFixedRows; i++)
  {
    int index = fixedRows[i] + 1 - oneIndexed;
    //if ((index > mFull) || (index < 1))
    assert((index <= mFull) && (index >= 1));
    index--;

    while (row<index)
    {
      xConstrained[numrows] = x[row];
      numrows++;
      row++;
    }

    row++; // skip the deselected row

    if (numrows > mFull)
    {
      printf("Error: too many rows specified.\n");
      exit(1);
    }
  }

  while (row < mFull)
  {
    xConstrained[numrows] = x[row];

    numrows++;
    row++;

    if (numrows > mFull)
    {
      printf("Error: too many rows specified.\n");
      exit(1);
    }
  }
}

void ConstrainedDOFs::FullDOFsToConstrainedDOFs(int mFull, int numDOFs, int * DOFsConstrained, const int * DOFs, int numFixedRows, const int * fixedRows, int oneIndexed)
{
  if (numDOFs == 0)
    return;

  int dof = 0;

  for(int i=0; i<numFixedRows; i++)
  {
    // each iteration processes one bucket of fixed vertices

    // correct for (optional) 1-indexing
    int index = fixedRows[i] + 1 - oneIndexed;
    assert((index <= mFull) && (index >= 1));
    //if ((index > mFull) || (index < 1))
    index--;

    while (DOFs[dof] < index)
    {
      DOFsConstrained[dof] = DOFs[dof] - i; 
      dof++;
      if (dof >= numDOFs)
        return;
    }

    // assign -1 to the fixed DOF
    if (DOFs[dof] == index)
    {
      DOFsConstrained[dof] = -1;
      dof++;
    }
  }

  while (dof < numDOFs)
  {
    DOFsConstrained[dof] = DOFs[dof] - numFixedRows; 
    dof++;
  }
}

void ConstrainedDOFs::FindFreeDOFs(int mFull, int * freeDOFs, int numFixedDOFs, const int * fixedDOFs, int oneIndexed)
{
  assert(mFull >= numFixedDOFs);
  int mFree = mFull - numFixedDOFs;
  int constrainedDOF = 0;
  int superDOF = 0;
  for(int i=0; i<numFixedDOFs; i++)
  { 
    int nextFixedDOF = fixedDOFs[i];
    nextFixedDOF -= oneIndexed;
    assert((nextFixedDOF < mFull) && (nextFixedDOF >= 0));

    while (superDOF < nextFixedDOF)
    {
      assert(constrainedDOF < mFree); // check if too many DOFs specified
      freeDOFs[constrainedDOF] = superDOF;
      constrainedDOF++;
      superDOF++;
    }

    superDOF++; // skip the deselected DOF
  }
  while (superDOF < mFull)
  { 
    assert(constrainedDOF < mFree);
    freeDOFs[constrainedDOF] = superDOF; // check if too many DOFs specified

    constrainedDOF++;
    superDOF++;
  }
}

bool ConstrainedDOFs::CheckValidSortedDOFs(int mFull, int numDOFs, const int * DOFs, int oneIndexed)
{
  int lowBound = 0 + oneIndexed;
  int highBound = mFull + oneIndexed;
  for(int i = 0; i < numDOFs; i++)
  {
    if (DOFs[i] < lowBound || DOFs[i] >= highBound)
      return false;
    if (i > 0 && DOFs[i] <= DOFs[i-1])
      return false;
  }
  return true;
}

