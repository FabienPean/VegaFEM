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

#ifndef _CONSTRAINED_DOFS_H_
#define _CONSTRAINED_DOFS_H_

/*
  Insert or remove given components from a linear array.
  The dynamic solver uses these routines to fix the specified vertices and remove
  rows from mass and stiffness matrices, as necessary.
*/
class ConstrainedDOFs
{
public:
  // The input fixedRows array must be pre-sorted

  // inserts zero entries into an array, at the specified locations
  // the locations must be given with respect to the full array 
  // input: xConstrained
  // output: x
  static void InsertDOFs(int mFull, const double * xConstrained, double * x, int numFixedRows, const int * fixedRows, int oneIndexed=0);   

  // removes entries at the specified locations from an array
  // the locations must be given with respect to the full array 
  // input: x
  // output: xConstrained
  static void RemoveDOFs(int mFull, double * xConstrained, const double * x, int numFixedRows, const int * fixedRows, int oneIndexed=0);   

  // translates the array indices from original indices to indices after removal of the specified entries
  // input: DOFs (must be sorted) (0-indexed)
  // output: DOFsConstrained (0-indexed)
  // oneIndexed applies only to fixedRows array, NOT to DOFsConstrained or DOFs
  static void FullDOFsToConstrainedDOFs(int mFull, int numDOFs, int * DOFsConstrained, const int * DOFs, int numFixedRows, const int * fixedRows, int oneIndexed=0);   

  // find the free DOFs that do not belong to fixedDOFs
  // input: fixedDOFs, pre-sorted
  // output freeDOFs, sorted
  // oneIndexed applies only to fixedDOFs array, not freeDOFs
  // example: mFull=7, fixedDOFs=(2, 3, 5), oneIndexed=0 => freeDOFs=(0,1,4,6)
  static void FindFreeDOFs(int mFull, int * freeDOFs, int numFixedDOFs, const int * fixedDOFs, int oneIndexed=0);  

  // return true if DOFs is sorted (from small to large) and in range of [0, mFull) if (oneIndexed == 0) or [1, mFull] if (oneIndex == 1)
  static bool CheckValidSortedDOFs(int mFull, int numDOFs, const int * DOFs, int oneIndexed=0);

};
#endif

