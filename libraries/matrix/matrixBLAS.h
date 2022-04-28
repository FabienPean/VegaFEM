/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "matrix" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC         *
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

/*
  Wrappers to matrix BLAS routines.
  See also matrix.h.
*/

#ifndef _MATRIX_BLAS_H_
#define _MATRIX_BLAS_H_

#include <stdlib.h>

template<class real>
real * SumMatrices(int m, int n, const real * mtx1, const real * mtx2, real * output = NULL);

template<class real>
real * SubtractMatrices(int m, int n, const real * mtx1, const real * mtx2, real * output = NULL);

// output = mtx1 * mtx2
// mtx1 is m x p
// mtx2 is p x n
template<class real>
real * MultiplyMatrices(int m, int p, int n, const real * mtx1, const real * mtx2, real * output = NULL);
// output += mtx1 * mtx2
template<class real>
void MultiplyMatricesAdd(int m, int p, int n, const real * mtx1, const real * mtx2, real * output);
// output -= mtx1 * mtx2
template<class real>
void MultiplyMatricesSub(int m, int p, int n, const real * mtx1, const real * mtx2, real * output);

// output = trans(mtx1) * mtx2
// trans(mtx1) is m x p
// mtx2 is p x n
template<class real>
real * MultiplyMatricesT(int m, int p, int n, const real * mtx1, const real * mtx2, real * output = NULL);
// output += trans(mtx1) * mtx2
template<class real>
void MultiplyMatricesTAdd(int m, int p, int n, const real * mtx1, const real * mtx2, real * output);
// output -= trans(mtx1) * mtx2
template<class real>
void MultiplyMatricesTSub(int m, int p, int n, const real * mtx1, const real * mtx2, real * output);

template<class real>
real * ScalarMultiplyMatrix(int m, int n, real alpha, const real * mtx, real * output = NULL);

// transposes the matrix (without making a separate copy)
template <class real>
void InPlaceTransposeMatrix(int m, int n, real * mtx);

// computes Euclidean norm of a vector
template <class real>
real VectorNorm(int m, const real * vec);

#endif

