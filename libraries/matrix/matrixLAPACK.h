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
  Wrapper to matrix LAPACK routines.
  See also matrix.h.
*/

#ifndef _MATRIX_LAPACK_H_
#define _MATRIX_LAPACK_H_

#include <stdlib.h>

template<class real>
real * InverseMatrix(int m, const real * mtx, real * output = NULL);

template<class real>
real * PseudoInverseMatrix(int m, int n, const real * mtx, real singularValueThreshold = 1e-6, int * rank = NULL, real * output = NULL);

template<class real>
real * MatrixLeastSquareSolve(int m, int n, int nRhs, const real * mtx, const real * b, real rcond, int * rank = NULL, real * output = NULL);

template<class real>
void MatrixLUSolve(int n, int nRhs, const real * mtx, real * x, const real * b);

// mtx * x = lambda * x
template<class real>
void SymmetricMatrixEigenDecomposition(int m, real * mtx, real * Q, real * Lambda);

// mtx * x = lambda * mtx2 * x
// Warning: mtx2 will be modified after this function call.
template<class real>
void SymmetricMatrixGeneralEigenDecomposition(int m, real * mtx, real * mtx2, real * Q, real * Lambda);

// mtx * x = lambda * x
template<class real> 
int MatrixEigenDecomposition(int m, real * mtx, real * EigenVectors, real * LambdaRe, real * LambdaIm);

// U is m x MIN(m,n)
// Sigma is a MIN(m,n) vector
// VT is MIN(m,n) x n
template<class real> 
void MatrixSVD(int m, int n, real * mtx, real * U, real * Sigma, real * VT);

#endif

