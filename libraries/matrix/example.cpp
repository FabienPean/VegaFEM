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

#include "matrix.h"

/*
  A simple example using the matrix class.
*/

int main()
{
  double Av[9] = {2.0, 1.0, 1.0,
                  1.0, 2.0, 1.0,
                  1.0, 1.0, 2.0, };

  Matrix<double> A(3, 3, Av);
  Matrix<double> B = 0.5 * Transpose(A) * A;
  printf("Matrix:\n");
  B.Print();
  char outputFilename[96] = "BMatrix";
  B.Save(outputFilename);

  Matrix<double> Q(3,3); // will hold eigenvectors
  Matrix<double> Lambda(3,1); // will hold eigenvalues
  B.SymmetricEigenDecomposition(Q, Lambda);

  printf("Eigenvalues:\n");
  Lambda.Print();

  return 0;
}

