/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Hongyi Xu                                *
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

#ifndef VEGAFEM_GENERATEGRADIENTMATRIX_H
#define VEGAFEM_GENERATEGRADIENTMATRIX_H

#include "sparseMatrix.h"
#include "tetMesh.h"

namespace vegafem
{

// creates the sparse matrix G that computes the gradient
// of a 3D vector field on a tet mesh
// optionally, computes G^T G as well

class GenerateGradientMatrix
{
public:

  // for a 3D vector field on a tet mesh
  // weights serve as GTw^2G
  static void Generate(const TetMesh * tetMesh, SparseMatrix ** G, SparseMatrix ** GTG = NULL, double * GTGElementWeights = NULL); // of the entire mesh 

  static void GenerateElementMatrix(const TetMesh * tetMesh, int el, double * G); // of a single element; G is 9 x 12

  // the version that operates on scalar fields
  static void GenerateForScalarField(const TetMesh * tetMesh, SparseMatrix ** GTG, double * GTGElementWeights); // of the entire mesh

protected:
  static void GenerateElementMatrixEntries(const TetMesh * tetMesh, int el, double * m); // of a single element, condensed array "m" of 12 doubles
};


}//namespace vegafem

#endif


