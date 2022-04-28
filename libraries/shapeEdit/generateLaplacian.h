/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "shapeEdit" library , Copyright (C) 2018 USC                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Hongyi Xu, Koki Nagano, Yijing Li, Jernej Barbic        *
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

#ifndef GENERATELAPLACIAN_H
#define GENERATELAPLACIAN_H

#include "sparseMatrix.h"
#include "tetMesh.h"
#include "objMesh.h"

// generate Laplacian matrix for tet mesh [Wang 2015 Linear Subspace Design for Real-Time Shape Deformation]
SparseMatrix * generateLaplacian(const TetMesh * tetMesh, bool linearlyPrecise);

// generate Laplacian matrix for obj mesh [Sorkine 2007 ARAP]
// if clampCotangent == true, clamp the cotangent weights in the Lapalcian matrix to be >= 0
SparseMatrix * generateLaplacian(const ObjMesh * objMesh, bool clampCotangent);

// return ||A 1_n||_F
double testConstantPrecision(const SparseMatrix * A);

// return ||A \bar{V}||_F, \bar{V} is a nx3 matrix storing rest positions of the mesh
double testLinearPrecision(const SparseMatrix * A, const VolumetricMesh * mesh);
//double testLinearPrecision(const SparseMatrix * A, const ObjMesh * mesh);

#endif
