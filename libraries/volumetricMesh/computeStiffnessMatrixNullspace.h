/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Large Modal Deformation Factory",                                    *
 * a pre-processing utility for model reduction of                       *
 * deformable objects undergoing large deformations.                     *
 *                                                                       *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2018 USC                           *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef VEGAFEM_COMPUTESTIFFNESSMATRIXNULLSPACE_H
#define VEGAFEM_COMPUTESTIFFNESSMATRIXNULLSPACE_H

#include <cstdlib>
#include <cstdio>
#include <cmath>

namespace vegafem
{

/*
  Computes the nullspace of the tangent stiffness matrix at the specified deformation. At the origin (u=0), stiffness matrix has nullspace of dimension 6, elsewhere, it has nullspace of dimension 3.

  Use "includeRotationalNullspace" to set the desired nullspace size (3 when deformed, or 6 when undeformed).
*/

class ComputeStiffnessMatrixNullspace
{
public:
  // n is the number of vertices
  // vertexPos is a 3n array giving vertex positions
  // (if undeformed, this will be the rest positions, otherwise the current world-coordinate deformed positions)
  static void ComputeNullspace(int n, const double * vertexPos, double * basis, int includeRotationalNullspace, int generateOrthogonalBasis=0);

  static void RemoveNullspaceComponent(int n, int nullspaceDimension, const double * nullspaceOrthonormalBasis, double *x);

};


}//namespace vegafem

#endif


