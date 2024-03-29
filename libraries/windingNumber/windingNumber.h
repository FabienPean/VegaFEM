/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "windingNumber" library , Copyright (C) 2018 USC                      *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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
  Computes the winding number of a 3D point with respect to a given 
  closed manifold triangle mesh, as described in:

  Alec Jacobson, Ladislav Kavan, Olga Sorkine-Hornung:
  Robust Inside-Outside Segmentation using Generalized Winding Numbers, SIGGRAPH 2013

  This implementation does not use the hierarchical acceleration presented
  in the paper above; instead, it loops over all the triangles.
*/

#ifndef VEGAFEM_WINDINGNUMBER_H
#define VEGAFEM_WINDINGNUMBER_H

#include "objMesh.h"
#include "minivector.h"

namespace vegafem
{

class WindingNumber 
{
public:

  // computes the winding number of point p with respect to the given mesh
  static double computeWindingNumber(ObjMesh * objMesh, const Vec3d & p);
};


}//namespace vegafem

#endif


