/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "distance field" library , Copyright (C) 2007 CMU, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Danyong Zhao, Jernej Barbic                             *
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
  This code implements marching cubes with topological guarantees, 
  similar to the following publications. Our code has been implemented
  from scratch by Danyong Zhao and Jernej Barbic.

  Thomas Lewiner, Hélio Lopes, Antônio Wilson Vieira, Geovan Tavares:
  Efficient implementation of marching cubes cases with topological guarantees
  Journal of Graphics Tools 8(2): pp. 1-15, 2003

  E. V. Chernyaev. Marching Cubes 33: construction of topologically correct isosurfaces. 
  Technical Report CERN CN 95–17, CERN, 1995.
*/

#ifndef _MARCHINGCUBES_H_
#define _MARCHINGCUBES_H_

#include <vector>
#include <map>

#include "triple.h"
#include "objMesh.h"
#include "distanceFieldBase.h"

class MarchingCubes
{
public:

  // computes the isosurface mesh, using marching cubes with topological guarantees
  // the input distance field can be complete, or narrow band
  // isoValue is the isosurface value to be meshed (given in absolute units, not grid units)
  // output: a triangle mesh corresponding to the isosurface
  static ObjMesh * compute(const DistanceFieldBase * distanceField, float isoValue = 0.0);

protected:
  MarchingCubes(const DistanceFieldBase * distanceField, float isoValue = 0.0);
  virtual ~MarchingCubes() {}

  const DistanceFieldBase * distanceFieldBase;
  float isoValue;

  int resolutionX, resolutionY, resolutionZ;

  // executes the marching cubes algorithm
  ObjMesh * compute();

  // void computeTriangleVertices(int i, int j, int k, bool center, int vtx[13]);
  bool faceTest(int face, float cube[8]);
  int interiorTest(int edge, float cube[8]);
  inline float getDistance(int i, int j, int k) const
  {
    float offset = distanceFieldBase->distance(i, j, k) - isoValue;
//    return offset + (offset == 0.0f) * FLT_EPSILON;
    if (offset == 0.0f) return FLT_EPSILON;
    return offset;
  }

  // computes all the tables needed for marching cubes
  static void createTable();
  static void printTable();
};

#endif

