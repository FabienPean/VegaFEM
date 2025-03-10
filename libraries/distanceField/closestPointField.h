/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "distance field" library , Copyright (C) 2007 CMU, 2018 USC           *
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
  ClosestPointField is similar to a distance field, 
  except it also stores the nearest point, not just the distance to it.
*/

#ifndef VEGAFEM_CLOSEST_POINT_FIELD_H
#define VEGAFEM_CLOSEST_POINT_FIELD_H

#include "distanceField.h"

namespace vegafem
{

class ClosestPointField : public DistanceField
{
public:
  ClosestPointField();
  ~ClosestPointField();

  int computeUnsignedField(ObjMesh * objMesh, int resolutionX, int resolutionY, int resolutionZ, int maxTriCount=15, int maxDepth=10, int zMin = -1, int zMax = -1);
  int computeSignedField(ObjMesh * objMesh, int resolutionX, int resolutionY, int resolutionZ, int maxTriCount=15, int maxDepth=10, int zMin = -1, int zMax = -1);

  virtual int load(const std::string& filename);
  virtual int save(const std::string& filename, bool doublePrecision=true);
  virtual void set(int resolutionX, int resolutionY, int resolutionZ, Vec3d bmin_, Vec3d bmax_, float * distanceData, float * closestPointData);

  void closestPoint(float pos[3], float result[3]);

  inline void closestPoint(int i, int j, int k, float result[3]) { int pos = 3*((k * (resolutionY+1) + j ) * (resolutionX+1) + i); result[0] = closestPointData[pos]; result[1] = closestPointData[pos+1]; result[2] = closestPointData[pos+2]; }

  virtual bool sanityCheck();

  void getClosestPointData(float ** floatBuffer);

  inline void setClosestPoint(int i, int j, int k, float value[3]) { int pos = 3*((k * (resolutionY+1) + j ) * (resolutionX+1) + i); closestPointData[pos] = value[0]; closestPointData[pos+1] = value[1]; closestPointData[pos+2] = value[2]; }

protected:
  float * closestPointData;

  virtual int ZigZagSigned(void * objMeshOctree, void * meshGraph);
  virtual int ZigZagUnsigned(void * objMeshOctree, void * meshGraph);
  virtual int ZigZagSigned(void * objMeshOctree, void * meshGraph, int zLo, int zHi, int asterisk=0);
  virtual int ZigZagUnsigned(void * objMeshOctree, void * meshGraph, int zLo, int zHi, int asterisk=0);
};


}//namespace vegafem

#endif


