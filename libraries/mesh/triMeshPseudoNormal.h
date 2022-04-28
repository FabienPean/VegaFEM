/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
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

#ifndef TRIMESHPSEUDONORMAL_H
#define TRIMESHPSEUDONORMAL_H

#include "triMeshGeo.h"
#include <map>

// compute pseudo-normal on a triangle mesh
// assuming TriMesh is edge-manifold, required by TriMeshNeighbor
class TriMeshPseudoNormal
{
public:
  TriMeshPseudoNormal() {} // empty
  TriMeshPseudoNormal(TriMeshRef triMesh, const Vec3d * triangleNormals = nullptr);

  int numVertices() const { return vtxNormals.size(); }
  int numTriangles() const { return triNormals.size(); }

  // return Vec3d(0.0) if this vtx has no nearby triangles
  const Vec3d & vtxNormal(int vtxID) const { return vtxNormals[vtxID]; }

  // assert (vtxID0, vtxID1) is a valid edge
  const Vec3d & edgeNormal(int vtxID0, int vtxID1) const;

  const Vec3d & triNormal(int triID) const { return triNormals[triID]; }

protected:
  std::vector<Vec3d> vtxNormals;
  std::vector<std::map<int, Vec3d>> edgeNormals;
  std::vector<Vec3d> triNormals;
};

#endif

