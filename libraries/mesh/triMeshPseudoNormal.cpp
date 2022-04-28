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

#include "triMeshPseudoNormal.h"
#include "triMeshNeighbor.h"
#include <cassert>
using namespace std;

TriMeshPseudoNormal::TriMeshPseudoNormal(TriMeshRef triMesh, const Vec3d * extTriNormals) :
  vtxNormals(triMesh.numVertices()), edgeNormals(triMesh.numVertices()), triNormals(triMesh.numTriangles())
{
  TriMeshNeighbor nbr(triMesh);
  for(int triID = 0; triID < triMesh.numTriangles(); triID++)
  {
    if (extTriNormals)
      triNormals[triID] = extTriNormals[triID];
    else
      triNormals[triID] = triMesh.computeTriangleNormal(triID);
    assert(triNormals[triID].hasNaN() == false);
  }

  for(int vtxID = 0; vtxID < triMesh.numVertices(); vtxID++)
  {
    auto & vtxNormal = vtxNormals[vtxID];
    vtxNormal = Vec3d(0.0);
    for(int triID : nbr.getVtxNearbyTriangles(vtxID))
    {
      double angle = triMesh.getTriangleAngleAtVertexRobust(triID, vtxID);
      vtxNormal += angle * triNormals[triID];
    }
    if (len2(vtxNormal) > 0)
    {
      vtxNormal.normalize();
      assert(vtxNormal.hasNaN() == false);
    }
  }

  for(int vtxID = 0; vtxID < triMesh.numVertices(); vtxID++)
  {
    for(const auto & p : nbr.getVertexEdgeTriMap(vtxID))
    {
      int vtxID2 = p.first;
      int triID = p.second;

      if (vtxID == vtxID2)
      {
        throw 1; // error, triMesh is not valid!
      }
      int v0 = vtxID, v1 = vtxID2;
      if (v0 > v1) swap(v0, v1);

      auto iter = edgeNormals[v0].find(v1);
      if (iter == edgeNormals[v0].end())
      {
        edgeNormals[v0][v1] = triNormals[triID];
      }
      else
      {
        iter->second += triNormals[triID];
        iter->second.normalize();
      }
    }
  }
}

const Vec3d & TriMeshPseudoNormal::edgeNormal(int vtxID0, int vtxID1) const
{
  if (vtxID0 > vtxID1) { swap(vtxID0, vtxID1); }
  const auto & edgeMap = edgeNormals[vtxID0];
  auto iter = edgeMap.find(vtxID1);
  assert(iter != edgeMap.end());
//  if (iter == edgeMap.end()) { return Vec3d(0.0); }
  return iter->second;
}
