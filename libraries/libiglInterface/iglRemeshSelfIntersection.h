/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "libiglInterface" library , Copyright (C) 2018 USC                    *
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

#ifndef IGLREMESHSELFINTERSECTION_H
#define IGLREMESHSELFINTERSECTION_H

#include "triMeshGeo.h"
#include "vec3ER.h"
namespace iglInterface
{

struct SelfCutMeshData
{
  TriMeshGeo cutMesh; // the new vtx at triangle intersections stitch the intersected triangles together
  std::vector<Vec3ER> cutPosExact; // exact cut positions
  std::vector<std::pair<int, int>> interTriPairs; // Indices of the uncut mesh
  std::vector<int> oldTriIDs;
  std::vector<int> triPatchIDs; // triID -> patchID it belongs to
  std::vector<int> vtxStitchIDs; // vtxID -> new ID if vtx on the same positions merge

  // cellID around each patch
  // patchID -> <cellID on front side of patch, cellID on back side of patch>
  // by default, cellID = 0 is the infinity cell (outer space of the mesh)
  std::vector<std::pair<int, int>> cellIDsAtPatch;
  // winding numbers around each triangle
  // triID -> <#winding on front side of triID, #winding on back side of triID>
  std::vector<std::pair<int, int>> windAroundTri;
};

SelfCutMeshData remeshSelfIntersection(TriMeshRef mesh, bool stitch = true, bool computeWindingNumbers = true,
    bool computeDoublePrecisionPos = true, bool computePatchAndCells = true);

}

#endif
