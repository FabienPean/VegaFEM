/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "virtualTets" library , Copyright (C) 2018 USC                        *
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

#ifndef VIRTUALTETS_H
#define VIRTUALTETS_H

#include "tetMeshGeo.h"
#include "triMeshGeo.h"
#include "barycentricCoordinates.h"
#include "vec3ER.h" // 3-vector of exact reals
#include "tetTriCuttingData.h"
#include "profiler.h"

/*
  This file implements our "virtual tets" method, as described in:

  Yijing Li, Jernej Barbic: Immersion of Self-Intersecting Solids and Surfaces,
  ACM Transactions on Graphics 37(4) (SIGGRAPH 2018), Vancouver, Canada, 2018.

  Creates a virtual tet mesh that embeds an intersection-free triangle mesh,
  given an input tet mesh that covers the space occupied by the triangle mesh.

  The triangle mesh is assumed to be manifold, free of self-intersections and have no degenerate triangles.
  The method duplicates the tets and tet vertices to follow the topology of the triangle mesh.

  Inputs: tetMesh, triMesh.

  Outputs: the virtual tet mesh (the return value)
           coord: barycentric coordinates embedding the input triangle mesh into the output tet mesh
           newTetVtxID2OldTetVtxID: the mapping: new tet vtxID -> old tet vtxID
           newTetID2OldTetID: the mapping: new tetID -> old tetID
           tetTris: the mapping: new tetID -> triangleIDs that this tet embeds

  If a pointer is provided as nullptr, then that output is not provided.
*/

// The simplest routine with minimal input.
TetMeshGeo createVirtualTetsMesh(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo & triMesh,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool verbose = false,
    Profiler * profiler = nullptr);

// If you have already previously computed the input triangle mesh vertex positions in exact arithmetic,
// use this routine to avoid re-computing them.
// exactTriPos: triangle positions in exact arithmetics
TetMeshGeo createVirtualTetsMesh(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo & triMesh,
    const std::vector<Vec3ER> & exactTriPos,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool verbose = false, Profiler * profiler = nullptr);

// This is an advanced routine. It is used as a subroutine by our immersion algorithm for self-intersecting inputs.
// cuttingData: Detailed information on tet mesh and triangle mesh intersection.
// removeAirTets: "Air tets" are tets that are completely outside the triangle mesh.
//                When set to true, the method will find all air tets and remove them.
//                When false, the method assumes that the input tet mesh is free of air tets.
TetMeshGeo createVirtualTetsMesh(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo & triMesh,
    const TetTriCuttingData & cuttingData,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool removeAirTets = true,
    bool verbose = false,
    Profiler * profiler = nullptr);

#endif

