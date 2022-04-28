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

#ifndef VIRTUALTETSVIACSG_H
#define VIRTUALTETSVIACSG_H

#include "tetMeshGeo.h"
#include "triMeshGeo.h"
#include "barycentricCoordinates.h"
#include "vec3ER.h"
#include "tetTriCuttingData.h"
#include "profiler.h"

/*
  This file implements "virtual tets" using constructive solid geometry, as described in:

  Eftychios Sifakis, Kevin Der, and Ronald Fedkiw. 2007.
  Arbitrary Cutting of Deformable Tetrahedralized Objects.
  In Symp. on Computer Animation (SCA 2017). 73–80.

  The usage of these routines is the same as in virtualTets.h . Please see that file for usage.
*/

TetMeshGeo createVirtualTetsMeshViaCSG(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo triMesh,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool verbose = false, Profiler * profiler = nullptr);

TetMeshGeo createVirtualTetsMeshViaCSG(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo triMesh,
    const std::vector<Vec3ER> & exactTriPos,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool verbose = false, Profiler * profiler = nullptr);

TetMeshGeo createVirtualTetsMeshViaCSG(
    const TetMeshGeo & tetMesh,
    const TriMeshGeo triMesh,
    const std::vector<Vec3ER> & exactTriPos,
    const TetTriIntersectingData & tetTriIntersectingData,
    BarycentricCoordinates * coord = nullptr,
    std::vector<int> * newTetVtxID2OldTetVtxID = nullptr,
    std::vector<int> * newTetID2OldTetID = nullptr,
    std::vector<std::vector<int>> * tetTris = nullptr,
    bool verbose = false, Profiler * profiler = nullptr);

#endif

