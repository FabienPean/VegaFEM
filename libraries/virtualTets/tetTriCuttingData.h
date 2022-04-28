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

#ifndef TETTRICUTTINGDATA_H
#define TETTRICUTTINGDATA_H

#include "tetTriMeshCutting.h"
#include <vector>
#include <mutex>

/*
  Datastructure to use our virtual tets in the immersion algorithm.
*/

struct TetTriCuttingData
{
  // Definition:
  //   cut triangle mesh: the new triangle mesh which is the result of the cut on the input triangle mesh by all the tets
  //     The newly geneated cuts are modeled by adding new vertices on the cut triangles.
  //     New vertices are shared by neighboring cut triangles. So after the cut, the topology of the input triangle mesh
  //     is preserved. Only its surface is subdivided by the tets' cut.
  //   cut feature: a pair of tetKey and triKey, used to store which geometry feature
  //     (tet: interior/face/edge/vtx, tri: interior, edge, vtx) a cut vertex lies on
  // TetTriMeshCutting::cutTriGroups: tetID -> cut triangle info in each tet
  //   cut triangle info includes: 1) the 3 cutVtxID of the triangle, 2) the original triangle ID where this cut triangle is from
  std::vector<Vec3d> cutVtxPos;           // cut triangle mesh vertex positions
  std::vector<Vec3ER> cutVtxPosER;        // cut triangle mesh vertex positions in exact arithmetic
  std::vector<TetTriCutFeature> features; // cut feature of each cutVtx
  std::vector<Vec3ER> tetPosER;           // tet mesh position in exact arithmetic
  std::vector<TetTriMeshCutting::CutTriGroup> cutTriGroups; // cut triangle info in each tet


  // used to protect ER because they are not thread-safe.
//  mutable bool multiThreading = false;
//  mutable std::mutex * cutVtxPosMutex = nullptr;
//  mutable std::mutex * tetPosMutex = nullptr;
//  mutable std::vector<std::mutex> * tetPosVecMutex = nullptr;
//  mutable std::vector<std::mutex> * cutVtxPosVecMutex = nullptr;
};

struct TetTriIntersectingData
{
  std::vector<std::vector<int>> triInTet;
  std::vector<Vec3ER> tetPosER;            // exact real for tet pos
};
#endif
