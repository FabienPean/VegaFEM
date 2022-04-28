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

#ifndef TETTRIPIECE_H
#define TETTRIPIECE_H

#include "tetTriMeshCutting.h"
#include "triMeshPseudoNormal.h"
#include "triMeshNeighbor.h"
//#include "triMeshOctreeK.h"
#include "tetKey.h"
#include "tetTriCuttingData.h"
#include "profiler.h"
#include <vector>
#include <map>

/*
  Represents one "piece" (a contiguous parts of a triangle belonging to one tet), as defined in our paper submission.
*/

struct TetTriPiece
{
  using CutTriGroup = TetTriMeshCutting::CutTriGroup;
  int tetID = -1;
  std::vector<int> groupTriID; // the triangle IDs in the CutTriGroup passed in the constructor
  std::vector<Vec3i> gnrCutTri;  // cut triangles (with cut vtxIDs) of this piece which are in general positions
                                 // general position means that the triangles don't completely lie on tet faces
  std::vector<int> vtx;    // cut vtxIDs inside tri
  std::vector<int> gnrOriTriID; // original, input, uncut triangle IDs where triangles in gnrCutTri are from
//  ExactTriMeshOctreeK octree;
  TriangleNeighbor gnrTriNbr;
  std::set<UTriKey> touchedTetFaces; // these tet faces are intersected by triangles in tri
  signed char inOutWhenAllTriOnTetFaces = 0; // the inOut value if all triangles are exactly on tet faces

  TetTriPiece(int tet); // constructor for a tet with no pieces
  // triIDs: the triangle IDs in the input CutTriGroup group
  TetTriPiece(int tet, const std::vector<int> & triIDs, const CutTriGroup & group, const TetTriCuttingData & cutting,
    const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape);

  //  void buildTree(const std::vector<Vec3ER> & cutPositions);
  // get the index of the closest triangle in the piece
  // also returns: feature: which site (edge/vtx/interior) the closest point on the triangle lies
  //               dist2: the squared distance between the query point and the closest site
  int getClosestTriangle(const TetTriCuttingData & cutting, const Vec3d & queryPoint, const Vec3ER & queryPointK, int & feature, ER & dist2);

  // build data for computing pseudo-normal on piece boundaries during inside/outside test
  // we don't need to use CSG to create a close manifold mesh to do inside/outside test
  // But we have to find the correct pseudo-normals on piece boundaries, which are the pseudo-normals of the meshes which
  // would be created by CSG
  void buildBoundaryData(const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape,
      Profiler * profiler = nullptr);
  void tryBuildingBoundaryData(const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape,
      Profiler * profiler = nullptr);

  int inOutTest(const Vec3d & queryPoint, const Vec3ER & queryPointK, const TetTriCuttingData & cuttingData,
      const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape, Profiler * profiler = nullptr);

  TriMeshRef mesh(const TetTriCuttingData & cut) const;

protected:
  //  bool treeBuilt = false;
  bool initializeBoundaryData = false; // flag to represent whether we have boundary data computed
  std::map<OEdgeKey, std::pair<double, UTriKey>> edgeBouAngle; // piece boundary edge -> angle between the tet face, tet face
  std::multimap<int, std::pair<double, Vec3d>> vtxBouNormal; // boundary vtx weighted normal on tet faces
};

TetTriPiece mergePiece(const TetTriPiece & comp0, const TetTriPiece & comp1, int tet, const TetTriMeshCutting::CutTriGroup & group,
    const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape);

#endif

