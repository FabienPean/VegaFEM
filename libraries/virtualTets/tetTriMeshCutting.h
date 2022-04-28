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

#ifndef TETTRIMESHCUTTING_H
#define TETTRIMESHCUTTING_H

#include "triMeshGeo.h"
#include "tetMeshGeo.h"
#include "tetKey.h"
#include "triKey.h"
#include "vec3ER.h"
#include "tetTriCutFeature.h"

/*
  Our implementation of the Sutherland-Hodgman algorithm for tet vs tri intersection.
  We assume input tet is in relatively good quality, and input triangle mesh has no degenerate triangles
*/

// exact intersection and cutting between a tet mesh and a triangle mesh
class TetTriMeshCutting
{
public:
  TetTriMeshCutting(const TetMeshRef & tetMesh, TriMeshRef triMesh);

  // compute which tet intersects which triangle, and which triangle intersects which tet
  // using the exact geometry predicates
  // call this function before using getTrianglesIntersectingTet() or computeCutTriangles()
  // you can get the result of this computation from getTrianglesIntersectingTet()
  // note that, if tet and triangle only intersect at boundary, i.e. only touch each other,
  // we still count them as intersecting
  // note 2: triangle completely inside a tet is considered as intersecting
  void computeIntersectingTriTets();
  const std::vector<int> & getTrianglesIntersectingTet(int tetID) const { return tet2tri[tetID]; }

  // compute the triangles cut by all the tets from the tet mesh
  // call this before calling saveCutTriMesh(), exportCutTriMesh and other functions below
  // note that, if tet and triangle only intersect at boundary, i.e. only touch each other,
  // we DON'T count the triangle as cut by the tet
  // optional: provide triangle mesh vtx positions in exact arithmetric
  // note, we use multi-threading to accelerate exact arithmetic operations. This REQUIRES
  // input Vec3ERs in exactTriVtxPos have no dependencies with each other
  // If all Vec3ERs are built directly from Vec3ds, then there is no dependecies
  // However, if one Vec3ER is built by, e.g. the summation of the other two Vec3ER in exactTriVtxPos, then
  // there is a dependency between Vec3ERs in exactTriVtxPos and multi-threading will crash!
  void computeCutTriangles(const std::vector<Vec3ER> * exactTriVtxPos = nullptr);

  // save the cut triangle mesh to disk
  void saveCutTriMesh(std::string filename) const;
  // export the cut triangle mesh
  TriMeshGeo exportCutTriMesh() const;

  const std::vector<Vec3d> & cutTriPositions() const { return cutTriVtxPos; }
  const std::vector<TetTriCutFeature> & cutVertexFeatures() const { return cutVtxFeatures; }

  const std::vector<Vec3ER> & cutTriPositionsER() const { return cutTriVtxPosER; }
  const std::vector<Vec3ER> & tetPositionsER() const { return tetVtxPosER; }

  struct CutTriGroup
  {
    std::vector<Vec3i> tri; // cut triangles with cut vtxID
    std::vector<int> oriID; // original uncut triangle ID for each cut triangle
    std::map<UTriKey,std::vector<int>> cutTriIDsOnFace; // record degenerate case where cut triangles are completely on a tet face
                                                        // store sorted triangle IDs on each tet face if they exist
  };

  const CutTriGroup & getCutTriGroup(int tetID) const { return cutTrisInTet[tetID]; }

protected:
  TetMeshRef tetMesh; // input tet mesh
  TriMeshRef triMesh; // input tri mesh

  std::vector<Vec3ER> tetVtxPosER; // exact arithmetic tet vtx pos

  std::vector<std::vector<int> > tri2tet; // tri ID -> intersecting tets ID
  std::vector<std::vector<int> > tet2tri; // tet ID -> sorted intersecting tri ID

  // positions and TetTriCutFeatures of cut tri vtx
  std::vector<Vec3d> cutTriVtxPos;
  std::vector<Vec3ER> cutTriVtxPosER;
  std::vector<TetTriCutFeature> cutVtxFeatures;

  std::vector<CutTriGroup> cutTrisInTet;

};

#endif
