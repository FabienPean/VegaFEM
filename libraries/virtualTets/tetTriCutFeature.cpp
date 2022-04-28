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

#include "tetTriCutFeature.h"
#include <cassert>
#include <array>
using namespace std;

TetShape::TetShape(const TetMeshRef & tetMesh, int tetID) : tetID(tetID)
{
  Vec4i tet = tetMesh.tet(tetID);
  for(int v : tet)
  {
    position[v] = tetMesh.pos(v);
  }
  assert(position.size() == 4);

  OTetKey tetkey(tet);
  for(int i = 0; i < 4; i++)
  {
    auto otrikey = tetkey.oFaceKey(i);
    int v0 = otrikey[0], v1 = otrikey[1], v2 = otrikey[2];
    Vec3d te1 = tetMesh.pos(v1) - tetMesh.pos(v0);
    Vec3d te2 = tetMesh.pos(v2) - tetMesh.pos(v0);
    Vec3d normal = cross(te1, te2);
    normal.normalize();
    assert(normal.hasNaN() == false);
    faceNormal[otrikey.uTriKey()] = normal;
  }
  assert(faceNormal.size() == 4);
}


bool TetShape::hasFace(const UTriKey & face) const
{
  return faceNormal.find(face) != faceNormal.end();
}

Vec3d TetShape::getNormal(const UTriKey & face) const
{
  auto iter = faceNormal.find(face);
  assert(iter != faceNormal.end());
  return iter->second;
}

Vec3d TetShape::getFaceCenter(const UTriKey & face) const
{
  auto iter = faceNormal.find(face);
  assert(iter != faceNormal.end());
  Vec3d center(0.0);
  for(int i = 0; i < 3; i++)
    center += getPos(face[i]);
  return center / 3.0;
}

Vec3d TetShape::getTetCenter() const
{
  Vec3d center(0.0);
  for(const auto & p : position)
    center += p.second;
  return center / 4.0;
}

Vec3d TetShape::getPos(int vtxID) const
{
  auto iter = position.find(vtxID);
  assert(iter != position.end());
  return iter->second;
}

array<UTriKey, 2> TetShape::getNeighboringFace(const UEdgeKey & edge) const
{
  array<UTriKey, 2> ret;
  int sz = 0;
  for(const auto & p : faceNormal)
  {
    if (p.first.hasUEdge(edge))
    {
      ret[sz++] = p.first;
    }
  }
  if (sz != 2)
  {
    cout << "Error in TetShape: ";
    for(auto p : faceNormal) { cout << p.first << " "; }
    cout << " querying edge: " << edge << endl;
  }
  assert(sz == 2);
  return ret;
}

array<UTriKey, 3> TetShape::getNeighboringFace(int vtxID) const
{
  array<UTriKey, 3> ret;
  int sz = 0;
  for(const auto & p : faceNormal)
  {
    if (p.first.hasIndex(vtxID)) { ret[sz++] = p.first; }
  }
  assert(sz == 3);
  return ret;
}

vector<UTriKey> TetTriCutFeature::getTouchingTetFaces(const TetShape & tetShape) const
{
  if (isInsideTet()) return {};
  if (isTetFace()) return { getTetFace() };
  if (isTetEdge())
  {
    auto faces = tetShape.getNeighboringFace(getTetEdge());
    assert(faces[0].isValidTriangle() && faces[1].isValidTriangle());
    return { faces[0], faces[1] };
  }
  // tet vtx
  auto faces = tetShape.getNeighboringFace(getTetVertex());
  return { faces[0], faces[1], faces[2] };
}
