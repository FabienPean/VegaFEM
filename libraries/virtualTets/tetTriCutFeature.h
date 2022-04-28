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

#ifndef TETTRICUTFEATURE_H
#define TETTRICUTFEATURE_H

#include "triKey.h"
#include "tetKey.h"
#include "hashHelper.h"
#include "tetMeshGeo.h"
#include <vector>
#include <map>


// to perform geometric queries on a tet
struct TetShape
{
  std::map<UTriKey, Vec3d> faceNormal;
  std::map<int, Vec3d> position; // vtxID in the mesh -> vtx position
  int tetID = -1;

  TetShape(const TetMeshRef & tetMesh, int tetID);

  bool hasFace(const UTriKey & face) const;

  Vec3d getNormal(const UTriKey & face) const;

  Vec3d getPos(int vtxID) const; // vtxID is the global vtx ID on the mesh, NOT local ID: [0,3]
  Vec3d getFaceCenter(const UTriKey & face) const;
  Vec3d getTetCenter() const;

  // get the two tet faces that share the tet edge
  std::array<UTriKey, 2> getNeighboringFace(const UEdgeKey & edge) const;
  // get the three tet faces that sahre the vtx with global ID: vtxID
  std::array<UTriKey, 3> getNeighboringFace(int vtxID) const;
};


// a unique pair of tet and tri, used to represent an intersection vertex between a tet and a triangle
// from a tet mesh and triangle mesh, respectively
struct TetTriCutFeature
{
  TetTriCutFeature() {} // initialize to be a null feature, with triFeature = {-1,-1,-1} and tetFeature = {-1,-1,-1,-1}
  TetTriCutFeature(const UTriKey & tri, const UTetKey & tet) : triFeature(tri), tetFeature(tet) {}
  UTriKey triFeature; // default initialization to {-1,-1,-1}
  UTetKey tetFeature; // default initialization to {-1,-1,-1,-1}
  // the indices stored in triFeature and tetFeature are sorted

  // tetFeature can be {-1,-1,-1, a}, {-1,-1, a, b}, {-1, a, b, c} or {a, b, c, d}, where a,b,c,d >= 0
  // {-1,-1,-1, a} means the represented vertex is a tet mesh vtx with tetVtxID = a
  // {-1,-1, a, b} means the represented vertex is on a tet edge between tetVtxID a and b
  // {-1, a, b, c} means the represented vertex is on the interior of a tet face, whose vertices are a, b and c
  // { a, b, c, d} means the represented vertex is on the interior of a tet, whose vertices are a,b,c,d
  bool isInsideTet() const { return tetFeature[0] >= 0; }
  bool isTetFace() const { return tetFeature[0] < 0 && tetFeature[1] >= 0; }
  bool isTetEdge() const { return tetFeature[1] < 0 && tetFeature[2] >= 0; }
  bool isTetVertex() const { return tetFeature[2] < 0 && tetFeature[3] >= 0; }

  // triFeature can be {-1,-1, a}, {-1, a, b} or {a, b, c}, where a,b,c >= 0
  // {-1,-1, a} means the represented vertex is a triangle vertex with triVtxID = a
  // {-1, a, b} means the represented vertex is on a triangle edge between triVtxID a and b
  // { a, b, c} means the represented vertex is on the interior of a triangle, whose vertices are a, b and c
  bool isInsideTri() const { return triFeature[0] >= 0; }
  bool isTriEdge() const { return triFeature[0] < 0 && triFeature[1] >= 0; }
  bool isTriVertex() const { return triFeature[1] < 0 && triFeature[2] >= 0; }

  UTetKey getTet() const { return tetFeature; }
  UTriKey getTetFace() const { return *(UTriKey*)(&tetFeature[1]); }
  UEdgeKey getTetEdge() const { return *(UEdgeKey*)(&tetFeature[2]); }
  int getTetVertex() const { return tetFeature[3]; }

  UTriKey getTri() const { return triFeature; }
  UEdgeKey getTriEdge() const { return *(UEdgeKey*)(&triFeature[1]); }
  int getTriVertex() const { return triFeature[2]; }

  bool operator == (const TetTriCutFeature & f2) const { return triFeature == f2.triFeature && tetFeature == f2.tetFeature; }
  bool operator != (const TetTriCutFeature & f2) const { return !((*this) == f2); }
  bool operator < (const TetTriCutFeature & f2) const
  {
    if (triFeature < f2.triFeature) return true;
    if (f2.triFeature < triFeature) return false;
    return tetFeature < f2.tetFeature;
  }

  // return the tet faces the represented vertex touches on tetShape
  // if the represented vertex is on the interior of the tet, return {}
  // if it is on the tet face, return only that face
  // if it is on the tet edge, return the two tet faces sharing this edge
  // else, it can only be a tet vtx, return the three tet faces sharing this vtx
  std::vector<UTriKey> getTouchingTetFaces(const TetShape & tetShape) const;
};

namespace std
{
  template <>
  struct hash<TetTriCutFeature>
  {
    size_t operator()(const TetTriCutFeature & k) const
    {
      size_t v = hash<int>()(k.triFeature[0]);
      v = hashCombine(v, hash<int>()(k.triFeature[1]));
      v = hashCombine(v, hash<int>()(k.triFeature[2]));
      v = hashCombine(v, hash<int>()(k.tetFeature[0]));
      v = hashCombine(v, hash<int>()(k.tetFeature[1]));
      v = hashCombine(v, hash<int>()(k.tetFeature[2]));
      v = hashCombine(v, hash<int>()(k.tetFeature[3]));
      return v;
    }
  };
}


#endif
