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

#include "createTriMesh.h"
#include "triMeshGeo.h"
#include "triKey.h"
#include "basicAlgorithms.h"
#include <iostream>
#include <cstring>
#include <cassert>
#include <fstream>
#include <cmath>
#include <functional>
#include <map>
using namespace std;

TriMeshGeo createBoxMesh(const Vec3d & bmin, const Vec3d & bmax)
{
// vtx order in box:
//
//     3 - - - 2
//    /|      /|
//   7 - - - 6 |       y
//   | |     | |       |
//   | 0 - - | 1       |_ _ _x
//   |/      |/       /
//   4 - - - 5       z
  vector<Vec3d> positions = { bmin, { bmax[0], bmin[1], bmin[2] }, { bmax[0], bmax[1], bmin[2] }, { bmin[0], bmax[1], bmin[2] },
                            { bmin[0], bmin[1], bmax[2] }, { bmax[0], bmin[1], bmax[2] }, bmax, { bmin[0], bmax[1], bmax[2] } };
  vector<Vec3i> triangles = { {0,3,2}, {0,2,1}, {4,5,6}, {4,6,7}, {0,1,5}, {0,5,4},
                              {3,7,6}, {3,6,2}, {1,2,6}, {1,6,5}, {0,4,7}, {0,7,3} };
  return TriMeshGeo(move(positions), move(triangles));
}

TriMeshGeo createCubeMesh(const Vec3d & center, double size)
{
  Vec3d side(size/2.0);
  return createBoxMesh(center - side, center + side);
}

TriMeshGeo createTetSurfaceMesh(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2, const Vec3d & v3)
{
  vector<Vec3d> positions = { v0, v1, v2, v3 };
  vector<Vec3i> triangles = { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
  return TriMeshGeo(move(positions), move(triangles));
}

TriMeshGeo createSingleTriangleMesh(const TriMeshRef mesh, int triID)
{
  vector<Vec3d> positions = { mesh.pos(triID, 0), mesh.pos(triID, 1), mesh.pos(triID, 2) };
  vector<Vec3i> triangles = { { 0, 1, 2 } };
  return TriMeshGeo(move(positions), move(triangles));
}

TriMeshGeo createSingleTriangleMesh(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2)
{
  vector<Vec3d> positions = { v0, v1, v2 };
  vector<Vec3i> triangles = { { 0, 1, 2 } };
  return TriMeshGeo(move(positions), move(triangles));
}

// create a cylinder without cap, centered at origin and axis at y-direction
TriMeshGeo createCylinderWallMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight)
{
  assert(subdivisionAxis >= 1);
  assert(subdivisionHeight >= 1);

  vector<Vec3d> positions(subdivisionAxis * (subdivisionHeight+1));
  vector<Vec3i> triangles;
  Vec3d x(1,0,0), nz(0,0,-1);
  double da = 2 * M_PI / subdivisionAxis;
  double dh = height / subdivisionHeight;
  for(int i = 0; i < subdivisionHeight+1; i++)
  {
    Vec3d base(0, -height/2.0 + dh*i, 0);
    for(int j = 0; j < subdivisionAxis; j++)
    {
      Vec3d p = base + radius * Vec3d(cos(j*da), 0, -sin(j*da));
      positions[i*subdivisionAxis + j] = p;
    }
  }

  for(int i = 0; i < subdivisionHeight; i++)
  {
    for(int j = 0; j < subdivisionAxis; j++)
    {
      int nj = (j+1 < subdivisionAxis ? j+1 : 0);
      int a = i*subdivisionAxis+j;
      int b = i*subdivisionAxis+nj;
      int c = (i+1)*subdivisionAxis+nj;
      int d = (i+1)*subdivisionAxis+j;
      triangles.emplace_back(a, b, c);
      triangles.emplace_back(a, c, d);
    }
  }

  return TriMeshGeo(move(positions), move(triangles));
}

// create a cylinder mesh centered at origin and axis at y-direction
TriMeshGeo createCylinderMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight)
{
  TriMeshGeo mesh = createCylinderWallMesh(radius, height, subdivisionAxis, subdivisionHeight);
  int lowerCenterID = mesh.numVertices();
  int upperCenterID = mesh.numVertices() + 1;
  mesh.addPos(Vec3d(0, -height/2.0, 0));
  mesh.addPos(Vec3d(0, height/2.0, 0));
  for(int i = 0; i < subdivisionAxis; i++)
  {
    int ni = (i+1)%subdivisionAxis;
    mesh.addTri(Vec3i(ni, i, lowerCenterID));
  }

  int offsetToUpper = subdivisionAxis * subdivisionHeight;
  for(int i = 0; i < subdivisionAxis; i++)
  {
    int ni = (i+1)%subdivisionAxis;
    mesh.addTri(Vec3i(offsetToUpper+i, offsetToUpper+ni, upperCenterID));
  }

  return mesh;
}
