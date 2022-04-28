/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "objMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC        *
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

#include "objMesh.h"
using namespace std;

// create a cylinder without cap, centered at origin and axis at y-direction
ObjMesh createCylinderWallObjMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight)
{
  assert(subdivisionAxis >= 1);
  assert(subdivisionHeight >= 1);

  ObjMesh mesh;
  Vec3d x(1,0,0), nz(0,0,-1);
  double da = 2 * M_PI / subdivisionAxis;
  double dh = height / subdivisionHeight;
  for(int i = 0; i < subdivisionHeight+1; i++)
  {
    Vec3d base(0, -height/2.0 + dh*i, 0);
    for(int j = 0; j < subdivisionAxis; j++)
    {
      Vec3d p = base + radius * Vec3d(cos(j*da), 0, -sin(j*da));
//      positions[i*subdivisionAxis + j] = p;
      mesh.addVertexPosition(p);
    }
  }

  ObjMesh::Group group;
  for(int i = 0; i < subdivisionHeight; i++)
  {
    for(int j = 0; j < subdivisionAxis; j++)
    {
      int nj = (j+1 < subdivisionAxis ? j+1 : 0);
      int a = i*subdivisionAxis+j;
      int b = i*subdivisionAxis+nj;
      int c = (i+1)*subdivisionAxis+nj;
      int d = (i+1)*subdivisionAxis+j;
      ObjMesh::Face face(a, b, c, d);
      group.addFace(move(face));
    }
  }

  mesh.addGroup(move(group));
  return mesh;
}

// create a cylinder mesh centered at origin and axis at y-direction
ObjMesh createCylinderObjMesh(double radius, double height, int subdivisionAxis, int subdivisionHeight)
{
  ObjMesh mesh = createCylinderWallObjMesh(radius, height, subdivisionAxis, subdivisionHeight);
  int lowerCenterID = mesh.getNumVertices();
  int upperCenterID = lowerCenterID + 1;
  mesh.addVertexPosition(Vec3d(0, -height/2.0, 0));
  mesh.addVertexPosition(Vec3d(0, height/2.0, 0));
  auto & group = mesh.getGroup(0);
  for(int i = 0; i < subdivisionAxis; i++)
  {
    int ni = (i+1)%subdivisionAxis;
    ObjMesh::Face face(ni, i, lowerCenterID);
    group.addFace((move(face)));
  }

  int offsetToUpper = subdivisionAxis * subdivisionHeight;
  for(int i = 0; i < subdivisionAxis; i++)
  {
    int ni = (i+1)%subdivisionAxis;
    ObjMesh::Face face(offsetToUpper+i, offsetToUpper+ni, upperCenterID);
    group.addFace((move(face)));
  }

  return mesh;
}

ObjMesh createBoundingBoxObjMesh(Vec3d bmin, Vec3d bmax)
{
  for(int i = 0; i < 3; i++)
    assert(bmin[i] <= bmax[i]);

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
  vector<Vec4i> faces = { {0,3,2,1}, {4,5,6,7}, {0,1,5,4}, {3,7,6,2}, {1,2,6,5}, {0,4,7,3} };

  ObjMesh mesh;
  mesh.addVertexPositions(positions);
  mesh.addGroup(ObjMesh::Group());
  for(size_t i = 0; i < faces.size(); i++)
    mesh.addFaceToGroup(ObjMesh::Face(faces[i]), 0);
  return mesh;
}
