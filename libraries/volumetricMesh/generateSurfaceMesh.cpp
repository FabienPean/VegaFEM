/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li                                *
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

#include <float.h>
#include "objMesh.h"
#include "generateSurfaceMesh.h"
#include "cubicMesh.h"
#include "triKey.h"
#include "rectKey.h"
using namespace std;

// the main routine
ObjMesh * GenerateSurfaceMesh::ComputeMesh(const VolumetricMesh * mesh, bool triangulate, bool allElementFaces)
{
  int numElementVertices = mesh->getNumElementVertices();
  int faceDegree = 0;

  if (numElementVertices == 4)
  {
    faceDegree = 3;
    triangulate = false;
  }

  if (numElementVertices == 8)
    faceDegree = 4;

  if (faceDegree == 0)
  {
    printf("Error: unsupported mesh type encountered.\n");
    return nullptr;
  }  

  // create an empty surface mesh
  ObjMesh * objMesh = new ObjMesh;
  objMesh->addGroup(ObjMesh::Group());

  // add all vertices
  for(int i=0; i<mesh->getNumVertices(); i++)
    objMesh->addVertexPosition(mesh->getVertex(i));

  // build unique list of all surface faces

  if (numElementVertices == 4) // tet mesh
  {
    map<OTriKey,int> surfaceFaces;
    for (int i=0; i<mesh->getNumElements(); i++)
    {
      // compute determinant to establish orientation
      double det = dot(mesh->getVertex(i, 1) - mesh->getVertex(i, 0), cross(mesh->getVertex(i, 2) - mesh->getVertex(i, 0), mesh->getVertex(i, 3) - mesh->getVertex(i, 0)));

      auto processFace = [&](int q0, int q1, int q2)
      {
        OTriKey key(mesh->getVertexIndex(i,q0),mesh->getVertexIndex(i,q1),mesh->getVertexIndex(i,q2));
        if (allElementFaces) // get all faces
        {
          objMesh->addFaceToGroup(ObjMesh::Face(key[0], key[1], key[2]), 0);
          return;
        }
        auto it = surfaceFaces.find(key);
        if (it != surfaceFaces.end())
          it->second++;
        else 
        {
          auto revKey = key.getReversedTriKey();
          it = surfaceFaces.find(revKey);  
          if (it != surfaceFaces.end()) 
            it->second--;
          else
            surfaceFaces.emplace(key, 1);
        }
      };
  
      if (det >= 0)
      {
        processFace(1,2,3);
        processFace(2,0,3);
        processFace(3,0,1);
        processFace(1,0,2);
      }
      else
      {
        processFace(3,2,1);
        processFace(3,0,2);
        processFace(1,0,3);
        processFace(2,0,1);
      }
    }

    if (allElementFaces == false) // we build surface mesh
    {
      for(const auto & p : surfaceFaces)
      {
        if (p.second == 0) continue; // inner face
        auto key = p.first;
        int numFaces = abs(p.second);
        if (p.second < 0) key.reverse();
        for(int i = 0; i < numFaces; i++)
        {
          objMesh->addFaceToGroup(ObjMesh::Face(key[0], key[1], key[2]), 0);
        }
      }
    }
  }
  else if (numElementVertices == 8) // cubic mesh
  {
    map<ORectKey,int> surfaceFaces;
    for (int i = 0; i < mesh->getNumElements(); i++)
    {
      auto processFace = [&](int q0, int q1, int q2, int q3)
      {
        ORectKey key(mesh->getVertexIndex(i,q0),mesh->getVertexIndex(i,q1),mesh->getVertexIndex(i,q2),mesh->getVertexIndex(i,q3));
        if (allElementFaces)
        {
          objMesh->addFaceToGroup(ObjMesh::Face(key[0], key[1], key[2], key[3]), 0);
          return;
        }
        auto it = surfaceFaces.find(key);
        if (it != surfaceFaces.end())
          it->second++;
        else 
        {
          auto revKey = key.getReversedRectKey();
          it = surfaceFaces.find(revKey);  
          if (it != surfaceFaces.end()) 
            it->second--;
          else
            surfaceFaces.emplace(key, 1);
        }
      };
  
      processFace(0,3,2,1);
      processFace(4,5,6,7);
      processFace(0,1,5,4);
      processFace(3,7,6,2);
      processFace(1,2,6,5);
      processFace(0,4,7,3);
    }
    if (allElementFaces == false) // we build surface mesh
    {
      for(const auto & p : surfaceFaces)
      {
        if (p.second == 0) continue; // inner face
        auto key = p.first;
        int numFaces = abs(p.second);
        if (p.second < 0) key.reverse();
        for(int i = 0; i < numFaces; i++)
        {
          if (triangulate)
          {
            objMesh->addFaceToGroup(ObjMesh::Face(key[0], key[1], key[2]), 0);
            objMesh->addFaceToGroup(ObjMesh::Face(key[2], key[3], key[0]), 0);
          }
          else
            objMesh->addFaceToGroup(ObjMesh::Face(key[0], key[1], key[2], key[3]), 0);
        }
      }
    }
  }

  if (mesh->getElementType() == CubicMesh::elementType())
  {
    // cubic mesh
    objMesh->setNormalsToFaceNormals();
  }
  else
  {
    // other types of meshes (e.g., tet)
    objMesh->computePseudoNormals();
    objMesh->setNormalsToPseudoNormals();
  }

  objMesh->setSingleMaterial(ObjMesh::Material());

  return objMesh;
}

// advanced routine, not used very often
ObjMesh * GenerateSurfaceMesh::ComputeMesh(const VolumetricMesh * mesh, const ObjMesh * superMesh, bool triangulate)
{
  ObjMesh * surfaceMesh = ComputeMesh(mesh, false);

  // for each volumetric mesh vertex, find the nearest obj file vertex
  vector<int> closestObjVertex(mesh->getNumVertices());
  for(int i=0; i<mesh->getNumVertices(); i++)
  {
    const Vec3d & pos = mesh->getVertex(i);
    double dist;
    closestObjVertex[i] = superMesh->getClosestVertex(pos, &dist);
  }

  // create a new objMesh
  ObjMesh * objMesh = new ObjMesh();
  objMesh->addGroup("Default");

  // build the list of triangles from superMesh
  set<UTriKey> superMeshFaces;
  for(unsigned int i=0; i < superMesh->getNumGroups(); i++)
  {
    const ObjMesh::Group * groupHandle = superMesh->getGroupHandle(i);
    for(unsigned int iFace = 0; iFace < groupHandle->getNumFaces(); iFace++)
    {
      const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(iFace);
      if (faceHandle->getNumVertices() != 3)
      {
        printf("Error: input superMesh is not triangulated.\n");
        return nullptr;
      }

      superMeshFaces.emplace(faceHandle->getVertexPositionIndex(0), faceHandle->getVertexPositionIndex(1), faceHandle->getVertexPositionIndex(2));
    }
  }

  for(unsigned int i=0; i < surfaceMesh->getNumGroups(); i++)
  {
    const ObjMesh::Group * groupHandle = surfaceMesh->getGroupHandle(i);
    for(unsigned int iFace = 0; iFace < groupHandle->getNumFaces(); iFace++)
    {
      const ObjMesh::Face * faceHandle = groupHandle->getFaceHandle(iFace);
      if (faceHandle->getNumVertices() != 3) continue;
      UTriKey key(closestObjVertex[faceHandle->getVertexPositionIndex(0)], 
                  closestObjVertex[faceHandle->getVertexPositionIndex(1)], 
                  closestObjVertex[faceHandle->getVertexPositionIndex(2)]);
      if (superMeshFaces.find(key) == superMeshFaces.end()) continue;
      // only include faces from superMesh
      objMesh->addFaceToGroup(ObjMesh::Face(faceHandle->getVertexPositionIndex(0), 
                                            faceHandle->getVertexPositionIndex(1), 
                                            faceHandle->getVertexPositionIndex(2)), 0);
    }
  }
  delete surfaceMesh;

  return objMesh;
}

