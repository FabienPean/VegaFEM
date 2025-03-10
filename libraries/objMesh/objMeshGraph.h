/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "objMesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC        *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic,                                          *
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

#ifndef VEGAFEM_OBJMESHGRAPH_H
#define VEGAFEM_OBJMESHGRAPH_H

/*
   Jernej Barbic, CMU, MIT, USC, 2007-2012

   A graph where the nodes are obj mesh vertices, edges and faces.
   Two nodes are connected if they are adjacent in the mesh.

   There is also a static function, "GenerateVertexGraph", which computes a graph where the nodes
   are obj mesh vertices, and two nodes are connected if they are adjacent in the mesh.
*/

#include "triple.h"
#include "graph.h"
#include "objMesh.h"

namespace vegafem
{

class ObjMeshGraph : public Graph
{
public:
  ObjMeshGraph(ObjMesh * objMesh); // will be triangulated if necessary

  int GetVertexID(int meshVertex);
  int GetEdgeID(int meshEdge);
  int GetFaceID(int meshFace);

  //void meshID(int graphVertex, int & siteType, int & siteID);
  // converts graph vertex to the integer indices of the particular site
  // sitetpye: 0 = vtx, 1 = edge, 2 = face
  void meshID(int graphVertex, int & siteType, int meshVtxData[3]);

  // converts distance-field style (face, siteIndex) pair to the graph vertex ID
  // siteIndex:
  //  0: vertex0
  //  1: vertex1
  //  2: vertex2
  //  3: edge among 01
  //  4: edge among 12
  //  5: edge among 20
  //  6: the face itself
  int graphID(int faceID, int siteIndex);

  // Creates a graph where the nodes are obj mesh vertices. Two nodes are connected if they are adjacent in the mesh.
  // Face clique:
  // 0: connect every vertex to the next and previous vertex along the face perimeter (but not to other vertices of that face)
  // 1: connect every vertex to all other face vertices
  static Graph * GenerateVertexGraph(const ObjMesh * objMesh, int faceClique=0);

protected:
  int nObj, eObj, fObj;
  std::vector<std::pair<int, int> > meshEdgesVector; // maps each integer identifier to the corresponding edge
  std::vector<triple<int, int, int> > meshFaceVerticesVector; // for every mesh triangle, gives indices of its vertices
  std::vector<triple<int, int, int> > meshFaceEdgesVector; // for every mesh triangle, gives indices of its edges
  ObjMesh * objMesh;
};

inline int ObjMeshGraph::GetVertexID(int vertex)
{
  return vertex;
}

inline int ObjMeshGraph::GetEdgeID(int edge)
{
  return nObj + edge;
}

inline int ObjMeshGraph::GetFaceID(int face)
{
  return nObj + eObj + face;
}


}//namespace vegafem

#endif


