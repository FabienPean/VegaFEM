/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Virtual tets" driver application,                                    *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2018 USC                           *
 *                                                                       *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*  
  This is an example driver for generating a virtualized tet mesh 
  based on the given input tet mesh and input manifold nearly self-intersecting triangle mesh.
  It calls routines from the virtualTets library.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <vector>
#include "virtualTets.h"
#include "virtualTets-via-csg.h"
#include "commandLineParser.h"
#include "filterIterator.h"
#include "predicates.h"
#include "exactOctree.h"
#include "tetMesh.h"
#include "objMesh.h"
#include "profiler.h"
using namespace std;

int main(int argc, char ** argv)
{
  initPredicates();

  int numFixedArgs = 4;
  if (argc < numFixedArgs) 
  {
    cout << "Usage: " << argv[0] << " <tet mesh> <tri mesh> <output virtualized tet mesh> -w <output barycentric weight file>" << endl;
    return 0;
  }

  char * tetMeshFilename = argv[1];
  char * triMeshFilename = argv[2];
  char * outputFilename = argv[3];
  string weightFilename;

  CommandLineParser parser;
  parser.addOption("w", weightFilename);

  int ret = parser.parse(argc, argv, 4);
  if (ret != argc) 
  {
    cout << "Failure parsing option: " << argv[ret] << endl;
    return 1;
  }

  // load the tet mesh and convert it to "TetMeshGeo" object
  TetMesh vegaTetMesh(tetMeshFilename);
  int numVertices, numTets, numEleVtx;
  double * vertices;
  int * tets;
  vegaTetMesh.exportMeshGeometry(&numVertices, &vertices, &numTets, &numEleVtx, &tets);
  TetMeshGeo tetMesh(numVertices, vertices, numTets, tets);
  free(vertices);
  free(tets);

  // load the triangle mesh and convert it to "TriMeshGeo" object
  ObjMesh objMesh(triMeshFilename);
  int numTriangles;
  int * triangles;
  objMesh.exportGeometry(&numVertices, &vertices, &numTriangles, &triangles);
  TriMeshGeo triMesh(numVertices, vertices, numTriangles, triangles);
  free(vertices);
  free(triangles);

  // check if the input mesh is self-intersecting
  ExactTriMeshOctree triMeshOctree;
  int maxDepth = 5;
  int maxNumTrianglesPerNode = 10;
  triMeshOctree.build(triMesh, maxDepth, maxNumTrianglesPerNode);
  vector<pair<int,int>> selfIintersectVector;
  triMeshOctree.selfIntersectionExact(triMesh, selfIintersectVector);
  if (selfIintersectVector.size() > 0) 
  {
    cout << "Error: the input mesh is self-intersecting." << endl;
    return 1;
  }

  // perform the virtualization
  BarycentricCoordinates bc;
  TetMeshGeo newTetMesh;
  try 
  {
    newTetMesh = createVirtualTetsMesh(tetMesh, triMesh, &bc);
  } catch(int) 
  {
    cout << "Failed to run virtual tets algorithm" << endl;
    return 1;
  }
  assert(newTetMesh.numTets() > 0);

  // save to disk
  newTetMesh.save(outputFilename);
  if (weightFilename.size() > 0) 
    bc.saveInterpolationWeights(weightFilename);

  return 0;
}

