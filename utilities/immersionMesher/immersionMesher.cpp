/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Immersion mesher" driver application,                                *
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
  This is an example driver for immersion meshing, i.e., generating a tet mesh with 
  properly duplicated tets to mesh the space occupied by a self-intersecting input triangle mesh.
  It calls routines from the immersionMeshing library.
*/

#include "containerHelper.h"
#include "sceneObjectDeformable.h"
#include "getopts.h"
#include "initPredicates.h"
#include "lighting.h"
#include "listIO.h"
#include "tetKey.h"
#include "basicAlgorithms.h"
#include "camera.h"
#include "performanceCounter.h"
#include "inputDevice.h"
#include "commandLineParser.h"
#include "saveScreenShot.h"
#include "vec4d.h"
#include "openGLHelper.h"
#include "immersionMesher.h"
#include "stringHelper.h"
#ifdef WIN32
  #include <windows.h>
#endif
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <climits>
#include <iostream>
#include <sstream>
#include <fstream>
#include <memory>
#include <math.h>
#include <time.h>
#include <GL/glui.h>

//#define USE_GPERF_TOOLS
#ifdef USE_GPERF_TOOLS
  #include "gperftools/profiler.h"
#endif

using namespace std;

#define PRINT_SWITCH(sw) printf(#sw " is %s\n", (sw) ? "ON" : "OFF")
#define PRINT_VALUE(value) cout << #value << " is now " << value << endl

const char * surfaceMeshFilename = nullptr;
const char * tetMeshFilename = nullptr;

static ObjMesh * uncutSurfaceMesh = nullptr;
static ObjMesh * cutSurfaceMesh = nullptr;

static ImmersionMesher im;
static PerformanceCounter counter;

static string finalTetMeshFilename = "finalTet.veg";
static string finalBarycentricWeightsFilename = "finalTet.interp";

static void runMesher()
{
  initPredicates(); // Initialize exact arithmetics.
  srand(time(nullptr)); // Initialize random number generator.

  // Load the input triangle mesh, and create the input "TriMeshGeo" object.
  uncutSurfaceMesh = new ObjMesh(surfaceMeshFilename);
  // Check if the input mesh has isolated vertices (there should be none).
  int numVertices = uncutSurfaceMesh->getNumVertices();
  uncutSurfaceMesh->removeIsolatedVertices(); 
  if ((int)uncutSurfaceMesh->getNumVertices() != numVertices)
  {
    cout << "Error: input mesh has isolated vertices" << endl;
    exit(1);
  }
  vector<Vec3i> uncutSurfaceTriangles;
  uncutSurfaceMesh->exportTriangles(uncutSurfaceTriangles);
  TriMeshGeo triMesh(uncutSurfaceMesh->getNumVertices(), &uncutSurfaceMesh->getPosition(0), move(uncutSurfaceTriangles));

  // Load the input tet mesh, and create the input "TetMeshGeo" object.
  TetMeshGeo tetMesh;
  {
    assert(tetMeshFilename);
    TetMesh inputTetMesh(tetMeshFilename);
    int numVertices, numTets, numEleVtx;
    double * vertices;
    int * tets;
    inputTetMesh.exportMeshGeometry(&numVertices, &vertices, &numTets, &numEleVtx, &tets);
    tetMesh = TetMeshGeo(numVertices, vertices, numTets, tets);
    free(vertices);
    free(tets);
  }

  PerformanceCounter pc;

  // Create the output datastructures.
  vector<TetMeshGeo> tetMeshes; // output tet meshes (see the immersionMesher library)
  vector<BarycentricCoordinates> weights; // output barycentric embedding weights (see the immersionMesher library)
  //vector<vector<TriMeshGeo>> allCellMeshes; // uncomment this if you also want to export the cell geometry (& need to then pass it into "run")
  //vector<vector<BarycentricCoordinates>> allCellWeights;
  // Run the mesher.
  im.run(triMesh, tetMesh, tetMeshes, weights);

  // cout << "Performance information:" << endl;
  // cout << im.getProfiler().toString() << endl;

  // Output the results to disk.
  for(int i = 0; i < sizei(tetMeshes); i++)
  {
    // Create filenames for the output meshes.
    string graphIDStr = "." + to_string(i);
    string finalTetFilename = finalTetMeshFilename;
    if (sizei(tetMeshes) > 1)
    {
      if (iendWith(finalTetFilename, ".veg"))
        finalTetFilename = finalTetFilename.substr(0, finalTetFilename.size()-4) + graphIDStr + ".veg";
      else
        finalTetFilename = finalTetFilename + graphIDStr;
    }
    string interpFilename = finalBarycentricWeightsFilename;
    if (sizei(tetMeshes) > 1)
    {
      if (iendWith(interpFilename, ".interp"))
        interpFilename = interpFilename.substr(0, interpFilename.size()-7) + graphIDStr + ".interp";
      else
        interpFilename = interpFilename + graphIDStr;
    }
    //string allCellFilename = "allCellMesh.obj";
    //if (sizei(tetMeshes) > 1) 
    //  allCellFilename = "allCellMesh" + graphIDStr + ".obj";
    //string allInterpFilename = "allCellMesh.interp";
    //if (sizei(tetMeshes) > 1) 
    //  allInterpFilename = "allCellMesh" + graphIDStr + ".interp";

    // Save the outputs to disk.
    tetMeshes[i].save(finalTetFilename);
    weights[i].saveInterpolationWeights(interpFilename);
    //allCellMeshes[i].save(allCellFilename);
    //allCellWeights[i].saveInterpolationWeights(allInterpFilename);
  }

  pc.StopCounter();
  cout << "Total  mesher execution time: " << pc.GetElapsedTime() << "s." << endl;

  return;
}

void exitHandler()
{
  delete uncutSurfaceMesh;
  delete cutSurfaceMesh;
}

int main (int argc, char ** argv)
{
  int numArgs = 3;
  if (argc < numArgs)
  {
    cout << "Usage: " << argv[0] << " <obj mesh> <tet mesh> [-o <output tet mesh>] [-S] [-v] [-w <output barycentric weights>]" << endl;
    cout << "Optional: " << endl;
    cout << "  -o: Specify output tet mesh filename. Default: finalTet.veg." << endl;
    cout << "  -S: Use CSG instead of the default implementation of [Li2018]. CSG is slower and equally precise, so [Li2018] is preferred." << endl;
    cout << "  -v: be verbose." << endl;
    cout << "  -w: Specify output filename to store the weights to embed the input obj mesh into the output tet mesh. Default: finalTet.interp." << endl;
    return 0;
  }

  // Parse command-line arguments.

  surfaceMeshFilename = argv[1];
  tetMeshFilename = argv[2];
  
  CommandLineParser parser;
  bool verbose = false;
  bool useCSGForVirtualTets = false;

  parser.addOption("v", verbose);
  parser.addOption("S", useCSGForVirtualTets);
  parser.addOption("o", finalTetMeshFilename);
  parser.addOption("w", finalBarycentricWeightsFilename);

  int ret = parser.parse(argc, argv, numArgs);
  if (ret != argc)
  {
    cout << "Error parsing option " << argv[ret] << endl;
    return 1;
  }
  im.setVerbose(verbose);
  im.useCSGForVirtualTets(useCSGForVirtualTets);

  cout << "Use CSG for virtual tets method ? " << useCSGForVirtualTets << endl;
  cout << "Verbose ? " << verbose << endl;

  #ifdef USE_GPERF_TOOLS
    ProfilerStart("out.prof");
  #endif
    // execute the actual meshing operation, and save the results
    runMesher();
  #ifdef USE_GPERF_TOOLS
    ProfilerStop();
  #endif

  atexit(exitHandler);
  return(0);
}

