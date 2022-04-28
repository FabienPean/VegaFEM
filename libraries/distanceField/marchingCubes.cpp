/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "distance field" library , Copyright (C) 2007 CMU, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Danyong Zhao, Jernej Barbic                             *
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

/*
  The labelings of the vertices, edges and faces of the cube follow the one in:

  Thomas Lewiner, Hélio Lopes, Antônio Wilson Vieira, Geovan Tavares:
  Efficient implementation of marching cubes cases with topological guarantees
  Journal of Graphics Tools 8(2): pp. 1-15, 2003
*/

// when a static table is not used, the table is generated from scratch using our code ("createTable()")
// otherwise, the marching cubes table is read from a static array (that was previously generated using "createTable()")
#define USE_STATIC_TABLE

#include <float.h>
#include <memory.h>
#include <string>
#include <algorithm>
#include <set>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <unordered_map>
using namespace std;
#include "vec3d.h"
#include "marchingCubes.h"
#include "performanceCounter.h"

#ifdef USE_TBB
  #include <tbb/tbb.h>
#else
  #include "range.h"
#endif

static const unsigned char edge_vertices[12][2] =
{
  { 0, 1 },
  { 1, 2 },
  { 2, 3 },
  { 0, 3 },
  { 4, 5 },
  { 5, 6 },
  { 6, 7 },
  { 4, 7 },
  { 0, 4 },
  { 1, 5 },
  { 2, 6 },
  { 3, 7 } };

static const unsigned char face_vertices[6][4] =
{
  { 0, 1, 4, 5 },
  { 1, 2, 5, 6 },
  { 2, 3, 6, 7 },
  { 0, 3, 4, 7 },
  { 0, 1, 2, 3 },
  { 4, 5, 6, 7 } };

/*
 The 24 symmetries are:
 Identity
 3 rotations (by +- pi/2 or pi) about centres of 3 pairs of opposite faces. (9)
 1 rotation (by pi) about centres of 6 pairs of opposite edges. (6)
 2 rotations (by +- 2pi/3) about 4 pairs of opposite vertices (diagonals). (8)
*/

// this table was generated manually, by physically manufacturing two cubes, 
// labeling the corners and applying each transformation
static const unsigned char cubeSymmetryVertices[25][8] =
{
    { 0, 1, 2, 3, 4, 5, 6, 7 }, // identity                                     1
    { 3, 2, 6, 7, 0, 1, 5, 4 }, // rotation by +pi/2 about faces 1,3            2
    { 1, 5, 6, 2, 0, 4, 7, 3 }, // rotation by +pi/2 about faces 0,2            3
    { 1, 2, 3, 0, 5, 6, 7, 4 }, // rotation by +pi/2 about faces 4,5            4
    { 4, 5, 1, 0, 7, 6, 2, 3 }, // rotation by -pi/2 about faces 1,3            5
    { 4, 0, 3, 7, 5, 1, 2, 6 }, // rotation by -pi/2 about faces 0,2            6
    { 3, 0, 1, 2, 7, 4, 5, 6 }, // rotation by -pi/2 about faces 4,5            7
    { 7, 6, 5, 4, 3, 2, 1, 0 }, // rotation by pi about faces 1,3               8
    { 5, 4, 7, 6, 1, 0, 3, 2 }, // rotation by pi about faces 0,2               9
    { 2, 3, 0, 1, 6, 7, 4, 5 }, // rotation by pi about faces 4,5               10
    { 1, 0, 4, 5, 2, 3, 7, 6 }, // rotation by pi about edges 0, 6              11
    { 6, 2, 1, 5, 7, 3, 0, 4 }, // rotation by pi about edges 1, 7              12
    { 6, 7, 3, 2, 5, 4, 0, 1 }, // rotation by pi about edges 2, 4              13
    { 3, 7, 4, 0, 2, 6, 5, 1 }, // rotation by pi about edges 3, 5              14
    { 4, 7, 6, 5, 0, 3, 2, 1 }, // rotation by pi about edges 8, 10             15
    { 6, 5, 4, 7, 2, 1, 0, 3 }, // rotation by pi about edges 9, 11             16
    { 0, 4, 5, 1, 3, 7, 6, 2 }, // rotation by +2pi/3 about vertices 0, 6       17
    { 5, 1, 0, 4, 6, 2, 3, 7 }, // rotation by +2pi/3 about vertices 1, 7       18
    { 5, 6, 2, 1, 4, 7, 3, 0 }, // rotation by +2pi/3 about vertices 2, 4
    { 7, 4, 0, 3, 6, 5, 1, 2 }, // rotation by +2pi/3 about vertices 3, 5
    { 0, 3, 7, 4, 1, 2, 6, 5 }, // rotation by -2pi/3 about vertices 0, 6
    { 2, 1, 5, 6, 3, 0, 4, 7 }, // rotation by -2pi/3 about vertices 1, 7
    { 7, 3, 2, 6, 4, 0, 1, 5 }, // rotation by -2pi/3 about vertices 2, 4
    { 2, 6, 7, 3, 1, 5, 4, 0 }, // rotation by -2pi/3 about vertices 3, 5
    { 4, 5, 6, 7, 0, 1, 2, 3 } };

// The positions of the eight vertices
static const unsigned char vertex_position[8][3] =
{
  { 0, 0, 0 },
  { 1, 0, 0 },
  { 1, 1, 0 },
  { 0, 1, 0 },
  { 0, 0, 1 },
  { 1, 0, 1 },
  { 1, 1, 1 },
  { 0, 1, 1 } 
};

// The four corners of the six faces 
static const unsigned char face_vertex[6][4] =
{
  { 0, 1, 5, 4 },
  { 1, 2, 6, 5 },
  { 3, 2, 6, 7 },
  { 0, 3, 7, 4 },
  { 0, 1, 2, 3 },
  { 4, 5, 6, 7 } 
};

// edge_map and edge_vertex: how vertices are placed around an edge
/*
  edge_map
  first index can be 0, 1 or 2. 0 means this edge is parallel to X axis. 1 means parallel to Y axis. 2 means parallel to Z axis.
  second index: which specific place is the edge in the X, Y or Z direction.
  This to indices point to the start position of edge_vertex
*/
static const unsigned char edge_map[12][2] =
{
  { 0, 0 },
  { 1, 0 },
  { 0, 2 },
  { 1, 2 },
  { 0, 6 },
  { 1, 6 },
  { 0, 4 },
  { 1, 4 },
  { 2, 0 },
  { 2, 6 },
  { 2, 4 },
  { 2, 2 } };

/*
  edge_vertex: how the 8 vertices are placed around an edge
  The start vertex is determined by edge_map, and then pick up 8 consecutive vertices.
  The first two vertices are the two endpoints of this edge.
  The next 6 vertices make up 3 pairs(two consecutive vertices make up one pair). Each pair is the two endpoints of an edge whose direction is the same of this edge.
  The four edges are in clockwise or counter clockwise (depending on how you look).
*/

static const unsigned char edge_vertex[3][14] =
{
  { 0, 1, 3, 2, 7, 6, 4, 5, 0, 1, 3, 2, 7, 6 },  // X direction
  { 1, 2, 0, 3, 4, 7, 5, 6, 1, 2, 0, 3, 4, 7 },  // Y direction
  { 0, 4, 3, 7, 2, 6, 1, 5, 0, 4, 3, 7, 2, 6 }   // Z direction
};

#ifdef USE_STATIC_TABLE
  #include "marchingCubesTable.h"
#else
  static char cubeSymmetryEdges[25][13]; // determined algorithmically from cubeSymmetryVertices

  static char cubeSymmetryFaces[24][6]; // determined algorithmically from cubeSymmetryVertices

  // if first entry is negative, this is a representative case, second entry has no meaning (it is always set to 0)
  // if first entry i is non-negative, this is a symmetry and/or complement of the representative case i; second entry specifies the symmetry and/or complement;
  // if second entry j is positive, this is a symmetry using the permutation j from cubeSymmetryVertices; 
  // if second entry j is zero or negative,  this is symmetry |j|, plus complement
  // can determine algorithmically from cubeSymmetryVertices
  static int marchingCubeSymmetries[256][3];

  //used to build the look up table for case 13 mentioned in Chernyaev’s paper
  const int oppositeFace[3][2] = {{1, 3}, {2, 4}, {5, 6}};
  const int vtxAdjacentFaces[8][3] = {{1,4,5}, {1,2,5}, {2,3,5}, {3,4,5}, {1,4,6}, {1,2,6}, {2,3,6}, {3,4,6}};

  static const unsigned char caseNumber[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 9, 14, 12, 10, 13 };

  static const char representativeCase1[3] = { 8, 3, 0 };
  static const char representativeCase2[6] = { 3, 1, 8, 9, 8, 1 };
  static const char representativeCase3_1[6] = { 1, 2, 10, 8, 3, 0 };
  static const char representativeCase3_2[12] = { 8, 3, 10, 10, 1, 0, 0, 8, 10, 2, 10, 3 };
  static const char representativeCase4_1[6] = { 6, 5, 10, 3, 0, 8 };
  static const char representativeCase4_2[18] = { 8, 6, 5, 8, 5, 0, 6, 3, 10, 0, 5, 10, 0, 10, 3, 3, 6, 8 };
  static const char representativeCase5[9] = { 8, 3, 2, 10, 8, 2, 10, 9, 8 };
  static const char representativeCase6_1_1[9] = { 6, 5, 10, 8, 1, 9, 8, 3, 1 };
  static const char representativeCase6_2[15] = { 8, 6, 5, 10, 3, 1, 8, 3, 6, 9, 8, 5, 6, 3, 10 };
  static const char representativeCase7_1[9] = { 1, 2, 10, 9, 5, 4, 8, 3, 0 };
  static const char representativeCase7_2_0[15] = { 1, 2, 10, 4, 8, 3, 5, 4, 3, 0, 5, 3, 5, 0, 9 };
  static char representativeCase7_2_1[15]; // determined algorithmically
  static char representativeCase7_2_2[15]; // determined algorithmically
  static const char representativeCase7_3_0[27] = { 12, 1, 2, 12, 9, 1, 5, 12, 10, 9, 12, 0, 12, 8, 3, 12, 5, 4, 0, 12, 3, 10, 12, 2, 4, 8, 12 };
  static char representativeCase7_3_1[27]; // determined algorithmically
  static char representativeCase7_3_2[27]; // determined algorithmically
  static const char representativeCase7_4_1[15] = { 5, 4, 8, 3, 2, 8, 9, 1, 0, 2, 5, 8, 5, 2, 10 };
  static const char representativeCase8[6] = { 10, 9, 8, 8, 11, 10 };
  static const char representativeCase9[12] = { 3, 5, 8, 3, 10, 5, 2, 10, 3, 8, 5, 4 };
  static const char representativeCase10_1_1[12] = { 11, 7, 10, 3, 1, 8, 1, 9, 8, 7, 5, 10 };
  static char representativeCase10_1_1_[12]; // determined algorithmically
  static const char representativeCase10_1_2[24] = { 11, 3, 10, 9, 7, 5, 9, 5, 10, 1, 10, 3, 7, 9, 8, 3, 7, 8, 9, 10, 1, 7, 3, 11 };
  static const char representativeCase10_2[24] = { 12, 3, 1, 7, 12, 11, 12, 7, 5, 12, 8, 3, 1, 10, 12, 11, 12, 10, 5, 9, 12, 12, 9, 8 };
  static char representativeCase10_2_[24]; // determined algorithmically
  static const char representativeCase11[12] = { 4, 7, 9, 2, 10, 9, 7, 3, 2, 7, 2, 9 };
  static const char representativeCase12_1_1[12] = { 3, 2, 10, 8, 3, 10, 6, 11, 7, 8, 10, 9 };
  static char representativeCase12_1_1_[12]; // determined algorithmically
  static const char representativeCase12_2[24] = { 12, 6, 11, 9, 12, 10, 8, 7, 12, 2, 12, 3, 3, 12, 11, 12, 9, 8, 12, 2, 10, 6, 12, 7 };
  static char representativeCase12_2_[24]; // determined algorithmically

  /* The subcases of Case 13 are not clearly defined in Lewiner's paper, they are described in Chernyaev’s paper */
  static const char representativeCase13_1[12] = { 10, 1, 2, 7, 6, 11, 0, 8, 3, 5, 4, 9 };
  static char representativeCase13_1_[12]; // determined algorithmically
  static const char representativeCase13_2_0[18] = { 1, 2, 10, 11, 7, 6, 5, 4, 3, 0, 9, 5, 4, 8, 3, 5, 3, 0 };
  static char representativeCase13_2[5][18]; // determined algorithmically
  static char representativeCase13_2_0_[18]; // determined algorithmically
  static char representativeCase13_2_[5][18]; // determined algorithmically
  static const char representativeCase13_3_0[30] = { 4, 12, 5, 12, 2, 10, 1, 12, 9, 12, 1, 2, 12, 8, 3, 5, 12, 10, 9, 12, 0, 6, 11, 7, 3, 0, 12, 8, 12, 4 };
  static char representativeCase13_3[11][30]; // determined algorithmically
  static char representativeCase13_3_0_[30]; // determined algorithmically
  static char representativeCase13_3_[11][30]; // determined algorithmically
  static const char representativeCase13_4_0[36] = { 12, 4, 8, 10, 5, 12, 0, 12, 3, 4, 12, 7, 2, 12, 1, 12, 0, 9, 1, 12, 9, 12, 11, 7, 8, 3, 12, 6, 12, 5, 2, 10, 12, 6, 11, 12 };
  static char representativeCase13_4[3][36]; // determined algorithmically
  static const char representativeCase13_5_1_0[18] = { 2, 10, 4, 4, 8, 2, 1, 0, 9, 2, 8, 3, 6, 11, 7, 4, 10, 5 };
  static char representativeCase13_5_1[3][18]; // determined algorithmically
  static const char representativeCase14[12] = { 2, 6, 3, 6, 5, 9, 6, 9, 8, 3, 6, 8 };

  // the numbers of faces should be used in face test for all 256 cases
  static vector<char> faceTest_num[256];
  // may be removed later 
  static char interiorTest_num[256];

  static vector<char> ambiguityTable[15];
   
  /*  
    "triangleTable" is a three-dimensional data structure.
    The first dimension is the 256 cases, depending on the signs of values of eight corners.
    The second dimension is for each case, there may exist some subcases depending on the results of face tests and interior tests.
    The third dimension is an array, the midpoints of each three consecutive edges form a triangle and a boolean value indicating whether this case has an point in the center of a cube, encoded as 12.
  */
  static vector<vector<unsigned char> > triangleTable[256];
  static vector<bool> centerVertexNeeded[256];

  /*  
    "rc" is a three-dimensional data structure.
    The first dimension is the 15 representive cases, depending on the signs of values of eight corners.
    The second dimension is for each representive case, there may exist some subcases depending on the results of face tests and interior tests.
    The third dimension is an array, the center of the cube, encoded as 12, the midpoints of each three consecutive edges form a triangle.
  */
  static std::vector<std::vector<char> > rc[15];

  static const char faceTest3[1] = { 5 };
  static const char faceTest6[1] = { 2 };
  static const char faceTest7[3] = { 1, 2, 5 };
  static const char faceTest10[2] = { 2, 4 };
  static const char faceTest12[2] = { 4, 3 };
  static const char faceTest13[6] = { 1, 2, 3, 4, 5, 6 };
  static const char * faceTest_number[15] = { 0, 0, 0, faceTest3, 0, 0, faceTest6, faceTest7, 0, 0, faceTest10, 0, faceTest12, faceTest13, 0 };
  static const char faceTest_size[15] = { 0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 2, 0, 2, 6, 0 };
  static const char interiorTest_number[15] = { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };

  // Whether the look up tables have been loaded. Only need to be loaded once. 
  static bool tableLoaded = false;
#endif

MarchingCubes::MarchingCubes(const DistanceFieldBase * field, float iv) : distanceFieldBase(field), isoValue(iv)
{
  resolutionX = distanceFieldBase->getResolutionX();
  resolutionY = distanceFieldBase->getResolutionY();
  resolutionZ = distanceFieldBase->getResolutionZ();
}

/*
  Perform the face test: whether the product of the two positive values is larger than the product of the two negative values.
  If the face is negative, the result needs to be reversed.
*/
bool MarchingCubes::faceTest(int face_, float cube[8])
{
  int face = abs(face_) - 1;
  return (((cube[face_vertex[face][0]] * cube[face_vertex[face][2]] > cube[face_vertex[face][1]] * cube[face_vertex[face][3]])
      == (face_ > 0)) == (cube[face_vertex[face][0]] > 0));
}

/*
  Perform the interior test.
  If the face is negative, the result needs to be reversed.
*/
int MarchingCubes::interiorTest(int edge, float cube[8])
{
  //printf("%d", edge);
  //for (int p = 0; p < 8; p++) printf(" %.5f", cube[p]); printf("\n");
  double aux1 = (cube[1] - cube[0]) * (cube[6] - cube[7]) - (cube[5] - cube[4]) * (cube[2] - cube[3]);
  double aux2 = cube[0] * (cube[6] - cube[7]) - cube[4] * (cube[2] - cube[3]) + cube[7] * (cube[1] - cube[0]) - cube[3] * (cube[5] - cube[4]);
  double s = -aux2 / (2.0 * aux1);
  if ((s < 0.0) || (s > 1.0))
    return edge > 0;

  double A = cube[0] + (cube[1] - cube[0]) * s;
  double B = cube[4] + (cube[5] - cube[4]) * s;
  double C = cube[7] + (cube[6] - cube[7]) * s;
  double D = cube[3] + (cube[2] - cube[3]) * s;

  int result = ((A >= 0)) | ((B >= 0) << 1) | ((C >= 0) << 2) | ((D >= 0) << 3);
  switch (result)
  {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 6:
    case 8:
    case 9:
    case 12:
      return edge > 0;
    case 7:
    case 11:
    case 13:
    case 14:
    case 15:
      return edge < 0;
    case 5:
      return (A * C < B * D) == (edge > 0);
    case 10:
      return (A * C >= B * D) == (edge > 0);
    default:
      return -1;
  }
}

ObjMesh * MarchingCubes::compute()
{
  #ifndef USE_STATIC_TABLE
    if (tableLoaded == false)
      createTable(); // data needs to be loaded
    //printTable();
  #endif
  
  int gpResX = distanceFieldBase->getResolutionX()+1; // grid point resolution
  int gpResY = distanceFieldBase->getResolutionY()+1;
  int gpResZ = distanceFieldBase->getResolutionZ()+1;
  if (gpResX <= 1 || gpResY <= 1 || gpResZ <= 1)
    return new ObjMesh();
  
  double gridSize[3];
  distanceFieldBase->getGridSpacing(gridSize, gridSize+1, gridSize+2);

//  cout << "Begin computing: " << endl;
  PerformanceCounter pc;

//
//  vector<Vec3d> allVertices;
//  vector<int> allTriangles;
//  // edgeInersectionID: an ID to each edge if it intersects with the iso-surface
//  unordered_map<int64_t, int> edgeIntersectionID2MeshVtxID;


  pc.StartCounter();

  vector<Vec3d> allVertices;
  vector<int> allTriangles;
  unordered_map<int64_t, int> edgeInterID2AllMeshVtxID;

#ifdef USE_TBB
  struct ThreadLocalData
  {
    vector<Vec3d> vertices;
    vector<int> triangles;
    unordered_map<int64_t, int> edgeInterID2LocalVtxID;
  };
  
  tbb::enumerable_thread_specific<ThreadLocalData> threadLocalData;
  // loop over all voxels again to compute triangles
  tbb::parallel_for(tbb::blocked_range<int>(0, distanceFieldBase->getResolutionZ()), [&](const tbb::blocked_range<int> & rng)
//  for (int k = 0; k < distanceFieldBase->getResolutionZ(); k++)
  {
    auto & localData = threadLocalData.local();
    auto & edgeInterID2MeshVtxID = localData.edgeInterID2LocalVtxID;
    auto & meshVertices = localData.vertices;
    auto & meshTriangles = localData.triangles;
#else
    Range<int> rng(0, distanceFieldBase->getResolutionZ());
    auto & edgeInterID2MeshVtxID = edgeInterID2AllMeshVtxID;
    auto & meshVertices = allVertices;
    auto & meshTriangles = allTriangles;
#endif
    for (int k = rng.begin(); k != rng.end(); ++k)
      for (int j = 0; j < distanceFieldBase->getResolutionY(); ++j)
        for (int i = 0; i < distanceFieldBase->getResolutionX(); ++i)
        {
          float cube[8];
          //get distance field values at 8 cube vertices
          for (int p = 0; p < 8; p++)
          {
            cube[p] = getDistance(i + vertex_position[p][0], j + vertex_position[p][1], k + vertex_position[p][2]);
          }

          // raw case [0,255] -> major case [0,14]

          int rawCaseNumber = 0; //raw case number, 0-255
          for (int p = 0; p < 8; p++)
            rawCaseNumber |= ((cube[p] > 0) << p);

          //get the major case number, 0-14
          int majorCaseNumber = marchingCubeSymmetries[rawCaseNumber][0];
          if ((majorCaseNumber == 0) || (majorCaseNumber == 255))
            continue;   //No triangles

          int faceAndInteriorTestResult = 0;

          int faceTestNumCase = 0;
          #ifdef USE_STATIC_TABLE
            // the first entry in the faceTest_num table stores number of entries
            faceTestNumCase = faceTest_num[rawCaseNumber][0];
            if (faceTestNumCase > 0)
              for (int p = faceTestNumCase; p > 0; p--)
              {
                int o = faceTest(faceTest_num[rawCaseNumber][p], cube);
                faceAndInteriorTestResult = (faceAndInteriorTestResult << 1) | o;
              }
          #else
            faceTestNumCase = int(faceTest_num[rawCaseNumber].size());
            if (!faceTest_num[rawCaseNumber].empty())
              for (vector<char>::const_iterator p = faceTest_num[rawCaseNumber].end() - 1; p >= faceTest_num[rawCaseNumber].begin(); p--)
              {
                int o = faceTest(*p, cube);
                faceAndInteriorTestResult = (faceAndInteriorTestResult << 1) | o;
              }
          #endif

          if (interiorTest_num[rawCaseNumber] != 0)
            faceAndInteriorTestResult = faceAndInteriorTestResult | (interiorTest(interiorTest_num[rawCaseNumber], cube) << faceTestNumCase);

          int ambiguityNumber = ambiguityTable[majorCaseNumber][faceAndInteriorTestResult];
          if (ambiguityNumber < 0)
          {
            cerr << "Internal error in computing marching cubes, ambiguityNumber is: " << ambiguityNumber << endl;
            continue;
          }

          //Get the list of triangle vertices
          #ifdef USE_STATIC_TABLE
            unsigned char * triangles = NULL;
            int numTri = 0;
            unsigned char ** tableEntry = triangleTable[rawCaseNumber];
            bool hasCenter = false;
            assert(tableEntry);
            {
              unsigned char * tableEntry2 = tableEntry[ambiguityNumber];
              numTri = tableEntry2[0] / 3;
              hasCenter = (tableEntry2[1] != 0);
              triangles = tableEntry2 + 2;
            }
          #else
            const vector<unsigned char> & triangles = triangleTable[rawCaseNumber][ambiguityNumber];
            bool hasCenter = centerVertexNeeded[rawCaseNumber][ambiguityNumber];
            int numTri = triangles.size() / 3;
          #endif

          int cubeVtxIndices[13] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
          Vec3d center = Vec3d(0, 0, 0); // an additional vertex in the middle of the cube may be needed
          int sum = 0;
          for (int edgeID = 0; edgeID < 12; edgeID++) // go to 12 edges of the cube
          {
            int node = edge_vertex[edge_map[edgeID][0]][edge_map[edgeID][1]];
            int direction = edge_map[edgeID][0];   //can be 0, 1 or 2. 0 means this edge is parallel to X axis. 1 means parallel to Y axis. 2 means parallel to Z axis.
            assert(direction >= 0 && direction <= 2);
            int _i = i + vertex_position[node][0];
            int _j = j + vertex_position[node][1];
            int _k = k + vertex_position[node][2];

            float edgeStartValue = getDistance(_i, _j, _k);
            float edgeEndValue = getDistance(_i+(direction==0), _j+(direction==1), _k+(direction==2));
            if (edgeStartValue * edgeEndValue <= 0 && edgeEndValue != 0) // TODO: check degenerate case here
            {
              // int64_t packedIndex = 3 * (_j * (resX+1) + _i) + direction;
              int64_t edgeInterID = 3 * ((int64_t)gpResX*gpResY*_k + _j*gpResX + _i) + direction;
              auto it = edgeInterID2MeshVtxID.find(edgeInterID);
              if (it == edgeInterID2MeshVtxID.end())
              {
                // compute pos
                Vec3d pos(0.0);
                Vec3d gridPoint = distanceFieldBase->getGridPosition(_i,_j,_k);
                if (direction == 0)
                {
                  pos = Vec3d(gridPoint[0] + gridSize[0] * edgeStartValue / (edgeStartValue - edgeEndValue), gridPoint[1], gridPoint[2]);
                }
                else if (direction == 1)
                {
                  pos = Vec3d(gridPoint[0], gridPoint[1] + gridSize[1] * edgeStartValue / (edgeStartValue - edgeEndValue), gridPoint[2]);
                }
                else // direction == 2
                {
                  pos = Vec3d(gridPoint[0], gridPoint[1], gridPoint[2] + gridSize[2] * edgeStartValue / (edgeStartValue - edgeEndValue));
                }
                cubeVtxIndices[edgeID] = meshVertices.size();
                meshVertices.push_back(pos);
                edgeInterID2MeshVtxID.emplace(edgeInterID, cubeVtxIndices[edgeID]);
              }
              else
              {
                cubeVtxIndices[edgeID] = it->second;
              }

              if (hasCenter)
              {
                center += meshVertices[cubeVtxIndices[edgeID]];
                sum++;
              }
            }
          }
          if (hasCenter)
          {
            center /= sum;
            // store this vertex globally in allVertices
            cubeVtxIndices[12] = meshVertices.size();
            meshVertices.push_back(center);
            int64_t centerID = - (int64_t)gpResX*gpResY*k + j*gpResX + i - 1; // define an ID to store center vtx into edgeInterID2LocalVtxID
            edgeInterID2MeshVtxID.emplace(centerID, cubeVtxIndices[12]);
          }

          for(int triID = 0; triID < numTri; triID++)
          {
            if ((triangles[3*triID+0] >= 13) || (triangles[3*triID+1] >= 13) || (triangles[3*triID+2] >= 13)) // this is for safety
              continue;
            int v1 = cubeVtxIndices[triangles[3*triID+0]];
            int v2 = cubeVtxIndices[triangles[3*triID+1]];
            int v3 = cubeVtxIndices[triangles[3*triID+2]];
            assert(v1 != -1 && v2 != -1 && v3 != -1);
            meshTriangles.push_back(v1);
            meshTriangles.push_back(v2);
            meshTriangles.push_back(v3);
          }
        } // end for (distance loop)
#ifdef USE_TBB
  }); //end for locations

//  pc.StopCounter();
//  cout << "Large loop: " << pc.GetElapsedTime() << endl;
//  pc.StartCounter();

  for(const auto & local : threadLocalData)
  {
    vector<int> localVtxID2AllVtxID(local.vertices.size());
    assert(local.edgeInterID2LocalVtxID.size() == local.vertices.size());
    for(const auto & localPair : local.edgeInterID2LocalVtxID)
    {
      int localID = localPair.second;
      assert(localID < (int)localVtxID2AllVtxID.size());
      auto allIt = edgeInterID2AllMeshVtxID.find(localPair.first);
      if (allIt == edgeInterID2AllMeshVtxID.end())
      {
        int globalID = allVertices.size();
        edgeInterID2AllMeshVtxID.emplace(localPair.first, globalID);
        localVtxID2AllVtxID[localID] = globalID;
        allVertices.push_back(local.vertices[localID]);
      }
      else
      {
        localVtxID2AllVtxID[localID] = allIt->second;
      }
    }

    for(size_t i = 0; i < local.triangles.size(); i++)
    {
      allTriangles.push_back(localVtxID2AllVtxID[local.triangles[i]]);
    }
  }
#endif

  ObjMesh * objMesh = new ObjMesh(allVertices.size(), (const double*)&allVertices[0], allTriangles.size() / 3, allTriangles.data());
//  pc.StopCounter();
//  cout << "build objMesh: " << pc.GetElapsedTime() << endl;
  return objMesh;
}

namespace
{

class ObjSorter
{
public:
  ObjSorter(const ObjMesh * objMesh) : objMesh(objMesh) {}
  bool operator() (const int a, const int b) 
  {
	 return objMesh->getPosition(a) < objMesh->getPosition(b);
  }
  bool operator () (const ObjMesh::Face a, ObjMesh::Face b) 
  {
	 return triple<Vec3d, Vec3d, Vec3d>(objMesh->getPosition(a.getVertex(0)), objMesh->getPosition(a.getVertex(1)), objMesh->getPosition(a.getVertex(2))) < triple<Vec3d, Vec3d, Vec3d>(objMesh->getPosition(b.getVertex(0)), objMesh->getPosition(b.getVertex(1)), objMesh->getPosition(b.getVertex(2))) ;
  }
protected:
  const ObjMesh * objMesh;
};

void sortMesh(ObjMesh * objMesh)
{
  ObjSorter sorter(objMesh);
  //sort vertices
  vector <int> a(objMesh->getNumVertices());
  for (size_t i = 0; i < objMesh->getNumVertices(); i++) 
    a[i] = i;
  sort(a.begin(), a.end(), sorter);
  vector <int> b(objMesh->getNumVertices());
  for (size_t i = 0; i < objMesh->getNumVertices(); i++) 
    b[a[i]] = i;
  objMesh->renumberVertices(b);

  //sort faces
  const ObjMesh::Group *group = objMesh->getGroupHandle(0);
  int n = group->getNumFaces();
  vector <ObjMesh::Face> faces(n);
  for (int i = 0; i < n; i++)
  {
    const ObjMesh::Face *face = group->getFaceHandle(i);
    int _min = 0;
    int m = face->getNumVertices();
    for(int j = 1; j < m; j++)
      if (objMesh->getPosition(face->getVertex(j)) < objMesh->getPosition(face->getVertex(_min))) _min = j;
    for (int j = 0; j < m; j++)
      faces[i].addVertex(face->getVertex((j + _min) % m));
  }
  sort(faces.begin(), faces.end(), sorter);
  objMesh->removeGroup(0);
  objMesh->addGroup("Sorted faces");
  for (int i = 0; i < n; i++) 
    objMesh->addFaceToGroup(faces[i], 0);
}

}

ObjMesh * MarchingCubes::compute(const DistanceFieldBase * distanceFieldBase, float isoValue)
{
  MarchingCubes m(distanceFieldBase, isoValue);
  ObjMesh * ret = m.compute();
  sortMesh(ret);
  return ret;
}

#ifdef USE_STATIC_TABLE
  //create dummy functions
  void MarchingCubes::printTable() {}
  void MarchingCubes::createTable() {}
#else
  // --- the functions below are used to create the case tables ---

  static void binaray(int num, int * array)
  {
    for (int i = 0; i < 8; i++)
    {
      array[i] = num & 0x1;
      num >>= 1;
    }
  }

  static int permute(int src_, int dst_, int transform)
  {
    int src[8], dst[8], permute[8];
    binaray(src_, src);
    binaray(dst_, dst);

    int i;
    for (i = 0; i < 8; i++)
      permute[i] = src[cubeSymmetryVertices[transform][i]];

    for (i = 0; i < 8; i++)
    {
      if (permute[i] != dst[i])
        break;
    }

    if (i == 8)
      return 1;

    for (i = 0; i < 8; i++)
    {
      if (permute[i] != 1 - dst[i])
        break;
    }

    if (i == 8)
      return -1;

    return 0;
  }

  static inline int f_hash_edge(int a, int b)
  {
    return a * a + b * b;
  }

  static inline int f_hash_face(int a, int b, int c, int d)
  {
    return a + b + c + d + (a > 0 && b > 0 && c > 0 && d > 0);
  }

  static void permute_face(const char *src, char *dst, int n, int permute)
  {
    int permute_ = abs(permute);
    for (int i = 0; i < n; i++)
      dst[i] = cubeSymmetryFaces[permute_][src[i] - 1] + 1;

    if (permute <= 0)
    {
      for (int i = 0; i < n; i++)
        dst[i] = -dst[i];
    }
  }

  static void permute_edge(const char *src, char *dst, int n, int permute)
  {
    int permute_ = abs(permute);
    for (int i = 0; i < n; i++)
      dst[i] = cubeSymmetryEdges[permute_][(unsigned char) src[i]];

    if (permute <= 0)
    {
      for (int i = 1; i < n; i += 3)
        swap(dst[i - 1], dst[i + 1]);
    }
  }

  static vector<char> make_vector_impl(const char *p, int n)
  {
    vector<char> ret;
    ret.resize(n);
    memcpy(&ret[0], p, sizeof(char) * n);
    return ret;
  }

  static bool violate(int k, vector<vector<int> >& binaryRepresentation) 
  {
    int p = 0, n = 0;
    for (int i = 0; i < 3; i++) 
    {
      bool r0 = find(binaryRepresentation[k].begin(), binaryRepresentation[k].end(), oppositeFace[i][0]) == binaryRepresentation[k].end();
      bool r1 = find(binaryRepresentation[k].begin(), binaryRepresentation[k].end(), oppositeFace[i][1]) == binaryRepresentation[k].end();
      n += r0 && r1;
      p += !(r0 || r1);
    }

    return (n > 0) && (p > 0);
  }

  static bool compare(const int *a, const int *b) 
  {
    for (int i = 0; i < 3; i++)
      if (a[i] != b[i]) return false;
    return true;
  }

  static bool cmp(const vector<int> a, const vector<int> b) 
  {
    bool interiora = ((a.size() == 3) && (compare(&a[0], vtxAdjacentFaces[1]) || compare(&a[0], vtxAdjacentFaces[3]) || compare(&a[0], vtxAdjacentFaces[4]) || compare(&a[0], vtxAdjacentFaces[6])));

    bool interiorb = ((b.size() == 3) && (compare(&b[0], vtxAdjacentFaces[1]) || compare(&b[0], vtxAdjacentFaces[3]) || compare(&b[0], vtxAdjacentFaces[4]) || compare(&b[0], vtxAdjacentFaces[6])));

    if (!interiora && interiorb) 
      return true;

    if (interiora && !interiorb) 
      return false;

    if (a.size() < b.size()) 
      return true;

    if (a.size() > b.size()) 
      return false;

    if (a.size() <= 3) 
    {
      for (unsigned int i = 0; i < a.size(); i++)
      {
        if (a[i] < b[i]) 
          return true;
        if (a[i] > b[i]) 
          return false;
      }
      return true;
    }
    else 
    {
      for (int i = 0; i < 6; i++) 
      {
        bool ra = find(a.begin(), a.end(), i) == a.end();
        bool rb = find(b.begin(), b.end(), i) == b.end();
        if (ra && !rb) 
          return true;
        if (!ra && rb) 
          return false;
      }
      return true;
    }
  }

  static vector<int> make_vector13(unsigned int p) 
  {
    vector <int> ret;
    for (int i = 1; i <= 6; i++) 
    {
      if ((p | ((1 << i) >> 1)) == p)
        ret.push_back(i);
    }
    return ret;
  }

  #define make(x, y, z) permute_edge(x, y, sizeof(x) / sizeof(x[0]), z)
  #define make_vector(x) make_vector_impl(x, sizeof(x) / sizeof(x[0]))

  void MarchingCubes::createTable()
  {
    // fill the ambiguity table
    ambiguityTable[0].push_back(0);
    ambiguityTable[1].push_back(0);
    ambiguityTable[2].push_back(0);
    ambiguityTable[3].push_back(0);
    ambiguityTable[3].push_back(1);
    ambiguityTable[4].push_back(1);
    ambiguityTable[4].push_back(0);
    ambiguityTable[5].push_back(0);

    ambiguityTable[6].push_back(0);
    ambiguityTable[6].push_back(1);

    ambiguityTable[7].push_back(0);
    ambiguityTable[7].push_back(1);
    ambiguityTable[7].push_back(2);
    ambiguityTable[7].push_back(4);
    ambiguityTable[7].push_back(3);
    ambiguityTable[7].push_back(5);
    ambiguityTable[7].push_back(6);
    ambiguityTable[7].push_back(7);

    ambiguityTable[8].push_back(0);
    ambiguityTable[9].push_back(0);

    ambiguityTable[10].push_back(2);
    ambiguityTable[10].push_back(3);
    ambiguityTable[10].push_back(4);
    ambiguityTable[10].push_back(1);
    ambiguityTable[10].push_back(0);
    ambiguityTable[10].push_back(3);
    ambiguityTable[10].push_back(4);
    ambiguityTable[10].push_back(1);

    ambiguityTable[11].push_back(0);

    ambiguityTable[12].push_back(0);
    ambiguityTable[12].push_back(2);
    ambiguityTable[12].push_back(3);
    ambiguityTable[12].push_back(1);
    ambiguityTable[14].push_back(0);

    ambiguityTable[13].resize(64);

    vector<vector<int> > binaryRepresentation;

    for (int i = 0; i < 64; i++)
      binaryRepresentation.push_back(make_vector13(i));

    for (int i = 63; i >= 0; i--)
    {
      if (violate(i, binaryRepresentation)) 
        binaryRepresentation.erase(binaryRepresentation.begin() + i);
    }

    sort(binaryRepresentation.begin(), binaryRepresentation.end(), cmp);
    memset(&ambiguityTable[13][0], 0xff, sizeof(ambiguityTable[13][0]) * 64);
    int pp = 0;
    for (unsigned int i = 0; i < binaryRepresentation.size(); i++)
    {
      int k = 0;
      for (unsigned int j = 0; j < binaryRepresentation[i].size(); j++)
        k += (1 << binaryRepresentation[i][j]) >> 1;
      ambiguityTable[13][k] = pp++;
    }

    // build marchingCubeSymmetries
    int major_case = 0;
    int sub_case[100];
    const unsigned char sequence[256] =
    { 0, 1, 2, 4, 8, 16, 32, 64, 128, 3, 5, 9, 17, 33, 65, 129, 6, 10, 18, 34, 66, 130, 12, 20, 36, 68, 132, 24, 40, 72, 136, 48, 80, 144, 96, 160, 192, 7, 11, 19, 35, 67, 131, 13, 21, 37, 69, 133, 25, 41, 73, 137, 49, 81, 145, 97, 161, 193, 14, 22, 38, 70, 134, 26, 42, 74, 138, 50, 82, 146, 98, 162, 194, 28, 44, 76, 140, 52, 84, 148, 100, 164, 196, 56, 88, 152, 104, 168, 200, 112, 176, 208, 224, 15, 23, 39, 71, 135, 27, 43, 75, 139, 51, 83, 147, 99, 163, 195, 29, 45, 77, 141, 53, 85, 149, 101, 165, 197, 57, 89, 153, 105, 169, 201, 113, 177, 209, 225, 30, 46, 78, 142, 54, 86, 150, 102, 166, 198, 58, 90, 154, 106, 170, 202, 114, 178, 210, 226, 60, 92, 156, 108, 172, 204, 116, 180, 212, 228, 120, 184, 216, 232, 240, 31, 47, 79, 143, 55, 87, 151, 103, 167, 199, 59, 91, 155, 107, 171, 203, 115, 179, 211, 227, 61, 93, 157, 109, 173, 205, 117, 181, 213, 229, 121, 185, 217, 233, 241, 62, 94, 158, 110, 174, 206, 118, 182, 214, 230, 122, 186, 218, 234, 242, 124, 188, 220, 236, 244, 248, 63, 95, 159, 111, 175, 207, 119, 183, 215, 231, 123, 187, 219, 235, 243, 125, 189, 221, 237, 245, 249, 126, 190, 222, 238, 246, 250, 252, 127, 191, 223, 239, 247, 251, 253, 254, 255 };

    vector<char> c[15];
    pair<int, int> candidate;
    set<int> am;
    for (int i_ = 0; i_ < 256; i_++)
    {
      int i = sequence[i_];
      int j = 0, k, p, j_;
      for (j_ = 0; j_ < i_; j_++)
      {
        j = sequence[j_];
        for (k = 0; k < 24; k++)
        {
          p = permute(i, j, k);
          if (p > 0) 
          {
            candidate = make_pair(k, p);
            break;
          }
        }

        if (k < 24)
          break;

        for (k = 0; k < 24; k++)
        {
          p = permute(i, j, k);
          if (p < 0) 
          {
            candidate = make_pair(k, p);
            break;
          }
        }

        if (k < 24) 
          break;
      }

      if (i_ == j_)
      {
        marchingCubeSymmetries[i][0] = caseNumber[major_case];
        marchingCubeSymmetries[i][1] = 0;
        marchingCubeSymmetries[i][2] = 0;
        sub_case[caseNumber[major_case++]] = 0;
      }
      else
      {
        k = candidate.first;
        p = candidate.second;
        am.insert(marchingCubeSymmetries[j][0]);
        marchingCubeSymmetries[i][0] = marchingCubeSymmetries[j][0];
        sub_case[marchingCubeSymmetries[j][0]]++;
        marchingCubeSymmetries[i][1] = sub_case[marchingCubeSymmetries[j][0]];
        marchingCubeSymmetries[i][2] = p * k;
      }
      c[marchingCubeSymmetries[i][0]].push_back(i);
    }

    unsigned char hash_edge[86];
    memset(hash_edge, 0xff, sizeof(hash_edge));
    for (unsigned int i = 0; i < 12; i++)
      hash_edge[f_hash_edge(edge_vertices[i][0], edge_vertices[i][1])] = i;

    unsigned char hash_face[23];
    memset(hash_face, 0xff, sizeof(hash_face));
    for (unsigned int i = 0; i < 6; i++)
      hash_face[f_hash_face(face_vertices[i][0], face_vertices[i][1], face_vertices[i][2], face_vertices[i][3])] = i;

    //build cubeSymmetryEdges
    for (unsigned int i = 0; i < 25; i++)
    {
      for (unsigned int j = 0; j < 12; j++)
      {
        unsigned int p0 = cubeSymmetryVertices[i][edge_vertices[j][0]];
        unsigned int p1 = cubeSymmetryVertices[i][edge_vertices[j][1]];
        cubeSymmetryEdges[i][j] = hash_edge[f_hash_edge(p0, p1)];
      }
      cubeSymmetryEdges[i][12] = 12;
    }

    // build cubeSymmetryFaces
    for (unsigned int i = 0; i < 24; i++)
    {
      for (unsigned int j = 0; j < 6; j++)
      {
        unsigned char p0 = cubeSymmetryVertices[i][face_vertices[j][0]];
        unsigned int p1 = cubeSymmetryVertices[i][face_vertices[j][1]];
        unsigned int p2 = cubeSymmetryVertices[i][face_vertices[j][2]];
        unsigned int p3 = cubeSymmetryVertices[i][face_vertices[j][3]];
        cubeSymmetryFaces[i][j] = hash_face[f_hash_face(p0, p1, p2, p3)];
      }
    }

    //fix 13 for marchingCubeSymmetries
    marchingCubeSymmetries[90][2] = 3;

    //build rc
    make(representativeCase7_2_0, representativeCase7_2_1, 17);
    make(representativeCase7_2_0, representativeCase7_2_2, 21);
    make(representativeCase7_3_0, representativeCase7_3_1, 21);
    make(representativeCase7_3_0, representativeCase7_3_2, 17);
    make(representativeCase10_1_1, representativeCase10_1_1_, -1);
    make(representativeCase10_2, representativeCase10_2_, 12);
    make(representativeCase12_1_1, representativeCase12_1_1_, 24);
    make(representativeCase12_2, representativeCase12_2_, 24);
    make(representativeCase13_1, representativeCase13_1_, -3);

    make(representativeCase13_2_0, representativeCase13_2[0], 23);
    make(representativeCase13_2_0, representativeCase13_2[1], 9);
    make(representativeCase13_2_0, representativeCase13_2[2], 16);
    make(representativeCase13_2_0, representativeCase13_2[3], 20);
    make(representativeCase13_2_0, representativeCase13_2[4], 19);

    make(representativeCase13_2_0, representativeCase13_2_0_, -2);
    make(representativeCase13_2_0_, representativeCase13_2_[0], 23);
    make(representativeCase13_2_0_, representativeCase13_2_[1], 9);
    make(representativeCase13_2_0_, representativeCase13_2_[2], 16);
    make(representativeCase13_2_0_, representativeCase13_2_[3], 20);
    make(representativeCase13_2_0_, representativeCase13_2_[4], 19);

    make(representativeCase13_3_0, representativeCase13_3_0_, 24);
    make(representativeCase13_3_0, representativeCase13_3[0], 8);
    make(representativeCase13_3_0_, representativeCase13_3_[0], 8);
    make(representativeCase13_3_0, representativeCase13_3_[1], -5);
    make(representativeCase13_3_0_, representativeCase13_3[1], -5);
    make(representativeCase13_3_0, representativeCase13_3_[2], -2);
    make(representativeCase13_3_0_, representativeCase13_3[2], -2);
    make(representativeCase13_3_0, representativeCase13_3_[3], -3);
    make(representativeCase13_3_0_, representativeCase13_3[3], -3);
    make(representativeCase13_3_0, representativeCase13_3[4], 17);
    make(representativeCase13_3_0_, representativeCase13_3_[4], 17);
    make(representativeCase13_3_0, representativeCase13_3_[5], -4);
    make(representativeCase13_3_0_, representativeCase13_3[5], -4);
    make(representativeCase13_3_0, representativeCase13_3[6], 9);
    make(representativeCase13_3_0_, representativeCase13_3_[6], 9);
    make(representativeCase13_3_0, representativeCase13_3_[7], -11);
    make(representativeCase13_3_0_, representativeCase13_3_[7], -11);
    make(representativeCase13_3_0, representativeCase13_3[8], 18);
    make(representativeCase13_3_0_, representativeCase13_3_[8], 18);
    make(representativeCase13_3_0, representativeCase13_3_[9], -10);
    make(representativeCase13_3_0_, representativeCase13_3[9], -10);
    make(representativeCase13_3_0, representativeCase13_3[10], 16);
    make(representativeCase13_3_0_, representativeCase13_3_[10], 16);

    make(representativeCase13_4_0, representativeCase13_4[0], 8);
    make(representativeCase13_4_0, representativeCase13_4[1], 17);
    make(representativeCase13_4_0, representativeCase13_4[2], 9);
    make(representativeCase13_5_1_0, representativeCase13_5_1[0], 8);
    make(representativeCase13_5_1_0, representativeCase13_5_1[1], 7);
    make(representativeCase13_5_1_0, representativeCase13_5_1[2], 9);
    rc[1].push_back(make_vector(representativeCase1));
    rc[2].push_back(make_vector(representativeCase2));
    rc[3].push_back(make_vector(representativeCase3_1));
    rc[3].push_back(make_vector(representativeCase3_2));
    rc[4].push_back(make_vector(representativeCase4_1));
    rc[4].push_back(make_vector(representativeCase4_2));
    rc[5].push_back(make_vector(representativeCase5));
    rc[6].push_back(make_vector(representativeCase6_1_1));
    rc[6].push_back(make_vector(representativeCase6_2));
    rc[7].push_back(make_vector(representativeCase7_1));
    rc[7].push_back(make_vector(representativeCase7_2_0));
    rc[7].push_back(make_vector(representativeCase7_2_1));
    rc[7].push_back(make_vector(representativeCase7_2_2));
    rc[7].push_back(make_vector(representativeCase7_3_0));
    rc[7].push_back(make_vector(representativeCase7_3_1));
    rc[7].push_back(make_vector(representativeCase7_3_2));
    rc[7].push_back(make_vector(representativeCase7_4_1));
    rc[8].push_back(make_vector(representativeCase8));
    rc[9].push_back(make_vector(representativeCase9));
    rc[10].push_back(make_vector(representativeCase10_1_1));
    rc[10].push_back(make_vector(representativeCase10_1_1_));
    rc[10].push_back(make_vector(representativeCase10_1_2));
    rc[10].push_back(make_vector(representativeCase10_2));
    rc[10].push_back(make_vector(representativeCase10_2_));
    rc[11].push_back(make_vector(representativeCase11));
    rc[12].push_back(make_vector(representativeCase12_1_1));
    rc[12].push_back(make_vector(representativeCase12_1_1_));
    rc[12].push_back(make_vector(representativeCase12_2));
    rc[12].push_back(make_vector(representativeCase12_2_));
    rc[13].push_back(make_vector(representativeCase13_1));

    rc[13].push_back(make_vector(representativeCase13_2_0));
    for (unsigned int i = 0; i < 5; i++)
      rc[13].push_back(make_vector(representativeCase13_2[i]));

    rc[13].push_back(make_vector(representativeCase13_3_0));
    for (unsigned int i = 0; i < 11; i++)
      rc[13].push_back(make_vector(representativeCase13_3[i]));

    rc[13].push_back(make_vector(representativeCase13_4_0));
    for (unsigned int i = 0; i < 3; i++)
      rc[13].push_back(make_vector(representativeCase13_4[i]));

    rc[13].push_back(make_vector(representativeCase13_3_0_));
    for (unsigned int i = 0; i < 11; i++)
      rc[13].push_back(make_vector(representativeCase13_3_[i]));

    rc[13].push_back(make_vector(representativeCase13_2_0_));
    for (unsigned int i = 0; i < 5; i++)
      rc[13].push_back(make_vector(representativeCase13_2_[i]));

    rc[13].push_back(make_vector(representativeCase13_1_));

    rc[13].push_back(make_vector(representativeCase13_5_1_0));
    for (unsigned int i = 0; i < 3; i++)
      rc[13].push_back(make_vector(representativeCase13_5_1[i]));

    rc[14].push_back(make_vector(representativeCase14));

    //build faceTest vector
    for (unsigned int i = 0; i < 15; i++)
    {
      for (unsigned int j = 0; j < c[i].size(); j++)
      {
        const unsigned char k = c[i][j];
        const char n = faceTest_size[i];
        faceTest_num[k].resize(n);
        if (n == 0)
          continue;
        if (j == 0)
          memcpy(&faceTest_num[k][0], faceTest_number[i], sizeof(faceTest_number[i][0]) * n);
        else
          permute_face(faceTest_number[i], &faceTest_num[k][0], n, marchingCubeSymmetries[k][2]);
      }
    }

    //build interiorTest vector
    for (unsigned int i = 0; i < 15; i++)
    {
      for (unsigned int j = 0; j < c[i].size(); j++)
      {
        const unsigned char k = c[i][j];
        interiorTest_num[k] = interiorTest_number[i];
        if (marchingCubeSymmetries[k][2] < 0 || (j > 0 && marchingCubeSymmetries[k][2] == 0))
          interiorTest_num[k] = -interiorTest_num[k];
      }
    }

    //build triangleTable
    for (unsigned int i = 0; i < 256; i++)
    {
      const int major_case = marchingCubeSymmetries[i][0];
      triangleTable[i].resize(rc[major_case].size());
      centerVertexNeeded[i].resize(rc[major_case].size());
      for (unsigned int j = 0; j < rc[major_case].size(); j++)
      {
        triangleTable[i][j].resize(rc[major_case][j].size());
        if (marchingCubeSymmetries[i][1] == 0)   //
        {
          memcpy(&triangleTable[i][j][0], &rc[major_case][j][0], rc[major_case][j].size() * sizeof(rc[major_case][j][0]));
        }
        else
        {
          permute_edge(&rc[major_case][j][0], (char*) &triangleTable[i][j][0], rc[major_case][j].size(), marchingCubeSymmetries[i][2]);
        }
        centerVertexNeeded[i][j] = false;
        for (unsigned int k = 0; k < triangleTable[i][j].size(); k++)
          centerVertexNeeded[i][j] = centerVertexNeeded[i][j] | (triangleTable[i][j][k] == 12);// |= (triangleTable[i][j][k] == 12);
      }
    }

    tableLoaded = true;
  }

  void MarchingCubes::printTable()
  {
    // once the file is OK, rename it to marchingCubesTable.h
    ofstream fout("marchingCubesTable_.h", ios::binary);
    if (!fout) 
      return;
    fout << endl;

    //static vector<char> ambiguityTable[15];
    unsigned int maxSize = 0;
    for(unsigned int i = 0; i < 15; i++) 
    {
      if (ambiguityTable[i].size() > maxSize) 
        maxSize = ambiguityTable[i].size();
    }
    fout << "static char ambiguityTable[15][" << maxSize << "] = " << endl << "{" << endl;
    //char buffer[100];
    //fout << 
    for(int i = 0; i < 15; i++) 
    {
      fout << "  { ";
      unsigned int j = 0;
      for(j = 0; j < ambiguityTable[i].size(); j++) 
      {
        fout << int(ambiguityTable[i][j]);
        if (j < maxSize-1)  
          fout << ", ";
      }
      for(; j < maxSize; j++) 
      {
        fout << "0";
        if (j < maxSize-1)  
          fout << ", ";
      }
      fout << " }";
      if (i < 15 - 1)
        fout << ", ";
      fout << endl;
    }
    fout << "};" << endl;

    fout << endl;

    // static int marchingCubeSymmetries[256][3];
    fout << "static int marchingCubeSymmetries[256][3] = {" << endl;
    for(int i = 0; i < 256; i++) 
    {
      fout << "  {";
      for(int j = 0; j < 3; j++) 
      {
        fout << setw(3) << marchingCubeSymmetries[i][j];
        if (j < 3 - 1) 
          fout << ",";
      }
      fout << " }";
      if (i < 256 -1)
        fout << ", ";
      if (i % 8 == 7) 
        fout << endl;
    }
    fout << "};" << endl;

    fout << endl;

    //static vector<char> faceTest_num[256];
    maxSize = 0;
    for(unsigned int i = 0; i < 256; i++) 
    {
      if (faceTest_num[i].size() > maxSize) 
        maxSize = faceTest_num[i].size();
    }
    fout << "static char faceTest_num[256][" << maxSize+1 << "] = " << endl << "{" << endl;
    //char buffer[100];
    for(int i = 0; i < 256; i++) 
    {
      fout << "  {";
      fout << setw(2) << faceTest_num[i].size();
      //if (faceTest_num[i].size() > 0)
      fout << ",";

      unsigned int j = 0;
      for(j = 0; j < faceTest_num[i].size(); j++) 
      {
        fout << setw(2) << int(faceTest_num[i][j]);
        if (j < maxSize-1)  
          fout << ",";
      }
      for(; j < maxSize; j++) 
      {
        fout << " 0";
        if (j < maxSize-1)  
          fout << ",";
      }
      fout << " }";
      if (i < 256 - 1)
        fout << ",";
      if (i % 8 == 7)
        fout << endl;
    }
    fout << "};" << endl;

    fout << endl;

    //static char interiorTest_num[256];
    fout << "static char interiorTest_num[256] = " << endl << "  { ";
    for(int i = 0; i < 256; i++) 
    {
      fout << int(interiorTest_num[i]);
      if (i < 256 - 1)
        fout << ", ";
    }
    fout << " };" << endl;

    fout << endl;

    // static vector< vector<unsigned char> > triangleTable[256];
    for(int i = 0; i < 256; i++) 
    {
      if (triangleTable[i].size() == 0) 
        continue;

      for(unsigned int j = 0; j < triangleTable[i].size(); j++) 
      {
        vector<unsigned char> & table = triangleTable[i][j];
        unsigned int num = table.size();
        fout << "static unsigned char triangleTable_" << i << "_" << j << "[" << num+2 << "] = " << " { ";
        fout << num << ", " << centerVertexNeeded[i][j];
        if (num > 0)
          fout << ", ";
        for(unsigned int k = 0; k < table.size(); k++) 
        {
          fout << int(table[k]);
          if (k < table.size() - 1)
            fout << ", ";
        }
        fout << " };" << endl;
      }

      unsigned int num = triangleTable[i].size();
      fout << "static unsigned char * triangleTable_" << i << "[" << num << "] = ";
      if (num <= 7) 
      {
        fout << "{ ";
        for(unsigned int k = 0; k < num; k++) 
        {
          fout << "triangleTable_" << i << "_" << k;
          if (k < num - 1)
            fout << ", ";
        }
        fout << " };" << endl;
      }
      else 
      {
        fout << endl << "{ " << endl;
        for(unsigned int k = 0; k < num; k++) 
        {
          if (k % 8 == 0)
            fout << "  ";
          fout << "triangleTable_" << i << "_" << k;
          if (k < num - 1)
            fout << ", ";
          if (k % 8 == 7)
          fout << endl;
        }
        if (num % 8 != 0)
          fout << endl;

        fout << "};" << endl;

      }
      fout << endl;
    } //end 256

    fout << "static unsigned char ** triangleTable[256] = " << endl << "{ " << endl;
    for(unsigned int k = 0; k < 256; k++) 
    {
      if (k % 8 == 0)
        fout << "  ";
      if (triangleTable[k].size() == 0)
        fout << "NULL";
      else
        fout << "triangleTable_" << k;

      if (k < 256 - 1)
        fout << ", ";
      if (k % 8 == 7)
        fout << endl;
    }
    fout << "};" << endl << endl; 
    fout.close();
    //exit(0);
  }
#endif
