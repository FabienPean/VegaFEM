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

#include "virtualTets.h"
#include "triMeshNeighbor.h"
#include "triMeshPseudoNormal.h"
#include "containerHelper.h"
#include "disjointSet.h"
#include "tetTriMeshCutting.h"
#include "basicAlgorithms.h"
#include "profiler.h"
#include "windingNumberTree.h"
#include "listIO.h"
#include "geometryQuery.h"
#include "tetTriCuttingData.h"
#include "valueIndex.h"
#include "createTriMesh.h"
#include "labelOuterTets.h"
#include "tetTriPiece.h"
#include <iostream>
#include <sstream>
#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif
using namespace std;

// The simplest routine with minimal input.
TetMeshGeo createVirtualTetsMesh(const TetMeshGeo & tetMesh, const TriMeshGeo & triMesh, BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris,
    bool verbose, Profiler * profiler)
{
  // compute triangle mesh postions in exact arithmetic
  vector<Vec3ER> triPosER;
  for(int i = 0; i < triMesh.numVertices(); i++)
  {
    Vec3d p = triMesh.pos(i);
    triPosER.emplace_back(p[0], p[1], p[2]);
  }
  // call the overloading function with exact triMesh positions as additional input
  return createVirtualTetsMesh(tetMesh, triMesh, triPosER, coord, newTetVtxID2OldTetVtxID, newTetID2OldTetID, tetTris,
      verbose, profiler);
}

// The routine with triangle mesh postions in exact arithmetic as additional input
TetMeshGeo createVirtualTetsMesh(const TetMeshGeo & tetMesh, const TriMeshGeo & triMesh,
    const std::vector<Vec3ER> & exactTriPos, BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris,
    bool verbose, Profiler * profiler)
{
  TetTriCuttingData cuttingData;
  ProfilerSection tetTriCutSec(profiler, "tetTriCutting");

  // use class TetTriMeshCutting to compute data in TetTriCuttingData
  tetTriCutSec.start();
  TetTriMeshCutting cutting(tetMesh, triMesh);

  cout << "Computing tet-triangle intersection..." << endl;
  // compute which tet intersect with triangle
  cutting.computeIntersectingTriTets();
  
  // Definition:
  //   cut triangle mesh: the new triangle mesh which is the result of the cut on the input triangle mesh by all the tets
  //     The newly geneated cuts are modeled by adding new vertices on the cut triangles.
  //     New vertices are shared by neighboring cut triangles. So after the cut, the topology of the input triangle mesh
  //     is preserved. Only its surface is subdivided by the tets' cut.

  // Here we compute the cut triangle mesh
  cout << "Performing Sutherland-Hodgman tet-triangle clipping algorithm..." << endl;
  cutting.computeCutTriangles(&exactTriPos);
  cout << "Finished clipping." << endl;
  tetTriCutSec.stop();

  // Definition:
  //   cut feature: a pair of tetKey and triKey, used to store which geometry feature
  //     (tet: interior/face/edge/vtx, tri: interior, edge, vtx) a cut vertex lies on
  // cuttingData.cutTriGroups: tetID -> cut triangle info in each tet
  //   cut triangle info includes: 1) the 3 cutVtxID of the triangle, 2) the original triangle ID where this cut triangle is from

  cuttingData.cutVtxPos = cutting.cutTriPositions();     // cut triangle mesh vertex positions
  cuttingData.cutVtxPosER = cutting.cutTriPositionsER();  // cut triangle mesh vertex positions in exact arithmetic
  cuttingData.features = cutting.cutVertexFeatures(); // cut feature of each cut triangle mesh vertex
  cuttingData.tetPosER = cutting.tetPositionsER();     // tet mesh position in exact arithmetic
  for(int i = 0; i < tetMesh.numTets(); i++)
    cuttingData.cutTriGroups.push_back(cutting.getCutTriGroup(i)); // cut triangle info in each tet

  // call the last, most complicated routine to do the work
  return createVirtualTetsMesh(tetMesh, triMesh, cuttingData, coord, newTetVtxID2OldTetVtxID, newTetID2OldTetID, tetTris,
      true, verbose, profiler);
}

// the last, most complicated routine where TetTriCuttingData and removeAirTets are additional input
TetMeshGeo createVirtualTetsMesh(const TetMeshGeo & tetMesh, const TriMeshGeo & triMesh,
    const TetTriCuttingData & cuttingData, BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris, bool removeAirTets,
    bool verbose, Profiler * profiler)
{
  ProfilerSection buildCopiesSec(profiler, "buildTetCopies");
  ProfilerSection connectCopiesSec(profiler, "connectTetCopies");
  ProfilerSection finalStuffSec(profiler, "finalStuffInVirtualTet");

  if (verbose)
    cout << "Begin virtual tets..." << endl;
  if (newTetVtxID2OldTetVtxID) newTetVtxID2OldTetVtxID->clear();
  if (newTetID2OldTetID) newTetID2OldTetID->clear();

  // check whether the input triangle mesh is manifold
  if (areTrianglesManifold(triMesh.triangles()) == false)
  {
    cout << "Error: input triangle mesh for virtual tets algorithm is not manifold" << endl;
    throw 1;
  }
  TriangleNeighbor triangleNeighbor(triMesh.triangles());

  // compute normals of each input triangle
  // we use exact arithmetic to compute normals only when the precision is not enough
  vector<Vec3d> triNormals(triMesh.numTriangles());
  for(int triID = 0; triID < triMesh.numTriangles(); triID++)
  {
    triNormals[triID] = triMesh.ref().computeTriangleNormal(triID);
    if (triNormals[triID].hasNaN())
    {
      Vec3ER pa = cuttingData.cutVtxPosER[triMesh.triVtxID(triID, 0)];
      Vec3ER pb = cuttingData.cutVtxPosER[triMesh.triVtxID(triID, 1)];
      Vec3ER pc = cuttingData.cutVtxPosER[triMesh.triVtxID(triID, 2)];
      Vec3ER pNormal = cross(pb - pa, pc- pa);
      Vec3d dNormal(0.0);
      for(int i = 0; i < 3; i++)
        dNormal[i] = ER_toDouble(pNormal[i]);
      dNormal.normalize();
      assert(dNormal.hasNaN() == false);
      triNormals[triID] = dNormal;
      // nanNormalTriIDs.push_back(triID);
    }
  }
  // TODO: we can input triMeshNormal from caller to reduce exact-arithmetic calls
  TriMeshPseudoNormal triMeshNormal(triMesh, triNormals.data());

  auto triangleBoundaryLoops = triangleNeighbor.findBoundaryLoops(triMesh.triangles());
  if (triangleBoundaryLoops.size() > 0)
  {
    cout << "Error: the input mesh has boundary" << endl;
    throw 1;
  }

  // #vtx in the input triangle mesh.
  // the cut mesh is assumed that the first numUncutVtx vertices are the same vertices from the input triangle mesh
  // the rest vertices are those added due to cuts
  int numUncutVtx = triMesh.numVertices();

  vector<bool> tetOutside; // if removeAirTets, stores the mapping: tetID -> whether the tet is outside the triangle mesh
  if (removeAirTets)
  {
    ProfilerSection airTetsSec(profiler, "airTets");
    airTetsSec.start();

    auto isTetBou = [&](int tetID) // return true if this tet intersects the triangle mesh
    {
      return cuttingData.cutTriGroups[tetID].tri.size() > 0;
    };

    // winding number tree from Jacobson 2013, Robust Inside-Outside Segmentation using Generalized Winding Numbers
    // it detects whether a point is inside/outside a triangle mesh
    WindingNumberTree wnTree;
    wnTree.build(triMesh);
    auto isTetOuter = [&](int tetID) // return true if this tet is outside the triangle mesh
    {
      double wn = wnTree.windingNumber(triMesh, tetMesh.ref().computeTetCenter(tetID));
      return wn < 0.5; // winding number = 0: outside, = 1: inside
    };

    // build the neighboring structure of the input tet mesh
    TetNeighbor tetNeighbor(tetMesh.tets());

    tetOutside = labelOuterTets(tetMesh, tetNeighbor, isTetBou, isTetOuter);
    airTetsSec.stop();
  }

  const vector<Vec3d> & cutPositions = cuttingData.cutVtxPos;
  const vector<Vec3ER> & cutPositionsER = cuttingData.cutVtxPosER;

#ifdef USE_VT_MULTITHREADING
  vector<mutex> cutVtxPosVecMutex(cuttingData.cutVtxPos.size());
  vector<mutex> tetPosVecMutex(cuttingData.tetPosER.size());
  cuttingData.multiThreading = true;
  cuttingData.cutVtxPosVecMutex = &cutVtxPosVecMutex;
  cuttingData.tetPosVecMutex = &tetPosVecMutex;
#endif

  // TetTriPiece: represent a piece in our paper
  // Pieces can be merged through the algorithm to represent all the boundaries of a '-' region in a tet
  // In the end, each such '-' region is assigned a unique tet
  vector<TetTriPiece> tetPieces;
  vector<vector<int>> piecesInTet(tetMesh.numTets()); // tetID -> pieces inside the tet, with pieceIDs in tetPieces

  buildCopiesSec.start();
  if (verbose)
    cout << "Begin to visit each tet for building comps" << endl;

  // go to each input tet, compute the pieces inside them, build the region graph as explained in our paper,
  // and stores the boundaries of each computed '-' region into tetPieces
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    if (removeAirTets && tetOutside[tetID]) continue; // if it's an outer tet, skip
    if (verbose && tetID % 1000 == 0)
      cout << "process comps in tetID " << tetID << "/" << tetMesh.numTets() << endl;

    TetShape tetShape(tetMesh, tetID); // class to do geometry query on the tet
    const auto & cutTriGroup = cuttingData.cutTriGroups[tetID];
    vector<vector<int>> triCCs = getConnectedComponentsByVertex(cutTriGroup.tri);
    for (const auto & cc : triCCs)
      assert(cc.size() > 0);

    // if no triangles inside this tet, then the entire volume of the tet is one region
    if (triCCs.size() == 0)
    {
      piecesInTet[tetID].push_back(tetPieces.size());
      tetPieces.emplace_back(tetID);
      continue;
    }
    // there is only one piece, then it is simple to process
    if (triCCs.size() == 1)
    {
      piecesInTet[tetID].push_back(tetPieces.size());
      tetPieces.emplace_back(tetID, triCCs[0], cutTriGroup, cuttingData, triMeshNormal, tetShape);
      continue;
    }

    // sample points on each piece
    // a sample point is used to do pseudo-normal test for building the region graph
    vector<Vec3d> samplePoint(triCCs.size()); // sample point position on each piece
    vector<Vec3ER> samplePointER(triCCs.size()); // exact arithmetric sample point position

    // find an appropriate sample point for each piece
    // to avoid degeneracy, sample points should be away from the tet boundary
    for(size_t i = 0; i < triCCs.size(); i++)
    {
      const vector<int> & cc = triCCs[i]; // local triangleIDs in cutTriGroup.tri
      bool found = false;
      for(int triID : cc)
      {
        Vec3i t = cutTriGroup.tri[triID];
        for(int j = 0; j < 3; j++)
        {
          int v = t[j];          // get cutVtxID
          if (v < numUncutVtx) // cutVtxID < numUncutVtx means this vtx is an original vertex in the input triangle mesh
          {
            samplePoint[i] = cutPositions[v];
            samplePointER[i] = cutPositionsER[v];
            found = true; // then we pick this vertex
            break;
          }
        }
        if (found) break;
      }
      // Here we choose the sample point as one original vertex because generally, original vertex is on the interior of the tet
      // We can make it more robust by finding an original vertex that's not on the tet faces exactly
      if (found) continue;
      // else, found == false, we then use the center of a triangle as the sample point
      // again, we don't consider the case where this triangle lies completely inside a tet face
      // we can change this to make the code more robust
      Vec3i t = cutTriGroup.tri[triCCs[i][0]];
      samplePoint[i] = (cutPositions[t[0]] + cutPositions[t[1]] + cutPositions[t[2]]) / 3.0;
      samplePointER[i] = (cutPositionsER[t[0]] + cutPositionsER[t[1]] + cutPositionsER[t[2]]) / ER(3.0);
    } // end for i on triCC

    vector<TetTriPiece> pieces; // pieces of this tet
    for(size_t i = 0; i < triCCs.size(); i++)
    {
      pieces.emplace_back(tetID, triCCs[i], cutTriGroup, cuttingData, triMeshNormal, tetShape);
    }
//    for (const auto & c : pieces)
//      assert(c.mesh(cuttingData).numTriangles() > 0);

    auto debugPieces = [&]()
    {
      auto tetSurfaceMesh = createTetSurfaceMesh(tetMesh.pos(tetID, 0), tetMesh.pos(tetID, 1), tetMesh.pos(tetID, 2), tetMesh.pos(tetID, 3));
      tetSurfaceMesh.save("tet.obj");
      for(size_t i = 0; i < pieces.size(); i++)
      {
        TriMeshGeo pieceMesh = removeIsolatedVertices(pieces[i].mesh(cuttingData));
        ostringstream os;
        os << "comp" << i << ".obj";
        pieceMesh.save(os.str());
      }
    };

    // if #pieces == 2, the region graph building algorithm can be simplied to
    // doing pseudo-normal inside/outside test only once:
    if (triCCs.size() == 2)
    {
      // check the sample point on pieces[1] against pieces[0]
      int ret = pieces[0].inOutTest(samplePoint[1], samplePointER[1], cuttingData, triMeshNormal, tetShape, profiler);
      assert(ret != 0);

      if (ret > 0) // pieces[1] is outside pieces[0]
      {
        // then because the input triangle mesh is closed manifold, free of self-intersection,
        // there must be two '-' regions in the tet
        piecesInTet[tetID].push_back(tetPieces.size());
        tetPieces.emplace_back(move(pieces[0]));
        piecesInTet[tetID].push_back(tetPieces.size());
        tetPieces.emplace_back(move(pieces[1]));
      }
      else // pieces[1] is inside pieces[0]
      {
        // then because the input triangle mesh is closed manifold, free of self-intersection,
        // there must be one '-' region which has both pieces as its boundaries
        piecesInTet[tetID].push_back(tetPieces.size());
        tetPieces.emplace_back(mergePiece(pieces[0], pieces[1], tetID, cutTriGroup, cuttingData, triMeshNormal, tetShape));
      }
      continue; // go to next iteration in the loop of tetID
    }

    int debugTetID = -1; // used to trigger debugging information
    if (tetID == debugTetID)
    {
      for(size_t i = 0; i < triCCs.size(); i++)
        cout << "sample pt " << samplePoint[i] << endl;
    }

    // more than three pieces in the tet, complicated situation
    // we have to run the full region graph algorithm

    // first initialize the region with one piece (pieces[0]) and two regions
    // it stores: regionID -> list of pairs: <pieceID as its boundary, whether the triangle normals on the piece pointing outward of this region>
    vector<vector<pair<int, bool>>> regions;
    regions.push_back({{0, true}});  // region 0 is inside the pieces[0]
    regions.push_back({{0, false}}); // region 1 is outside the pieces[0]


    auto debugRegion = [&]()
    {
      debugPieces();
      cout << "Dump data for debugging tetID " << tetID << "..." << endl;
      for(size_t i = 0; i < regions.size(); i++)
      {
        cout << "region " << i << ": ";
        for(auto p : regions[i])
          cout << "(" << p.first << ", " << p.second << ") ";
        cout << endl;
      }
    };

    // do the iteration of the region graph algorithm
    // at each iteration, a new piece is added to the graph, splitting an old region in the graph
    for(size_t newPieceID = 1; newPieceID < triCCs.size(); newPieceID++)
    {
      const auto & samplePt = samplePoint[newPieceID];
      const auto & samplePtER = samplePointER[newPieceID];
      const int UNINIT = -999; // lable for uninitialized result
      // vector inOutTest: pieceID of the existing pieces in the graph in the iteration so far -> in/out test result from pieces[newPieceID]
      vector<int> inOutTest(newPieceID, UNINIT);
      int regIDToCut = -1; // the regionID that will be cut by the new piece
      for(size_t regID = 0; regID < regions.size(); regID++)
      {
        bool sampleIsOutside = false;
        for(auto p : regions[regID]) // for each piece serving as the boundary of this region
        {
          int pieceID = p.first;
          if (inOutTest[pieceID] == UNINIT) // check whether the new piece is in/out of the pieces[pieceID]
            inOutTest[pieceID] = pieces[pieceID].inOutTest(samplePt, samplePtER, cuttingData, triMeshNormal, tetShape, profiler);
          // if the newPiece is outside pieces[pieceID] and the normal of pieces[pieceID] is pointing outward the region
          // or the newPiece is inside pieces[pieceID] and the normal of pieces[pieceID] is pointing inward the region
          // then newPiece is outside the region
          if ((inOutTest[pieceID] > 0) == p.second)
          {
            sampleIsOutside = true;
            break;
          }
        }
        if (sampleIsOutside == false) // in this case, the new piece is inside the region
        {
          // the newCCID is inside this region
          regIDToCut = regID; // we shall split this region
          break;
        }
      }
      assert(regIDToCut >= 0);
      if (tetID == debugTetID) // debugging
        cout << "Found cutting CC " << newPieceID << " is in region " << regIDToCut << endl;
      // now we split this region: regions[regIDToCut]
      regions.emplace_back(); // add a new empty region to regions
      vector<pair<int, bool>> regInsideNewPiece; // the boundary piece data for the newly generated region that is inside the new piece
      vector<pair<int, bool>> & regOutsideNewPiece = regions.back(); // the new region that is outside the new piece
      for(auto p : regions[regIDToCut]) // go through all the boundary pieces of the split region
      {
        int pieceID = p.first;
        const auto & samplePt = samplePoint[pieceID];
        const auto & samplePtER = samplePointER[pieceID];

        // check whether pieces[pieceID] is in/out of the new piece
        int test = pieces[newPieceID].inOutTest(samplePt, samplePtER, cuttingData, triMeshNormal, tetShape, profiler);
        assert(test != 0);
        // place pieces[pieceID] to one of the new region based on the in/out result
        if (test > 0) regOutsideNewPiece.emplace_back(pieceID, p.second);
        else regInsideNewPiece.emplace_back(pieceID, p.second);
      }
      // finally, add the new piece to the two new regions
      regInsideNewPiece.emplace_back(newPieceID, true);
      regOutsideNewPiece.emplace_back(newPieceID, false);
      // assign regInsideNewPiece back to the vector: regions
      regions[regIDToCut] = move(regInsideNewPiece);

      if (tetID == debugTetID) // debugging
      {
        cout << "After cut, regions are" << endl;
        for(size_t i = 0; i < regions.size(); i++)
        {
          cout << "region " << i << ": ";
          for(auto p : regions[i])
            cout << "(" << p.first << ", " << p.second << ") ";
          cout << endl;
        }
      }
    } // end each newPieceID, end of the region graph iteraion

    // now regions have been created, we will visit each region to check whether everything is fine
    for(const auto & region : regions)
    {
      assert(region.size() > 0);
      // after the region graph is built, every piece connecting the same region should have
      // the same normal orientation relative to the region
      // if this is broken, sth. is wrong
      bool insideTriMesh = region[0].second;
      for(auto p : region)
      {
        if (p.second != insideTriMesh) { debugRegion(); }
        assert(p.second == insideTriMesh);
      }
    }

    // in the final region graph, each piece belongs to only one '-' region
    // the vector below: piecesVisited is used to check this
    vector<bool> piecesVisited(pieces.size(), false);
    for(const auto & region : regions)
    {
      bool insideTriMesh = region[0].second;

      if (insideTriMesh == false) continue; // if this region is '+', we continue

      // now this region is a '-' region, we should merge all its regions to create the
      // final boundary for this region and store it in tetPieces

      vector<int> groupTriIDs; // store all cut triangles of the merged boundary
      for(auto p : region)
      {
        int pieceID = p.first;
        assert(piecesVisited[pieceID] == false);
        piecesVisited[pieceID] = true;
        if (region.size() == 1)
        {
          piecesInTet[tetID].push_back(tetPieces.size());
          tetPieces.emplace_back(move(pieces[pieceID]));
          break;
        }
        const auto & IDs = pieces[pieceID].groupTriID;
        groupTriIDs.insert(groupTriIDs.end(), IDs.begin(), IDs.end());
      }
      if (region.size() == 1) { continue; }

      sort(groupTriIDs.begin(), groupTriIDs.end());
//      assert(unique(groupTriIDs.begin(), groupTriIDs.end()) == groupTriIDs.end());
      piecesInTet[tetID].push_back(tetPieces.size());
      // TODO: Here we create a merged boundary by first find the union of all cut triangles from several regions,
      // then build a new TetTriPiece from the union.
      // But in this way, internal data built in individual regions is lost and will be recomputed again
      // we can write a better implementation by reusing the data
      tetPieces.emplace_back(tetID, groupTriIDs, cutTriGroup, cuttingData, triMeshNormal, tetShape);
    } // end each region

    // assert all pieces have been added or merged
    assert(find(piecesVisited.begin(), piecesVisited.end(), false) == piecesVisited.end());
  } // end each tetID
  if (verbose)
    cout << "Done." << endl;
  buildCopiesSec.stop();

  connectCopiesSec.start();

  // now each piece stored in tetPieces corresponds to a virtualized tet to appear in the final tet mesh
  // we will begin connect those tets

  // first, we build a structure to find neighboring tets and their shared tet faces
  vector<map<UTriKey, int>> tetNbrs(tetMesh.numTets()); // input tetID -> a tet face -> the neigbor tetID sharing this face
  map<UTriKey, vector<int>> tetFaceMap; // tet face -> all tets sharing this tet face
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    Vec4i tet = tetMesh.tet(tetID);
    UTetKey tetkey(&tet[0]);
    for(int i = 0; i < 4; i++)
    {
      UTriKey face = tetkey.uFaceKey(i);
      tetFaceMap[face].push_back(tetID);
    }
  }

  for(const auto & p : tetFaceMap)
  {
    assert(p.second.size() <= 2);
    const auto & tetIDs = p.second;
    if (tetIDs.size() == 2)
    {
      // we only store tetIDs[0] -> tet face -> tetIDs[1] and don't store tetIDs[1] -> tet face -> tetIDs[0]
      // because later we need to visit each tet connection only once
      tetNbrs[tetIDs[0]].emplace(p.first, tetIDs[1]);
    }
  }

  vector<Vec3d> tetCenters(tetMesh.numTets(), Vec3d(0.0));
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    tetCenters[tetID] = tetMesh.ref().computeTetCenter(tetID);
  }

  // We build a DisjointSet for each input tet mesh vertex
  // each tet copy is first assigned with four unique tet copy vertex IDs
  // these unique tet copy vertex IDs are not explicitly built and stored on memory, because they are just:
  // [tetCopyID * 4, tetCopyID * 4 + 3] for tetCopyID
  // Then the unique tet copy vertex IDs belonging to the same input tet vertex will be thrown into the disjointSet
  // structure of that input tet vertex
  // if any tet copy face is merged at the input tet vertex, corresponding unique tet copy vertex IDs are merged in
  // the disjoint set
  vector<DisjointSetDynamic> dset(tetMesh.numVertices());

  // merge the vertices of the two tet copies: copyID and nbrCopyID on the shared face
  auto mergeVtxOnTetFace = [&](int copyID, int nbrCopyID, UTriKey face)
  {
    assert(copyID >= 0 && copyID < (int)tetPieces.size());
    assert(nbrCopyID >= 0 && nbrCopyID < (int)tetPieces.size());
    int tetID = tetPieces[copyID].tetID;
    int nbrtetID = tetPieces[nbrCopyID].tetID;
    Vec4i tet = tetMesh.tet(tetID);
    Vec4i nbrtet = tetMesh.tet(nbrtetID);

    for(int vtxID = 0; vtxID < 3; vtxID++)
    {
      int v = face[vtxID];
      int i = tet.getInvertedIndex(v);
      assert(i >= 0);
      int nbri = nbrtet.getInvertedIndex(v);
      assert(nbri >= 0);
      // now the two tet copy vertex to be merged are: 4*copyID + i, 4*nbrCopyID + nbri
      dset[v].unionSet(4*copyID + i, 4*nbrCopyID + nbri); // merge them in the disjoint set
    }
  };

  // use to store the potential pairs of tet copy IDs and their shared faces for later processing
  vector<tuple<int,int,UTriKey>> mergeCandList;
  for(size_t copyID = 0; copyID < tetPieces.size(); copyID++) // for each virtualied tet: tet copy
  {
    int tetID = tetPieces[copyID].tetID; // tetID on the input tet mesh

    UTetKey tetkey(tetMesh.tet(tetID));
    for(int faceID = 0; faceID < 4; faceID++)
    {
      UTriKey face = tetkey.uFaceKey(faceID);
      auto it = tetNbrs[tetID].find(face);
      if (it == tetNbrs[tetID].end()) continue;
      int nbrTetID = it->second; // find a neighboring input tet
//      TetShape nbrTetShape(tetMesh, nbrTetID);

      for(int nbrCopyID : piecesInTet[nbrTetID]) // for each tet copy on this nbrTetID
      {
        if (intersectSorted(tetPieces[copyID].vtx, tetPieces[nbrCopyID].vtx)) // if the tet copies share a cut vtx
        {
          // we shall merge tet copy vertices at this face
          mergeVtxOnTetFace(copyID, nbrCopyID, face);
          continue;
        }
        // although the two copies are not connected through cut triMesh vtx,
        // they can still be connected through "inner space"

        // if any piece vtx is on this face, then this face is not the part of the inner space of this region
        // touchedTetFaces stores all the tet faces that the piece vtx lie on
        if (setFind(tetPieces[copyID].touchedTetFaces, face)) continue;
        if (setFind(tetPieces[nbrCopyID].touchedTetFaces, face)) continue;

        // store them for later processing
        mergeCandList.emplace_back(copyID, nbrCopyID, face);
      } // end nbrCopyID
    }
  }

  // stores whether the candidate nbring pair should be connected
  // we use char as elements of the vector because it is thread-safe
  // Note that vector<bool> is not thread-safe when modifying individual elements concurrently
  vector<char> mergeCandResult(mergeCandList.size(), false);

#ifdef USE_VT_MULTITHREADING
  mutex inOutTestMutex;
  vector<mutex> copyInitMutex(tetPieces.size());
  // use tbb here
#endif
  for(size_t i = 0; i < mergeCandList.size(); i++)
  {
    int copyID = get<0>(mergeCandList[i]);
    int nbrCopyID = get<1>(mergeCandList[i]);
    UTriKey face = get<2>(mergeCandList[i]);
    int tetID = tetPieces[copyID].tetID;
    int nbrTetID = tetPieces[nbrCopyID].tetID;
    TetShape tetShape(tetMesh, tetID);
    TetShape nbrTetShape(tetMesh, nbrTetID);

    // here we check whether the two tet copies can be connected via inner space
    ProfilerExtraSection InOutTestInConnectionSec(profiler, "inOutTestInConnection");
    InOutTestInConnectionSec.start();
    // if this tet copy has no pieces, then its an inner tet, each of whose face is inner space
    // otherwise, we use the pieces to do in/out check on the face center to determine
    // whether this face is inside/outside the pieces, thus the region this tet copy represents
    if (tetPieces[copyID].vtx.size() > 0)
    {
      Vec3d faceCenter = (tetMesh.pos(face[0]) + tetMesh.pos(face[1]) + tetMesh.pos(face[2])) / 3.0;
      Vec3d queryPt = 1e-6 * norm(faceCenter - tetCenters[tetID]) + faceCenter;
      // we offset queryPoint a little bit away from faceCenter, to account for numerical errors in inside/outside test
      Vec3ER queryPtK(queryPt[0], queryPt[1], queryPt[2]);
#ifdef USE_VT_MULTITHREADING
      lock_guard<mutex> lock(copyInitMutex[copyID]);
#endif
      int test = tetPieces[copyID].inOutTest(queryPt, queryPtK, cuttingData, triMeshNormal, tetShape, profiler);
      assert(test != 0);
      if (test > 0) continue; // the tet face is outside, therefore not inner space
    }
    if (tetPieces[nbrCopyID].vtx.size() > 0)
    {
      Vec3d faceCenter = (tetMesh.pos(face[0]) + tetMesh.pos(face[1]) + tetMesh.pos(face[2])) / 3.0;
      Vec3d queryPt = 1e-6 * norm(faceCenter - tetCenters[nbrTetID]) + faceCenter;
      // we offset queryPoint a little bit away from faceCenter, to account for numerical errors in inside/outside test
      Vec3ER queryPtK(queryPt[0], queryPt[1], queryPt[2]);
#ifdef USE_VT_MULTITHREADING
      lock_guard<mutex> lock(copyInitMutex[nbrCopyID]);
#endif
      int test = tetPieces[nbrCopyID].inOutTest(queryPt, queryPtK, cuttingData, triMeshNormal, nbrTetShape, profiler);
      assert(test != 0);
      if (test > 0) continue; // outside
    }
    InOutTestInConnectionSec.stop();
    mergeCandResult[i] = true;
  }

  for(size_t i = 0; i < mergeCandList.size(); i++)
  {
    if (mergeCandResult[i])
    {
      int copyID = get<0>(mergeCandList[i]);
      int nbrCopyID = get<1>(mergeCandList[i]);
      UTriKey face = get<2>(mergeCandList[i]);
      mergeVtxOnTetFace(copyID, nbrCopyID, face);
    }
  }


  // now we will assign the new tet copy vertex ID to the merged unique tet copy vertex IDs
  vector<Vec4i> outputTets(tetPieces.size()); // stores the new tet copy vertex IDs of each tet copy
  vector<Vec3d> outputTetPositions;           // new tet copy vertex ID -> its position
  // input tet mesh vtx ID -> the representive setID in the disjiont set ->  new tet copy vtx ID
  vector<map<int, int>> repSetID2NewID(tetMesh.numVertices());
  for(size_t tetCopyID = 0; tetCopyID < tetPieces.size(); tetCopyID++)
  {
    for(int j = 0; j < 4; j++)
    {
      int uniqueID = (int)tetCopyID * 4 + j; // the original unique tet copy vtx ID
      int oldVtxID = tetMesh.tet(tetPieces[tetCopyID].tetID)[j]; // the input tet mesh vtx ID
      DisjointSetDynamic & vdset = dset[oldVtxID];
      int repSetID = vdset.findSet(uniqueID); // get the representive setID of this uniqueID
      auto it = repSetID2NewID[oldVtxID].find(repSetID);
      int newID = -1;
      if (it == repSetID2NewID[oldVtxID].end()) // this is a new representive setID not met before
      {
        newID = outputTetPositions.size();      // we then create a new tet copy vtx ID
        repSetID2NewID[oldVtxID].emplace(repSetID, newID);
        outputTetPositions.push_back(tetMesh.pos(oldVtxID));
      }
      else
      {
        newID = it->second;
      }
      outputTets[tetCopyID][j] = newID;        // store the new tet copy vtx ID to outputTets

      if (newTetVtxID2OldTetVtxID) // create the mapping: new tet copy vtx ID -> input tet vtx ID for optional output of this function
      {
        if ((*newTetVtxID2OldTetVtxID).size() <= (size_t)newID)
        {
          (*newTetVtxID2OldTetVtxID).resize(newID+1);
        }
        (*newTetVtxID2OldTetVtxID)[newID] = tetMesh.tet(tetPieces[tetCopyID].tetID)[j];
      }
    }
  }

  // After connecting neighboring tet copies, some tet copies might have the same set of tet copy vertices.
  // This means that they should be merged.
  // We will do this merge here.
  map<Vec4i, int> tetCopyVtx2CopyID; // tet copy vertices -> tet copy IDs
  for(size_t copyID = 0; copyID < tetPieces.size(); copyID++)
  {
    Vec4i tet = outputTets[copyID];
    auto it = tetCopyVtx2CopyID.find(tet);
    if (it == tetCopyVtx2CopyID.end())
      tetCopyVtx2CopyID.emplace(tet, copyID);
    else
    { // merge the tet copy stored in tetCopyVtx2CopyID and this tet copy with index of copyID
      int prevCopyID = it->second;
      // naturally, the two copies to be merged should come from the same input tet
      if (!(tetPieces[prevCopyID].tetID == tetPieces[copyID].tetID))
      {
        cout << "copyID " << copyID << " preCopyID " << prevCopyID << " tet " << tet << " tetID" << tetPieces[copyID].tetID
            << " pre tetID " << tetPieces[prevCopyID].tetID << endl;
      }
      assert(tetPieces[prevCopyID].tetID == tetPieces[copyID].tetID);
      int tetID = tetPieces[copyID].tetID;
      TetShape tetShape(tetMesh, tetID);
      const auto & cutTriGroup = cuttingData.cutTriGroups[tetID];
      // merge the pieces of the two tet copies
      tetPieces[prevCopyID] = mergePiece(tetPieces[prevCopyID], tetPieces[copyID], tetID, cutTriGroup,
          cuttingData, triMeshNormal, tetShape);
    }
  }

  // now we have merged pieces of those copies with same tet copy vtx,
  // we should remove those redundant tet copies in outputTets and tetPieces
  outputTets.clear();
  vector<TetTriPiece> newTetPieces;
  for(auto p : tetCopyVtx2CopyID)
  {
    outputTets.push_back(p.first);
    if (newTetID2OldTetID) newTetID2OldTetID->push_back(tetPieces[p.second].tetID);
    newTetPieces.emplace_back(move(tetPieces[p.second]));
  }
  // now newTetPieces stores the final tetPieces and outputTets stores the final tet copies

  connectCopiesSec.stop();
 //  cout << if (profiler) profiler->toString() << endl;

  finalStuffSec.start();

  // build the final tet mesh
  TetMeshGeo newTetMesh(move(outputTetPositions), move(outputTets));

  if (coord) // build optional data, coord
  {
    // compute barycentric weights for embedding the input traingle mesh into the new virtualized tet mesh
    vector<double> weights(triMesh.numVertices() * 4);
    vector<int> tetVtxIndices(triMesh.numVertices() * 4, -1);
    vector<int> embeddingTetIndices(triMesh.numVertices());
    for(size_t copyID = 0; copyID < newTetPieces.size(); copyID++)
    {
      int tetID = newTetPieces[copyID].tetID;
      for(int v : newTetPieces[copyID].vtx)
      {
        if (v >= triMesh.numVertices()) continue; // this tri vtx is a cut vtx
        if (tetVtxIndices[4*v] >= 0) continue; // this tri vtx has been assigned a weight before
        getTetBarycentricWeights(triMesh.pos(v), tetMesh.pos(tetID, 0), tetMesh.pos(tetID, 1), tetMesh.pos(tetID, 2), tetMesh.pos(tetID, 3),
            &weights[4*v]);
        newTetMesh.tet(copyID).convertToArray(&tetVtxIndices[4*v]);
        embeddingTetIndices[v] = copyID;
      }
    }

    // if some input triangle vtx are outside the tet mesh,
    // we will find the closet tets to embed them
    vector<int> outerTriVtxIDs;
    for(int triVtxID = 0; triVtxID < triMesh.numVertices(); triVtxID++)
    {
      if (tetVtxIndices[4*triVtxID] >= 0) continue;
      outerTriVtxIDs.push_back(triVtxID);
    }

    if (outerTriVtxIDs.size() > 0)
    {
      TetMeshRef newTetMeshRef = newTetMesh.ref();
  #ifdef USE_TBB
      tbb::parallel_for(0, sizei(outerTriVtxIDs), [&](const int & i)
  #else
      for(int i = 0; i < sizei(outerTriVtxIDs); i++)
  #endif
      {
        int triVtxID = outerTriVtxIDs[i];
        // now this vtx is not embedded because it is out of tet mesh
        Vec3d pos = triMesh.pos(triVtxID);

        MinValueIndex vi;
        for(int newTetID = 0; newTetID < newTetMesh.numTets(); newTetID++)
        {
          vi.update(newTetMeshRef.computeSquaredDistanceToTet(newTetID, pos), newTetID);
        }

        assert(vi.index >= 0);
        newTetMeshRef.computeTetBarycentricWeights(vi.index, pos, &weights[4*triVtxID]);
        memcpy(&tetVtxIndices[4*triVtxID], newTetMesh.tet(vi.index).data(), sizeof(int) * 4);
        embeddingTetIndices[triVtxID] = vi.index;
  #ifdef USE_TBB
      });
  #else
      }
  #endif
    }

    BarycentricCoordinates bary(triMesh.numVertices(), 4, tetVtxIndices.data(), weights.data(), embeddingTetIndices.data());
    *coord = move(bary);
  }

  if (tetTris) // build optional data, tetTris
  {
    tetTris->resize(newTetMesh.numTets());
    for(int tetCopyID = 0; tetCopyID < newTetMesh.numTets(); tetCopyID++)
    {
      auto tris = newTetPieces[tetCopyID].gnrOriTriID;
      sortAndDeduplicate(tris);
      (*tetTris)[tetCopyID] = move(tris);
    }
  }

  finalStuffSec.stop();

  return newTetMesh;
}
