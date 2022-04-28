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

#include "virtualTets-via-csg.h"
#include "triMeshNeighbor.h"
#include "triMeshPseudoNormal.h"
#include "containerHelper.h"
#include "disjointSet.h"
#include "tetTriMeshCutting.h"
#include "exactOctree.h"
#include "basicAlgorithms.h"
#include "profiler.h"
#include "windingNumberTree.h"
#include "tetTriCuttingData.h"
#include "createTriMesh.h"
#include "geometryQuery.h"
#include "planeER.h"
#include "iglRemeshSelfIntersection.h"
#include "labelOuterTets.h"
#include "tetTriPiece.h"
#include "valueIndex.h"
#include <iostream>
#include <sstream>
#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif

using namespace std;

class TetER
{
public:
  TetER(const Vec3ER & v0, const Vec3ER & v1, const Vec3ER & v2, const Vec3ER & v3);
  int outside(const Vec3ER & v) const;

protected:
  FastPlaneER plane[4];
};

TetER::TetER(const Vec3ER & v0, const Vec3ER & v1, const Vec3ER & v2, const Vec3ER & v3)
: plane { {v1, v2, v3}, {v0, v3, v2}, {v0, v1, v3}, {v0, v2, v1} }
{}

int TetER::outside(const Vec3ER & v) const
{
  for(int i = 0; i < 4; i++)
  {
    int o = plane[i].outside(v);
    if (o >= 0) return o;
  }
  return -1;
}

TetMeshGeo createVirtualTetsMeshViaCSG(const TetMeshGeo & tetMesh, const TriMeshGeo triMesh, BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris,
    bool verbose, Profiler * profiler)
{
  vector<Vec3ER> triPosK;
  for(int i = 0; i < triMesh.numVertices(); i++)
  {
    Vec3d p = triMesh.pos(i);
    triPosK.emplace_back(p[0], p[1], p[2]);
  }
  return createVirtualTetsMeshViaCSG(tetMesh, triMesh, triPosK, coord, newTetVtxID2OldTetVtxID, newTetID2OldTetID, tetTris,
      verbose, profiler);
}

struct TetTriCompData
{
  int tetID = -1;
  vector<int> vtx; // original tri mesh vtxIDs in this tet copy
  vector<int> tri; // original tri mesh triIDs intersected with this tet copy
  map<UTriKey, int> tetFaceRel; // tet face relation: tet face key -> 1: face is outside this comp, 0: touching, -1: inside
  void check() const
  {
    assert(tetID >= 0);
    for(int v : vtx) { assert(v >= 0); }
    for(int t : tri) { assert(t >= 0); }
  }
};


TetMeshGeo createVirtualTetsMeshViaCSG(const TetMeshGeo & tetMesh, const TriMeshGeo triMesh,
    const std::vector<Vec3ER> & exactTriPos, BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris,
    bool verbose, Profiler * profiler)
{
  if (profiler) profiler->startTimer("tetTriEmbedding");
  TetTriMeshCutting cutting(tetMesh, triMesh);
  cout << "computing embedding..." << endl;
  cutting.computeIntersectingTriTets();
  cout << "computing cut..." << endl;
  if (profiler) profiler->stopLastTimer();

  TetTriIntersectingData data;
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    data.triInTet.push_back(cutting.getTrianglesIntersectingTet(tetID));
  }
  for(int vtxID = 0; vtxID < tetMesh.numVertices(); vtxID++)
  {
    Vec3d p = tetMesh.pos(vtxID);
    Vec3ER pER(p[0], p[1], p[2]);
    data.tetPosER.emplace_back(move(pER));
  }

  return createVirtualTetsMeshViaCSG(tetMesh, triMesh, exactTriPos, data, coord, newTetVtxID2OldTetVtxID, newTetID2OldTetID, tetTris,
      verbose, profiler);
}


TetMeshGeo  createVirtualTetsMeshViaCSG(const TetMeshGeo & tetMesh, const TriMeshGeo triMesh,
    const std::vector<Vec3ER> & exactTriPos, const TetTriIntersectingData & tetTriIntersectingData,
    BarycentricCoordinates * coord,
    vector<int> * newTetVtxID2OldTetVtxID, vector<int> * newTetID2OldTetID, vector<vector<int>> * tetTris,
    bool verbose, Profiler * profiler)
    {

  cout << "Begin virtual tets..." << endl;
  if (newTetVtxID2OldTetVtxID) newTetVtxID2OldTetVtxID->clear();
  if (newTetID2OldTetID) newTetID2OldTetID->clear();

  TriMeshNeighbor triMeshNeighbor;
  TriangleNeighbor triangleNeighbor;
  try
  {
    triMeshNeighbor = TriMeshNeighbor(triMesh);
    triangleNeighbor = TriangleNeighbor(triMesh.triangles());
  } catch(int)
  {
    cout << "Error: the input mesh is not edge-manifold" << endl;
    throw 1;
  }

  auto triangleBoundaryLoops = triangleNeighbor.findBoundaryLoops(triMesh.triangles());
  if (triangleBoundaryLoops.size() > 0)
  {
    cout << "Error: the input mesh has boundary" << endl;
    throw 1;
  }

  if (profiler) profiler->startTimer("airTets");

  auto isTetBou = [&](int tetID)
  {
    return tetTriIntersectingData.triInTet[tetID].size() > 0;
  };

  WindingNumberTree wnTree;
  wnTree.build(triMesh);
  auto isTetOuter = [&](int tetID)
  {
    double wn = wnTree.windingNumber(triMesh, tetMesh.ref().computeTetCenter(tetID));
    return wn < 0.5;
  };

  TetNeighbor tetNeighbor(tetMesh.tets());
  auto tetOutside = labelOuterTets(tetMesh, tetNeighbor, isTetBou, isTetOuter);
  if (profiler) profiler->stopLastTimer();

  vector<TetTriCompData> tetCopies;
  if (profiler) profiler->startTimer("buildComps");

  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    bool debug = false;
    if (tetOutside[tetID]) continue;
    if (tetID % 100 == 0)
      cout << "process comps in tetID " << tetID << "/" << tetMesh.numTets() << endl;

//    if (tetID == 828) debug = true;

    const vector<int> & interTriSet = tetTriIntersectingData.triInTet[tetID];
    if (debug)
    {
      cout << "intersecting tris: " << streamRange(interTriSet) << endl;
    }
    if (interTriSet.size() == 0)
    {
      // inner tet:
      TetTriCompData data;
      data.tetID = tetID;
      TetShape tetShape(tetMesh, tetID);
      for(auto p : tetShape.faceNormal)
      {
//        cout << "inner tet " << tetID << " face " << p.first << endl;
        data.tetFaceRel[p.first] = -1;
      }
      tetCopies.emplace_back(move(data));
      continue;
    }

    TriMeshGeo tetSurface = createTetSurfaceMesh(tetMesh.pos(tetID, 0), tetMesh.pos(tetID, 1), tetMesh.pos(tetID, 2), tetMesh.pos(tetID, 3));
    vector<Vec3i> subTris;
    vector<int> inputTriIDs(interTriSet.begin(), interTriSet.end());
    getSubMesh(triMesh, inputTriIDs, subTris);
    TriMeshRef subMesh(triMesh.positions(), subTris);
    vector<int> sub2triVtx;
    map<int,int> tri2subVtx;
    TriMeshGeo rivSubMesh = removeIsolatedVerticesWithMap(subMesh, &sub2triVtx, &tri2subVtx);
    TriMeshGeo mergedMesh = mergeMesh(tetSurface, rivSubMesh);

    ProfilerExtraSection cellSelfInterTimer(profiler, "cellSelfInter");
    cellSelfInterTimer.start();
    auto selfInterRet = iglInterface::remeshSelfIntersection(mergedMesh, false, false, debug);
    cellSelfInterTimer.stop();

    const auto & cutMesh = selfInterRet.cutMesh;

    if (debug)
    {
      cutMesh.save("debugViaCSGCutMesh.obj");
    }

    map<int, vector<int>> cellID2cutTriIDs;
    map<int, set<UTriKey>> cellID2InnerTetFace;
    map<int, set<UTriKey>> cellID2TouchingTetFace;
    for(int cutTriID = 0; cutTriID < cutMesh.numTriangles(); cutTriID++)
    {
      int patchID = selfInterRet.triPatchIDs[cutTriID];
      int cellID = selfInterRet.cellIDsAtPatch[patchID].second; // we only get the cell inside this patch
      int mergedTriID = selfInterRet.oldTriIDs[cutTriID];
      if (mergedTriID < 4) // these are from tet faces
      {
        if (debug)
        {
          cout << "found tet surface on cut mesh: " << cutTriID << endl;
        }
        bool allFromTetFace = true;
        for(int i = 0; i < 3; i++)
        {
          if (cutMesh.tri(cutTriID)[i] >= 4) { allFromTetFace = false; break; }
        }
        Vec3i tetSurfaceTri = tetSurface.tri(mergedTriID);  // [0, 4)
        Vec3i tetFace;
        for(int i = 0; i < 3; i++)
          tetFace[i] = tetMesh.tet(tetID)[tetSurfaceTri[i]];
        if (allFromTetFace)
        {
          if (debug && setNotFind(cellID2InnerTetFace[cellID], UTriKey(tetFace)))
            cout << "cellID " << cellID << " has inner face " << UTriKey(tetFace) << endl;
          cellID2InnerTetFace[cellID].emplace(tetFace);
        }
        else
        {
          if (debug && setNotFind(cellID2TouchingTetFace[cellID], UTriKey(tetFace)))
            cout << "cellID " << cellID << " has touching face " << UTriKey(tetFace) << endl;
          cellID2TouchingTetFace[cellID].emplace(tetFace);
        }
      }
      else // these are from tri mesh
      {
        if (cellID == 0) continue; // these are triangles outside the tet
        cellID2cutTriIDs[cellID].push_back(cutTriID);
      }
    }

    const auto & tetPosER = tetTriIntersectingData.tetPosER;
    TetER tetER(tetPosER[tetMesh.tetVtxID(tetID, 0)], tetPosER[tetMesh.tetVtxID(tetID, 1)],
                tetPosER[tetMesh.tetVtxID(tetID, 2)], tetPosER[tetMesh.tetVtxID(tetID, 3)]);
    map<int, int> cutVtxID2Outside;
    if (debug)
    {
      for(const auto & p : cellID2cutTriIDs)
      {
        int cellID = p.first;
        cout << "At cell " << cellID << " tris " << streamRange(p.second) << ", vtx inout: ";
        set<int> cutTriVtx;
        for(int cutTriID : p.second)
        {
          for(int cutVtxID : cutMesh.tri(cutTriID))
          {
            cutTriVtx.insert(cutVtxID);
          }
        }

        for(int cutVtxID : cutTriVtx)
        {
          auto it = cutVtxID2Outside.find(cutVtxID);
          if (it == cutVtxID2Outside.end())
          {
            int o = tetER.outside(selfInterRet.cutPosExact[cutVtxID]);
            cutVtxID2Outside[cutVtxID] = o;
          }
          cout << cutVtxID << "->" << cutVtxID2Outside[cutVtxID] << " ";
        }
        cout << endl;
      }
    }
    for(const auto & p : cellID2cutTriIDs)
    {
      int cellID = p.first;

      // first check whether it's outside or inside the tet
      int cellInOut = 0;
      for(int cutTriID : p.second)
      {
        for(int cutVtxID : cutMesh.tri(cutTriID))
        {
          auto it = cutVtxID2Outside.find(cutVtxID);
          if (it != cutVtxID2Outside.end())
          {
            int p = it->second;
            if (p != 0) // outside or inside
            {
              cellInOut = p;
              break;
            }
          }
          // compute inOut of this cutVtxID
          int o = tetER.outside(selfInterRet.cutPosExact[cutVtxID]);
          cutVtxID2Outside[cutVtxID] = o;
          if (o != 0)
          {
            cellInOut = o;
            break;
          }
        }
        if (cellInOut != 0) break;
      }

      if (cellInOut > 0) continue; // outer patch

      if (debug)
      {
        cout << "CellID " << cellID << " " << streamRange(p.second) << endl;
      }

      TetTriCompData data;
      data.tetID = tetID;

      for(int cutTriID : p.second)
      {
        for(int j = 0; j < 3; j++)
        {
          int cutVtx = cutMesh.tri(cutTriID)[j];
          // if cutVtx is tet surface vtx or vtx from cut itnersections
          if (cutVtx < 4 || cutVtx >= mergedMesh.numVertices()) continue;
          int subVtxID = cutVtx - 4;
          assert(subVtxID < rivSubMesh.numVertices());
          int triVtxID = sub2triVtx[subVtxID];
          assert(inVectorRange(triVtxID, triMesh.positions()));
          data.vtx.push_back(triVtxID);
        }

        int mergedTriID = selfInterRet.oldTriIDs[cutTriID];
        assert(mergedTriID >= 4);
        int subMeshTriID = mergedTriID - 4;
        assert(inVectorRange(subMeshTriID, inputTriIDs));
        int inputTriID = inputTriIDs[subMeshTriID];
        data.tri.push_back(inputTriID);
      }

      sortAndDeduplicate(data.vtx);
      sortAndDeduplicate(data.tri);

      for(UTriKey key : cellID2InnerTetFace[cellID])
      {
        assert(mapNotFind(data.tetFaceRel, key));
        data.tetFaceRel[key] = -1;
      }
      for(UTriKey key : cellID2TouchingTetFace[cellID])
      {
        if (mapNotFind(data.tetFaceRel, key) == false)
        {
          cout << "Error: found tet face " << key << " at both inner and touching at tetID " << tetID << endl;
          TetShape tetShape(tetMesh, tetID);
          cout << "All tet faces: " ;
          for(auto p : tetShape.faceNormal) cout << p.first << " ";
          cout << endl;
          cout << "cellID2InnerTetFace: ";
          for(auto key : cellID2InnerTetFace[cellID]) cout << key << " ";
          cout << endl;
          cout << "cellID2TouchingTetFace: ";
          for(auto key : cellID2TouchingTetFace[cellID]) cout << key << " ";
          cout << endl;
        }
        assert(mapNotFind(data.tetFaceRel, key));
        data.tetFaceRel[key] = 0;
      }

      data.check();
      tetCopies.emplace_back(move(data));
    } // end cellID2cutTriIDs
  } // end for each tetID

//  for(const auto & c : tetCopies) c.check();
  if (profiler) profiler->stopLastTimer();

  if (profiler) profiler->startTimer("connectComps");
  vector<vector<int>> compsInTet(tetMesh.numTets());
  for(size_t i = 0; i < tetCopies.size(); i++)
  {
    compsInTet[tetCopies[i].tetID].push_back(i);
  }

  // now let's connect those tet copies
  vector<map<UTriKey, int>> tetNbrs(tetMesh.numTets());
  map<UTriKey, vector<int>> tetFaceMap;
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
    const auto & nbr = p.second;
    if (nbr.size() == 2)
    {
      tetNbrs[nbr[0]].emplace(p.first, nbr[1]); // we need to visit each tet connection only once
    }
  }

  vector<Vec4i> outputTets(tetCopies.size());

  vector<DisjointSetDynamic> dset(tetMesh.numVertices());

  auto mergeVtxOnTetFace = [&](int copyID, int nbrCopyID, UTriKey face)
  {
    assert(copyID >= 0 && copyID < (int)tetCopies.size());
    assert(nbrCopyID >= 0 && nbrCopyID < (int)tetCopies.size());
    int tetID = tetCopies[copyID].tetID;
    int nbrtetID = tetCopies[nbrCopyID].tetID;
    Vec4i tet = tetMesh.tet(tetID);
    Vec4i nbrtet = tetMesh.tet(nbrtetID);

    for(int vtxID = 0; vtxID < 3; vtxID++)
    {
      int v = face[vtxID];
      int i = tet.getInvertedIndex(v);
      assert(i >= 0);
      int nbri = nbrtet.getInvertedIndex(v);
      assert(nbri >= 0);
      dset[v].unionSet(4*copyID + i, 4*nbrCopyID + nbri);
    }
  };

  for(size_t copyID = 0; copyID < tetCopies.size(); copyID++)
  {
    int tetID = tetCopies[copyID].tetID;
    TetShape tetShape(tetMesh, tetID);
    UTetKey tetkey(tetMesh.tet(tetID));
    for(int faceID = 0; faceID < 4; faceID++)
    {
      UTriKey face = tetkey.uFaceKey(faceID);
      auto it = tetNbrs[tetID].find(face);
      if (it == tetNbrs[tetID].end()) continue;
      int nbrTetID = it->second;
      TetShape nbrTetShape(tetMesh, nbrTetID);

      for(int nbrCopyID : compsInTet[nbrTetID])
      {
        // TODO: this is not safe, if one triangle is exactly touch the tet face at a triangle edge,
        // then this could miss
        if (intersectSorted(tetCopies[copyID].tri, tetCopies[nbrCopyID].tri)
            || intersectSorted(tetCopies[copyID].vtx, tetCopies[nbrCopyID].vtx))
            {
          // connect them
          mergeVtxOnTetFace(copyID, nbrCopyID, face);
          continue;
        }
        // although the two copies are not connected through cut triMesh vtx,
        // they can still be connected through inner space
        auto faceIt = tetCopies[copyID].tetFaceRel.find(face);
        if (faceIt == tetCopies[copyID].tetFaceRel.end() || faceIt->second != -1) continue;
        auto nbrFaceIt = tetCopies[nbrCopyID].tetFaceRel.find(face);
        if (nbrFaceIt == tetCopies[nbrCopyID].tetFaceRel.end() || nbrFaceIt->second != -1) continue;

//        cout << "Found tet merging..." << endl;
        mergeVtxOnTetFace(copyID, nbrCopyID, face);
      } // end nbrCopyID
    }
  }

  vector<Vec3d> outputTetPositions;
  vector<map<int, int>> repSetID2NewID(tetMesh.numVertices());
  for(size_t i = 0; i < tetCopies.size(); i++)
  {
    for(int j = 0; j < 4; j++)
    {
      int uniqueID = (int)i * 4 + j;
      int oriVtxID = tetMesh.tet(tetCopies[i].tetID)[j];
      DisjointSetDynamic & vdset = dset[oriVtxID];
      int repSetID = vdset.findSet(uniqueID);
      auto it = repSetID2NewID[oriVtxID].find(repSetID);
      int newID = -1;
      if (it == repSetID2NewID[oriVtxID].end())
      {
        newID = outputTetPositions.size();
        repSetID2NewID[oriVtxID].emplace(repSetID, newID);
        outputTetPositions.push_back(tetMesh.pos(oriVtxID));
      }
      else
      {
        newID = it->second;
      }
      outputTets[i][j] = newID;

      if (newTetVtxID2OldTetVtxID)
      {
        if ((*newTetVtxID2OldTetVtxID).size() <= (size_t)newID)
        {
          (*newTetVtxID2OldTetVtxID).resize(newID+1);
        }
        (*newTetVtxID2OldTetVtxID)[newID] = tetMesh.tet(tetCopies[i].tetID)[j];
      }
    }
  }

  // remove redundant tet copies
  map<Vec4i, int> removeMap;
  for(size_t copyID = 0; copyID < tetCopies.size(); copyID++)
  {
    Vec4i tet = outputTets[copyID];
    auto it = removeMap.find(tet);
    if (it == removeMap.end())
      removeMap.emplace(tet, copyID);
    else // merge the stored tetCopy and this tetCopy
    {
      int prevCopyID = it->second;
      if (!(tetCopies[prevCopyID].tetID == tetCopies[copyID].tetID))
      {
        cout << "copyID " << copyID << " preCopyID " << prevCopyID << " tet " << tet << " tetID" << tetCopies[copyID].tetID
            << " pre tetID " << tetCopies[prevCopyID].tetID << endl;
      }
      assert(tetCopies[prevCopyID].tetID == tetCopies[copyID].tetID);

      vectorInsertRangeBack(tetCopies[prevCopyID].vtx, tetCopies[copyID].vtx);
      vectorInsertRangeBack(tetCopies[prevCopyID].tri, tetCopies[copyID].tri);
      sortAndDeduplicate(tetCopies[prevCopyID].vtx);
      sortAndDeduplicate(tetCopies[prevCopyID].tri);

      tetCopies[prevCopyID].check();
    }
  }

  outputTets.clear();
  vector<TetTriCompData> newTetCopies;
  for(auto p : removeMap)
  {
    outputTets.push_back(p.first);
    if (newTetID2OldTetID) newTetID2OldTetID->push_back(tetCopies[p.second].tetID);
    newTetCopies.emplace_back(move(tetCopies[p.second]));
  }
  for(const auto & c : newTetCopies)
  {
    c.check();
  }

  if (profiler) profiler->stopLastTimer();
//  cout << if (profiler) profiler->toString() << endl;

  TetMeshGeo newTetMesh(move(outputTetPositions), move(outputTets));

  if (coord)
  {
    // compute barycentric weights
    vector<double> weights(triMesh.numVertices() * 4);
    vector<int> tetVtxIndices(triMesh.numVertices() * 4, -1);
    vector<int> embeddingTetIndices(triMesh.numVertices());
    for(size_t copyID = 0; copyID < newTetCopies.size(); copyID++)
    {
      int tetID = newTetCopies[copyID].tetID;
      for(int v : newTetCopies[copyID].vtx)
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
    *coord = bary;
  }

  if (tetTris)
  {
    tetTris->resize(newTetMesh.numTets());
    for(int tetCopyID = 0; tetCopyID < newTetMesh.numTets(); tetCopyID++)
    {
      auto tris = newTetCopies[tetCopyID].tri;
      sortAndDeduplicate(tris);
      (*tetTris)[tetCopyID] = move(tris);
    }
  }

  // XXX
//  exit(0);

  return newTetMesh;
}
