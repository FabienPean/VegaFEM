/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "immersionMesher" library , Copyright (C) 2018 USC                    *
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

#include "disjointSet.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "tetKey.h"
#include "vec4d.h"
#include "triMeshGeo.h"
#include "virtualTets.h"
#include "virtualTets-via-csg.h"
#include "immersionMesher.h"
#include "windingNumberTree.h"
#include "labelOuterTets.h"
using namespace std;

namespace
{

struct PointInTet
{
  Vec4i tetVtxIndices{-1};
  Vec4d weights{0.0};
  int tetID = -1;
};

PointInTet getPointInTet(const BarycentricCoordinates & bc, int pointID)
{
  PointInTet pit;
  assert(pointID >= 0 && pointID < (int)bc.getNumLocations());
  pit.tetID = bc.getEmbeddingElement(pointID);
  assert(pit.tetID >= 0);
  pit.tetVtxIndices = bc.getEmbeddingVertexIndices(pointID);
  pit.weights = bc.getEmbeddingWeights(pointID);
  return pit;
}

}

void ImmersionMesher::buildCellSurfaceMeshes()
{
  const auto & cutMesh = selfCutMesh.cutMesh;

  // first, prepare buffers
  cellPatchStart.resize(numCells);
  manifoldCellMesh.resize(numCells);
  manifoldCellMeshOriVtxIDs.resize(numCells);
  manifoldCellMeshPatchOri2NewVtxIDMap.resize(numCells);
  manifoldCellOriTriIDs.resize(numCells);
  manifoldCellOriTris.resize(numCells);
  manifoldCellOri2NewTriIDMap.resize(numCells);

  for(int cellID = 1; cellID < numCells; cellID++)
  {
    vector<Vec3i> cellTris; // records the triangles of the cell. the normal of each cell triangle points outward (relative to the cell)
                            // cell triangle vtx are in cutVtxIDs
    vector<int> patchIDs; // B-patches of the cell
    map<int, int> patchTriIDStart; // patchID -> cell triID start in cellTris
    map<int, int> patchTriIDEnd; // patchID -> cell triID end in cellTris
    map<int, pair<int, bool>> cutTriID2CellTriIDMap; // cutTriID -> <cellTriID, orientation>
    vector<int> cellTri2cutTriIDs; // cell tri ID -> cutTriID

    cellPatchStart[cellID].clear();
    for(auto p : cellPatches[cellID])
    {
      int patchID = p.first;
      patchIDs.push_back(patchID);
      patchTriIDStart[patchID] = (cellTris.size());
      cellPatchStart[cellID][patchID] = cellTris.size();
      int IDStart = cellTris.size();
      if (p.second) // triangle orientation agrees with the cell, a.k.a, triangle normals point outward
      {
        cellTris.insert(cellTris.end(), patchTris[patchID].begin(), patchTris[patchID].end());
      }
      else
      {
        for(Vec3i t : patchTris[patchID])
          cellTris.emplace_back(t[0], t[2], t[1]); // reverse orientation
      }
      for(size_t i = 0; i < patchTris[patchID].size(); i++)
      {
        int cutTriID = patchTriIDs[patchID][i];
        cutTriID2CellTriIDMap[cutTriID] = make_pair(IDStart+i, p.second);
        cellTri2cutTriIDs.push_back(cutTriID);
      }
      patchTriIDEnd[patchID] = (cellTris.size());
    }

    // to create a manifold surface mesh for the cell, we should merge vertices on arcs shared by B-patches

    // we properly connect triangles by first separate all triangles, assuming each triangle has three unique
    // vtx, with IDs to be called unqID (unique ID), for cellTriID, its vtx have IDs: 3*cellTriID+1, ... , 3*cellTriID+3
    // each triangle has three unqID that are not shared by any other triangles
    // later we will merge them according to the triangle connection

    // we use DisjointSet to merge unqIDs
    DisjointSet dset(cellTris.size() * 3);

    // we first merge unqIDs if they are from the same cutVtxID on the same patch
    for(size_t lpID = 0; lpID < patchIDs.size(); lpID++) // local patch ID
    {
      int patchID = patchIDs[lpID];
      // first, let's reconnect triangles that are in the same patch
      map<int, set<int>> oldVtxID2unqID;
      for(int cellTriID = patchTriIDStart[patchID]; cellTriID < patchTriIDEnd[patchID]; cellTriID++)
      {
        // this cellTriID is the ID of triangles in cellTris
        for(int i = 0; i < 3; i++)
        {
          int unqID = cellTriID * 3 + i;
          int cutVtxID = cellTris[cellTriID][i]; // old cut vtx ID
          oldVtxID2unqID[cutVtxID].insert(unqID);
        }
      }
      // merge unqIDs that belong to the same oldVtxID in this patch
      for(const auto & p : oldVtxID2unqID) { dset.unionRange(p.second); }
    }

    // now let's merge unqiue IDs based on B-patch geometric neighbor information,
    // which means, if two B-patches are geometric neighbors along a particular arc on the cell,
    // then the vertex pairs from both patches on the arc should be merged
    // cellPatchBouNbrs: cellID -> patchID -> bouID -> <nbr bouID, nbr patchID >
    for(const auto & p : cellPatchBouNbrs[cellID])
    {
      int patchID = p.first;
      for(const auto & p2 : p.second) // p.second: bouID -> <nbr bouID, nbr patchID >
      {
        int bouID = p2.first;
        int nbrBouID = p2.second.first;
        int nbrPatchID = p2.second.second;
        // cout << cellID << " " << patchID << " " << bouID << " " << nbrBouID << " " << nbrPatchID << endl;
        assert(mapFind(patchBouNbrs[nbrPatchID],nbrBouID));
        assert(mapFind(patchBouNbrs[patchID],bouID));
        if (patchID > nbrPatchID) continue; // visit each patchID pair only once

        // we first build a map from cutVtxIDs on the bouID from patchID to unqID
        map<int, int> bouVtxID2unqID; // cutVtxOnBouID -> unqID
        for(int cellTriID = patchTriIDStart[patchID]; cellTriID < patchTriIDEnd[patchID]; cellTriID++)
        {
          for(int i = 0; i < 3; i++)
          {
            int unqID = cellTriID * 3 + i;
            int cutVtxID = cellTris[cellTriID][i]; // old vtx ID
            if (bouVertices[bouID].find(cutVtxID) == bouVertices[bouID].end()) continue;
            // in this loop, bouVtxID2unqID[cutVtxID] may be assigned more than once
            // which means, one cutVtxID on the bou corresponds to more than one unqID
            // but this is fine, because in the previous code we already merge those unqIDs on the same patch
            // so those unqIDs are already in the same set in DisjointSet
            bouVtxID2unqID[cutVtxID] = unqID;
          }
        }
        assert(bouVtxID2unqID.size() == bouVertices[bouID].size());

        // then we build a map from cutVtxIDs on the nbrBouID from nbrPatchID to unqID
        map<int, int> nbrBouVtxID2unqID; // cutVtxOnNbrBouID -> unqID
        for(int cellTriID = patchTriIDStart[nbrPatchID]; cellTriID < patchTriIDEnd[nbrPatchID]; cellTriID++)
        {
          for(int i = 0; i < 3; i++)
          {
            int unqID = cellTriID * 3 + i;
            int cutVtxID = cellTris[cellTriID][i]; // old vtx ID
            if (bouVertices[nbrBouID].find(cutVtxID) == bouVertices[nbrBouID].end()) continue;
            nbrBouVtxID2unqID[cutVtxID] = unqID; // we only need to record one unqID since they are already in the same disjoint set
          }
        }
        assert(nbrBouVtxID2unqID.size() == bouVertices[nbrBouID].size());

        // finally, we use unqID to merge cutVtxIDs on bouID and nbrBouID
        UEdgeKey bouEdge(bouID, nbrBouID);
        for(auto p : cellBouCorrespondence[cellID][bouEdge]) // for each vtx pair on the bou-pair
        {
          int bouVtx = p.first, nbrBouVtx = p.second;
          if (bouEdge[0] == nbrBouID)
          {
            swap(bouVtx, nbrBouVtx); // make sure bouVtx is from bouID and nbrBouVtx is from nbrBouID
          }
          assert(setFind(bouVertices[bouID], bouVtx));
          assert(setFind(bouVertices[nbrBouID], nbrBouVtx));
          assert(mapFind(bouVtxID2unqID,bouVtx));
          assert(mapFind(nbrBouVtxID2unqID,nbrBouVtx));
          int unqID = bouVtxID2unqID[bouVtx];
          int nbrUnqID = nbrBouVtxID2unqID[nbrBouVtx];
          assert(cellTris[unqID/3][unqID%3] == bouVtx);
          assert(cellTris[nbrUnqID/3][nbrUnqID%3] == nbrBouVtx);
          dset.unionSet(unqID, nbrUnqID); // merge unqID and nbrUnqID belonging to the same cutVtxID on the arc
        }
      } // end each neighboring patch
    } // end each patch on the cell

    // now we have finished merging unqIDs
    auto unqID2mergedID = dset.createOldToNewIDMapping(); // unqID -> mergedID
    // the mergedID will become the ID for the final manifold cell mesh vertices

    // now we build newPos, mergedID2cutVtxIDs and manifoldCellTris
    vector<set<int>> mergedID2cutVtxIDs; // mergedID -> original IDs in cutMesh
    vector<Vec3d> mergedPos;
    vector<Vec3i> mergedCellTris(cellTris.size());
    for(size_t cellTriID = 0; cellTriID < cellTris.size(); cellTriID++)
    {
      for(int i = 0; i < 3; i++)
      {
        int unqID = cellTriID * 3 + i;
        int mergedID = unqID2mergedID[unqID];

        if (mergedID >= (int)mergedPos.size())
        {
          mergedPos.resize(mergedID + 1);
          mergedID2cutVtxIDs.resize(mergedID + 1);
        }
        mergedPos[mergedID] = cutMesh.pos(cellTris[cellTriID][i]);
        mergedID2cutVtxIDs[mergedID].insert(cellTris[cellTriID][i]);
        mergedCellTris[cellTriID][i] = mergedID;
      }
    }

    // next, build a mapping: patchCutVtx2MergedVtxID: patchID -> cutVtxID -> mergedVtxID
    map<int, map<int,int>> patchCutVtx2MergedVtxID;
    for(size_t lpID = 0; lpID < cellPatches[cellID].size(); lpID++) // for each local patch ID
    {
      int patchID = patchIDs[lpID];
      // first, let's reconnect triangles that are in the same patch
      for(int cellTriID = patchTriIDStart[patchID]; cellTriID < patchTriIDEnd[patchID]; cellTriID++)
      {
        for(int i = 0; i < 3; i++)
        {
          int mergedID = mergedCellTris[cellTriID][i];
          int cutVtxID = cellTris[cellTriID][i];
          patchCutVtx2MergedVtxID[patchID][cutVtxID] = mergedID;
        }
      }
    }

    TriMeshGeo cellMesh(move(mergedPos), move(mergedCellTris));

    if (verbose)
      cellMesh.save("cell" + to_string(cellID) + ".obj");

    if (areTrianglesManifold(cellMesh.triangles()) == false) // if the cell mesh is not manifold
    {
      // Ideally, if the input is not degenerate, all cell meshes should be manifold.
      unordered_map<OEdgeKey, int> oedgeTri; 
      bool isEdgeManifold = getOEdgeTriMap(cellMesh.triangles(), oedgeTri); // then check if it is edge-manifold
      if (isEdgeManifold == false)
      {
        cout << "Error: created cell mesh for cellID " << cellID << " is not edge-manifold" << endl;
        throw 1;
      }

      // if it is edge-manifold, then we can still fix this!
      map<int, int> fixedNewVtxID2mergedVtxID;
      fixNonManifoldVerticesOnEdgeManifoldTriangles(cellMesh, oedgeTri, &fixedNewVtxID2mergedVtxID);
      // We fix this mesh by duplicating the non-manifold vtx to make entire mesh manifold.
      // Now we have more vertices, we should modify mergedID2cutVtxIDs as well.
      mergedID2cutVtxIDs.resize(cellMesh.numVertices());
      for(auto p : fixedNewVtxID2mergedVtxID)
      {
        int newVtxID = p.first, oriMergedVtxID = p.second;
        mergedID2cutVtxIDs[newVtxID] = mergedID2cutVtxIDs[oriMergedVtxID];
      }
      // cellMesh.save("cell" + to_string(cellID) + ".obj");
      // exit(1);
    }

    manifoldCellMesh[cellID] = (move(cellMesh));
    manifoldCellMeshOriVtxIDs[cellID] = (move(mergedID2cutVtxIDs));
    manifoldCellMeshPatchOri2NewVtxIDMap[cellID] = (move(patchCutVtx2MergedVtxID));
    manifoldCellOriTriIDs[cellID] = move(cellTri2cutTriIDs);
    manifoldCellOriTris[cellID] = move(cellTris);
    manifoldCellOri2NewTriIDMap[cellID] = move(cutTriID2CellTriIDMap);
  }
}

// build a tet mesh for each cell, using virtual tets algorithm
void ImmersionMesher::buildCellTetMeshes()
{
  const auto & cutMesh = selfCutMesh.cutMesh;

  ProfilerSection cellSurfaceMeshSec(&profiler, "buildCellSurfaceMesh");
  ProfilerSection entireCuttingSec(&profiler, "entireCutting");
  ProfilerSection searchTetsSec(&profiler, "searchTetsInEachCell");
  ProfilerSection virtualTetsSec(&profiler, "prepareVirtualTets");
  ProfilerExtraSection tetTriInterExtraSec(&profiler, "tetTriIntersecting");

  // first, let's build a manifold surface mesh for each cell, which will be used
  // as input to the virtual tets algorithm
  cout << "Building manifold cell surface mesh..." << endl;

  cellSurfaceMeshSec.start();
  buildCellSurfaceMeshes();
  cellSurfaceMeshSec.stop();

  cellTetMeshInterps.resize(numCells);
  cellTetMeshTri2TetIDs.resize(numCells);
  cellTetMeshes.resize(numCells);
  cellTetMeshOriTetIDs.resize(numCells);

  // build tet vs tri cutting
  // cout << "cutting on the entire mesh" << endl;
  assert(sizei(selfCutMesh.cutPosExact) == cutMesh.numVertices());
  entireCuttingSec.start();
  TetTriMeshCutting entireCutting(inputTetMeshGeo, cutMesh);
  cout << "Computing tet-triangle intersection..." << endl;

  // compute the mapping: tetID -> triangleIDs intersect with the tet
  tetTriInterExtraSec.start();
  entireCutting.computeIntersectingTriTets();
  tetTriInterExtraSec.stop();

  // note: that in other places of the code, we use "cut trinagle or cutTri(ID)" to refer to the triangles in the cutMesh
  // from the libigl::selfIntersection code, where cells and patches are created
  // but here in this function we also do some cutting (tet cut tri), so to avoid confusion,
  // from now on in this function we still use cut triangle to refer to the original cut tri,
  // and we use ecTri (entire cutting tri) to refer to the triangles after the cut from tet tets, and
  // ecVtx to refer to the vertices after the tet cut

  // cutTriID -> the ecVtxIDs inside cutTriID due to entireCutting
  vector<vector<int>> cutTri2innerECVtxIDs(cutMesh.numTriangles());
  // UEdgeKey in cutVtxIDs -> tet ecVtxIDs on the interior of this edge
  map<UEdgeKey, vector<int>> cutEdge2ECVtxIDs;
  vector<Vec3ER> tetPosER;
  if (useCSGInVT)  // The CSG method does not need cut, but it needs exact tet positions
  {
    for(int vtxID = 0; vtxID < inputTetMeshGeo.numVertices(); vtxID++)
    {
      Vec3d p = inputTetMeshGeo.pos(vtxID);
      Vec3ER pER(p[0], p[1], p[2]);
      tetPosER.emplace_back(move(pER));
    }
  }
  else // we perform Sutherland-Hodgman tet-triangle clipping algorithm
  {
    cout << "Performing Sutherland-Hodgman tet-triangle clipping algorithm..." << endl;
    entireCutting.computeCutTriangles(&selfCutMesh.cutPosExact);
    cout << "Finished clipping." << endl;

    // we build cutTri2innerECVtxIDs and cutEdge2ECVtxIDs

    // for each new vtx generated by the tet cut:
    for(size_t ecVtxID = cutMesh.numVertices(); ecVtxID < entireCutting.cutTriPositions().size(); ecVtxID++) // ecVtx: entire cutting vtx
    {
      const auto & f = entireCutting.cutVertexFeatures()[ecVtxID]; // get the feature the vtx lies on

      if (f.isInsideTri()) // ecVtxID is on the interior of a cut triangle
      {
        assert(mapFind(cutMeshTriKeyIDs, f.triFeature));
        cutTri2innerECVtxIDs[cutMeshTriKeyIDs[f.triFeature]].push_back(ecVtxID);
      }
      else if (f.isTriEdge())
      {
        cutEdge2ECVtxIDs[f.getTriEdge()].push_back(ecVtxID);
      }
    }
  }
  entireCuttingSec.stop();

  // for each cell, let's build the tet mesh serving as the input to the virtual tets algorithm
  // we call this tet mesh, bvTetMesh (before virtual tet mesh)
  searchTetsSec.start();

  // remove the outer tets (tets embedding empty space outside the cutTriMesh)
  function<bool(int tetID)> isTetBou1 = [&](int tetID) ->bool { return entireCutting.getTrianglesIntersectingTet(tetID).size() > 0; };
  function<bool(int tetID)> isTetBou2 = [&](int tetID) ->bool { return entireCutting.getCutTriGroup(tetID).tri.size() > 0; };

  WindingNumberTree wnTree; // we use winding number to check inside/outside of the cutMesh
  wnTree.build(cutMesh);
  auto isTetOuter = [&](int tetID)
  {
    double wn = wnTree.windingNumber(cutMesh, inputTetMeshGeo.ref().computeTetCenter(tetID));
    return wn < 0.5;
  };

  TetNeighbor tetNeighbor(inputTetMeshGeo.tets());
  // inputTetOutside: inputTetID -> whether tet is outside of cutMesh
  auto inputTetOutside = labelOuterTets(inputTetMeshGeo, tetNeighbor, (useCSGInVT ? isTetBou1 : isTetBou2), isTetOuter);

  // prepare buffers
  vector<TetMeshGeo> cellBVTetMesh(numCells); // cellID -> bvTetMesh for the virtual tets algorithm on this cell
  vector<vector<int>> cellBVInputTetIDs(numCells);  // input tet IDs for each cell's bvTetMesh: cellID -> bvTetID -> inputTetIDs
  vector<vector<int>> cellBVTetVtxID2InputTetVtxID(numCells);
  vector<map<int,int>> cellInputTetVtxID2BVTetVtxID(numCells);

  // now, let's find all the tets that belong to each bvTetMesh

  int numTets = inputTetMeshGeo.numTets();

  vector<vector<int>> cutTri2cutTets(cutMesh.numTriangles()); // cutTriID -> tetIDs that intersect the triangle
  vector<set<int>> tetCells(numTets); // tetID -> cells need this tetID for BVTetMesh

  const int TET_OUTER = -100; // tet outside cutMesh
  const int TET_CUR = INT_MAX;
  const int TET_UNINIT = -1;
  const int TET_INTER_TRI = 0; // tet intersecting cutTris
  vector<int> tetLabel(numTets, TET_UNINIT);

  // first, build cutTri2cutTets and set TET_OUTER and TET_INTET_TRI to tetLabel
  for(int tetID = 0; tetID < numTets; tetID++)
  {
    vector<int> cutTriIDs;
    if (useCSGInVT)
    {
      cutTriIDs = entireCutting.getTrianglesIntersectingTet(tetID);
    }
    else
    {
      cutTriIDs = entireCutting.getCutTriGroup(tetID).oriID;
      sortAndDeduplicate(cutTriIDs);
    }

    for(int cutTriID : cutTriIDs) { cutTri2cutTets[cutTriID].push_back(tetID); }
    if (inputTetOutside[tetID]) { tetLabel[tetID] = TET_OUTER; }
    else if (cutTriIDs.size() > 0) { tetLabel[tetID] = TET_INTER_TRI; }
  }

  // then, put tets intersecting the surface cell meshes into cellBVInputTetIDs, and tetCells
  for(int cellID = 1; cellID < numCells; cellID++)
  {
    for(int cutTriID : manifoldCellOriTriIDs[cellID])
    {
      for(int tetID : cutTri2cutTets[cutTriID])
      {
        cellBVInputTetIDs[cellID].push_back(tetID);
        tetCells[tetID].insert(cellID);
      }
    }
  }

  vector<WindingNumberTree> cellWinTrees(numCells); // build winding number tree for each cell mesh
  vector<bool> cellWinTreesBuilt(numCells, false);

  // we then go to every tet not labeld yet (which are not outer tets and does not intersect cutTris),
  // determing which cell they belong to
  // for those tets, each of them belongs to one and only one cell
  for(int tetID = 0; tetID < numTets; tetID++)
  {
    if (tetLabel[tetID] != TET_UNINIT) continue; // skip if has a label already
    // we don't query one tet to all cell meshes every time, which is too costly
    // instead, we use a flooding method, a basic implementaion is: once a seed tet is checked against all cell meshes
    // and found which cell it belongs to, we flood-fill its neighboring non-labeled tets to belong to that cell
    // then, we have a better method to reduce inside/outside queries: we don't need to check a seed tet against
    // all cell meshes, we only check those cell meshes whose triangles intersect a tet neighboring the group of
    // tets found by this flood-fill
    int seed = tetID;
    vector<int> ffCand = {seed}; // ffCand: flood-fill candidates
    tetLabel[seed] = TET_CUR;
    int ffBegin = 0, ffEnd = 1;
    int targetCellID = -1; // which cell this seed belongs to
    set<int> candCellIDs;
    while(ffBegin < ffEnd)
    {
      for(int i = ffBegin; i < ffEnd; i++)
      {
        int ffTetID = ffCand[i];
        for(int ffnbrTetID: tetNeighbor.getTetNeighbors(ffTetID))
        {
          if (ffnbrTetID < 0) continue;

          if(targetCellID < 0 && tetLabel[ffnbrTetID] == TET_INTER_TRI)
          {
            assert(tetCells[ffnbrTetID].size() > 0);
            const set<int> & rest = tetCells[ffnbrTetID]; // for each cell whose triangles intersect this nbrTet
            if (rest.size() == 1) // this can only be the cell this seed tet belongs to
            {
              // we found this cell
              int cellID = *rest.begin();
              targetCellID = cellID;
            }
            else
            {
              candCellIDs.insert(rest.begin(), rest.end());
            }
          }

          if (tetLabel[ffnbrTetID] != TET_UNINIT) { continue; }
          tetLabel[ffnbrTetID] = TET_CUR;
          ffCand.push_back(ffnbrTetID);
        }
      }
      ffBegin = ffEnd;
      ffEnd = ffCand.size();
    }
    // end of floodFill
    // now ffCand stores all the tets found by this flood-fill

    if(targetCellID < 0) // let's find out which cell this seed tet belongs to
    {
      assert(candCellIDs.size() > 0);
      for(int cellID : candCellIDs)
      {
        if (cellWinTreesBuilt[cellID] == false)
        {
          cellWinTreesBuilt[cellID] = true;
          cellWinTrees[cellID].build(manifoldCellMesh[cellID]);
        }
        double wn = cellWinTrees[cellID].windingNumber(manifoldCellMesh[cellID], inputTetMeshGeo.ref().computeTetCenter(seed));
        if (wn > 0.5) // seed is inside, found the cell!
        {
          targetCellID = cellID;
          break;
        }
      }
      if (targetCellID < 0)
      {
        cout << "error on tet " << seed << ", vtx: " << inputTetMeshGeo.tet(seed) << endl;
        cout << "cand cells are " << streamRange(candCellIDs) << endl;
        for(int cellID : candCellIDs)
        {
          double wn = cellWinTrees[cellID].windingNumber(manifoldCellMesh[cellID], inputTetMeshGeo.ref().computeTetCenter(seed));
          cout << "cell " << cellID << " wn is " << wn << endl;
        }
      }
    }
    assert(targetCellID >= 0);

    for(int fftetID : ffCand)
    {
      // since cellID == 0 is the outer cell (empty space) that we will never used to run virtual tets algorithm
      // we don't need to offset targetCellID to assign to tetLabel
      tetLabel[fftetID] = targetCellID;
      cellBVInputTetIDs[targetCellID].push_back(fftetID);
    }
  } // end for tetID, end of building tetLabel and cellBVInputTetIDs

  // create cellBVTetMesh
  for(int cellID = 1; cellID < numCells; cellID++)
  {
    sortAndDeduplicate(cellBVInputTetIDs[cellID]);

    cellBVTetMesh[cellID] = getSubTetMesh(inputTetMeshGeo, cellBVInputTetIDs[cellID], &cellBVTetVtxID2InputTetVtxID[cellID],
        &cellInputTetVtxID2BVTetVtxID[cellID]);

    if (verbose)
    {
      cout << "tet at cell " << cellID << " #tet: " << cellBVInputTetIDs[cellID].size() << endl;
      string filename = "cellBVTetMesh" + to_string(cellID) + ".veg";
      cellBVTetMesh[cellID].save(filename);
    }
  }

  searchTetsSec.stop();

  // now let's build the tet-cut-tri data for virtual tets algorithm
  for(int cellID = 1; cellID < numCells; cellID++)
  {
    virtualTetsSec.start();
    const TriMeshGeo & cellMesh = manifoldCellMesh[cellID];
    const vector<set<int>> & cellMeshVtxID2CutVtxIDs = manifoldCellMeshOriVtxIDs[cellID];
    const auto & cellTri2cutTriIDs = manifoldCellOriTriIDs[cellID];
    const auto & cellCutTris = manifoldCellOriTris[cellID];
    const auto & cutTri2CellTriIDMap = manifoldCellOri2NewTriIDMap[cellID];
    cout << "Using virtual tets method to create tet mesh for cell mesh ID " << cellID << "." << endl;
    vector<Vec3ER> triPosER; // build exact positions for the manifold cell mesh
    for(int i = 0; i < cellMesh.numVertices(); i++)
    {
      triPosER.emplace_back(selfCutMesh.cutPosExact[*cellMeshVtxID2CutVtxIDs[i].begin()]);
    }

    // the vars below are the output of virtual tets algorithm
    TetMeshGeo cellTetMesh;                               // output virtualized cell tet mesh
    BarycentricCoordinates bc;                            // barycentric coord for embedding cellMesh into cellTetMesh
    vector<vector<int>> tetTris;                          // virtualized cellTetMesh tetID -> triangleIDs that this tet embeds
    vector<vector<int>> triTets(cellMesh.numTriangles()); // cellTriID -> virtualized tetIDs that intersect this tri
    vector<int> tetMeshOriTetID;                          // cellTetID -> bvTetID

    if (useCSGInVT)
    {
      TetTriIntersectingData data; // build this data:

      for(int bvTetVtxID = 0; bvTetVtxID < cellBVTetMesh[cellID].numVertices(); bvTetVtxID++)
        data.tetPosER.push_back(tetPosER[cellBVTetVtxID2InputTetVtxID[cellID][bvTetVtxID]]);

      // for each tet in the BVTetMesh for this cell
      for(int bvTetID = 0; bvTetID < cellBVTetMesh[cellID].numTets(); bvTetID++)
      {
        vector<int> group; // store cellTriIDs intersecting this tet
        // we have entireCutting, which stores the intersection between inputTetIDs and cutTriIDs
        // we need to transfer between bvTetIDs and inputTetIDs, and between cutTriIDs and cellTriIDs
        int inputTetID = cellBVInputTetIDs[cellID][bvTetID];
        const auto & cutTriGroup = entireCutting.getTrianglesIntersectingTet(inputTetID);
        for(size_t j = 0; j < cutTriGroup.size(); j++)
        {
          int cutTriID = cutTriGroup[j];
          if (mapNotFind(cutTri2CellTriIDMap, cutTriID)) continue;
          auto it = cutTri2CellTriIDMap.find(cutTriID);
          int cellTriID = it->second.first;
          group.push_back(cellTriID);
        }
        data.triInTet.emplace_back(move(group));
      }

      virtualTetsSec.stop();
      cellTetMesh = createVirtualTetsMeshViaCSG(cellBVTetMesh[cellID], cellMesh, triPosER, data, &bc, nullptr, &tetMeshOriTetID, &tetTris,
          verbose, &profiler);
    }
    else
    {
      TetTriCuttingData cellCutting;

      // The code below in this context has very complicated logic.
      // So we will go through the meshes and indices we used in this algorithm to make sure people don't get confused
      // when reading the following lengthy code:
      //
      // the input to the WHOLE class is:
      //   inputTriMesh and inputTetMesh, they have indices defined as: inputTriVtxID, inputTriID, inputTetVtxID, inputTetID
      // after we call libigl::remeshSelfIntersection, inputTriMesh is cut by itself to become: cutMesh, which has indices: cutVtxID, cutTriID
      // note that all the inputTriVtxIDs are inherited by cutMesh, so the first #inputTriVtxID cutVtx are exactly the inputTriVtx with same IDs
      // the new cutVtx created are on the self-intersection cuts. Since two triangles intersect to form a cut,
      // we can create two groups of cutVtx on the two triangles, respectively, to represent this cut.
      // But since they intersect, the two groups also co-position in space.
      // we have another index: st(Vtx)ID, stitched vtx (named from libigl interface) to represent the unique spatial location
      // of cutVtx. Therefore, an "arc" is always expressed in stVtxIDs, where the "bou"s are represented in cutVtxIDs
      // so that we know which two bous lie on the same arc.
      //
      // We also created (manifold) cell meshes: cellMeshes, whose IDs are: cellVtxID and cellTriID
      // cellVtx are subset of cutVtx and cellTris are subset of cutTris and we have mapping between the IDs
      // note that more than one cutVtx can map to the same cellVtx, because cellMesh is created by merging bous on its B-patches
      // so those cutVtx on the shared bous may get merged and mapped to the same cellVtx
      // Things get complicated on a self-touching cell, where the cell "touches" itself at a particular arc (or more arcs).
      // At that arc, originally there are two bous which means two groups of cutVtxIDs. After building the manifold cell mesh,
      // there are also two groups of cellVtxIDs on that arc, but they are not the copy of the old groups.
      // An illustration is as follows:
      //            b3                             c3
      //            |                               |
      //  b0------b1-------b2            c0------c1 c4-------c2
      //            b4                           |
      //            |                            |
      //            b5                           c5
      // on the self-touching arc on the cutMesh (left), b1 and b4 are on the same location, mapping to the same stVtxID.
      // on the cell mesh (right), we have two copies of the vtx with same stVtxIDs, c1 and c4, but they connect to
      // neighboring vtx in a different pattern than on the cutMesh.
      // luckly the cellTri mapping is simpler, each cellTri maps to one and only one cutTri
      // however, note that some cellTris have opposite orientations as the corresponding cellTris, because we want cellMesh to
      // be consistently oriented. That's why we also need to store a bool variable in manifoldCellOri2NewTriIDMap.

      // then, we do a tet vs tri cut between inputTetMesh and cutMesh, the result is: ecMesh (entire(Tet)CuttingMesh)
      // the IDs of ecMesh are: ecVtxID and ecTriID
      // note that all the cutVtxIDs are inherited by ecMesh, so the first #cutVtxID ecVtx are exactly the cutVtx with same IDs

      // finally, we will create a ccMesh (cell(Tet)CuttingMesh), the IDs of which are: ccVtxID and ccTriID
      // and a bvTetMesh (beforeVirtual(Tets)TetMesh), whose IDs are bvVtxID and bvTetID
      // bvTetMesh is just a subset of inputTetMesh, serving as the input for the virtual tet algorithm
      // ccMesh is the result of the tet vs tri cut between bvTetMesh and cellMesh
      // note that all the cellVtxIDs are inherited by ccMesh, so the first #cellVtxID ccVtx are exactly the cellVtx with same IDs

      // well, that's the thing...
      // I know an easier implementation is just to run the tet vs tri cutting algorithm for each (cellMesh, bvTetMesh) pair,
      // but it requires more exact arithmetic operations because a patch shared by two cells will be cut twice
      // to avoid this redundancy, we come at this lengthy and complicated algorithm, which firt do tet vs tri cut on inputTetMesh and cutMesh,
      // then use this cuttind data to construct cellCutting, the data for cellMesh vs bvTetMesh.

      // OK, let's look into what cellCutting is made of:
      // vector<Vec3d> cutVtxPos;              // ccVtx pos
      // vector<Vec3ER> cutVtxPosER;           // exact ccVtx pos
      // vector<TetTriCutFeature> features;    // cut feature for each ccVtx, made by triFeature in cellTriIDs and tetFeature in bvTetIDs
      // vector<Vec3ER> tetPosER;              // exact bvTetVtx pos
      // vector<CutTriGroup> cutTriGroups;     // for each bvTet, stores the ccTris indexed in ccVtxIDs and the cellTriIDs those ccTris are from

      // the following code will build the above data

      // add cell mesh vtx positions and features to cellCutting
      // we first add those vertices of the cell surface mesh
      for(int i = 0; i < cellMesh.numVertices(); i++)
      {
        int cutVtxID = *cellMeshVtxID2CutVtxIDs[i].begin();
        cellCutting.cutVtxPos.push_back(cutMesh.pos(cutVtxID));
        cellCutting.cutVtxPosER.push_back(selfCutMesh.cutPosExact[cutVtxID]);
        const auto & f = entireCutting.cutVertexFeatures()[cutVtxID]; // get feature
        // the feature is for the vtx that belongs to the vtx group of the ecMesh which is inherited from cutMesh
        assert(f.isTriVertex()); // f.triFeature must be (-1,-1, cutVtxID), we should convert it into cellVtxIDs
        UTriKey newTriKey(-1,-1, i);
        // we should also convert from the inputTetIDs used in f.tetFeature into bvTetIDs
        // but we will do this later, so for now let's assign cellCutting.features in this way:
        cellCutting.features.emplace_back(newTriKey, f.tetFeature);
      }

      // next, let's add those vertices generated by cutting of the tets to cellMesh
      // those vertices are found in entireCutting, we should convert their indices before adding them to cellCutting
      map<TetTriCutFeature, int> edgeFeatureMap; // <triKey in cellVtxIDs, tetKey in cutVtxIDs> -> ecVtxID
      // This edgeFeatureMap is used because we are adding those vertices generated by tet cuts on ecMesh to cellCutting
      // and we are searching for those vertices by looping over each cell triangle.
      // Since some cell edges are shared by two cell triangles, to avoid adding those vertices on the cell edge interior
      // multiple times, we have to use this edgeFeatureMap

      vector<map<int,int>> cellTriEC2CC(cellMesh.numTriangles()); // cellTriID -> ecVtxID -> ccVtxID, ccVtxID: cell cutting vtx
      // We also build this cellTriEC2CC for creating a local map: ecVtxID -> ccVtxID on each cellTri
      for(int cellTriID = 0; cellTriID < cellMesh.numTriangles(); cellTriID++) // for each cellMesh triangles
      {
        const auto & cutTri = cellCutTris[cellTriID];
        for(int i = 0; i < 3; i++)
        {
          int cellVtxID = cellMesh.triVtxID(cellTriID, i);
          int cutVtxID = cutTri[i];
          cellTriEC2CC[cellTriID][cutVtxID] = cellVtxID;
        }

        // On ecMesh, some vertices are inherited from cutMesh. The part of them which belongs to this ceMesh is added to entireCutting
        // in the last loop. The other part is the new vertices generated by the tet cuts.
        // In the following loop, we will add those new vertices that are on the interior of a cutTriangle in entireCutting
        int cutTriID = cellTri2cutTriIDs[cellTriID];
        for(int ecVtxID : cutTri2innerECVtxIDs[cutTriID]) // the ecVtxIDs on the interior of cutTriID
        {
          int ccVtxID = cellCutting.cutVtxPos.size(); // create a new ccVtxID to represent this ecVtxID
          cellCutting.cutVtxPos.push_back(entireCutting.cutTriPositions()[ecVtxID]);
          cellCutting.cutVtxPosER.push_back(entireCutting.cutTriPositionsER()[ecVtxID]);
          const auto & f = entireCutting.cutVertexFeatures()[ecVtxID]; // get feature
          // we will convert the triFeature in f from cutVtxID to cellVtxID
          // again we leave the conversion of tetFeature for now
          cellCutting.features.emplace_back(cellMesh.tri(cellTriID), f.tetFeature);
          assert(mapNotFind(cellTriEC2CC[cellTriID], ecVtxID));
          cellTriEC2CC[cellTriID][ecVtxID] = ccVtxID;
        }

        // In the following loop, we will add those new vertices that are on the edge interior of a cutTriangle in entireCutting
        for(int i = 0; i < 3; i++)
        {
          int cutVtx0 = cutTri[i], cutVtx1 = cutTri[(i+1)%3]; // give the edge in cutVtxIDs
          UEdgeKey cutEdge(cutVtx0, cutVtx1);
          if (mapNotFind(cutEdge2ECVtxIDs, cutEdge)) continue;
          for(int ecVtxID : cutEdge2ECVtxIDs[cutEdge]) // the ecVtxIDs on the edge interior
          {
            const auto & f = entireCutting.cutVertexFeatures()[ecVtxID]; // get feature
            assert(f.isTriEdge());
            UTriKey newTriKey(-1, cellMesh.triVtxID(cellTriID, i), cellMesh.triVtxID(cellTriID, (i+1)%3)); // convert from cutVtxID to cellVtxID
            TetTriCutFeature nf(newTriKey, f.tetFeature);
            // again leave tetFeature for now
            if (mapFind(edgeFeatureMap, nf)) // if the feature has been visited before, then we won't create a new ccVtx for it
            {
              // note that here is the usage of edgeFeatureMap
              // we don't use ecVtx as the key for this mapping because, on the cell mesh,
              // more than one cutVtx will map to one cellVtx (due to merging of bous among cel B-patches),
              // so, more than one ecVtx will map to one ccVtx, therefore, ecVtx shouldn't be the key
              cellTriEC2CC[cellTriID][ecVtxID] = edgeFeatureMap[nf];
              continue;
            }

            // otherwise, this is a new ccVtx
            int ccVtxID = cellCutting.cutVtxPos.size();
            cellCutting.cutVtxPos.push_back(entireCutting.cutTriPositions()[ecVtxID]);
            cellCutting.cutVtxPosER.push_back(entireCutting.cutTriPositionsER()[ecVtxID]);
            cellCutting.features.emplace_back(nf);
            edgeFeatureMap[nf] = ccVtxID;
            assert(mapNotFind(cellTriEC2CC[cellTriID], ecVtxID));
            cellTriEC2CC[cellTriID][ecVtxID] = ccVtxID;
          }
        }
      }

      // now we finished building cutVtxPos, cutVtxPosER and features (minus tetFeature indices)
      // we will convert inputTetIDs to bvTetIDs in this loop
      for(auto & f : cellCutting.features)
      {
        int t[4];
        for(int i = 0; i < 4; i++)
        {
          if (f.tetFeature[i] < 0) t[i] = -1;
          else t[i] = cellInputTetVtxID2BVTetVtxID[cellID][f.tetFeature[i]];
        }
        f = TetTriCutFeature(f.triFeature, UTetKey(t));
      }

      // next, build cellCutting.cutTriGroups
      // we will convert CutTriGrop.oriID from cutTriID to cellTriID, and
      // convert CutTriGroup.tri from ecVtxID to ccVtxID
      for(int bvTetID = 0; bvTetID < cellBVTetMesh[cellID].numTets(); bvTetID++)
      {
        int inputTetID = cellBVInputTetIDs[cellID][bvTetID];
        TetTriMeshCutting::CutTriGroup group; // store the group in cellCutting
        const auto & oldGroup = entireCutting.getCutTriGroup(inputTetID);
        for(const auto & p : oldGroup.cutTriIDsOnFace)
        {
          const UTriKey & inputTriKey = p.first;
          int t[3];
          for(int i = 0; i < 3; i++)
            t[i] = cellInputTetVtxID2BVTetVtxID[cellID][inputTriKey[i]]; // convert from inputTetVtxID to BVTetVtxID
          // p.second stores the groupTriIDs of the triangles exactly on the tet face (the UTriKey)
          // this groupTriIDs are for indices of oldGroup.tri, so it is the same for group
          group.cutTriIDsOnFace.emplace(UTriKey(t), p.second);
        }
        for(size_t j = 0; j < oldGroup.tri.size(); j++)
        {
          int cutTriID = oldGroup.oriID[j];
          // one input tet may embed cutTris from more than one cell
          // so we skip if this cutTriID is not in cellMesh
          if (mapNotFind(cutTri2CellTriIDMap, cutTriID)) continue;
          auto it = cutTri2CellTriIDMap.find(cutTriID);
          int cellTriID = it->second.first;
          group.oriID.push_back(cellTriID);
          Vec3i ccTri;
          for(int k = 0; k < 3; k++)
          {
            // here we use cellTriEC2CC to find the mapping: ecVtxID -> ccVtxID
            // the reason why we build cellTriEC2CC as a per-cellTri mapping, instead of a cell-level ecVtxID->ccVtxID mapping,
            // is because that on a self-touching cell, along that self-touching arc, the mapping between
            // ecVtxID and ccVtxID is very complicated (see illustration above before building cellCutting).
            // So we use a per-cellTri mapping to be safe.
            int ecVtxID = oldGroup.tri[j][k];
            assert(mapFind(cellTriEC2CC[cellTriID], ecVtxID));
            int ccVtxID = cellTriEC2CC[cellTriID][ecVtxID];
            ccTri[k] = ccVtxID;
          }
          if (it->second.second == false) // this cell triangle has opposite orientation from the cutTri
            swap(ccTri[1], ccTri[2]);
          group.tri.push_back(ccTri);
        }
        cellCutting.cutTriGroups.emplace_back(move(group));
      }

      // finally, we build tetPosER
      for(int bvTetVtxID = 0; bvTetVtxID < cellBVTetMesh[cellID].numVertices(); bvTetVtxID++)
        cellCutting.tetPosER.push_back(entireCutting.tetPositionsER()[cellBVTetVtxID2InputTetVtxID[cellID][bvTetVtxID]]);

      virtualTetsSec.stop();
      // calling virtual tets algorithm
      cellTetMesh = createVirtualTetsMesh(cellBVTetMesh[cellID], cellMesh, cellCutting, &bc, nullptr,
          &tetMeshOriTetID, &tetTris, false, verbose, &profiler);
    }

    assert(tetTris.size() == tetMeshOriTetID.size());
    for(size_t tetID = 0; tetID < tetTris.size(); tetID++)
    {
      for(int triID : tetTris[tetID]) { triTets[triID].push_back(tetID); }
    }

    if(verbose)
    {
      string filename = "cellTetMesh" + to_string(cellID) + ".veg";
      cellTetMesh.save(filename.c_str());
    }

    cellTetMeshInterps[cellID] = (move(bc));
    cellTetMeshTri2TetIDs[cellID] = move(triTets);
    cellTetMeshes[cellID] = move(cellTetMesh);
    for(int & oldID : tetMeshOriTetID)          // here we use celTetID to refer to virtualized tet ID
      oldID = cellBVInputTetIDs[cellID][oldID]; // convert the mapping from  cellTetID -> bvTetID to cellTetID -> inputTetID
    cellTetMeshOriTetIDs[cellID] = (move(tetMeshOriTetID));
  } // end for cellID
}

void ImmersionMesher::produceFinalTetMesh(vector<ImmersionGraphNode> & cellNodes, int graphID,
    TetMeshGeo & outputTetMesh, BarycentricCoordinates & outputInterpWeight,
    std::vector<TriMeshGeo> * allCellMesh, std::vector<BarycentricCoordinates> * allCellInerpWeight)
{
  // first, creates a tet mesh without stitching cell tetMeshes
  // This means, this tet mesh is only the union of the tet meshes from each nodes, which don't connect any one
  // we call this tet mesh: allTetMesh, with ID: allTetVtxID and allTetID
  vector<Vec4i> allTets;            // allTetID -> tets in allTetVtxIDs
  vector<Vec3d> allTetVertices;     // allTetVtxID -> positions
  vector<int> nodeTetStart;         // nodeID -> allTetID start in allTets
  vector<int> allTetID2InputTetID;  // allTetID -> inputTetID

  // in this loop, we build the above vars
  for(size_t nodeID = 0; nodeID < cellNodes.size(); nodeID++)
  {
    const ImmersionGraphNode & node = cellNodes[nodeID];
    int cellID = node.getCellID();
    nodeTetStart.push_back(allTets.size());

    const auto & cellTetMesh = cellTetMeshes[cellID]; // virtualized tet mesh for cellID
    int tetVtxIDStart = allTetVertices.size();

    for(int i = 0; i < cellTetMesh.numVertices(); i++)
    {
      allTetVertices.push_back(cellTetMesh.pos(i));
    }
    assert(cellTetMeshOriTetIDs[cellID].size() == (size_t)cellTetMesh.numTets());
    for(int i = 0; i < cellTetMesh.numTets(); i++)
    {
      Vec4i vec;
      for(int j = 0; j < 4; j++)
      {
        vec[j] = cellTetMesh.tetVtxID(i, j) + tetVtxIDStart;
      }
      allTets.push_back(vec);
      allTetID2InputTetID.push_back(cellTetMeshOriTetIDs[cellID][i]);
    }
  }

  // now let's do merging!
  // we build a disjoint set structure on allTet vertices, and another one on allTets
  DisjointSet allTetVtxDSet(allTetVertices.size()), allTetsDSet(allTets.size());

  // now merge tets that share the patches connecting two copies
  for(int nodeID = 0; nodeID < sizei(cellNodes); nodeID++)
  {
    int tetStart = nodeTetStart[nodeID];
    ImmersionGraphNode & node = cellNodes[nodeID];
    int cellID = node.getCellID();
    for(auto nbr : node.getNbrs())
    {
      int patchID = nbr.first;
      int nbrNodeID = nbr.second;
      if (nbrNodeID < 0) continue;
      if (nodeID > nbrNodeID) continue; // we visit unordered each node pair only once

      ImmersionGraphNode & nbrNode = cellNodes[nbrNodeID];
      int nbrTetStart = nodeTetStart[nbrNodeID];
      int nbrCellID = nbrNode.getCellID();
      assert(mapFind(cellPatches[cellID],patchID));
      assert(mapFind(cellPatches[nbrCellID], patchID));
      //      map<int, int> patchIDVtx2CellTetID;
      for(size_t i = 0; i < patchTris[patchID].size(); i++) // for each triangle on the patchID shared by nodeID and nbrNodeID
      {
        int cellMeshTriID = cellPatchStart[cellID][patchID] + i; // from cutTri to cellTriID
        auto cellTetMeshIDs = cellTetMeshTri2TetIDs[cellID][cellMeshTriID]; // all the cellTetIDs intersecting this triangle
        int nbrCellMeshTriID = cellPatchStart[nbrCellID][patchID] + i;
        auto nbrCellTetMeshIDs = cellTetMeshTri2TetIDs[nbrCellID][nbrCellMeshTriID];

        // we don't assume cellTetMeshIDs and nbrCellTetMeshIDs having the same size because
        // on degenerate tet vs tri cases where one tri only touches at a tet, this tet may be counted only once
        // in cellTetMeshIDs and nbrCellTetMeshIDs
        // assert(cellTetMeshIDs.size() == nbrCellTetMeshIDs.size());
        map<int, int> localIDPair; // inputTetID -> cellTetID
        for(int cellTetID : cellTetMeshIDs)
        {
          int inputTetID = cellTetMeshOriTetIDs[cellID][cellTetID];
          localIDPair[inputTetID] = cellTetID;
        }
        for(int nbrCellTetID : nbrCellTetMeshIDs)
        {
          int nbrInputTetID = cellTetMeshOriTetIDs[nbrCellID][nbrCellTetID];
          if (mapNotFind(localIDPair, nbrInputTetID)) continue;
//          assert(localIDPair.find(oldID) != localIDPair.end());
          int allTetID = localIDPair[nbrInputTetID] + tetStart; // the tetID in allTetMesh
          int nbrAllTetID = nbrCellTetID + nbrTetStart;    // the nbr tetID in allTetMesh

          allTetsDSet.unionSet(allTetID, nbrAllTetID); // merge the two allTetIDs
          for(int k = 0; k < 4; k++)
          {
            int tetVtx = allTets[allTetID][k];
            int nbrTetVtx = allTets[nbrAllTetID][k];
            assert(allTetVertices[tetVtx] == allTetVertices[nbrTetVtx]); // two tet vtx should be on the same pos
            allTetVtxDSet.unionSet(tetVtx, nbrTetVtx); // merge the two allTetVtxIDs
          }
        }
      } // end for each shared patch
    } // end nbr in node
  } // end loop over node

  // now tets are merged, we will create this merged final tet mesh

  auto allTetVtxID2Final = allTetVtxDSet.createOldToNewIDMapping(); // allTetVtxID -> finalTetVtxID
  int numFinalTetVtx = 1 + *max_element(allTetVtxID2Final.begin(), allTetVtxID2Final.end());

  // build final tet vertex positions
  vector<Vec3d> finalTetVertices(numFinalTetVtx);
  for(size_t i = 0; i < allTetVertices.size(); i++)
  {
    int newID = allTetVtxID2Final[i];
    finalTetVertices[newID] = allTetVertices[i];
  }

  map<UTetKey, vector<int>> finalTetVtxIndex2allTetIDs; // finalTetKey -> allTetIDs
  for(size_t i = 0; i < allTets.size(); i++)
  {
    Vec4i ni;
    for(int j = 0; j < 4; j++)
    {
      ni[j] = allTetVtxID2Final[allTets[i][j]];
      assert(ni[j] >= 0);
    }
    UTetKey key(ni); // the tet key in finalTetVtxIDs
    assert(key.isValidTet());
    finalTetVtxIndex2allTetIDs[key].push_back(i);
  }

  // we use disjoint set to merge tets to get the final tets
  for(const auto & p : finalTetVtxIndex2allTetIDs)
  {
    if (p.second.size() > 1)
    {
      allTetsDSet.unionRange(p.second);
    }
  }

  auto allTetID2Final = allTetsDSet.createOldToNewIDMapping();
  int numFinalTets = 1 + *max_element(allTetID2Final.begin(), allTetID2Final.end());

  // build the finalTets
  vector<Vec4i> finalTetVtxIndices(numFinalTets);
  for(size_t i = 0; i < allTets.size(); i++)
  {
    Vec4i ni;
    for(int j = 0; j < 4; j++)
    {
      ni[j] = allTetVtxID2Final[allTets[i][j]];
      assert(ni[j] >= 0);
    }
    int newID = allTetID2Final[i];
    finalTetVtxIndices[newID] = ni;
  }

  // ========================================================================
  // Then we will try to find the correct embedding for inputTriMesh as well
  // ========================================================================
  int numInputTriVtx = inputTriMeshGeo.numVertices();
  vector<bool> patchUsed(numPatches, false); // whether all patch will be visited, used as integrity check
  vector<bool> embeddedVtxVisited(numInputTriVtx, false); // whether all embedding inputTriVtx are visited, used as integrity check
  vector<Vec4i> finalEmbeddingTetVtxID(numInputTriVtx);
  vector<Vec4d> finalEmbeddingWeights(numInputTriVtx);

  for(size_t nodeID = 0; nodeID < cellNodes.size(); nodeID++) // for each node
  {
    ImmersionGraphNode & node = cellNodes[nodeID];
    int cellID = node.getCellID();
    int tetStart = nodeTetStart[nodeID];

    for(auto cellPatchOri : cellPatches[cellID])  // find which patch this node owns
    {
      int patchID = cellPatchOri.first;
      if (node.getPatchOwnership(patchID) != OWNED) continue;
//      if (cellPatchOri.second == false) continue;
//      assert(mapFind(node.nbrs,patchID));
//      if (node.getNbrIDAtPatch(patchID) >= 0) continue;
      assert(patchUsed[patchID] == false);
      patchUsed[patchID] = true;
      // this node embeds this patch
      for(size_t i = 0; i < patchTris[patchID].size(); i++)
      {
        for(int j = 0; j < 3; j++)
        {
          int inputTriVtxID = patchTris[patchID][i][j];
          if (inputTriVtxID >= numInputTriVtx) continue; // we only need to compute embedding for inputTriVtxIDs
          int cellMeshVtxID = manifoldCellMeshPatchOri2NewVtxIDMap[cellID][patchID][inputTriVtxID];

          PointInTet pit = getPointInTet(cellTetMeshInterps[cellID], cellMeshVtxID);
          int allTetID = pit.tetID;
          assert(allTetID >= 0);
          int finalTetID = allTetID2Final[allTetID + tetStart];

          Vec4i finalTetVtxIndex = finalTetVtxIndices[finalTetID];
          finalEmbeddingTetVtxID[inputTriVtxID] = finalTetVtxIndex;
          finalEmbeddingWeights[inputTriVtxID] = pit.weights;

          if (verbose)
          {
            Vec4d w(0.0);
            TetMesh::computeBarycentricWeights(
                finalTetVertices[finalTetVtxIndex[0]],
                finalTetVertices[finalTetVtxIndex[1]],
                finalTetVertices[finalTetVtxIndex[2]],
                finalTetVertices[finalTetVtxIndex[3]],
                inputTriMeshGeo.pos(inputTriVtxID), &w[0]);
            if (len(w - pit.weights) > 1e-6)
            {
              cout << "Error: weights don't match: " << w << pit.weights << endl;
            }
          }
          embeddedVtxVisited[inputTriVtxID] = true;
        }
      }
    }
  }

  assert(allOf(patchUsed.begin(), patchUsed.end(), true));
  assert(allOf(embeddedVtxVisited.begin(), embeddedVtxVisited.end(), true));


  string graphIDStr = to_string(graphID);
  string interpFilename = "finalTet.interp";
  if (graphID >= 0) interpFilename = "finalTet" + graphIDStr + ".interp";

  // produce output data: outputInterpWeight & outputTetMesh
  outputInterpWeight = BarycentricCoordinates(finalEmbeddingTetVtxID, finalEmbeddingWeights);
  outputTetMesh = TetMeshGeo(finalTetVertices, finalTetVtxIndices);

  // ===========================================================
  // build additional allCellMesh for algorithm visualization
  // ===========================================================

  if (allCellMesh)
  {
    allCellMesh->clear();
    if (allCellInerpWeight) allCellInerpWeight->clear();

    for(size_t nodeID = 0; nodeID < cellNodes.size(); nodeID++) // for each node
    {
      vector<Vec4i> allCellMeshEmbTetVtxIDs;
      vector<Vec4d> allCellMeshEmbWeights;

      ImmersionGraphNode & node = cellNodes[nodeID];
      int cellID = node.getCellID();
      // create an objMesh::Group for each node
      allCellMesh->push_back(manifoldCellMesh[cellID]);

      if (allCellInerpWeight)
      {
        int tetStart = nodeTetStart[nodeID];
        for(int cellMeshVtxID = 0; cellMeshVtxID < manifoldCellMesh[cellID].numVertices(); cellMeshVtxID++)
        {
          PointInTet pit = getPointInTet(cellTetMeshInterps[cellID], cellMeshVtxID);
          int cellTetID = pit.tetID;
          assert(cellTetID >= 0);
          int finalTetID = allTetID2Final[cellTetID + tetStart];
          Vec4i finalTetVtxIndex = finalTetVtxIndices[finalTetID];
          allCellMeshEmbTetVtxIDs.push_back(finalTetVtxIndex);
          allCellMeshEmbWeights.push_back(pit.weights);
        }
      }

      if (allCellInerpWeight)
      {
        allCellInerpWeight->emplace_back(move(allCellMeshEmbTetVtxIDs), move(allCellMeshEmbWeights));
      }
    }
  }
}

void ImmersionMesher::generateImmersedTetMesh(vector<TetMeshGeo> & tetMeshes, vector<BarycentricCoordinates> & embeddingWeights,
    vector<std::vector<TriMeshGeo>> * allCellMeshes, vector<std::vector<BarycentricCoordinates>> * allCellMeshWeights)
{
  buildCellTetMeshes();

  cout << "Given " << computedCellNodes.size() << " graph" << (computedCellNodes.size() <= 1 ? "" : "s") << " to produce..." << endl;
  int count = 0;
  // now we go through each immersion and build a tet mesh for it
  for(int i = 0; i < sizei(computedCellNodes); i++)
  {
    bool same = false;
    // we will check whether two immersions are the same and avoid creating duplicate immersions
    for(int j = 0; j < i; j++)
    {
      if (computedCellNodes[i] == computedCellNodes[j])
      {
        cout << "Graph " << i << " and " << j << " are the same" << endl;
        same = true;
        break;
      }
    }

    if (same == false)
    {
      TetMeshGeo tetMesh;
      BarycentricCoordinates weight;
      vector<TriMeshGeo> allCell;
      vector<BarycentricCoordinates> allCellWeight;
      // produce the final tet mesh for the immersion
      produceFinalTetMesh(computedCellNodes[i], count++, tetMesh, weight,
          (allCellMeshes ? &allCell : nullptr), (allCellMeshWeights ? &allCellWeight : nullptr));
      tetMeshes.emplace_back(move(tetMesh));
      embeddingWeights.emplace_back(move(weight));
      if (allCellMeshes)
        allCellMeshes->emplace_back(move(allCell));
      if (allCellMeshWeights)
        allCellMeshWeights->emplace_back(move(allCellWeight));
    }
  }
}

