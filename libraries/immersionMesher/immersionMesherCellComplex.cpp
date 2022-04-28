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

#include "immersionMesher.h"
#include "predicates.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "performanceCounter.h"
#include "geometryQuery.h"
#include <iostream>
#include <fstream>
using namespace std;

ImmersionMesher::ImmersionMesher()
{
}

void ImmersionMesher::run(
           const TriMeshGeo & triMesh,
           const TetMeshGeo & tetMesh,
           std::vector<TetMeshGeo> & tetMeshes,
           std::vector<BarycentricCoordinates> & embeddingWeights,
           std::vector<std::vector<TriMeshGeo>> * allCellMeshes,
           std::vector<std::vector<BarycentricCoordinates>> * allCellMeshWeights)
{
  buildTopologyData(triMesh, tetMesh);
  runImmersionAlgorithm();
  generateImmersedTetMesh(tetMeshes, embeddingWeights, allCellMeshes, allCellMeshWeights);
}

void ImmersionMesher::buildTopologyData(const TriMeshGeo & triMesh, const TetMeshGeo & tetMesh)
{
  bool tmpDumpInfo = verbose;
  bool tmpCSG = useCSGInVT;
  *this = ImmersionMesher(); // clear all internal data
  verbose = tmpDumpInfo;
  useCSGInVT = tmpCSG;

  inputTetMeshGeo = tetMesh;
  inputTriMeshGeo = triMesh;
  cout << "Input mesh for ImmersionMesher:, #v: " << inputTriMeshGeo.numVertices() << " #t: " << inputTriMeshGeo.numTriangles() << "." << endl;
  for(int triID = 0; triID < inputTriMeshGeo.numTriangles(); triID++)
  {
    if (isTriangleDegenerate(inputTriMeshGeo.pos(triID, 0),inputTriMeshGeo.pos(triID, 1),inputTriMeshGeo.pos(triID, 2)))
    {
      cout << "Error: input triMesh has degenerate triangles" << endl;
      throw 1;
    }
  }
  if (areTrianglesManifold(inputTriMeshGeo.triangles()) == false)
  {
    cout << "Error: Input triMesh is not manifold" << endl;
    throw 1;
  }
  {
    set<Vec3d> inputPosSet;
    for(Vec3d p : inputTriMeshGeo.positions()) inputPosSet.insert(p);
    if (sizei(inputPosSet) < inputTriMeshGeo.numVertices())
    {
      cout << "Error: input triMesh has duplicating positions: #ori vtx " << inputTriMeshGeo.numVertices() << " #unique vtx " << sizei(inputPosSet) << endl;
      throw 1;
    }
    set<UTriKey> inputUTriKeys;
    for(Vec3i tri : inputTriMeshGeo.triangles()) inputUTriKeys.emplace(tri);
    if(sizei(inputUTriKeys) < inputTriMeshGeo.numTriangles())
    {
      cout << "Error: input triMesh has duplicating faces: #ori tri: " << inputTriMeshGeo.numTriangles() << " #unique tri: " << inputUTriKeys.size() << endl;
      throw 1;
    }
  }

  ProfilerSection basicDataSec(&profiler, "buildCells&Patches");
  basicDataSec.start();
  buildBasicData();
  basicDataSec.stop();

  ProfilerSection moreDataSec(&profiler, "moreDataPreparation");
  moreDataSec.start();
  prepareDataForSolve();
  moreDataSec.stop();
}

void ImmersionMesher::runImmersionAlgorithm()
{
  PerformanceCounter graphTime;
  computedCellNodes.clear();
  ProfilerSection searchSec(&profiler, "nodeSearch");
  searchSec.start();
  graphTime.StartCounter();
  runNodeSearchMethod(computedCellNodes);
  graphTime.StopCounter();
  cout << "Time to build immersion graph: " << graphTime.GetElapsedTime() << "s." << endl;
  searchSec.stop();
}

void ImmersionMesher::buildBasicData()
{
  ProfilerExtraSection iglRemeshSection(&profiler, "iglRemeshOnInput");
  iglRemeshSection.start();
  cout << "Run igl::remeshSelfIntersection..." << endl;
  selfCutMesh = iglInterface::remeshSelfIntersection(inputTriMeshGeo, false);
  iglRemeshSection.stop();

  assert(sizei(selfCutMesh.cutPosExact) == selfCutMesh.cutMesh.numVertices());

//  fixCutMeshConnection(selfCutMesh);
  const auto & cutMesh = selfCutMesh.cutMesh;
  cout << "Finished remeshing, the result cut mesh: #v " << cutMesh.numVertices() << " #t: " << cutMesh.numTriangles() << "." << endl;

  if (verbose)
  {
    cutMesh.save("cutMesh.obj");
    cutMeshSaved = true;
  }

  const auto & vtxStitchIDs = selfCutMesh.vtxStitchIDs; // cutVtxID -> stitchID
  assert(vtxStitchIDs.size() == (size_t)cutMesh.numVertices());
  for(int i = 0; i < cutMesh.numVertices(); i++) // build stitchPositions: stitchID -> position
  {
    int stID = vtxStitchIDs[i];
    if (stitchPositions.size() <= (size_t)stID) stitchPositions.resize(stID+1);
    stitchPositions[stID] = cutMesh.pos(i);
  }

  if (areTrianglesManifold(cutMesh.triangles()) == false)
  {
    cout << "Error: cut mesh is not manifold!!!!" << endl;
    cutMesh.save("cutMesh.obj");
    cutMeshSaved = true;
    exit(1);
  }

  cutMeshNbr = TriMeshNeighbor(cutMesh);
  auto cutBou = cutMeshNbr.findBoundaryTriangles(cutMesh.triangles());
  if (cutBou.size() > 0)
  {
    cout << "Error: cut mesh has boundary" << endl;
    for(auto p : cutBou) { cout << "(" << p.first << ", " << p.second << ") "; }
    cout << endl;
    exit(1);
  }

  const auto & triPatchIDs = selfCutMesh.triPatchIDs;       // triID -> patchID it belongs to
  const auto & cellIDsAtPatch = selfCutMesh.cellIDsAtPatch; // patchID -> <cellID on front side of patch, cellID on back side of patch>
  numPatches = cellIDsAtPatch.size();

  // build patchTris patchTriIDs cutMeshTriKeyIDs
  // patchTris: patchID -> cut triangles (Vec3i) of the patch
  // patchTriIDs: patchID -> cutTriIDs of the patch
  // cutMeshTriKeyIDs: cut UTriKey -> cutTriID
  for(int triID = 0; triID < cutMesh.numTriangles(); triID++)
  {
    int patchID = triPatchIDs[triID];
    if (sizei(patchTris) <= patchID) { patchTris.resize(patchID+1); patchTriIDs.resize(patchID+1); }
    patchTris[patchID].push_back(cutMesh.tri(triID));
    patchTriIDs[patchID].push_back(triID);
    cutMeshTriKeyIDs[cutMesh.tri(triID)] = triID;
  }

  if (verbose)
  {
    ObjMesh patchMesh;
    patchMesh.addVertexPositions(cutMesh.positions());
    for(int patchID = 0; patchID < numPatches; patchID++)
    {
      ObjMesh::Group g("p" + to_string(patchID));
      for(Vec3i tri : patchTris[patchID])
      {
        ObjMesh::Face f(tri);
        g.addFace(move(f));
      }
      patchMesh.addGroup(move(g));
    }
    patchMesh.save("patchMesh.obj");
  }

  // build cellTris, cellPatches, cellWindingNumbers and cellNeighobrsAtPatch
  for(size_t patchID = 0; patchID < cellIDsAtPatch.size(); patchID++)
  {
    int outCellID = cellIDsAtPatch[patchID].first;
    int inCellID = cellIDsAtPatch[patchID].second;

    int maxCellID = max(inCellID, outCellID);
    if (sizei(cellTris) <= maxCellID)
    {
      cellTris.resize(maxCellID+1);
      cellPatches.resize(maxCellID+1);
      cellWindingNumbers.resize(maxCellID+1, -1);
      cellNeighborsAtPatch.resize(maxCellID+1);
    }
    int triID = patchTriIDs[patchID][0]; // get the first triangle at patchID

    int wnOut = selfCutMesh.windAroundTri[triID].first; // windingNumber outside triID
    int wnIn = selfCutMesh.windAroundTri[triID].second; // windingNumber inside triID
    for(int tID : patchTriIDs[patchID])
    { // check all triangles on this patch has the same inner/outer winding numbers as the first one
      assert(selfCutMesh.windAroundTri[tID] == make_pair(wnOut, wnIn));
    }
    // cellWindingNumber is initialized to be -1
    // if they are not initialized, then their value must match the current value
    if (cellWindingNumbers[outCellID] != -1) assert(cellWindingNumbers[outCellID] == wnOut);
    if (cellWindingNumbers[inCellID] != -1) assert(cellWindingNumbers[inCellID] == wnIn);
    cellWindingNumbers[outCellID] = wnOut;
    cellWindingNumbers[inCellID] = wnIn;

    cellNeighborsAtPatch[inCellID][patchID] = outCellID;
    cellNeighborsAtPatch[outCellID][patchID] = inCellID;

    for(auto t : patchTris[patchID])
    {
      cellTris[inCellID].push_back(t);
      Vec3i tr(t[0], t[2], t[1]);        // get the reversed triangle
      cellTris[outCellID].push_back(tr); // add it to the outer cell
      cellPatches[inCellID][patchID] = true;
      cellPatches[outCellID][patchID] =  false;
    }
  }
  numCells = cellTris.size();
  // libigl selfIntersection code assumes the first cell is always the outer cell (empty space) and its winding number is 0
  assert(cellWindingNumbers[0] == 0);

  cout << "Found #cells = " << numCells << " and #Patches = " << numPatches << "." << endl;
  if (verbose)
  {
    ofstream fout("cellPatches.txt");
    for(int cellID = 0; cellID < numCells; cellID++)
    {
      fout << cellID << ": ";
      for(auto p : cellPatches[cellID]) { fout << "p" << p.first << ", "; }
      fout << endl;
    }
    fout.close();

    fout.open("cellWN.txt");
    for(int cellID = 0; cellID < numCells; cellID++)
      fout << cellID << ": " << cellWindingNumbers[cellID] << endl;
    fout.close();

    for(int i = 0; i < numPatches; i++)
    {
      TriMeshRef mesh(cutMesh.positions(), patchTris[i]);
      auto newMesh = removeIsolatedVertices(mesh);
      newMesh.save("patch" + to_string(i) + ".obj");
    }
  }
}


void ImmersionMesher::prepareDataForSolve()
{
  const auto & cutMesh = selfCutMesh.cutMesh;
  const auto & vtxStitchIDs = selfCutMesh.vtxStitchIDs;

  // build patchMeshRef and patchTriNbr
  // patchMeshRef: patchID -> patch mesh (cutVtxPos, patchTris[patchID])
  // patchTriNbr:  patchID -> TriangleNeighbor for the patch
  vector<TriMeshRef> patchMeshRef(numPatches);
  vector<TriangleNeighbor> patchTriNbr(patchTris.size());
  for(size_t patchID = 0; patchID < patchTris.size(); patchID++)
  {
    // cout << "patchID " << patchID << endl;
    patchMeshRef[patchID] = TriMeshRef(cutMesh.positions(), patchTris[patchID]);
    const auto & patchTri = patchTris[patchID];
    if (areTrianglesEdgeManifold(patchTri) == false)
    {
      cout << "Error: patch " << patchID << " is not edge manifold" << endl;
      throw 1;
    }
    else
    {
      if (areTrianglesManifold(patchTri) == false)
      {
        cout << "Warning: patch " << patchID << " is not vtx manifold" << endl;
      }
    }
//    assert(areTrianglesManifold(patchTri) == true);
    patchTriNbr[patchID] = TriangleNeighbor(patchTri);
  }


  // ===========================================================
  //               Compute data about bous
  // ===========================================================

  {
    map<UEdgeKey, set<pair<int, int>>> uedgesPatch; // uedge -> the <patchID, loopID>s this uedge belongs to
    for(int patchID = 0; patchID < numPatches; patchID++)
    {
      const auto & tris = patchTris[patchID];
      auto loops = patchTriNbr[patchID].findBoundaryLoops(tris); // get the boundary loops for each patch
      for(size_t loopID = 0; loopID < loops.size(); loopID++)    // for each loop
      {
        const auto & loop = loops[loopID];
        for(size_t i = 0; i < loop.size(); i++)
        {
          UEdgeKey ue(loop[i], loop[(i+1)%loop.size()]);
          uedgesPatch[ue].emplace(patchID, loopID);
        }
      }
    }
    // each boundary uedge should be shared by two and only two <patch, loop> pairs
    for(const auto & p : uedgesPatch)
    {
      if (p.second.size() != 2) // this uedge is not shared by exactly two <patch, loop> pairs, sth. wrong
      { // print debugging information
        cout << "Error: one uedge has " << p.second.size() << " patches around" << endl;
        cout << "patchIDs this uedge lies: ";
        for(auto & p2 : p.second) { cout << p2.first << " "; }
        cout << endl;

        UEdgeKey key = p.first; // uedge key
        cout << "UEdge is " << key << endl;
        cout << "Search for all patches to find this key: " << endl;
        for(int pID = 0; pID < numPatches; pID++)
        {
          bool found = false;
          const auto & tris = patchTris[pID];

          for(Vec3i t : tris)
          {
            for(int i = 0; i < 3; i++)
            {
              UEdgeKey e(t[i], t[(i+1)%3]);
              if (e == key)
              {
                cout << "Find this key on patch " << pID << " tri " << t << endl;
                found = true;
              }
            }
          }
          if (found)
          {
            cout << "now lets see the boundary loops of this patch: " << endl;
            auto loops = patchTriNbr[pID].findBoundaryLoops(tris);
            for(size_t loopID = 0; loopID < loops.size(); loopID++)
            {
              const auto & loop = loops[loopID];
              cout << "loop " << loopID << ": ";
              for(auto l : loop) cout << l << " ";
              cout << endl;
            }
          }
        }
        throw 1;
      }
    }

    map<set<pair<int,int>>, set<UEdgeKey>> patch2bous; // <patchID,loopID> pairs -> uedges on the pairs
    for(const auto & p : uedgesPatch)
    {
      patch2bous[p.second].insert(p.first);
    }

    // build patchBouNbrs, bous, bouVertices, bouPatchIDs
    patchBouNbrs.resize(numPatches);
    for(const auto & p : patch2bous)
    {
      // we need to find the connected components from uedgeSet
      // each of the CC is a bou
      map<int, bool> allVtxVisited;
      map<int, set<UEdgeKey>> allVtxNbr; // vtx -> nbr uedges
      for(auto ue : p.second) // for all uedges on the pairs
      {
        allVtxVisited[ue[0]] = false;
        allVtxVisited[ue[1]] = false;
        allVtxNbr[ue[0]].insert(ue);
        allVtxNbr[ue[1]].insert(ue);
      }

      for(auto & allVtxVisitedPair : allVtxVisited)
      {
        // if this vtx (allVtxVisitedPair,first) is visited, continue
        if (allVtxVisitedPair.second == true) { continue; }
        // else, label this vtx as visited
        allVtxVisitedPair.second = true;
        set<int> bouVtx = { allVtxVisitedPair.first }; // assign seed vtx
        // use BFS to find all the other vtx in the same connected component on uedges as the seed
        // also find the edges BFS visited and stored in set<UEdgeKey> bouEdges
        set<UEdgeKey> bouEdges;
        set<int> candidates = bouVtx, nextCandidates; // buffers for BFS
        while(candidates.size() > 0)
        {
          for(int vtx : candidates)
          {
            for(const UEdgeKey & nbrVtxPair : allVtxNbr[vtx])
            {
              bouEdges.insert(nbrVtxPair); // store the edge we visited
              for(int i = 0; i < 2; i++)
              {
                int nbrVtx = nbrVtxPair[i];
                if (allVtxVisited[nbrVtx] == true) { continue; }
                bouVtx.insert(nbrVtx); // bouVtx stores the vtx visited on the connected component on uedges
                nextCandidates.insert(nbrVtx);
                allVtxVisited[nbrVtx] = true;
              }
            }
          }
          candidates.swap(nextCandidates);
          nextCandidates.clear();
        }

        // record this bou
        int bouID = bous.size();
        bous.emplace_back(move(bouEdges));
        bouVertices.emplace_back(move(bouVtx));


        vector<pair<int, int>> pIDs(p.first.begin(), p.first.end()); // <patchID,loopID> pairs
        // get the two patchIDs to build patchBouNbrs:
        int pID0 = pIDs[0].first, pID1 = pIDs[1].first;
        if (pID0 > pID1) swap(pID0, pID1);
        patchBouNbrs[pID0][bouID] = pID1;
        patchBouNbrs[pID1][bouID] = pID0;
        bouPatchIDs.emplace_back(std::array<int,2>({pID0, pID1}));
      } // end for one BFS
    } // end for <patchID, loopID> pairs
    numBous = bous.size();
  }

  // ===========================================================
  //               Now Bous have been computed
  //               Lets build arcs
  // ===========================================================

  // build arcs from bous
  std::vector<std::map<UEdgeKey, UEdgeKey>> arcEdge2bouEdge(numBous); // bouID -> arcEdgeKey -> bouEdgeKey
  map<set<UEdgeKey>, set<int>> arcEdges2Bous; // arcEdgeKey group representing one arc -> bouIDs lie on this arc
  for(int bouID = 0; bouID < numBous; bouID++)
  {
    set<UEdgeKey> stueSet; // stores UEdgeKeys in stitchIDs, which are actually arcEdgeKeys
    for(auto ue : bous[bouID])
    {
      UEdgeKey stue (vtxStitchIDs[ue[0]], vtxStitchIDs[ue[1]]); // stitchID UEdge
      stueSet.insert(stue);
      if(arcEdge2bouEdge[bouID].find(stue) != arcEdge2bouEdge[bouID].end()) // if this stUEdge has been visited on this bouID
      {
        cout << "Error: one bou has edges collide with each other, bou " << bouID <<
            " edge is " << ue << " and " << arcEdge2bouEdge[bouID][stue] <<  " st edge " << stue << endl;
        if (cutMeshSaved == false && verbose) { cutMesh.save("cutMesh.obj"); cutMeshSaved = true; }
        cout << "saving this bou: bouID " << bouID << endl;
        for(auto ue: bous[bouID]) { cout << ue[0] << " " << ue[1] << " "; }
        cout << endl;
        cout << "Done" << endl;
        throw 1;
      }
      arcEdge2bouEdge[bouID][stue] = ue;
    }
    arcEdges2Bous[stueSet].insert(bouID);
  }

  // build arc2Bous, bou2Arcs and arcs
  numArcs = arcEdges2Bous.size();
  cout << "Found #pre-arcs = " << numBous << " and #arcs = " << numArcs << "." << endl;
  
  arc2Bous.resize(numArcs);
  bou2Arcs.resize(numBous, -1);
  for(const auto & p : arcEdges2Bous)
  {
    if (p.second.size() != 2) // exactly two bous lie on one arc
    {
      cout << "Error: one arc has " << p.second.size() << " corresponding bous" << endl;
      throw 1;
    }
    int arcID = arcs.size();
    arcs.push_back(p.first); // push the set of arcEdgeKeys into arcs
    for(int bouID : p.second) // for the bous lie on this arc
    {
      arc2Bous[arcID].insert(bouID);
      assert(bou2Arcs[bouID] < 0 || bou2Arcs[bouID] == arcID);
      bou2Arcs[bouID] = arcID;
    }
  }
  // check all bous have corresponding arcs
  assert(find(bou2Arcs.begin(), bou2Arcs.end(), -1) == bou2Arcs.end());
  if (verbose)
  {
    for(int arcID = 0; arcID < numArcs; arcID++)
      cout << "arcID " << arcID << " has bou: " << streamRange(arc2Bous[arcID]) << endl;
  }

  // build patchNbrArcs, arcNbrPatches, patchNbrArcBous
  patchNbrArcs.resize(numPatches);
  arcNbrPatches.resize(numArcs);
  patchNbrArcBous.resize(numPatches);
  for(int patchID = 0; patchID < numPatches; patchID++)
  {
    for(const auto & p : patchBouNbrs[patchID])
    {
      int bouID = p.first;
      int arcID = bou2Arcs[bouID];
      patchNbrArcs[patchID].insert(arcID);
      arcNbrPatches[arcID].insert(patchID);
      patchNbrArcBous[patchID][arcID].insert(bouID);
    }
  }

  if (verbose)
  {
    for(int patchID = 0; patchID < numPatches; patchID++)
    {
      cout << "patch " << patchID << " GeoNbr: ";
      for(auto nbr : patchBouNbrs[patchID])
      {
        cout << "(" << nbr.first << ", " << nbr.second << ") ";
      }
      cout << endl;
      for(auto p : patchNbrArcBous[patchID])
      {
        if (p.second.size() > 1)
          cout << "Note! patch " << patchID << " at arc " << p.first << " has " << p.second.size() << " bous!" << endl;
      }
    }
    for(int bouID = 0; bouID < numBous; bouID++)
    {
      string filename = "bou" + to_string(bouID) + ".txt";
      ofstream fout(filename.c_str());
      for(auto edge : bous[bouID])
        fout << edge[0] << " " << edge[1] << endl;
      fout.close();
    }
  }


  // ===========================================================
  //      Finally, build B-patch neigboring data on a cell
  // ===========================================================

  // build cellPatchNeighbors, cellBouCorrespondence
  cellPatchBouNbrs.resize(numCells);
  cellBouCorrespondence.resize(numCells);
  for(int cellID = 1; cellID < numCells; cellID++)
  {
    if (verbose)
      cout << "Build cellPatchBouNbrs: " << cellID << endl;
    // build cell patch neighbors for pID0 and pID1. The bou of pID0 at this arc is b0. The bou of pID1 at the arc is b1.
    auto buildCellPatchBouNbr = [&](int pID0, int pID1, int b0, int b1)
    {
      assert(b0 != b1);
      assert(bou2Arcs[b0] == bou2Arcs[b1]); // b0 and b1 must lie on the same arc

      if (verbose)
        cout << "add cell " << cellID << " patch nbr: " << pID0 << '-' << b0 << " " << pID1 << '-' << b1 << endl;
      // patchBouNbrs: patchID -> nbring bouID -> topological nbring patchID sharing the bou
      assert(mapFind(patchBouNbrs[pID0], b0));
      assert(mapFind(patchBouNbrs[pID1], b1));
      cellPatchBouNbrs[cellID][pID0][b0] = pair<int, int>(b1, pID1);
      cellPatchBouNbrs[cellID][pID1][b1] = pair<int, int>(b0, pID0);

      UEdgeKey bouPair(b0, b1);
      // the two bous lie on the same arc, so their #vertices should be the same
      assert(bouVertices[b0].size() == bouVertices[b1].size());
      if (mapFind(cellBouCorrespondence[cellID], bouPair)) return; // if this bou pair has been built, then return

      // now let's build the correspondence on b0 and b1
      vector<pair<int,int>> vtxCor;
      map<int, int> stvtx2b0vtx; // stitchID -> cutVtxID in b0
      for(int vID : bouVertices[b0])
        stvtx2b0vtx[vtxStitchIDs[vID]] = vID;
      if(stvtx2b0vtx.size() == bouVertices[b0].size()) // one to one mapping of stitchID and b0VtxID
      {
        for(int vID : bouVertices[b1])
        {
          int stID = vtxStitchIDs[vID];
          assert(mapFind(stvtx2b0vtx, stID));
          vtxCor.emplace_back(stvtx2b0vtx[stID], vID);
        }
      }
      else // more than one v0 vtx match one stID (stID = stitchID)
      {    // this can happen when the bou is a ring shape, then the start and end vtx of the ring have the same stID
        map<int, set<int>> stvtx2b0vtxSet; // stitchID -> set of b0VtxIDs
        for(int vID : bouVertices[b0])
          stvtx2b0vtxSet[vtxStitchIDs[vID]].insert(vID);

        if (verbose)
        {
          int nonSingleCount = 0;
          for(const auto & p : stvtx2b0vtxSet)
            if (p.second.size() > 1)
            {
              assert(p.second.size() == 2);
              nonSingleCount++;
            }
          cout << "found non-trivial boundary at b0 : " << b0 << " with nonSingleCount: " << nonSingleCount << endl;
        }

        // first pass, find correspondance on those unambiguous vertices
        for(int vID : bouVertices[b1])
        {
          int stID = vtxStitchIDs[vID];
          assert(mapFind(stvtx2b0vtxSet, stID));
          if (stvtx2b0vtxSet[stID].size() == 1) // the stID corresponds to one bou0 vtx, then no ambiguity
            vtxCor.emplace_back(stvtx2b0vtx[stID], vID);
        }
        // second pass, we use the neighboring information on the vertices to determine the matches
        // if b0a, b0b, b1a, b1b all have the same stID, then
        // we use the topological neighboring vtx of b0a to find whether it should match b1a or b1b
        // since no way that the neighboring vtx of b0a: b0c( and b0d, if any) are ambiguous as well,
        // otherwise, two UEdgeKeys on bou0 map to the same arcEdgeKey, which is forbiden and checked beforehand
        for(int b1vID : bouVertices[b1])
        {
          int stID = vtxStitchIDs[b1vID];
          // in this loop, we deal with those b1vID with ambiguious correspondence
          if (stvtx2b0vtxSet[stID].size() == 1) continue;

          map<int, int> targetCount; // candidate of b0vtx that maps to b1vID -> confidence count
          // find its neighbor
          auto b1vtxNbr = cutMeshNbr.getVtxNearbyVertices(b1vID, cutMesh);
          for(int nv : b1vtxNbr)  // visit b1 nbrs
          {
            if (setNotFind(bouVertices[b1], nv)) continue;
            // for a neighboring b1 vtx: nv
            int nbrstID = vtxStitchIDs[nv]; // nv's stID
            if (mapNotFind(stvtx2b0vtxSet, nbrstID)) continue;
            if (stvtx2b0vtxSet[nbrstID].size() > 1) continue;
            if(mapNotFind(stvtx2b0vtx, nbrstID)) continue;
            int nbrb0v = stvtx2b0vtx[nbrstID]; // the b0 vtx mapping to this nv
            for(int b0Cand : stvtx2b0vtxSet[stID]) // b0 vtx that map to the b1 vtx with stID
            {
              // if b0Cand is a neighbor of nbrb0v, then we treat it as a candidate for mapping to b1vID
              if (cutMeshNbr.areVerticesNeighbors(nbrb0v, b0Cand)
                  && setFind(bous[b0], UEdgeKey(nbrb0v, b0Cand)))
              {
                if (mapFind(targetCount, b0Cand)) targetCount[b0Cand]++;
                else targetCount[b0Cand] = 1;
              }
            }
          }
          assert(targetCount.size() == 1); // should only have one candidate
          int b0v = targetCount.begin()->first;
          vtxCor.emplace_back(b0v, b1vID);
        }
      }

      if (bouPair[0] == b1)  // switch order in pairs
      {
        for(auto & p : vtxCor) swap(p.first, p.second);
      }
      cellBouCorrespondence[cellID][bouPair] = move(vtxCor); // finally, build the correspondence on b0 and b1
    };

    // first, build arc2patches: arcID -> <patchIDs, bouIDs> on this cell
    // here we define <patchID, bouID> as "patch end"
    // around one arc, there are 4 and only 4 patch ends
    map<int, set<pair<int, int>>> arc2patches;
    for(auto p : cellPatches[cellID]) // patchID -> whether this patch points outward for this cell
    {
      int patchID = p.first;
      for(const auto & patchBou : patchBouNbrs[patchID])
      {
        int bouID = patchBou.first;
        int arcID = bou2Arcs[bouID];
        arc2patches[arcID].emplace(patchID, bouID);
      }
    }

    if (verbose)
      cout << arc2patches.size() << " arcs in cell " << cellID << endl;

    // then, try to build the B-patch correspondence around each arcID
    for(const auto & p : arc2patches)
    {
      int arcID = p.first;
      const auto & patchIDs = p.second;
      if(verbose)
      {
        cout << "in cell " << cellID <<", " << patchIDs.size() << " patch ends" << endl;
        for(auto p : patchIDs) cout << p.first << '-' << p.second << " ";
        cout << endl;
      }
      if (patchIDs.size() != 2 && patchIDs.size() != 4) // around one arc there can only be 2 or 4 patch ends
      {
        cout << "Error: in cell " << cellID << " an arc " << arcID << " has strange patch-bou pairs: ";
        for(auto pb : patchIDs) cout << "(" << pb.first << " " << pb.second << ") ";
        cout << endl;
        throw 1;
      }
      if (patchIDs.size() == 2) // simple case
      {
        auto it = patchIDs.begin();
        int pID0 = it->first, b0 = it->second;
        it++;
        int pID1 = it->first, b1 = it->second;
        // normal situation
        if (b0 == b1)
        {
          cout << "Error: in cell " << cellID << " patch " << pID0 << " and " << pID1 << " has an arc " << arcID;
          cout <<" but the bou is the same: " << b0 << endl;
        }
        assert(b0 != b1);
        buildCellPatchBouNbr(pID0, pID1, b0, b1);
      }
      else if (patchIDs.size() == 4) // difficult case, happens on the "self-touching" cell
      {
        // we first get one patch end: pID0, bou0
        auto patchIt = patchIDs.begin();
        int pID0 = patchIt->first, bou0 = patchIt->second;
        patchIt++;
        // then we try to find which of the rest three patch end should match <pID0, bou0>

        auto getArcEdgeLen2 = [&](UEdgeKey edge)
        {
          return len2(stitchPositions[edge[0]] - stitchPositions[edge[1]]);
        };
        // to improve robustness, we find the UEdge with the longest length
        auto maxLenEdge = *maxFunctionValue(arcs[arcID].begin(), arcs[arcID].end(), getArcEdgeLen2);

        auto bouEdge0 = arcEdge2bouEdge[bou0][maxLenEdge]; // the longest bou edge
        assert(bouEdge0[0] >= 0 && bouEdge0[1] >= 0);
        // get the triangleID on the patch pID0 at bouEdge0
        int triID0 = patchTriNbr[pID0].getTriangleAtEdge({bouEdge0[0], bouEdge0[1]});
        if (triID0 < 0) triID0 = patchTriNbr[pID0].getTriangleAtEdge({bouEdge0[1], bouEdge0[0]});
        assert(triID0 >= 0);
        Vec3i t0 = patchTris[pID0][triID0];
        int t0vtxID = getTriangleVertexOppositeEdge(t0, bouEdge0);
        assert(t0vtxID >= 0);
        OEdgeKey triEdge0 = getTriangleOEdge(t0, bouEdge0); // get the OEdgeKey on the triangle triID0
        assert(triEdge0[0] >= 0 && triEdge0[1] >= 0);
        OEdgeKey stTriEdge0(vtxStitchIDs[triEdge0[0]], vtxStitchIDs[triEdge0[1]]); // get the OEdgeKey in stitchID

        // get the exact pos of the triangle triID0
        const Vec3ER & t0Pos0 = selfCutMesh.cutPosExact[t0[0]];
        const Vec3ER & t0Pos1 = selfCutMesh.cutPosExact[t0[1]];
        const Vec3ER & t0Pos2 = selfCutMesh.cutPosExact[t0[2]];
        Vec3ER t0ScaledNormal = cross(t0Pos1 - t0Pos0, t0Pos2 - t0Pos0);

        Vec3d t0normal = patchMeshRef[pID0].computeTriangleNormal(triID0);
        if (t0normal.hasNaN()) // improve the normal of triID0 with exact computation
        {
          for(int i = 0; i < 3; i++)
            t0normal[i] = ER_toDouble(t0ScaledNormal[i]);
          assert(t0normal.hasNaN() == false);
        }

        if (cellPatches[cellID][pID0] == false) // if the triangle normal orientation of the patchID0 is different from the cell
        {
          t0ScaledNormal = (-1) * t0ScaledNormal; // then we reverse relevant data
          t0normal *= (-1);
          stTriEdge0.reverse();
        }

        if (verbose)
        {
          cout << "bouEdge0 " << bouEdge0 << endl;
          cout << "triEdge0 " << triEdge0 << endl;
          cout << "stTriEdge0 " << OEdgeKey(vtxStitchIDs[triEdge0[0]], vtxStitchIDs[triEdge0[1]]) << endl;
          cout << "ori " << cellPatches[cellID][pID0] << " " << stTriEdge0 << endl;
        }

        const Vec3ER & t0OpPos = selfCutMesh.cutPosExact[t0vtxID];
        const Vec3ER & edgePos = selfCutMesh.cutPosExact[bouEdge0[0]];
        const Vec3ER t0Vec = t0OpPos - edgePos; // a vector from one vtx on the edgeKey to the opposite vtx on triangle triID0
        double closestAngle = 2*M_PI;      // store the dihedral angle between the patch end <pID0,bou0> and the final found patch end at the edgeKey
        auto closestIt = patchIDs.begin(); // the iterator to point to the final patch end we will find
        for(; patchIt != patchIDs.end(); patchIt++) // go through the other three patch ends
        {
          int pID1 = patchIt->first;
          int bou1 = patchIt->second;
          auto bouEdge1 = arcEdge2bouEdge[bou1][maxLenEdge];

          assert(bouEdge1[0] >= 0 && bouEdge1[1] >= 0);
          // get the triangle sharing the same longest edgeKey
          int triID1 = patchTriNbr[pID1].getTriangleAtEdge({bouEdge1[0], bouEdge1[1]});
          if (triID1 < 0) triID1 = patchTriNbr[pID1].getTriangleAtEdge({bouEdge1[1], bouEdge1[0]});
          assert(triID1 >= 0);
          Vec3i t1 = patchTris[pID1][triID1];
          OEdgeKey triEdge1 = getTriangleOEdge(t1, bouEdge1);
          OEdgeKey stTriEdge1(vtxStitchIDs[triEdge1[0]], vtxStitchIDs[triEdge1[1]]);
          if (cellPatches[cellID][pID1] == false)
          {
            stTriEdge1.reverse();
          }

          if (verbose)
          {
            cout << pID1 << "-" << bou1 << " bouEdge1 " << bouEdge1 << endl;
            cout << "triEdge1 " << triEdge1 << endl;
            cout << "stTriEdge " << OEdgeKey(vtxStitchIDs[triEdge1[0]], vtxStitchIDs[triEdge1[1]]) << endl;
            cout << "ori " << cellPatches[cellID][pID1] << " " << stTriEdge1 << endl;
          }

          // we check whether this patch end can form a locally manifold shape with <pID0, bou0> on this cell
          // if the orientation of the stitchID edgeKey is the same, then the orientation will not be consistent
          // if this patch end is mapped to <pID0, bou0>, so we skip this patch end
          if (stTriEdge0 == stTriEdge1) // same orientation on the edge key
          {
            if (verbose)
            {
              cout << "skip " << pID1 << "-" << bou1 << endl;
            }
            continue;
          }

          const Vec3ER & t1Pos0 = selfCutMesh.cutPosExact[t1[0]];
          const Vec3ER & t1Pos1 = selfCutMesh.cutPosExact[t1[1]];
          const Vec3ER & t1Pos2 = selfCutMesh.cutPosExact[t1[2]];
          Vec3ER t1ScaledNormal = cross(t1Pos1 - t1Pos0, t1Pos2 - t1Pos0);
          Vec3d t1normal = patchMeshRef[pID1].computeTriangleNormal(triID1);
          if (cellPatches[cellID][pID1] == false) // again, reverse normals if patch orientation not agree with the cell
          {
            t1ScaledNormal = (-1) * t1ScaledNormal;
            t1normal *= (-1);
          }

          // compute the dihedral angle between the two triangle t0 and t1 from the two patch ends
          double angle = getVectorAngle(t0normal, t1normal);

          if (isnan(angle))
          {
            cout << "Error: angle is NAN: at cell" << cellID;
            for(auto p : patchIDs)
              cout << " patch " << p.first << " bou " << p.second;
            cout << endl;
            assert(isnan(angle) == false);
          }

          // the actual dihedral angle can be angle or 2*M_PI-angle
          // we choose the correct one from the two based on the computation of the following code
          if (ER_sign(dot(t0Vec, t1ScaledNormal)) > 0)
          {
            angle = 2*M_PI - angle;
          }

          if (verbose)
          {
            cout << "at cell " << cellID << " find an angle: patch " << pID0 << "-" << bou0 << " and " << pID1 << "-" << bou1 <<
                " is " << angle * 180.0 / M_PI << endl;
          }

          // store the patch end with the smallest dihedral angle
          if (angle < closestAngle)
          {
            closestIt = patchIt;
            closestAngle = angle;
          }
        } // end of patchIt

        if (closestAngle >= 2*M_PI) // fail to find a valid patch end to match <pID0, bou0>
        {
          cout << "Error at 4 patch case: " << cellID << " pID0 " << pID0 << " the rest: ";
          for(auto p : patchIDs)
            cout << p.first << " ";
          cout << endl;
          throw 1;
        }

        assert(closestAngle < 2 * M_PI);
        int pID1 = closestIt->first;
        int bou1 = closestIt->second;
        if (verbose)
          cout << "the closest patch to " << pID0 << "-" << bou0 << " is " << pID1 << "-" << bou1 << endl;

        buildCellPatchBouNbr(pID0, pID1, bou0, bou1); // build the mapping

        // now, build the mapping for the rest two patch endings
        vector<pair<int,int>> restPatches;
        for(auto p : patchIDs)
        {
          if (p != make_pair(pID0, bou0) && p != make_pair(pID1, bou1))
            restPatches.push_back(p);
        }
        assert(restPatches.size() == 2);
        buildCellPatchBouNbr(restPatches[0].first, restPatches[1].first, restPatches[0].second, restPatches[1].second);
      }
    } // end loop on arcID

    // we comment out the check below because if there can be a cell with only one B-patch and this B-patch will have no geometric nbrs
    // so the size of cellPatchBouNbrs[cellID] will be zero, while cellPatches[cellID] is one
    // assert(cellPatchBouNbrs[cellID].size() == cellPatches[cellID].size());
    for(const auto & patchPair : cellPatchBouNbrs[cellID])
    {
      int patchID = patchPair.first;
      assert(patchPair.second.size() == patchBouNbrs[patchID].size());
      for(const auto & p2 : patchPair.second)
      {
        int bouID = p2.first;
        int nbrBouID = p2.second.first;
        int nbrPatchID = p2.second.second;
        assert(patchBouNbrs[nbrPatchID].find(nbrBouID) != patchBouNbrs[nbrPatchID].end());
        assert(patchBouNbrs[patchID].find(bouID) != patchBouNbrs[patchID].end());
      }
    }
  } // end loop on cellID
}

