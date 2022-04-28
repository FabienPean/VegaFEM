/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
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

#include "triMeshNeighbor.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include <unordered_map>
#include <cassert>
#include <set>
#include <iostream>
using namespace std;

namespace
{

vector<vector<int>> findBoundaryLoops(const std::vector<Vec3i> & triangles, const std::vector<Vec3i> & triNbrs)
{
  assert(triNbrs.size() == triangles.size());
  vector<vector<int>> ret;

  vector<bool> triVisited(triangles.size()*3, false);
  for(size_t triID = 0; triID < triNbrs.size(); triID++)
  {
    const auto & n = triNbrs[triID];
    const auto & t = triangles[triID];
//    cout << "visit triID " << triID << endl;
    for(int i = 0; i < 3; i++)
    {
      if (triVisited[triID*3+i]) continue;
      triVisited[triID*3+i] = true;
      if (n[i] >= 0) continue;
      // find a boundary, let's go through the loop and find other edges in this loop
      int t1 = t[(i+1)%3];
      vector<int> loop = { t[i], t1 };

      int start = t[i];
      int end = t1;
      int curTriID = triID;
//      cout << "start a loop at " << t[i] << " " << t1 << endl;

      int nextI = (i+1) % 3;
      while(true)
      {
        // end = tri[curTriID][nextI+1]
        // check whether the edge <tri[curTriID][nextI] = end, tri[curTriID][nextI+1]> is boundary
        triVisited[curTriID*3+nextI] = true;
        if (triNbrs[curTriID][nextI] < 0) // the next edge on cur triID is a boundary
        {
          end = triangles[curTriID][(nextI+1)%3];
          if (end == start) { break; } // end this loop
          loop.push_back(end);
          nextI = (nextI+1)%3;
//          cout << "find next vtx on triangle: " << end << endl;
          continue;
        }

        curTriID = triNbrs[curTriID][nextI];
//        assert(curTriID != triID || triVisited[curTriID] == false);
        // next cutTriID has the edge <tri[old curTriID][nextI+1], tri[old curTriID][nextI] = end>
        nextI = triangles[curTriID].getInvertedIndex(end);
        assert(nextI >= 0);
//        cout << "move to another triangle triID " <<  curTriID << endl;
      }

      ret.emplace_back(move(loop));
      break; // because one triangle can only be in one loop
    } // end for i
  }

  return ret;
}

vector<pair<int, OEdgeKey>> findBoundaryTriangles(const std::vector<Vec3i> & triangles, const std::vector<Vec3i> & triNbrs)
{
  assert(triNbrs.size() == triangles.size());
  vector<pair<int, OEdgeKey>> ret;
  for(size_t triID = 0; triID < triNbrs.size(); triID++)
  {
    const auto & n = triNbrs[triID];
    const auto & t = triangles[triID];
    for(int i = 0; i < 3; i++)
    {
      if (n[i] < 0) // find a boundary
      {
        ret.push_back(make_pair(triID, OEdgeKey(t[i], t[(i+1)%3])));
      }
    }
  }
  return ret;
}

} // end anonymous namespace



TriangleNeighbor::TriangleNeighbor(const vector<Vec3i> & triangles) : numTriangles(triangles.size()),
    triNbrs(triangles.size(), Vec3i(-1))
{
  if (getOEdgeTriMap(triangles, oedgeTri) == false) throw 1;

  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int j = 0; j < 3; j++)
    {
      if (triNbrs[triID][j] < 0)
      {
        int v0 = triangles[triID][j];
        int v1 = triangles[triID][(j+1)%3];
        auto iter = oedgeTri.find(OEdgeKey(v1,v0)); // get the revserse key
        if (iter != oedgeTri.end())
        {
          int triID2 = iter->second;
          triNbrs[triID][j] = triID2;
          int j2 = triangles[triID2].getInvertedIndex(v1);
          assert(j2 >= 0);
          triNbrs[triID2][j2] = triID;
        }
      }
    }
  }
}

int TriangleNeighbor::getTriangleAtEdge(const OEdgeKey & edge)
{
  auto it = oedgeTri.find(edge);
  if (it == oedgeTri.end()) return -1;
  return it->second;
}

vector<pair<int, OEdgeKey>> TriangleNeighbor::findBoundaryTriangles(const std::vector<Vec3i> & triangles) const
{
  return ::findBoundaryTriangles(triangles, triNbrs);
}

vector<vector<int>> TriangleNeighbor::findBoundaryLoops(const std::vector<Vec3i> & triangles) const
{
  return ::findBoundaryLoops(triangles, triNbrs);
}

vector<int> TriangleNeighbor::findTrianglesArroundBoundaryVertex(int prevVtxID, int vtxID, const vector<Vec3i> & triangles) const
{
  assert(triNbrs.size() == triangles.size());
  vector<int> ret;
  OEdgeKey prevEdge(prevVtxID, vtxID);
  auto it = oedgeTri.find(prevEdge);
  if (it == oedgeTri.end()) return ret;
  int triID = it->second;
  do
  {
    ret.push_back(triID);
    int i = triangles[triID].getInvertedIndex(vtxID);
    assert(i >= 0);
    triID = triNbrs[triID][i];
  } while(triID >= 0);
  return ret;
}

// ==============================================================================
//                              TriMeshNeighbor
// ==============================================================================

TriMeshNeighbor::TriMeshNeighbor(TriMeshRef triMesh) : numVertices(triMesh.numVertices()), numTriangles(triMesh.numTriangles()),
    vtxNbrTri(triMesh.numVertices()), triNbrs(triMesh.numTriangles(), Vec3i(-1)), vtxEdgeTri(triMesh.numVertices())
{
  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int j = 0; j < 3; j++)
    {
      int v0 = triMesh.tri(triID)[j];
      int v1 = triMesh.tri(triID)[(j+1)%3];
      if (v0 == v1) // error, TriMesh is not valid!
      {
        throw 1;
      }
      vtxNbrTri[v0].push_back(triID);
      if (vtxEdgeTri[v0].find(v1) != vtxEdgeTri[v0].end())
      {
        throw 1; // error, TriMesh is not manifold!
      }
      vtxEdgeTri[v0][v1] = triID;
    }
  }

  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int j = 0; j < 3; j++)
    {
      if (triNbrs[triID][j] < 0)
      {
        int v0 = triMesh.tri(triID)[j];
        int v1 = triMesh.tri(triID)[(j+1)%3];
        auto iter = vtxEdgeTri[v1].find(v0);
        if (iter != vtxEdgeTri[v1].end())
        {
          int triID2 = iter->second;
          triNbrs[triID][j] = triID2;
          int j2 = triMesh.tri(triID2).getInvertedIndex(v1);
          assert(j2 >= 0);
          triNbrs[triID2][j2] = triID;
        }
      }
    }
  }
}

vector<pair<int, OEdgeKey>> TriMeshNeighbor::findBoundaryTriangles(const std::vector<Vec3i> & triangles) const
{
  return ::findBoundaryTriangles(triangles, triNbrs);
}

vector<vector<int>> TriMeshNeighbor::findBoundaryLoops(const std::vector<Vec3i> & triangles) const
{
  return ::findBoundaryLoops(triangles, triNbrs);
}

vector<int> TriMeshNeighbor::getVtxNearbyVertices(int vtxID, const TriMeshGeo & mesh) const
{
  assert(numVertices == mesh.numVertices() && numTriangles == mesh.numTriangles());
  vector<int> ret;
  for(int triID : vtxNbrTri[vtxID])
  {
    const auto & t = mesh.tri(triID);
    ret.insert(ret.end(), t.begin(), t.end());
  }
  sortAndDeduplicate(ret);
  return ret;
}

bool TriMeshNeighbor::areVerticesNeighbors(int vtxID0, int vtxID1) const
{
  return mapFind(vtxEdgeTri[vtxID0], vtxID1) || mapFind(vtxEdgeTri[vtxID1], vtxID0);
}

// ==============================================================================
//                          NonManifoldTriangleNeighbor
// ==============================================================================

NonManifoldTriangleNeighbor::NonManifoldTriangleNeighbor(const vector<Vec3i> & triangles) :
    numTriangles(triangles.size())
{
  if (getNonManifoldOEdgeTrisMap(triangles, oedgeTris) == false) throw 1;
}

vector<int> NonManifoldTriangleNeighbor::getTriangleAtOEdge(const OEdgeKey & oedge) const
{
  auto it = oedgeTris.find(oedge);
  if (it == oedgeTris.end()) return {};
  return it->second;
}

vector<OEdgeKey> NonManifoldTriangleNeighbor::findNonManifoldOEdges() const
{
  vector<OEdgeKey> ret;
  for(const auto & p : oedgeTris)
  {
    if (p.second.size() > 1) ret.push_back(p.first);
  }
  return ret;
}
/////////////////////////////////////////////////////////////////////////
//                       Other Functions
/////////////////////////////////////////////////////////////////////////

namespace
{

using Neighbor = std::function<const std::vector<int> &(int nodeA)>;
vector<vector<int>> getConnectedComponents(int numNodes, Neighbor getNeighbor)
{
  vector<vector<int>> ret;
  if (numNodes == 0) return ret;

  vector<bool> visited(numNodes, false);

  for(int seedID = 0; seedID < numNodes; seedID++)
  {
    if (visited[seedID]) continue;

    vector<int> component;  // record traversed triangles

    visited[seedID] = true;
    component.push_back(seedID);
    size_t candidateBegin = 0, candidateEnd = 1;

    while(candidateBegin != candidateEnd)
    {
      for(size_t i = candidateBegin; i < candidateEnd; i++)
      {
        int node = component[i];
        const auto & nbrs = getNeighbor(node);
        for(int nbr : nbrs)
        {
          if (nbr < 0) continue;
          if (visited[nbr]) continue;
          visited[nbr] = true;
          component.push_back(nbr);
        }
      }
      candidateBegin = candidateEnd;
      candidateEnd = component.size();
    }

    sort(component.begin(), component.end());
    assert(unique(component.begin(), component.end()) == component.end());
    ret.emplace_back(move(component));
  }
  return ret;
}

} // anonymous namespace

// build oedgeTri: <v0,v1> -> tri with ordered edge <v0, v1>
bool getOEdgeTriMap(const std::vector<Vec3i> & triangles, unordered_map<OEdgeKey, int> & oedgeTri)
{
  int numTriangles = triangles.size();
  oedgeTri.clear();

  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int j = 0; j < 3; j++)
    {
      int v0 = triangles[triID][j];
      int v1 = triangles[triID][(j+1)%3];
      if (v0 == v1 || v0 < 0) // error, TriMesh is not valid!
      {
        return false;
      }
      OEdgeKey edge(v0, v1);
      if (oedgeTri.find(edge) != oedgeTri.end())
      {
        return false; // this signed edge appears twice, not manifold!
      }
      oedgeTri.insert(make_pair(edge, triID));
    }
  }
  return true;
}

// go through a triangle fan around vtxID on an edge-manifold mesh
// input includes vtxID, a startingTriID to start the search,
// and the localVtxID == triangles[startingTriID].getInvertedIndex(vtxID), localVtxID: [0,3)
// once a new triangle is visited, its triID and the local vtxID [0,3) of the input vtxID is passed to
// processTriangle for custom processing
static void visitTriangleFanNeighboringVtx(const std::vector<Vec3i> & triangles,
    const unordered_map<OEdgeKey, int> & oedgeTri, int vtxID, int startingTriID, int localVtxID,
    function<void(int newTriID, int newLocalVtxID)> processTriangle)
{
  // visit the neighboring triangles to search for the fan shape.
  // ((lvID + edgeVtxOffset)%3, (lvID + edgeVtxOffset+1)%3) represents an OEdge on triangle triID
  // return true if we go 360 degrees around the vtx, finishing one full loop.
  // if it returns false, then we hit a boundary edge on the triangle mesh
  auto visitNbr = [&](int edgeVtxOffset)
  {
    //        cout << "begin search" << endl;
    // loop over triangles around this vtx, starting at the edge (vEdgeStart, vEdgeEnd)
    int vEdgeStart = triangles[startingTriID][(localVtxID + edgeVtxOffset)%3];
    int vEdgeEnd = triangles[startingTriID][(localVtxID + edgeVtxOffset+1)%3];
    // edge<vEdgeStart, vEdgeEnd> is an OEdge of triID
    // since we want to find the nbring tri on this edge, we should search the reversed edge in oedgeTri
    auto it = oedgeTri.find({vEdgeEnd, vEdgeStart});
    while(it != oedgeTri.end())
    {
      int nextTriID = it->second;
      if(nextTriID == startingTriID) return true; // we finished one full loop around the vtx v
      int nextlvID = triangles[nextTriID].getInvertedIndex(vtxID);
      assert(nextlvID >= 0);

      processTriangle(nextTriID, nextlvID);

      vEdgeStart = triangles[nextTriID][(nextlvID + edgeVtxOffset)%3];
      vEdgeEnd = triangles[nextTriID][(nextlvID + edgeVtxOffset+1)%3];
      it = oedgeTri.find({vEdgeEnd, vEdgeStart});
    }
    return false;
  };

  if (visitNbr(0) == false) // visitNbr(0) return true only if it finishes one complete loop
    visitNbr(2);
}

// The algorithm here is that for each vtx, we visit each neighboring triangles and record this "fan shape" we visited.
// If the mesh is vtx-manifold, then there will only be one fan shape for each vtx.
// But on a non-vtx-manifold but edge-manifold mesh, one vtx can have more than one fan shape.
vector<int> getNonManifoldVerticesOnEdgeManifoldTriangles(const std::vector<Vec3i> & triangles,
    const unordered_map<OEdgeKey, int> & oedgeTri)
{
  vector<int> ret;
  int numTriangles = triangles.size();
  set<int> vtxVisited;
  vector<bool> triVtxVisited(numTriangles * 3, false);

  for(int triID = 0;  triID < numTriangles; triID++)
  {
    for(int lvID = 0; lvID < 3; lvID++) // local vtx ID
    {
      if (triVtxVisited[triID*3+lvID]) continue;
      int curVtxID = triangles[triID][lvID];
      if (vtxVisited.find(curVtxID) != vtxVisited.end()) // if this v has been visited
      {
        ret.push_back(curVtxID); // not vtx-manifold
//        cout << "non-manifold: " << triID << " " << v << endl;
        continue;
      }
      vtxVisited.insert(curVtxID);
      triVtxVisited[triID*3+lvID] = true;

//      cout << "visit triID " << triID << " v " << v << endl;
      visitTriangleFanNeighboringVtx(triangles, oedgeTri, curVtxID, triID, lvID, [&](int nextTriID, int nextlvID)
      {
        assert(triVtxVisited[nextTriID*3+nextlvID] == false);
        triVtxVisited[nextTriID*3+nextlvID] = true;
      });
    }
  }
  sortAndDeduplicate(ret);
  return ret;
}

void fixNonManifoldVerticesOnEdgeManifoldTriangles(TriMeshGeo & triMesh, const std::unordered_map<OEdgeKey, int> & oedgeTri,
    std::map<int, int> * newVtx2OldVtxMap)
{
  vector<int> nmVtxIDs = getNonManifoldVerticesOnEdgeManifoldTriangles(triMesh.triangles(), oedgeTri);
  if (nmVtxIDs.size() == 0) return;

  map<int, vector<int>> nmVtxNbringTris;
  for(int triID = 0; triID < triMesh.numTriangles(); triID++)
    for(int i = 0; i < 3; i++)
    {
      int vtxID = triMesh.triVtxID(triID, i);
      if (binarySearchFound(nmVtxIDs, vtxID))
      {
        nmVtxNbringTris[vtxID].push_back(triID);
      }
    }

  for(int nmVtxID : nmVtxIDs)
  {
    const vector<int> & nbringTriIDs = nmVtxNbringTris[nmVtxID];
    assert(nbringTriIDs.size() > 1);
    set<int> nbringTriIDVisited;
    vector<vector<int>> triFans; // stores the triIDs in each triangle fan around nmVtxID
    for(int i = 0; i < sizei(nbringTriIDs); i++)
    {
      int triID = nbringTriIDs[i];
      if (setFind(nbringTriIDVisited, triID)) continue;
      nbringTriIDVisited.insert(triID);

      vector<int> newFan = { triID };
      int lvID = triMesh.tri(triID).getInvertedIndex(nmVtxID); // local vtx ID, [0,3)
      visitTriangleFanNeighboringVtx(triMesh.triangles(), oedgeTri, nmVtxID, triID, lvID, [&](int newTriID, int newLvID)
      {
        assert(binarySearchFound(nbringTriIDs, newTriID));
        assert(setNotFind(nbringTriIDVisited, newTriID));
        nbringTriIDVisited.insert(newTriID);
        newFan.push_back(newTriID);
      });
//      cout << "Found a new fan, size " << newFan.size() << endl;
      triFans.emplace_back(move(newFan));
    }

    // sanity check
    size_t sum = 0;
    for(const auto & fan : triFans)
      sum += fan.size();
//    cout << "Fans:  " << triFans.size() << endl;
//    cout << "sum " << sum << " " << nbringTriIDs.size() << endl;
    assert(sum == nbringTriIDs.size());
    assert(triFans.size() > 1);

    for(int i = 1; i < sizei(triFans); i++) // for each addtional fan
    {
      int newVtxID = triMesh.numVertices();
      if (newVtx2OldVtxMap)
        newVtx2OldVtxMap->emplace(newVtxID, nmVtxID);
      triMesh.addPos(triMesh.pos(nmVtxID));
      for(int triID : triFans[i])
      {
        int lvID = triMesh.tri(triID).getInvertedIndex(nmVtxID);
        triMesh.tri(triID)[lvID] = newVtxID;
      }
    }
  } // for each nmVtxID
}


bool areTrianglesEdgeManifold(const std::vector<Vec3i> & triangles)
{
  std::unordered_map<OEdgeKey, int> oedgeTri; // <v0,v1> -> tri with ordered edge <v0, v1>
  return getOEdgeTriMap(triangles, oedgeTri);
}

bool areTrianglesManifold(const std::vector<Vec3i> & triangles)
{
//  int numTriangles = triangles.size();
  unordered_map<OEdgeKey, int> oedgeTri; // <v0,v1> -> tri with ordered edge <v0, v1>

  if (getOEdgeTriMap(triangles, oedgeTri) == false) return false;
  // now the mesh is at least edge-manifold

  return (getNonManifoldVerticesOnEdgeManifoldTriangles(triangles, oedgeTri).size() == 0);
}

std::vector<std::vector<int>> getTriangleNeighborsByEdge(const std::vector<Vec3i> & triangles)
{
  vector<vector<int>> nbrs(triangles.size());
  map<UEdgeKey, vector<int>> uedgeMap;
  for(size_t triID = 0; triID < triangles.size(); triID++)
  {
    const Vec3i & tri = triangles[triID];
    for(int j = 0; j < 3; j++)
    {
      UEdgeKey edge(tri[j], tri[(j+1)%3]);
      uedgeMap[edge].push_back(triID);
    }
  }
  for(const auto & p : uedgeMap)
  {
    const auto & tris = p.second;
    for(size_t i = 0; i < tris.size(); i++)
      for(size_t j = i+1; j < tris.size(); j++)
      {
        nbrs[tris[i]].push_back(tris[j]);
        nbrs[tris[j]].push_back(tris[i]);
      }
  }
  for(auto & vec : nbrs) // sort and remove duplicates in each nbr vector
  {
    sortAndDeduplicate(vec);
  }

  return nbrs;
}

vector<vector<int>> getConnectedComponentsByEdge(const std::vector<Vec3i> & triangles)
{
  if (triangles.size() == 0) return {};

  auto triNbrs = getTriangleNeighborsByEdge(triangles);

  auto getNeighbor = [&](int triID) -> const vector<int> & { return triNbrs[triID]; };

  return getConnectedComponents(triangles.size(), getNeighbor);
}

vector<std::vector<int>> getConnectedComponentsByVertex(const std::vector<Vec3i> & triangles)
{
  if (triangles.size() == 0) return {};


  map<int, vector<int>> vtxNbringTri;
  for(size_t i = 0; i < triangles.size(); i++)
  {
    for(int j = 0; j < 3; j++)
      vtxNbringTri[triangles[i][j]].push_back(i);
  }
  vector<vector<int>> triNbrTri(triangles.size());
  for(const auto & p : vtxNbringTri)
  {
    const auto & t = p.second; // neighborhood
    for(size_t i = 0; i < t.size(); i++)
    {
      for(size_t j = i+1; j < t.size(); j++)
      {
        triNbrTri[t[i]].push_back(t[j]);
        triNbrTri[t[j]].push_back(t[i]);
      }
    }
  }
  for(auto & t : triNbrTri)
  {
    sortAndDeduplicate(t);
  }

  auto getNeighbor = [&](int triID) -> const vector<int> & { return triNbrTri[triID]; };

  return getConnectedComponents(triangles.size(), getNeighbor);
}

bool getNonManifoldOEdgeTrisMap(const vector<Vec3i> & triangles, unordered_map<OEdgeKey, vector<int>> & oedgeTris)
{
  int numTriangles = triangles.size();
  oedgeTris.clear();

  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int j = 0; j < 3; j++)
    {
      int v0 = triangles[triID][j];
      int v1 = triangles[triID][(j+1)%3];
      if (v0 == v1 || v0 < 0) // error, TriMesh is not valid!
      {
        return false;
      }
      OEdgeKey edge(v0, v1);
      oedgeTris[edge].push_back(triID);
    }
  }
  return true;
}

vector<OEdgeKey> getExteriorEdges(int numTriangles, const Vec3i * triangles)
{
  map<OEdgeKey, int> count;
  for(int triID = 0; triID < numTriangles; triID++)
  {
    for(int i = 0; i < 3; i++)
    {
      OEdgeKey e(triangles[triID][i], triangles[triID][(i+1)%3]);
      auto it = count.find(e);
      if (it != count.end()) { it->second++; }
      else
      {
        it = count.find(e.getReversedEdgeKey());
        if (it != count.end()) { it->second--; }
        else
        {
          count.emplace(e, 1);
        }
      }
    }
  }
  vector<OEdgeKey> ret;
  for(const auto & p : count)
  {
    if (p.second == 0) continue;
    if (p.second > 0)
    {
      for(int i = 0; i < p.second; i++)
        ret.push_back(p.first);
    }
    else
    {
      OEdgeKey e = p.first.getReversedEdgeKey();
      for(int i = 0; i < (-p.second); i++)
      {
        ret.push_back(e);
      }
    }
  }
  return ret;
}
