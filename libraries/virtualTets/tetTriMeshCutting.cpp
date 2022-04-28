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

#include "tetTriMeshCutting.h"
#include "tetKey.h"
#include "containerHelper.h"
#include "basicAlgorithms.h"
#include "exactOctree.h"
#include "predicates.h"
#include "profiler.h"
#include "performanceCounter.h"
#include "planeER.h"
#include "meshIntersection.h"
#include "plane.h"
#include <cassert>
#include <array>
#include <fstream>
#include <unordered_map>
#include <mutex>
#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif

using namespace std;

//static int debugInt = 0;

///////////////////////////////////////////////////////////////
//                    TetTriMeshCutting
///////////////////////////////////////////////////////////////


static Profiler profiler;
static Profiler cutProfiler;
//static StopWatch wnTimer;

// namespace {
// class ExitRegister {
// public:
//   ~ExitRegister() {
//     cout << "Time for tet tri cut: " << endl;
//     cout << profiler.toString() << endl;
//     cout << "Time inside the cut:" << endl;
//     cout << cutProfiler.toString() << endl;
//   }
// };
// static ExitRegister exitRegister_;
// }

TetTriMeshCutting::TetTriMeshCutting(const TetMeshRef & tetMesh, TriMeshRef triMesh)
: tetMesh(tetMesh), triMesh(triMesh)
{}

void TetTriMeshCutting::computeIntersectingTriTets()
{
  profiler.startTimer("embedding");

  tri2tet.clear();
  tri2tet.resize(triMesh.numTriangles());

  // tet2tri: tetID -> intersecting triIDs
  tet2tri = computeTrianglesIntersectingEachTetExact(tetMesh, triMesh);

  // build tri2tet from tet2tri
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    for(int tri : tet2tri[tetID]) {  tri2tet[tri].push_back(tetID); }
  }

  vector<int> outsideTriangles; // triangleIDs not intersected by any tet
  for(int i = 0; i < triMesh.numTriangles(); i++)
  {
    if (tri2tet[i].size() == 0) { outsideTriangles.push_back(i); }
  }

  if (outsideTriangles.size() > 0)
  {
    cout << "warning: " << outsideTriangles.size() << " triangle(s) not embedded in the tet mesh" << endl;
//    throw 1;
  }

  profiler.stopLastTimer();
}

void TetTriMeshCutting::computeCutTriangles(const std::vector<Vec3ER> * triVtxPosER)
{
  assert(tet2tri.size() > 0);

  profiler.startTimer("prepare to cut");

  ER ZERO = 0.0;

  // build tetVtxPosER
  tetVtxPosER.clear();
  cutTriVtxPosER.clear(); // stores vtx positions on the cut triangle mesh
  for(int i = 0; i < tetMesh.numVertices(); i++)
  {
    Vec3d v = tetMesh.pos(i);
    tetVtxPosER.emplace_back(v[0], v[1], v[2]);
  }

  // the cut tri mesh inherits the vtx from the input tri mesh
  // those vtx have the same vtxIDs in cut tri mesh as the input tri mesh
  if (triVtxPosER) // if triVtxPosER has been provided
  {
    assert(sizei(*triVtxPosER) == triMesh.numVertices());
    cutTriVtxPosER = (*triVtxPosER); // initialize the inherited vertices from the input tri mesh
  }
  else // compute triVtxPosER
  {
    for(int i = 0; i < triMesh.numVertices(); i++)
    {
      Vec3d v = triMesh.pos(i);
      cutTriVtxPosER.emplace_back(v[0], v[1], v[2]);
    }
  }

  // The basic algorithm works as follows: we go through each triangle/tet pair and compute vertices and triangles in the
  // cut triangle mesh, and store them in a global buffer.
  // We denote those new vertices as cutTriVtx.
  // We parallelize over each tet for multi-threading.
  // Because we are using CGAL's lazy evaluation exact kernel, we need to be careful as exact numbers are dependent on each other.
  // One computed cutTriVtx exact position is dependent on the triangle and tet pair where the intersection is from.
  // So this cutTriVtx posER is dependent on three triangle vtx posER and four tet vtx posER.
  // Therefore, we lock seven mutexes for one triangle/tet pair and after computing the cutTriVtxes, we lock one more mutex
  // to access the global buffer to store the cutTriVtxes.

  // stores new cutTriVtx on the cut triangle mesh
  // we use a list instead of a vector here because we want to avoid thread contest
  // If it is a vector, then vector will reallocate its memory when its capacity is not enough.
  // At this time, all the previous data stored in newCutTriVtxPosER will be touched.
  // However, since we parallelize over triangle/tet pairs, the previous data is still dependent on tet/tri values touched
  // by another thread. This might be dangerous.
  // So to be safe, we use a list.
  list<Vec3ER> newCutTriVtxPosER;

  vector<UTriKey> allTetFaces; // find all tet faces, which are used as cutting plane
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    OTetKey tetkey(tetMesh.tet(tetID));
    for(int j = 0; j < 4; j++)
    {
      UTriKey utri = tetkey.uFaceKey(j);
      allTetFaces.push_back(utri);
    }
  }
  sortAndDeduplicate(allTetFaces);

  // build the mapping from tet face key to tetFaceID
  // this mapping is useful because we want to save exact arithmetic operations
  // on converting tet faces to planes for cutting
  // TODO: since allTetFaces is already sorted, we can use a std::lower_bound to find the mapping from face key to faceID,
  // thus we don't need this map<UTriKey, int>
  map<UTriKey, int> tetFaceID;
  for(size_t i = 0; i < allTetFaces.size(); i++) { tetFaceID[allTetFaces[i]] = i; }

  // build planes from tet faces, those planes will be used in the Sutherland-Hodgman algorithm
  vector<FastPlane> allTetPlanes;
  for(size_t i = 0; i < allTetFaces.size(); i++)
  {
    UTriKey t = allTetFaces[i];
    allTetPlanes.emplace_back(tetMesh.pos(t[0]), tetMesh.pos(t[1]), tetMesh.pos(t[2]));
  }

  // initialize the double-precision buffer for cut triMesh vtx pos
  cutTriVtxPos = triMesh.exportPositions();
  cutVtxFeatures.clear();
  cutVtxFeatures.resize(cutTriVtxPos.size());

  typedef std::unordered_map<TetTriCutFeature, int> FeatureIDMap; // feature -> vtxID
  FeatureIDMap featureIDMap;

  const TetTriCutFeature featurePairNull;

  profiler.stopLastTimer();

  profiler.startTimer("cut");
  // Now we will loop over each tet, cut the triangles intersecting this tet with the four tet faces
  // We use parallelize this loop. The result of the cut, the cut vtx, feature and cut triangles are then
  // stored in cutTriVtxPos, cutVtxFeatures and cutTrisInTet

  cutTrisInTet.clear();
  cutTrisInTet.resize(tetMesh.numTets());
#ifdef USE_TBB
  mutex cutVtxMutex; // used to lock the access to the buffer storing the cut result
  vector<mutex> tetVtxMutex(tetMesh.numVertices()); // used to lock the access to each exact tet vtx pos
  vector<mutex> triVtxMutex(triMesh.numVertices()); // used to lock the access to each exact input tri vtx pos
  tbb::parallel_for(tbb::blocked_range<int>(0, tetMesh.numTets()), [&](const tbb::blocked_range<int> & rng)
  {
    for (int tetID = rng.begin(); tetID != rng.end(); ++tetID)
#else
    for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
#endif
    {
      if (tet2tri[tetID].size() == 0) { continue; } // no tris in this tet

      const auto & tris = tet2tri[tetID]; // triangles intersecting this tet

      OTetKey tetkey(tetMesh.tet(tetID));
      UTetKey utetkey(tetMesh.tet(tetID));
      CutTriGroup triGroupInTet;

      const int UNINIT = -100;

      for(int triID : tris) // for each triangle, cut it with the four tet faces
      {
        Vec3i t = triMesh.tri(triID);
        // struct CutVertex: temporary cut vertex data after cutting by one plane
        // There are four tet faces, so we do four iterations of the Sutherland-Hodgman algorithm
        // CutVertex stores the result after one iteration
        struct CutVertex
        {
          Vec3d pos;    // vtx pos
          Vec3ER posER; // vtx pos in exact arithmetic
          Vec3ER bary;  // barycentric coords of this vtx with respect to the three triangle vertices
          int ID = -1;  // cutVtxID
          UTriKey triFeature; // the triangle interior / edge / vertex this vertex lies in
          UTetKey tetFeature; // the tet interior / face / edge / vertex this vertex lies in
          int inOut[4] = { UNINIT, UNINIT, UNINIT, UNINIT }; // local face ID [0,3] -> whether this vtx is outside (+1) or inside (-1) the face
        };


#ifdef USE_TBB
        // guard to all the traingle vtx
        // this will not cause dead locks, because we acquire the lock of each tri vtx according to
        // the order of the tri vtx
        // UTriKey store its vtx indices strictly sorted
        // So utrikey[0] < utrikey[1] < utrikey[2]
        // Since it follows an order of tri vtx indices, deadlock will not occur.
        UTriKey utrikey(t);
        lock_guard<mutex> triVtxGd0(triVtxMutex[utrikey[0]]);
        lock_guard<mutex> triVtxGd1(triVtxMutex[utrikey[1]]);
        lock_guard<mutex> triVtxGd2(triVtxMutex[utrikey[2]]);

        // guard to all the tet faces
        // this will not cause dead locks, because we acquire the lock of each tet vtx according to
        // the order of the tet vtx
        // utetkey store its vtx indices strictly sorted
        // So utetkey[0] < utetkey[1] < utetkey[2] < utetkey[3]
        // Since it follows an order of tet vtx indices, deadlock will not occur.
        lock_guard<mutex> guard0(tetVtxMutex[utetkey[0]]);
        lock_guard<mutex> guard1(tetVtxMutex[utetkey[1]]);
        lock_guard<mutex> guard2(tetVtxMutex[utetkey[2]]);
        lock_guard<mutex> guard3(tetVtxMutex[utetkey[3]]);
#endif

        vector<CutVertex> cutVertices(3); // initialize the triangle before any cut


        for(int i = 0; i < 3; i++)
        {
          int vtxID = t[i];
          cutVertices[i].ID = vtxID;
          cutVertices[i].pos = triMesh.pos(vtxID);
          cutVertices[i].posER = cutTriVtxPosER[vtxID];
          cutVertices[i].triFeature = UTriKey(-1, -1, vtxID);
          Vec3i b(0);
          b[i] = 1.0;
          cutVertices[i].bary = Vec3ER(b[0], b[1], b[2]);
        }

        for(int fID = 0; fID < 4; fID++)
        {
          OTriKey otri = tetkey.oFaceKey(fID);
          UTriKey utri = tetkey.uFaceKey(fID);
          int planeID = tetFaceID[utri];
          assert(planeID >= 0);
          // planeOutward: 1: if vtx indices in otri is the same as utri, which means that the plane stored
          //                  in allTetPlanes(ER) has the orientation of facing outward the tet
          //              -1: otherwise
          int planeOutward = int(otri[1] == utri[1]) * 2 - 1;
          assert(planeOutward == 1 || planeOutward == -1);

  //        const auto & plane = allTetPlanesER[planeID];

          FastPlaneER plane(tetVtxPosER[utri[0]], tetVtxPosER[utri[1]], tetVtxPosER[utri[2]]);

          vector<CutVertex> outputCutVertices;

          for(size_t i = 0; i < cutVertices.size(); i++)
          {
            cutVertices[i].inOut[fID] = planeOutward * plane.outside(cutVertices[i].posER);
          }

  //        cutProfiler.startTimer("cutVtx");
          // loop over each edge of the current polygon stored in cutVertices
          assert(cutVertices.size() > 0);
          for(size_t cutVtxID = 0; cutVtxID < cutVertices.size(); cutVtxID++)
          {
            // get edge <v0, v1>
            const CutVertex & v0 = cutVertices[cutVtxID];
            const CutVertex & v1 = cutVertices[(cutVtxID+1)%cutVertices.size()];
            const int & l0 = v0.inOut[fID];
            const int & l1 = v1.inOut[fID];
            assert(l0 != UNINIT && l1 != UNINIT);
            // core Sutherland–Hodgman algorithm code: cut edge <v0, v1> with the plane

            // v0 is on the plane, or if both vtx inside or on clippling plane
            if (l0 == 0 || (l0 < 0 && l1 <= 0))
            {
              outputCutVertices.push_back(v0); // we retain v0 in the cut polygon
              continue;
            }

            // both vtx outside clippling plane, or if v0 outside and v1 on the plane
            else if (l0 > 0 && l1 >= 0)
            {
              continue;
            }

            // now, either: v0 inside and v1 outside
            //          or: v0 outside and v1 inside
            // in either case, we will generate a new cut vertex different than v0 or v1
            if (l0 < 0) // if v0 inside
            {
              outputCutVertices.push_back(v0); // we retain v0 in the cut polygon
            }

            // generate the new vertex
            CutVertex vi;

            // compute weight for where on the edge this new vtx lies
            ER a0 = plane.scaledDistance(v0.posER); // we use scaledDistance to avoid normalization operations for exact arithmetic
            ER a1 = plane.scaledDistance(v1.posER);
            ER asum = a0 + a1;

            // we compute both exact and double-precision values because we need both exact and double-precision positions
            double a0d = allTetPlanes[planeID].scaledDistance(v0.pos);
            double a1d = allTetPlanes[planeID].scaledDistance(v1.pos);

            vi.pos = (a1d * v0.pos + a0d * v1.pos) / (a0d + a1d);
            vi.posER = (a1 * v0.posER + a0 * v1.posER) / asum;
            vi.bary = (a1 * v0.bary + a0 * v1.bary) / asum;
            for(int k = 0; k < fID; k++) // we add back the missing inOut data for the new cut vtx
            {
              // integrity check: for all the previous faces, v0 and v1 should be considered inside or exactly on the face
              assert(v0.inOut[k] <= 0 && v1.inOut[k] <= 0);
              vi.inOut[k] = ((v0.inOut[k] == 0 && v1.inOut[k] == 0) ? 0 : -1); // assign vi.inOut
            }
            vi.inOut[fID] = 0;
  //          for(int k = fID+1; k < 4; k++)
  //            vi.inOut[k] = ((v0.inOut[k] == 0 && v1.inOut[k] == 0) ? 0 : -1);

            outputCutVertices.push_back(vi); // finally, store the new vtx vi
          } // end for(size_t cutVtxID = 0; cutVtxID < cutVertices.size(); cutVtxID++)
  //        cutProfiler.stopLastTimer();
          if (outputCutVertices.size() < 3) // if the cut resulting a degenerate triangle,
          {                                 // then we should skip this cut
            cutVertices.clear();            // set cutVertices to be empty so that we can get out of the outer loop
            break;
          }
          cutVertices.swap(outputCutVertices); // store the cut polygon to cutVertices for the next iteration
        } // end fID loop on tet faces

        if (cutVertices.size() < 3) continue; // we don't store degenerate triangles

        // now we build vtx cut feature
        // TODO: CutVertex::bary's only usage is to find the correct triFeature.
        // But we can achieve this without exact arithmetic by tracking the location of each cut vtx during cutting
        for(size_t i = 0; i < cutVertices.size(); i++)
        {
          CutVertex & v = cutVertices[i];
          bool zeroBary[3];
          for(int j = 0; j < 3; j++)
            zeroBary[j] = (v.bary[j] == ZERO);

          if (zeroBary[0] == true && zeroBary[1] == true) // v is a triangle vertex t[2]
          {
            assert(v.ID == t[2]);
            v.triFeature = UTriKey(-1,-1, t[2]);
          }
          else if (zeroBary[0] == true && zeroBary[2] == true)  // v is a triangle vertex t[1]
          {
            assert(v.ID == t[1]);
            v.triFeature = UTriKey(-1,-1, t[1]);
          }
          else if (zeroBary[1] == true && zeroBary[2] == true)  // v is a triangle vertex t[0]
          {
            assert(v.ID == t[0]);
            v.triFeature = UTriKey(-1,-1, t[0]);
          }
          else if (zeroBary[0] == true) // v is on triangle edge <t[1], t[2]>
          {
            v.triFeature = UTriKey(-1, t[1], t[2]);
          }
          else if (zeroBary[1] == true) // v is on triangle edge <t[0], t[2]>
          {
            v.triFeature = UTriKey(-1, t[0], t[2]);
          }
          else if (zeroBary[2] == true) // v is on triangle edge <t[0], t[1]>
          {
            v.triFeature = UTriKey(-1, t[0], t[1]);
          }
          else // v is on the interior of the triangle
          {
            v.triFeature = UTriKey(t);
          }

          for(int j = 0; j < 4; j++)
          {
            assert(v.inOut[j] <= 0);
          }
          int inOutSum = v.inOut[0] + v.inOut[1] + v.inOut[2] + v.inOut[3];
          assert(inOutSum <= -1);
          if (inOutSum == -1) // it's a tet vtx
          {
            int fID = 0;
            for(; fID < 4; fID++) { if (v.inOut[fID] < 0) break; }
            assert(fID < 4);
            v.tetFeature = UTetKey(-1,-1,-1,tetkey[fID]);
          }
          else if (inOutSum == -2)  // it's a tet edge
          {
            UTriKey f[2]; // first get the two tet faces sharing this tet edge
            for(int fID = 0; fID < 4; fID++)
            {
              if (v.inOut[fID] == 0)
              {
                if (f[0][0] < 0) { f[0] = tetkey.uFaceKey(fID); }
                else { assert(f[1][0] < 0); f[1] = tetkey.uFaceKey(fID); }
              }
            }
            assert(f[1][0] >= 0);
            UEdgeKey e = f[0].getSharedUEdge(f[1]);
            assert(e[0] >= 0);
            v.tetFeature = UTetKey(-1, -1, e[0], e[1]);
          }
          else if (inOutSum == -3)  // it's a tet face
          {
            int fID = 0;
            for(; fID < 4; fID++)
            {
              if (v.inOut[fID] == 0)  { break; }
            }
            assert(fID < 4);
            UTriKey f = tetkey.uFaceKey(fID);
            v.tetFeature = UTetKey(-1, f[0], f[1], f[2]);
          }
          else  // it's inside the tet
          {
            v.tetFeature = tetkey.uTetKey();
          }
        } // end for(size_t i = 0; i < cutVertices.size(); i++) to build vtx cut feature

        bool isOnTetFace = false; // whether this cut polygon is on a tet face
        UTriKey tetFaceThisTriIsOn;
        for(int fID = 0; fID < 4; fID++)
        {
          isOnTetFace = all_of(cutVertices.begin(), cutVertices.end(), [&](const CutVertex & v) { return v.inOut[fID] == 0; });
          if (isOnTetFace) { tetFaceThisTriIsOn = tetkey.uFaceKey(fID); break; }
        }
        // this triangle cannot be exactly on more than two tet faces, unless this triangle is degenerate and lies on one tet edge
        // but we require input triangles are free of degeneracy

        {
#ifdef USE_TBB
            lock_guard<mutex> cutVtxlock(cutVtxMutex);
#endif
          for(size_t i = 0; i < cutVertices.size(); i++)
          {
            CutVertex & v = cutVertices[i];
            TetTriCutFeature fp(v.triFeature, v.tetFeature);

            // now we have triFeature and tetFeature ready, we can store each vtx
            if (v.ID >= 0)  // they are input triangle mesh vtx
            {
              if (cutVtxFeatures[v.ID] != featurePairNull) { assert(cutVtxFeatures[v.ID] == fp); }
              else { cutVtxFeatures[v.ID] = fp; }
              continue;
            }

            // it's a new cut vtx, we search for its feature
            if (featureIDMap.find(fp) == featureIDMap.end())
            {
              // new feature, we will create a new cut vtx ID for it
              int newID = cutTriVtxPos.size();
              featureIDMap[fp] = newID;
              cutTriVtxPos.push_back(v.pos);
              newCutTriVtxPosER.push_back(v.posER);
              cutVtxFeatures.push_back(fp);
              assert(featureIDMap.size() + triMesh.numVertices() == cutTriVtxPos.size());
              v.ID = newID;
            }
            else
            {
              v.ID = featureIDMap[fp];
            }
          } // end loop on cutVertices
        } // end threading-safe block

        // now all cut vertices have been processed, we can triangulate the cut polygon into triangles
        vector<int> IDList(cutVertices.size());
        for(size_t i = 0; i < cutVertices.size(); i++)
          IDList[i] = cutVertices[i].ID;

        vector<Vec3i> triangleList = triangulatePolygon(IDList);
        assert(hasInvalidTriangles(triangleList) == false);

        int triIDStart = triGroupInTet.tri.size();
        triGroupInTet.tri.insert(triGroupInTet.tri.end(), triangleList.begin(), triangleList.end());
        for(size_t i = 0; i < triangleList.size(); i++)
          triGroupInTet.oriID.push_back(triID);
        if (isOnTetFace)
        {
          for(size_t i = 0; i < triangleList.size(); i++)
            triGroupInTet.cutTriIDsOnFace[tetFaceThisTriIsOn].push_back(triIDStart+i);
  //        cout << "has tri on tet face: " << triID << " " << tetID << endl;
        }
      } // end triID

      cutTrisInTet[tetID] = move(triGroupInTet);
    } // end each tetID
#ifdef USE_TBB
  }, tbb::auto_partitioner()); //end for locations
#endif

  cutTriVtxPosER.insert(cutTriVtxPosER.end(), newCutTriVtxPosER.begin(), newCutTriVtxPosER.end());
  profiler.stopLastTimer();
}

void TetTriMeshCutting::saveCutTriMesh(string filename) const
{
  vector<Vec3i> allTris;
  for(const auto & triGroup: cutTrisInTet)
  {
    const auto & tri = triGroup.tri;
    allTris.insert(allTris.end(), tri.begin(), tri.end());
  }
  TriMeshGeo mesh(cutTriVtxPos, move(allTris));
  mesh.save(filename);
}

TriMeshGeo TetTriMeshCutting::exportCutTriMesh() const
{
  vector<Vec3i> allTris;
  for(const auto & triGroup: cutTrisInTet)
  {
    const auto & tri = triGroup.tri;
    allTris.insert(allTris.end(), tri.begin(), tri.end());
  }
  return TriMeshGeo(cutTriVtxPos, allTris);
}

