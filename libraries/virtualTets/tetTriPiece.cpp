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

#include "tetTriPiece.h"

#include "minivector.h"
#include "basicAlgorithms.h"
#include "containerHelper.h"
#include "valueIndex.h"
#include "profiler.h"
//#include "triMeshOctreeK.h"
#include "geometryQuery.h"
#include "geometryQueryER.h"
#include "predicates.h"
#include <cassert>
using namespace std;

//extern Profiler profiler;
//#define USE_MULTITHREADING

// given two cut vtx on the piece boundary, try to find the tet face the piece boundary edge between them lies on
// we need the tet face because when we do pseudo-normal test on pieces, if the closest site is on the boundary of a piece,
// then we need to tet face to form a locally consistently-orientated shape to compute the correct pseudo-normal
// Note: if the piece boundary edge is on a tet edge, then there are two tet faces we can choose from
// we should choose the one that can form a local consistently-orientated shape around the closest site
static UTriKey getTetFaceFromFeatureOnEdge(TetTriCutFeature f0, TetTriCutFeature f1, const Vec3d & triNormal, const TetShape & tetShape)
{
  if (f0.isTetFace() && f1.isTetFace()) // if both two vtx are on the interior of tet faces
  {
    assert(f0.tetFeature == f1.tetFeature); // assert they are on the same tet face
    return f0.getTetFace(); // then return this face
  }
  // since the two vtx are on the piece boundary, they should be inside the tet
  assert(f0.isInsideTet() == false && f1.isInsideTet() == false);
  if (f0.isTetFace() || f1.isTetFace()) // if either one is on the tet face
  {
    if (f1.isTetFace()) swap(f0, f1); // assuming v0 with feature f0 is on the tet face
    UTriKey face = f0.getTetFace();
    if (f1.isTetEdge()) // if v1 is on a tet edge
    {
      assert(face.hasUEdge(f1.getTetEdge())); // then v1 must be on the edge of the same tet face as v0 is on
    }
    else if (f1.isTetVertex()) // if v1 is on a tet vtx
    {
      assert(face.hasIndex(f1.getTetVertex())); // then v1 must be on the vtx of the same tet face as v0 is on
    }
    return face;
  }

  // in the case that v0 and v1 are on the same tet edge, we have to use piece normal information to
  // find out which tet face they lie on
  auto processTetEdgeCase = [&](UEdgeKey e0)
  {
      auto faces = tetShape.getNeighboringFace(e0); // get the two tet faces sharing this tet edge
      assert(faces[1].isValidTriangle());

      // let's check which face is what we want
      // by checking which tet face lies inside/outside the plane of the triangle this input boundary seg lies
      int v0 = faces[0].getVertexOppositeEdge(e0);
      int v1 = faces[1].getVertexOppositeEdge(e0);

      Vec3d diff0 = tetShape.getPos(v0) - tetShape.getPos(e0[0]);
      Vec3d diff1 = tetShape.getPos(v1) - tetShape.getPos(e0[0]);
      double dot0 = dot(diff0, triNormal), dot1 = dot(diff1, triNormal);
      int minFace = 0; // face local index [0,1] with the smallest dot product
      if (dot0 > dot1) { swap(dot0, dot1); minFace = 1; }
      // now dot0 <= dot1
      // it should be assumed that the triangle of this edge e0 shoule be inside the dihedral angle of faces[0] and faces[1]
      // so dot1 should >= 0 and dot0 should <= 0
      // we don't use exact aritmethic here because unless the tet is nearly degenerate, the angle between
      // the two tet faces should be large enough so that dot0 and dot1 should be reasonably different
      assert(dot1 >= -1e-10);
      assert(dot0 <= 1e-10);
      return faces[minFace];
  };

  if (f0.isTetEdge() && f1.isTetEdge()) // if both are on tet edges
  {
    auto e0 = f0.getTetEdge(), e1 = f1.getTetEdge();
    assert(e0.shareIndex(e1)); // e0 and e1 must have common vtx, otherwise the piece edge <v0,v1> does ont lie on the tet boundary
    if (e0 != e1)
    {
      UTriKey face;
      // if one vtx of e0 is the same as e1[0], then we know the tet face the two tet edges lie on
      if (e0[0] == e1[0] || e0[1] == e1[0]) face = UTriKey(e0[0], e0[1], e1[1]);
      // else if one vtx of e0 is the same as e1[1], then we also know the tet face
      else if (e0[0] == e1[1] || e0[1] == e1[1]) face = UTriKey(e0[0], e0[1], e1[0]);
      assert(face.isValidTriangle());
      return face;
    }
    else // e0 == e1, two vtx lie on the same edge, we have to use piece normal information to find out the tet face
    {
      return processTetEdgeCase(e0);
    }
  }

  if (f0.isTetEdge() || f1.isTetEdge()) // if only one the cut vtx is on a tet edge
  {
    if (f1.isTetEdge()) swap(f0, f1); // make sure v0 is on the tet edge
    // now f0 is tet edge, then v1 is on a tet vtx
    auto e0 = f0.getTetEdge(); // get the tet edge
    int v1 = f1.getTetVertex();
    if (e0[0] == v1 || e0[1] == v1) // if the tet vtx v1 is on is also the vtx of the tet edge e0
    {
      return processTetEdgeCase(e0); // then we have two candidate tet faces
    }
    UTriKey face(e0[0], e0[1], v1); // otherwise, we can locate the tet face
    assert(tetShape.hasFace(face));
    return face;
  }
  // now f0 and f1 are both tet vtx
  // the edge <v0, v1> is identical to a tet edge, again we have two candidate tet faces
  assert(f0.isTetVertex() && f1.isTetVertex());
  int v0 = f0.getTetVertex(), v1 = f1.getTetVertex();
//  if (v0 == v1) {
//    cout << "Error in " << __func__ << ", v0 == v1, f0: " <<
//        f0.triFeature << " " << f0.tetFeature << ", f1: " <<
//        f1.triFeature << " " << f1.tetFeature << endl;
//  }
  assert(v0 != v1);
  return processTetEdgeCase(UEdgeKey(v0, v1));
}

TetTriPiece::TetTriPiece(int tet) : tetID(tet) {}

TetTriPiece::TetTriPiece(int tet, const std::vector<int> & triIDs, const CutTriGroup & group, const TetTriCuttingData & cutting,
  const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape) : tetID(tet)
{
  const vector<Vec3d> & cutPos = cutting.cutVtxPos;
  const vector<TetTriCutFeature> & features = cutting.features;

  if (triIDs.size() == 0) return;
  groupTriID = triIDs;
  int exactOnTetFaceTriID = -1;
  UTriKey tetFaceExactlyOn;
  for(int triID : triIDs)
  {
    for(int vID: group.tri[triID]) vtx.push_back(vID);
    // remove cut triangles completely on the tet face because they don't change the inOutTest result
    bool exactOnTetFace = false;
    for(const auto & p : group.cutTriIDsOnFace)
      if (binary_search(p.second.begin(), p.second.end(), triID))
      {
        exactOnTetFace = true;
        exactOnTetFaceTriID = triID; // we only record one such degenerate-pos triID because only one is needed afterwards
        tetFaceExactlyOn = p.first;
        break;
      }
    if (exactOnTetFace) { continue; }
    gnrCutTri.push_back(group.tri[triID]);
    gnrOriTriID.push_back(group.oriID[triID]);
  }
  sortAndDeduplicate(vtx);
  if (gnrCutTri.size() == 0) // this piece only contains triangles completely on tet faces
  {
    // we will compute the in/out result in this case here
    Vec3d oriTriNormal = oriMeshNormal.triNormal(group.oriID[exactOnTetFaceTriID]);
    double d = dot(oriTriNormal, tetShape.getNormal(tetFaceExactlyOn));
    assert(d > 1 - 1e-6 || d < -1 + 1e-6); // d should be either +1 or -1
    // if this degenerate-pos triangle has the same orientation as the tet face it is on, then this all points inside this tet
    // is considered as IN for this triangle; otherwise, OUT
    inOutWhenAllTriOnTetFaces = (d > 0 ? -1 : +1);
    return;
  }
//  getVerticesInTriangles(tri, vtx);
  assert(hasInvalidTriangles(gnrCutTri) == false);

  // build touchedTetFaces: tet faces which are intersected by triangles in tri
  for(int v : vtx)
  {
    auto faces = features[v].getTouchingTetFaces(tetShape);
    touchedTetFaces.insert(faces.begin(), faces.end());
  }

  try
  {
    gnrTriNbr = TriangleNeighbor(gnrCutTri);
  } catch(int)
  {
    cout << "Error, one triangle group in the tet is not edge-manifold" << endl;
    TriMeshRef triMesh(cutPos, gnrCutTri);
    triMesh.save("Debug.obj");
    throw 1;
  }
}

void TetTriPiece::buildBoundaryData(const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal,
    const TetShape & tetShape, Profiler * profiler)
{
  initializeBoundaryData = true;
  ProfilerExtraSection preprocesTime(profiler, "triCompPreProcess");
  preprocesTime.start();

  const vector<Vec3d> & cutPos = cutting.cutVtxPos;
  const vector<TetTriCutFeature> & features = cutting.features;
  const vector<Vec3ER> & cutPosER = cutting.cutVtxPosER;

  TriMeshRef triMesh(cutPos, gnrCutTri); // the triangle mesh for this piece

  auto triBou = gnrTriNbr.findBoundaryTriangles(gnrCutTri); // get the boundaries of the piece
  bool debug = false;

  // For each boundary edge in tri, compute its missing tet face and angle between them
  for(auto p : triBou)
  {
    OEdgeKey edge = p.second; // get one boundary edge
    int ev0 = edge[0], ev1 = edge[1];
    int triID = p.first; // the triangle ID where this boundary edge belongs to
    Vec3d triNormal = oriMeshNormal.triNormal(gnrOriTriID[triID]); // get the normal of this triangle

    // get the tet face that this tri boundary edge lies
    auto tetFace = getTetFaceFromFeatureOnEdge(features[ev0], features[ev1], triNormal, tetShape);
    assert(tetShape.hasFace(tetFace));

    // compute the angle between the tet face and the tri
    // the angle is [0, M_PI]
    double angle = M_PI - getVectorAngle(triNormal, tetShape.getNormal(tetFace));
    angle = clamp(angle, 0.0, M_PI);

    assert(edgeBouAngle.find(edge) == edgeBouAngle.end());
    edgeBouAngle[edge] = make_pair(angle, tetFace); // store the angle and the tet face
  }

  // next, compute pseudo-normal for piece boundary vertices
  auto loops = gnrTriNbr.findBoundaryLoops(gnrCutTri); // get the boundary edge loops
  for(const auto & loop : loops) // for one loop
  {
    int nl = loop.size(); // numVtx on this loop

    // first, find boundary triangles of this loop
    vector<int> bouTri(nl); // loop edgeID -> boundary triangleID in tri
    for(size_t i = 0; i < loop.size(); i++)
    {
      OEdgeKey edge(loop[i], loop[(i+1)%nl]);
      int tri = gnrTriNbr.getTriangleAtEdge(edge);
      assert(tri >= 0);
      bouTri[i] = tri;
    }

    // then, compute the pseudo-normal at vertex in the loop
    for(size_t i = 0; i < loop.size(); i++)
    {
      int v = loop[i]; // v: vtx ID on loop[i]
      if (debug)
      {
        cout << "compute correct pseudo-normal on cut vtx: " << v << endl;
      }
      // pv: previous vtxID, nv: next vtxID
      int pv = loop[(i+nl-1)%nl], nv = loop[(i+1)%nl];
      // pe: previous edge, ne: next edge
      OEdgeKey pe(pv, v), ne(v, nv);
      // ptriID: the triangle ID of the previous edge, ntriID: the triangle ID of the next edge
      int ptriID = bouTri[(i+nl-1)%nl], ntriID = bouTri[i];
      OTriKey pt(gnrCutTri[ptriID]), nt(gnrCutTri[ntriID]);
      assert(edgeBouAngle.find(pe) != edgeBouAngle.end());
      assert(edgeBouAngle.find(ne) != edgeBouAngle.end());
      // get the previous tet face for pe, and the next tet face for ne
      UTriKey ptetFace = edgeBouAngle[pe].second, ntetFace = edgeBouAngle[ne].second;

      // compute weighted normal on tri, without tet faces
      Vec3d normal(0.0);
      double normalLen = 0.0; // length of the un-normalized normal, used to determine whether we are in a degeneracy case
      // find the fan of triangles around v, starting at the boundary edge of <pv, v>
      auto nbrTris = gnrTriNbr.findTrianglesArroundBoundaryVertex(pv, v, gnrCutTri);
      assert(find(nbrTris.begin(), nbrTris.end(), ptriID) != nbrTris.end());
      assert(find(nbrTris.begin(), nbrTris.end(), ntriID) != nbrTris.end());

      auto getVectorAngleER = [&](Vec3ER & ej, Vec3ER & ek) -> double
      {
        double l2ej = ER_toDouble(len2(ej)), l2ek = ER_toDouble(len2(ek));
        double ejdotek = ER_toDouble(dot(ej, ek));
        assert(l2ej > 0 && l2ek > 0);
        double lejlek = sqrt(l2ej * l2ek);
        assert(lejlek > 0);
        double cosAngle = ejdotek / lejlek;
        cosAngle = clamp(cosAngle, -1.0, 1.0);
        double angle = acos(cosAngle);
        return angle;
      };

      auto getRobustAngle = [&](int vi, int vj, int vk) -> double
      {
        Vec3d pi = triMesh.pos(vi), pj = triMesh.pos(vj), pk = triMesh.pos(vk);
        Vec3d ej = pj - pi, ek = pk - pi;
        if (len2(ej) < 1e-10 || len2(ek) < 1e-10) // triangle too small, not enough precision
        {                                         // switch to exact arithmetic
          Vec3ER ej = cutPosER[vj] - cutPosER[v], ek = cutPosER[vk] - cutPosER[v];
          return getVectorAngleER(ej, ek);
        }

        double cosAngle = dot(ej, ek) / sqrt(len2(ej) * len2(ek));
        cosAngle = clamp(cosAngle, -1.0, 1.0);
        double angle = acos(cosAngle);
        return angle;
      };

      for(int triID : nbrTris) // for each triangle in the fan
      {
        int oriID = gnrOriTriID[triID];
        // we use the triangle normal from the original mesh because it's more robust than computing it
        // directly from the cut triangles
        Vec3d oriNormal = oriMeshNormal.triNormal(oriID);
        Vec3i tri = triMesh.tri(triID);
        int li = tri.getInvertedIndex(v); // li: [0,3), local index of the vtx in the triangle "tri"
        int lj = (li+1)%3, lk = (li+2)%3; // lj, lk: [0,3) local indices of the other two vtx in tri
        int vj = triMesh.tri(triID)[lj], vk = triMesh.tri(triID)[lk];
        double angle = getRobustAngle(v, vj, vk);
//        double angle = triMesh.getTriangleAngleAtVertexRobust(triID, v);

        normal += angle * oriNormal;
        if (debug) // debug code
        {
          cout << "triangle triID " << triID << " is " << triMesh.tri(triID) << endl;
          int i = triMesh.tri(triID).getInvertedIndex(v);
          int j = (i+1)%3, k = (i+2)%3;
          cout << "edge vector: " << triMesh.pos(triID, j) - triMesh.pos(triID, i) << " " << triMesh.pos(triID, k) - triMesh.pos(triID, i) << endl;

          int vj = triMesh.tri(triID)[j];
          int vk = triMesh.tri(triID)[k];
          Vec3ER ej = cutPosER[vj] - cutPosER[v];
          Vec3ER ek = cutPosER[vk] - cutPosER[v];
          cout << "edge vector ER: " << ej << " " << ek << endl;
          cout << "angle of comp tri is " << angle * 180.0 / M_PI << " triID " << triID << " tri normal " << oriNormal << endl;
        }
      }

      if (debug) // debug code
      {
        Vec3d n2 = normal;
        n2.normalize();
        cout << "triangle only normal is " << n2 << endl;
      }

      // try to find the tet faces for this vertex pseudo-normal
      if (ptetFace == ntetFace)
      {
        // this vtx lies inside one tet face
        // easy case:

        Vec3d tetFaceNormal = tetShape.getNormal(ptetFace);

        // the cut triangles on the piece craete an angle "/_" on the tet face
        // now we have to determine which value (angle or 2*M_PI - angle) to use for computing the angle-weighted pseudo-normal
        // we need to check this because locally around the vtx v, faces should be consistently oriented
        // we achieve this by doing following operations:
        double angle = getRobustAngle(v, pv, nv);

        if (debug) // debug code
        {
          cout << "angle between two boundary edges are " << angle * 180.0 / M_PI << endl;
        }

        // to determine whether we should use angle or 2M_PI - angle, we check the cross product of pe and ne
        Vec3ER crossResult = cross(cutPosER[v] - cutPosER[pv], cutPosER[nv] - cutPosER[v]);
        Vec3ER tetFaceNormalER(tetFaceNormal[0], tetFaceNormal[1], tetFaceNormal[2]);
//        Vec3d crossResult = cross(triMesh.pos(v) - triMesh.pos(pv), triMesh.pos(nv) - triMesh.pos(v));
        if (dot(crossResult, tetFaceNormalER) > 0)
        {
          angle = 2*M_PI - angle;
        }
        normal += angle * tetFaceNormal;

        if (debug) // debug
        {
          cout << "the angle is on the same tet face: " << ptetFace << endl;
          cout << "angle " << angle * 180.0 / M_PI << " tetFaceNormal: " << tetFaceNormal << endl;
        }
      }
      else // in this case, the piece boundary vtx is on a tet edge, this means two tet faces will take part in the pseudo-normal computation
      {
        assert(ptetFace.shareUEdge(ntetFace)); // the two tet faces must share a tet edge
        UEdgeKey tetEdge = ptetFace.getSharedUEdge(ntetFace);
        // the piece create a "cut" on the tet edge
        // now again we have a choice of picking which angles to form a consistently oriented local shape
        // The angles are determined by a direction, which is either the direction of the tet edge, or its reverse
        // after we have this direction, the angles between this direction and -pe / ne will be the two angles for the two faces
        int dsign = 0;

        const vector<Vec3ER> & tetPosER = cutting.tetPosER;
#ifdef USE_MULTITHREADING
        assert(cutting.multiThreading);
        assert(cutting.cutVtxPosVecMutex);
        assert(cutting.cutVtxPosVecMutex->size() == cutting.cutVtxPos.size());
        assert(cutting.tetPosVecMutex);
        assert(cutting.tetPosVecMutex->size() == cutting.tetPosER.size());
//          lock_guard<mutex> cutPosLock(*cutting.cutVtxPosMutex);
//          lock_guard<mutex> tetPosLock(*cutting.tetPosMutex);
        UTriKey utriKey(v, nv, pv);
        lock_guard<mutex> lock0(cutting.cutVtxPosVecMutex->at(utriKey[0]));
        lock_guard<mutex> lock1(cutting.cutVtxPosVecMutex->at(utriKey[1]));
        lock_guard<mutex> lock2(cutting.cutVtxPosVecMutex->at(utriKey[2]));
        lock_guard<mutex> lock3(cutting.tetPosVecMutex->at(tetEdge[0]));
        lock_guard<mutex> lock4(cutting.tetPosVecMutex->at(tetEdge[1]));
#endif
        // get triangle normal of triangle <pv, v, nv>
        Vec3ER triNormal = cross(cutPosER[v] - cutPosER[pv], cutPosER[nv] - cutPosER[v]);
        // get tet edge vector <te0, te1>
        Vec3ER tetEdgeVec = tetPosER[tetEdge[1]] - tetPosER[tetEdge[0]]; // te0 -> te1
        // check which tet vtx on tetEdge is under the triangle <pv, v, nv>
        ER d = dot(tetEdgeVec, triNormal);
        dsign = ER_sign(d); // we use exact-arithmetic here for robustness

        if (dsign == 0)
        {
          // degenerate case, this should not happen when
          cout << "Error: dsign == 0 when finding missing tet face normals and angles on a boundary piece vtx on a tet edge" << endl;
          cout << "tetID " << tetID << endl;
          cout << "tri: " << streamRange(gnrCutTri) << endl;
          cout << "v: " << triMesh.pos(v) << endl;
          cout << "trinormal: " << triNormal << endl;
          cout << "tetEdgeVec: " << tetEdgeVec << endl;
        }

        assert(dsign != 0); // assert that this triMesh plane is not parallel with tetEdge, otherwise this triMesh will be degenerate into a line, not a triangle
        if (dsign > 0) { tetEdgeVec = - tetEdgeVec; }

        // now we find the direction, we can compute the two angles, and finally, the correct pseudo-normal
        Vec3ER v2pv = cutPosER[pv] - cutPosER[v];
        Vec3ER v2nv = cutPosER[nv] - cutPosER[v];
        double alpha0 = getVectorAngleER(v2pv, tetEdgeVec);
        normal += alpha0 * tetShape.getNormal(ptetFace);
        double alpha1 = getVectorAngleER(v2nv, tetEdgeVec);
        normal += alpha1 * tetShape.getNormal(ntetFace);

        if (debug) // debug code
        {
          cout << "the angle is on two tet faces: " << ptetFace << " " << ntetFace << endl;
          cout << "angle " << alpha0 * 180.0 / M_PI << " tetFaceNormal: " << tetShape.getNormal(ptetFace) << endl;
          cout << "angle " << alpha1 * 180.0 / M_PI << " tetFaceNormal: " << tetShape.getNormal(ntetFace) << endl;
        }
      }

      normalLen = len2(normal);
      if (normalLen > 0) normal.normalize();
      // we store the len2 of normal because we need it to measure the degeneracy of this pseudo-normal
      vtxBouNormal.insert(make_pair(v, make_pair(normalLen, normal)));

      if (debug) // debug code
      {
        cout << "Finally, normalLen2: " << normalLen << " normal " << normal << endl;
      }
    } // end i in loop
  } // end loop in loops
}

void TetTriPiece::tryBuildingBoundaryData(const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal,
    const TetShape & tetShape, Profiler * profiler)
{
  if (initializeBoundaryData == false)
    buildBoundaryData(cutting, oriMeshNormal, tetShape, profiler);
}

//void TetTriComp::buildTree(const vector<Vec3ER> & cutPosK) {
//  if (tri.size() > 0) {
//    octree.build(cutPosK, tri, 5, 10);
//  }
//  treeBuilt = true;
//}

int TetTriPiece::getClosestTriangle(const TetTriCuttingData & cutting, const Vec3d & queryPoint,
    const Vec3ER & queryPointER, int & feature, ER & dist2)
{
  const auto & cutPosER = cutting.cutVtxPosER;
  // try using brute-force:

  auto getDist = [&](int triID, int & minFeature)
  {
    Vec3i t = gnrCutTri[triID];
#ifdef USE_MULTITHREADING
    assert(cutting.multiThreading);
    assert(cutting.cutVtxPosVecMutex);
    assert(cutting.cutVtxPosVecMutex->size() == cutting.cutVtxPos.size());
    UTriKey utri(t);
    lock_guard<mutex> lock0(cutting.cutVtxPosVecMutex->at(utri[0]));
    lock_guard<mutex> lock1(cutting.cutVtxPosVecMutex->at(utri[1]));
    lock_guard<mutex> lock2(cutting.cutVtxPosVecMutex->at(utri[2]));
#endif
    return squaredDistanceToTriangle(queryPointER, cutPosER[t[0]], cutPosER[t[1]], cutPosER[t[2]], minFeature);
  };

  // we use brute-force to search for the closest triangle and site
  // we also tried building an octree tree, but it was slower
  int minTriID = 0, minFeature = 0;
  dist2 = getDist(0, minFeature);
  for(size_t triID = 1; triID < gnrCutTri.size(); triID++)
  {
    ER d2 = getDist(triID, feature);
    if (d2 < dist2)
    {
      dist2 = d2;
      minTriID = triID;
      minFeature = feature;
    }
  }
  feature = minFeature;
  return minTriID;

//  if (treeBuilt == false) { buildTree(cutPosK); }
//  return octree.getClosestTriangle(cutPosK, tri, queryPointK, feature);
}

int TetTriPiece::inOutTest(const Vec3d & queryPoint, const Vec3ER & queryPointER, const TetTriCuttingData & cuttingData,
    const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape, Profiler * profiler)
{
  bool debug = false;

  if (debug)
  {
    cout << "vtx " << streamRange(vtx) << endl;
    cout << "#gnrCutTri: " << gnrCutTri.size() << " " << inOutWhenAllTriOnTetFaces << endl;
  }

  if (vtx.size() == 0) { return 1; } // interior tet, all query points are considered outside
  if (gnrCutTri.size() == 0) { return inOutWhenAllTriOnTetFaces; } // only triangles completely on tet faces

  const vector<Vec3d> & cutPositions = cuttingData.cutVtxPos;
  const vector<Vec3ER> & cutPosER = cuttingData.cutVtxPosER;

  ProfilerExtraSection testTime(profiler, "inOutTest");
  testTime.start();

  TriMeshRef mesh(cutPositions, gnrCutTri); // triangle mesh of this piece
  int feature = -1;
  ER closestDist;
  // the closest triangle on this piece to the query point
  int cTriID = getClosestTriangle(cuttingData, queryPoint, queryPointER, feature, closestDist);
  assert(cTriID >= 0 && feature >= 0);
  if (debug)
  {
    cout << "closest triID is " << cTriID << endl;
  }
  if (ER_toDouble(closestDist) == 0.0)
  {
    return 1; // mesh are touching, we treate it as separate
  }

  Vec3i closestTri = gnrCutTri[cTriID];

  Vec3d projedPt = getClosestPointToTriangleWithNormalAndFeature(queryPoint, mesh.pos(cTriID, 0), mesh.pos(cTriID, 1), mesh.pos(cTriID, 2),
      oriMeshNormal.triNormal(gnrOriTriID[cTriID]), feature);
  Vec3d diff = queryPoint - projedPt;
  Vec3d normal(0.0);


  if (debug)
  {
    cout << "query pt is " << queryPoint << endl;
    cout << "projected pt is " << projedPt << " f " << feature << endl;
  }

  auto computeFinalDotProduct = [&]() -> int
  {
    if (len2(diff) < 1e-10)
    {
//      cout << "warning: diff is too small in inout test: " << diff << endl;
      // to avoid errors, we use exact-arithmetic to compute the in/out result
      Vec3ER projectedPtER = getClosestPointToTriangleWithFeature(queryPointER,
          cutPosER[closestTri[0]], cutPosER[closestTri[1]], cutPosER[closestTri[2]], feature);
      Vec3ER diffER = queryPointER - projectedPtER;
      Vec3ER normalER(normal[0], normal[1], normal[2]);
      int sign = ER_sign(dot(normalER, diffER));
      assert(sign != 0);
      assert(sign == 1 || sign == -1);
      return sign;
    }

    if (dot(normal, diff) > 0) { return 1; }
    return -1;
  };

  if (feature < 3) // closest feature is a vtx
  {
    int vtxID = closestTri[feature];
    if (vtxID < oriMeshNormal.numVertices()) // it's a vtx from original input mesh
    {
      normal = oriMeshNormal.vtxNormal(vtxID);
      return computeFinalDotProduct();
    }
    // else, it's a cut vtx, which is on the piece boundary
    if (debug)
      cout << "projected pt is on boundary vtx: " << vtxID << ", has boundary normals: " << endl;

    // try building boundary data related to pseudo-normals on the boundary, if not yet
    tryBuildingBoundaryData(cuttingData, oriMeshNormal, tetShape, profiler);

    // get the vtx boundary pseudo-normal
    // vtxBouNormal is a multimap, because one boundary vtx can have more than one pseudo-normal
    // this may seem wrong, but actually since we allow the piece to be only edge-manifold, two
    // manifold triangle groups in the piece can be connected via only one sharing vtx. In this case,
    // the piece triangle mesh is edge-manifold, but not manifold.
    // but we can still run our method on this case.
    // each returned value in this multimap represents one local edge-manifold on this vtx: vtxID
    // we check in-out by: query point is out only if all local edge-manifolds return out
    auto iterPair = vtxBouNormal.equal_range(vtxID);
    assert(iterPair.first != vtxBouNormal.end());

    for(auto iter = iterPair.first; iter != iterPair.second; iter++)
    {
      const auto & p = iter->second;
      if (debug)
        cout << "scaled normal len2 " << p.first << " normal " << p.second << endl;
      if (p.first < 1e-6)
      {           // the scaled normal is too short, close to a degenerte feature
        continue; // we assume it's outside, so we do nothing in this loop
      }
      normal = p.second;
      if (computeFinalDotProduct() > 0) { continue; } // it's out, do nothing
      return -1; // it's inside one local edge-manifold, so it must return inside
    }
    // now no local edge-manifolds hold this query point. It's outside
    return 1;
  }
  else if (feature < 6) // closest feature is an edge
  {
    int nbrTriID = gnrTriNbr.getTriangleNeighbors(cTriID)[feature-3]; // nbr local triangleID sharing this edge
    int oriID = gnrOriTriID[cTriID];
    int nbrOriID = (nbrTriID >= 0 ? gnrOriTriID[nbrTriID] : oriID);
    // if cTriID and its nbr, nbrTriID are from the same input triangle,
    // or if cTriID has no nbr sharing this edge
    // then the normal is assigned to the normal of the original input triangle oriID
    if (oriID == nbrOriID) { normal = oriMeshNormal.triNormal(oriID); }
    else // else, TriID and its nbr, nbrTriID are from the two different input triangles,
    {    // then normal is assigned as the average of the two normals
      normal = oriMeshNormal.triNormal(oriID) + oriMeshNormal.triNormal(nbrOriID);
      normal.normalize();
      assert(normal.hasNaN() == false);
    }
    if (debug)
      cout << "clost pt is on an edge" << endl;
    if (nbrTriID >= 0) // this edge is not a boundary edge in this comp
    {
      if (debug)
        cout << "its not a boundary edge, the normal is " << normal << endl;
      return computeFinalDotProduct();
    }

    // else, this edge is on the boundary
    tryBuildingBoundaryData(cuttingData, oriMeshNormal, tetShape, profiler);

    // we can still use this normal but we have to know the dihedral angle as well
    Vec3i ctri = gnrCutTri[cTriID];
    OEdgeKey cedge(ctri[feature-3], ctri[(feature-2)%3]);
    auto it = edgeBouAngle.find(cedge);
    assert(it != edgeBouAngle.end());
    double angle = it->second.first; // the angle between the triangle cTriID and the tet face that should also contribute to the pseudo-normal computation
    UTriKey tetFace = it->second.second; // the tet face that should also contribute to the pseudo-normal computation
    if (debug)
      cout << "the edge is on a boundary, tri normal " << normal << ", cedge: " << cedge << endl;
    if (angle < 1e-2) // the angle between cTriID and tetFace is too small, we think the query point can only be outside the piece
    {
      if (debug)
      {
        cout << "the angle between the tet face & tri face " << rad2deg(angle) << " is too small, we think the query pt is outside" << endl;
        cout << "the tet face involved is " << tetFace << " normal " << tetShape.getNormal(tetFace) << endl;
      }
      return 1; // if the dihedral angle is too small, we think it cannot be inside (theoretically it cannot be inside if angle < 90)
    }
    if (angle > M_PI - 1e-2) // again, if the angle is too large, then we think the query point can only be inside the piece
    {
      if (debug)
      {
        cout << "the angle between the tet face & tri face " << rad2deg(angle) << " is too large, we think the query pt is inside" << endl;
        cout << "the tet face involved is " << tetFace << " normal " << tetShape.getNormal(tetFace) << endl;
      }
      return -1;
    }
    int oriTID = gnrOriTriID[cTriID];
    assert(oriTID >= 0 && oriTID < oriMeshNormal.numTriangles());
    // then the final normal is the averge of the triangle normal and the tet normal
    normal = oriMeshNormal.triNormal(oriTID) + tetShape.getNormal(tetFace);
    normal.normalize();
    assert(normal.hasNaN() == false);
    if (debug)
      cout << "the normal is " << normal << endl;
    return computeFinalDotProduct();
  }
  // else, closest feature is triangle interior, simplest case
  int oriTID = gnrOriTriID[cTriID];
  assert(oriTID >= 0 && oriTID < oriMeshNormal.numTriangles());
  normal = oriMeshNormal.triNormal(oriTID);
  return computeFinalDotProduct();
}

// merge the two pieces
TetTriPiece mergePiece(const TetTriPiece & comp0, const TetTriPiece & comp1, int tet, const TetTriMeshCutting::CutTriGroup & group,
    const TetTriCuttingData & cutting, const TriMeshPseudoNormal & oriMeshNormal, const TetShape & tetShape)
{
  vector<int> groupTriID = comp0.groupTriID;
  groupTriID.insert(groupTriID.end(), comp1.groupTriID.begin(), comp1.groupTriID.end());
  sort(groupTriID.begin(), groupTriID.end());
  assert(unique(groupTriID.begin(), groupTriID.end()) == groupTriID.end());

  return TetTriPiece(tet, groupTriID, group, cutting, oriMeshNormal, tetShape);
}

TriMeshRef TetTriPiece::mesh(const TetTriCuttingData & cut) const
{
  assert(gnrCutTri.size() > 0);
  return { cut.cutVtxPos, gnrCutTri };
}
