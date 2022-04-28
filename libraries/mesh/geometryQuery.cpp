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

#include "geometryQuery.h"
#include "basicAlgorithms.h"
#include <cassert>
using namespace std;

double getTriangleAngle(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2)
{
  Vec3d e1 = v1 - v0, e2 = v2 - v0;
  double cosAngle = dot(e1, e2) / (len(e1) * len(e2));
  cosAngle = clamp(cosAngle, -1.0, 1.0);
  return acos(cosAngle);
}

double getVectorAngle(const Vec3d & vec1, const Vec3d vec2)
{
  double cosAngle = dot(vec1, vec2) / sqrt(len2(vec1) * len2(vec2));
  cosAngle = clamp(cosAngle, -1.0, 1.0);
  return acos(cosAngle);
}

double getTriangleAngleRobust(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2)
{
  Vec3d e1 = v1 - v0;
  Vec3d e2 = v2 - v0;
  double l2e1 = len2(e1), l2e2 = len2(e2);
  double alpha = 0.0;
  if (l2e1 > 0 && l2e2 > 0)
  {
    double cosAlpha = dot(e1,e2) / sqrt(l2e1 * l2e2);
    cosAlpha = clamp(cosAlpha, -1.0, 1.0);
    alpha = acos(cosAlpha);
  }
  else if (l2e1 == 0 && l2e2 == 0)
    alpha = M_PI / 3;
  else
    alpha = M_PI / 2;
  return alpha;
}

double getTwoTriangleDihedralAngle(const Vec3d & e0, const Vec3d & e1, const Vec3d & t0, const Vec3d & t1)
{
  Vec3d n0 = getTriangleScaledNormal(e0, e1, t0);
  Vec3d n1 = getTriangleScaledNormal(e0, e1, t1);
  return getVectorAngle(n0, n1);
}

//double getTwoTriangleInnerDihedralAngle(const IndexedTriangle & it0, const IndexedTriangle & it1)
//{
//  Vec3d n0 = getTriangleScaledNormal(it0.pos[0], it0.pos[1], it0.pos[2]);
//  Vec3d n1 = getTriangleScaledNormal(it1.pos[0], it1.pos[1], it1.pos[2]);
//  Vec3i e = getSharedVertices(it0.vtxID, it1.vtxID);
//  assert(e[0] >= 0 && e[1] >= 0 && e[2] < 0); // assert shared is an edge
//  Vec3d t0, t1, edgePos;
//  for(int i = 0; i < 3; i++)
//    if (it0.vtxID[i] != e[0] && it0.vtxID[i] != e[1]) {
//      t0 = it0.pos[i];
//      edgePos = it0.pos[(i+1)%3]; // one pos on the joint edge
//      break;
//    }
//
//  for(int i = 0; i < 3; i++)
//    if (it1.vtxID[i] != e[0] && it1.vtxID[i] != e[1]) {
//      t1 = it1.pos[i];
//      break;
//    }
//
//  double normalAngle = getVectorAngle(n0, n1);
//  double out0 = dot(t0 - edgePos, n1); // whether t0 is outside n1
//  double out1 = dot(t1 - edgePos, n0); // whether t1 is outside n0
//  if (out0 > 0 && out1 > 0) { return M_PI + normalAngle; }
//  else if (out0 < 0 && out1 < 0) { return M_PI - normalAngle; }
//}


// assert len2(normal) > 0
Vec3d getClosestPointToPlaneWithScaledNormal(const Vec3d & queryPoint, const Vec3d & scaledNormal, const Vec3d & start)
{
  Vec3d diff = queryPoint - start;
  double d = dot(scaledNormal, diff);
  double normalL2 = len2(scaledNormal);
  assert(normalL2 > 0);
  return queryPoint - d * scaledNormal / normalL2;
}

Vec3d getClosestPointToPlaneWithNormal(const Vec3d & queryPoint, const Vec3d & normal, const Vec3d & start)
{
  Vec3d diff = queryPoint - start;
  double d = dot(normal, diff);
  return queryPoint - d * normal;
}

Vec3d getClosestPointToLineSegment(const Vec3d & queryPoint, const Vec3d & lineStart, const Vec3d & lineEnd)
{
  Vec3d lineVec = lineEnd - lineStart;

  double d = dot(lineVec, queryPoint - lineStart);

  // the closest point is lineStart
  if (d <= 0) return lineStart;

  double lineLen2 = len2(lineVec);
  // the closest point is lineEnd
  if (d >= lineLen2) return lineEnd;

  return lineStart + (d / lineLen2) * lineVec;
}

Vec3d getClosestPointToTriangleWithFeature(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2, int feature)
{
  if (feature == 0) return t0;
  else if (feature == 1) return t1;
  else if (feature == 2) return t2;
  else if (feature == 3) // edge 01
  {
    return getClosestPointToLineSegment(queryPoint, t0, t1);
  }
  else if (feature == 4) // edge 12
  {
    return getClosestPointToLineSegment(queryPoint, t1, t2);
  }
  else if (feature == 5) // edge 20
  {
    return getClosestPointToLineSegment(queryPoint, t2, t0);
  }
  // else, feature == 6, on triangle
  Vec3d scaledNormal = getTriangleScaledNormal(t0, t1, t2);
  return getClosestPointToPlaneWithScaledNormal(queryPoint, scaledNormal, t0);
}

Vec3d getClosestPointToTriangleWithNormalAndFeature(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2, const Vec3d & normal, int feature)
{
  if (feature == 0) return t0;
  else if (feature == 1) return t1;
  else if (feature == 2) return t2;
  else if (feature == 3) // edge 01
  {
    return getClosestPointToLineSegment(queryPoint, t0, t1);
  }
  else if (feature == 4) // edge 12
  {
    return getClosestPointToLineSegment(queryPoint, t1, t2);
  }
  else if (feature == 5) // edge 20
  {
    return getClosestPointToLineSegment(queryPoint, t2, t0);
  }
  // else, feature == 6, on triangle
  return getClosestPointToPlaneWithNormal(queryPoint, normal, t0);
}

namespace
{

// scaledTriangleNormal must not be zero-length, but need not be unit length
// is queryPoint to the left of the edge (edgeStart -> edgeEnd)
bool isToLeftOfTriangleEdge(const Vec3d & queryPoint, const Vec3d & scaledTriangleNormal, const Vec3d & edgeStart, const Vec3d & edgeEnd)
{
  double d = dot(cross(edgeEnd - edgeStart, queryPoint - edgeStart), scaledTriangleNormal);
  return d > 0;
}

// delta: vector from an arbitrary point on the plane to the query point
double getSquaredDistanceToPlaneWithScaledNormalAndDelta(const Vec3d & scaledPlaneNormal, const Vec3d & delta)
{
  double d = dot(scaledPlaneNormal, delta);
  double normalLen2 = len2(scaledPlaneNormal);
  return (d*d) / normalLen2;
}

} // anonymous namespace


double getSquaredDistanceToLineSegment(const Vec3d & queryPoint, const Vec3d & lineStart, const Vec3d & lineEnd)
{
  Vec3d lineVec = lineEnd - lineStart;

  Vec3d segStart2Query = queryPoint - lineStart;
  double d = dot(lineVec, segStart2Query);

  // return distance from segStart to query point
  if (d <= 0) return len2(segStart2Query);

  double lineLen2 = len2(lineVec);
  // return distance from segEnd to query point
  if (d > lineLen2) return len2(queryPoint - lineEnd);

  // return distance from the infinite line to query point
  return len2(cross(lineVec, segStart2Query)) / lineLen2;
}

// also returns the closest feature to the query point:
//  0: vertex0
//  1: vertex1
//  2: vertex2
//  3: edge among 01
//  4: edge among 12
//  5: edge among 20
//  6: the face itself
double getSquaredDistanceToTriangle(const Vec3d & queryPoint, const Vec3d & vertex0, const Vec3d & vertex1, const Vec3d & vertex2, int & feature)
{
  Vec3d scaledNormal = getTriangleScaledNormal(vertex0, vertex1, vertex2);

  if(scaledNormal != Vec3d(0.0) &&
      isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex0, vertex1) &&
      isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex1, vertex2) &&
      isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex2, vertex0))
  {
    // the closest point on triangle to queryPoint is the same closet point on the plane where the triangle lies to queryPoint
    feature = 6;
    return getSquaredDistanceToPlaneWithScaledNormalAndDelta(scaledNormal, queryPoint-vertex0);
  }
  else // the projection of the queryPoint onto the triangle plane is outside the triangle, or the triangle is degenerate
    {    // then we query the closest distance from the query point to all the three edges
      double d0 = getSquaredDistanceToLineSegment(queryPoint, vertex1, vertex0);
      double d1 = getSquaredDistanceToLineSegment(queryPoint, vertex2, vertex1);
      double d2 = getSquaredDistanceToLineSegment(queryPoint, vertex0, vertex2);

      pair<double, int> sortBuffer[3] = { { d0, 0 }, { d1, 1 }, { d2, 2 } };
      sort(sortBuffer, sortBuffer+3);
      if (sortBuffer[0].first == sortBuffer[1].first)
      {
        // closest feature is a vertex
        int edgeIDToFeatureMap[3][3] = { {-1, 1, 0}, {1, -1, 2}, {0, 2, -1} };
        feature = edgeIDToFeatureMap[sortBuffer[0].second][sortBuffer[1].second];
        assert(feature >= 0);
      }
      else
      { // closest feature is an edge
        feature = sortBuffer[0].second + 3;
      }
      return sortBuffer[0].first;
    }
}

void getTetBarycentricWeights(const Vec3d & queryPoint, const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, double weight[4])
{
  //       |x1 y1 z1 1|         |x  y  z  1|        |x1 y1 z1 1|        |x1 y1 z1 1|        |x1 y1 z1 1|
  //  D0 = |x2 y2 z2 1|   D1 =  |x2 y2 z2 1|   D2 = |x  y  z  1|   D3 = |x2 y2 z2 1|   D4 = |x2 y2 z2 1|
  //       |x3 y3 z3 1|         |x3 y3 z3 1|        |x3 y3 z3 1|        |x  y  z  1|        |x3 y3 z3 1|
  //       |x4 y4 z4 1|         |x4 y4 z4 1|        |x4 y4 z4 1|        |x4 y4 z4 1|        |x  y  z  1|
  //  wi = Di / D0

  double tetDet = getTetDeterminant(a, b, c, d);

  for(int i=0; i<4; i++)
  {
    // compute D[i+1]
    Vec3d buf[4] = { a, b, c, d };
    buf[i] = queryPoint;
    double D = getTetDeterminant(buf[0], buf[1], buf[2], buf[3]);
    weight[i] = D / tetDet;
  }
}

double getSquaredDistanceToTet(const Vec3d & queryPoint, const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d)
{
//  const int OTetKey::tetFaceIndex[4][3] =  { { 1, 2, 3 }, { 0, 3, 2 }, { 0, 1, 3 }, { 0, 2, 1 } };
  if (getScaledSignedDistanceToTrianglePlane(queryPoint, b, c, d) <= 0.0 ||
      getScaledSignedDistanceToTrianglePlane(queryPoint, a, d, c) <= 0.0 ||
      getScaledSignedDistanceToTrianglePlane(queryPoint, a, b, d) <= 0.0 ||
      getScaledSignedDistanceToTrianglePlane(queryPoint, a, c, b) <= 0.0)
    return 0.0; // inside the tet

  double d2 = getSquaredDistanceToTriangle(queryPoint, b, c, d);
  d2 = min(d2, getSquaredDistanceToTriangle(queryPoint, a, d, c));
  d2 = min(d2, getSquaredDistanceToTriangle(queryPoint, a, b, d));
  d2 = min(d2, getSquaredDistanceToTriangle(queryPoint, a, c, b));
  return d2;
}
