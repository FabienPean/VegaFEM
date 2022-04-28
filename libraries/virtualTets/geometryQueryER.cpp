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

#include "geometryQueryER.h"
using namespace std;

namespace
{

// scaledTriangleNormal must not be zero-length, but need not be unit length
// is queryPoint to the left of the edge (edgeStart -> edgeEnd)
bool isToLeftOfTriangleEdge(const Vec3ER & queryPoint, const Vec3ER & scaledTriangleNormal, const Vec3ER & edgeStart, const Vec3ER & edgeEnd)
{
  ER d = dot(cross(edgeEnd - edgeStart, queryPoint - edgeStart), scaledTriangleNormal);
  return ER_sign(d) > 0;
}

ER squaredDistanceToPlane(const Vec3ER & scaledPlaneNormal, const Vec3ER & delta)
{
  ER d = dot(scaledPlaneNormal, delta);
  ER normalLen2 = len2(scaledPlaneNormal);
  return (d*d) / normalLen2;
}

// distance from query point to line segment (segStart -> segEnd)
ER squaredDistanceToLineSegment(const Vec3ER & queryPoint, const Vec3ER & lineStart, const Vec3ER & lineEnd)
{
  Vec3ER lineVec = lineEnd - lineStart;

  Vec3ER segStart2Query = queryPoint - lineStart;
  ER d = dot(lineVec, segStart2Query);

  // return distance from segStart to query point
  if (ER_sign(d) <= 0) return len2(segStart2Query);

  ER lineLen2 = len2(lineVec);
  // return distance from segEnd to query point
  if (d > lineLen2) return len2(queryPoint - lineEnd);

  // return distance from the infinite line to query point
  return len2(cross(lineVec, segStart2Query)) / lineLen2;
}

Vec3ER getClosestPointToPlane(const Vec3ER & queryPoint, const Vec3ER & scaledPlaneNormal, const Vec3ER & planeStart)
{
  ER d = dot(scaledPlaneNormal, queryPoint - planeStart);
  ER normalLen2 = len2(scaledPlaneNormal);
  assert(normalLen2 > 0);
  return queryPoint - d * scaledPlaneNormal / normalLen2;
}

Vec3ER getClosestPointToLineSegment(const Vec3ER & queryPoint, const Vec3ER & lineStart, const Vec3ER & lineEnd)
{
  Vec3ER lineVec = lineEnd - lineStart;

  ER d = dot(lineVec, queryPoint - lineStart);

  // the closest point is lineStart
  if (ER_sign(d) <= 0) return lineStart;

  ER lineLen2 = len2(lineVec);
  // the closest point is lineEnd
  if (d >= lineLen2) return lineEnd;

  return lineStart + (d / lineLen2) * lineVec;
}

} // anonymous namespace

// also returns the closest feature to the query point:
//  0: vertex0
//  1: vertex1
//  2: vertex2
//  3: edge among 01
//  4: edge among 12
//  5: edge among 20
//  6: the face itself
ER squaredDistanceToTriangle(const Vec3ER & queryPoint, const Vec3ER & vertex0, const Vec3ER & vertex1, const Vec3ER & vertex2, int & feature)
{
  Vec3ER scaledNormal = cross(vertex1 - vertex0, vertex2 - vertex0);

  if(scaledNormal != VEC3ER_NULL
     && isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex0, vertex1)
     && isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex1, vertex2)
     && isToLeftOfTriangleEdge(queryPoint, scaledNormal, vertex2, vertex0))
  {
    // the closest point on triangle to queryPoint is the same closet point on the plane where the triangle lies to queryPoint
    feature = 6;
    return squaredDistanceToPlane(scaledNormal, queryPoint-vertex0);
  }
  else // the projection of the queryPoint onto the triangle plane is outside the triangle, or the triangle is degenerate
  {    // then we query the closest distance from the query point to all the three edges
    ER d0 = squaredDistanceToLineSegment(queryPoint, vertex1, vertex0);
    ER d1 = squaredDistanceToLineSegment(queryPoint, vertex2, vertex1);
    ER d2 = squaredDistanceToLineSegment(queryPoint, vertex0, vertex2);

    pair<ER, int> sortBuffer[3] = { { d0, 0 }, { d1, 1 }, { d2, 2 } };
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

Vec3ER getClosestPointToTriangleWithFeature(const Vec3ER & queryPoint, const Vec3ER & vertex0, const Vec3ER & vertex1, const Vec3ER & vertex2, int feature)
{
  if (feature == 0) return vertex0;
  else if (feature == 1) return vertex1;
  else if (feature == 2) return vertex2;
  else if (feature == 3) // edge 01
    return getClosestPointToLineSegment(queryPoint, vertex0, vertex1);
  else if (feature == 4) // edge 12
    return getClosestPointToLineSegment(queryPoint, vertex1, vertex2);
  else if (feature == 5) // edge 20
    return getClosestPointToLineSegment(queryPoint, vertex2, vertex0);
  // else, feature == 6, on triangle
  const Vec3ER e1 = vertex1 - vertex0;
  const Vec3ER oe3 = vertex2 - vertex0;
  const Vec3ER scaledNormal = cross(e1, oe3);
  return getClosestPointToPlane(queryPoint, scaledNormal, vertex0);
}

