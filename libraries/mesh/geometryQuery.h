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

#ifndef GEOMETRYQUERY_H
#define GEOMETRYQUERY_H

#include "minivector.h"
#include <cmath>

inline double rad2deg(double x) { return x * (180.0 / M_PI); }
inline double deg2rad(double x) { return x * (M_PI / 180.0); }

// get angle at v0: /_v1v0v2
double getTriangleAngle(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2);

// get angle between vec0 and vec1
double getVectorAngle(const Vec3d & vec0, const Vec3d vec1);

// get the scaled normal of triangle <v0, v1, v2>: cross(v1-v0, v2-v0) = 2 * tri_area * n
// its length is twice the triangle area
inline Vec3d getTriangleScaledNormal(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2) { return cross(v1 - v0, v2 - v0); }

// get the normal of triangle <v0, v1, v2>
inline Vec3d getTriangleNormal(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2) { return norm(getTriangleScaledNormal(v0, v1, v2)); }

inline double getTriangleArea(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2) { return 0.5 * len(getTriangleScaledNormal(v0, v1, v2)); }

// return dot (queryPoint - v0, cross(v1-v0, v2-v0) )
// the returned value is positive if queryPoint is above triangle, negative if under and zero if on the plane where the triangle lies
// the returned value is 2 * tri_area * signed_distance
inline double getScaledSignedDistanceToTrianglePlane(const Vec3d & queryPoint, const Vec3d & v0, const Vec3d & v1, const Vec3d & v2) { return dot(queryPoint - v0, getTriangleScaledNormal(v0,v1,v2)); }

// robust computation of the angle at v0 in a triangle v0, v1, v2
// if any two points have the same position, it treats the triangle as a degenerate one with angle 0, 90 and 90
// if all three points have the same position, it treats the triangle angles as 60, 60, 60
double getTriangleAngleRobust(const Vec3d & v0, const Vec3d & v1, const Vec3d & v2);

// get the dihedral angle between triangle (e0, e1, t0) and (e0, e1, t1), triangle ordering does not matter
// return value of [0, PI]
double getTwoTriangleDihedralAngle(const Vec3d & e0, const Vec3d & e1, const Vec3d & t0, const Vec3d & t1);

// assert len(scaledNormal) > 0.0
Vec3d getClosestPointToPlaneWithScaledNormal(const Vec3d & queryPoint, const Vec3d & scaledNormal, const Vec3d & planeStart);
// asssert len(normal) == 1.0
Vec3d getClosestPointToPlaneWithNormal(const Vec3d & queryPoint, const Vec3d & normal, const Vec3d & planeStart);

Vec3d getClosestPointToLineSegment(const Vec3d & queryPoint, const Vec3d & segStart, const Vec3d & segEnd);

// also input the closest feature to the query point:
//  0: vertex0
//  1: vertex1
//  2: vertex2
//  3: edge among 01
//  4: edge among 12
//  5: edge among 20
//  6: the face itself
Vec3d getClosestPointToTriangleWithFeature(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2, int feature);
Vec3d getClosestPointToTriangleWithNormalAndFeature(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2, const Vec3d & normal, int feature);

double getSquaredDistanceToTriangle(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2, int & feature);
inline double getSquaredDistanceToTriangle(const Vec3d & queryPoint, const Vec3d & t0, const Vec3d & t1, const Vec3d & t2) { int f = 0; return getSquaredDistanceToTriangle(queryPoint,t0,t1,t2, f); }

double getSquaredDistanceToLineSegment(const Vec3d & queryPoint, const Vec3d & segStart, const Vec3d & segEnd);

// compute the inner dihedral angle of the two triangles
// triangle normals should agree with each other:
// n0 /                    /                           n0 /
//  \/                    /\n0                          \/
//  /            or      /     | n1          but not    /   | n1
// ----------           ----------                     ----------
//     | n1
// IndexedTriangle it0 and it1 should have two shared vtx indices which form the joint edge
// return [0, 2PI), return nan if normals don't agree
//double getTwoTriangleInnerDihedralAngle(const IndexedTriangle & it0, const IndexedTriangle & it1);


// computes det(A), for the 4x4 matrix A
//     [ 1 a ]
// A = [ 1 b ]
//     [ 1 c ]
//     [ 1 d ]
// It can also be computed as det(A) = dot(d - a, cross(b - a, c - a))
// When det(A) > 0, the tet has positive orientation.
// When det(A) = 0, the tet is degenerate.
// When det(A) < 0, the tet has negative orientation.
// The orientation can also be determined as:
// if a is under the plane of the triangle formed by <b, c, d>, then it has positive orientation
inline double getTetDeterminant(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d)
{
  return dot(d - a, cross(b - a, c - a));
}

// compute barycentric weights of tet <a,b,c,d> for queryPoint
void getTetBarycentricWeights(const Vec3d & queryPoint, const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, double weight[4]);

double getSquaredDistanceToTet(const Vec3d & queryPoint, const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d);

#endif
