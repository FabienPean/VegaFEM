/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Christopher Twigg, Daniel Schroeder      *
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

/*
  A triangle
  Jernej Barbic, CMU
*/

#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include <math.h>
#include <vector>
#include <iostream>
#include "boundingBox.h"
#include "minivector.h"

class TriangleBasic
{
public:

  TriangleBasic(const Vec3d & first, const Vec3d & second, const Vec3d & third): index_(0) { vertex[0] = first; vertex[1] = second; vertex[2] = third; }

  // accessors
  inline const Vec3d & first() const {return vertex[0] ;}
  inline const Vec3d & second() const {return vertex[1] ;}
  inline const Vec3d & third() const {return vertex[2];}
  inline const Vec3d & operator[] (size_t i) const { return vertex[i]; }
  inline int index() const { return index_; }

  // index can be used to keep track of what triangle this is
  inline void setIndex(int index_) { this->index_ = index_; }

  // squared 3d distance to a point
  double distanceToPoint2(const Vec3d & point) const { std::cout << "Unimplemented..." << std::endl; return 1;} // unimplemented (it is implemented in class "TriangleWithCollisionInfo" below)
  double distanceToPoint(const Vec3d & point) const { return sqrt(distanceToPoint2(point));}

  // return the triangle normal
  Vec3d triangleNormal() const { return norm(cross(vertex[1] - vertex[0], vertex[2] - vertex[0])); }

  // whether a point is inside/outside/on the triangle
  // return: 1 = point is outside the triangle (assuming normal points outward)
  //        -1 = point is inside
  //         0 = point is on the triangle
  int pointSide(const Vec3d & point) const;

  bool doesIntersectBox(const BoundingBox & bbox) const;

  // int lineSegmentIntersection(Vec3d segmentStart, Vec3d segmentEnd, Vec3d * intersectionPoint);
  int lineSegmentIntersection(const Vec3d & segmentStart, const Vec3d & segmentEnd, Vec3d * intersectionPoint = NULL, double * t = NULL, double * barycentricWeight = NULL) const;
    //    Output: intersection point (when it exists)
    //    Return: -1 = triangle is degenerate (a segment or point)
    //             0 = disjoint (no intersect)
    //             1 = intersect in unique point I1
    //             2 = are in the same plane


  void render() const;
  void renderEdges() const;

  Vec3d getBarycentricLocation(double alpha, double beta, double gamma) const { return alpha * vertex[0] + beta * vertex[1] + gamma * vertex[2]; }

protected:
  Vec3d vertex[3];
  int index_;
};

// makes the triangle list unique, using the "index_" field of the <TriangleClass>
template<class TriangleClass>
void makeUniqueList(const std::vector<TriangleClass*> & triangleList, std::vector<TriangleClass*> & uniqueList);

template<class TriangleClass>
void makeUniqueList(std::vector<TriangleClass*> & triangleList); // overwrites triangleList with the unique list

class TriangleWithCollisionInfo : public TriangleBasic
{
public:
  TriangleWithCollisionInfo(const Vec3d & first_g, const Vec3d & second_g, const Vec3d & third_g):
        TriangleBasic(first_g,second_g,third_g) { ComputeCollisionInfo(); }

  // squared 3d distance to a point
  // also returns the closest feature to the query point:
  //  0: vertex0
  //  1: vertex1
  //  2: vertex2
  //  3: edge among 01
  //  4: edge among 12
  //  5: edge among 20
  //  6: the face itself
  double distanceToPoint2(const Vec3d & point, int * closestFeature) const;
  double distanceToPoint2(const Vec3d & point, int * closestFeature, double * alpha, double * beta, double * gamma) const; // also returns the barycentric coordinates of the closest point

protected:

  // note: the following collision detection parameters are pre-computed with respect to a permuted set of triangle indices (a cycle)
    Mat3d Q;
    Vec3d x0;

    double sidea,sideb,sidec,area;
    Vec3d S1,S2,N11,N12,N21,N22;

  int permutation[6]; // used internally so that first triangle edge is always the longest
  // transformation(x) equals Q * x + x0;

  void ComputeCollisionInfo();
};

class TriangleWithCollisionInfoAndPseudoNormals : public TriangleWithCollisionInfo
{
public:
  // pseudoNormals is a caller-provided list of 7 normals, indexed the same as closest features below
  TriangleWithCollisionInfoAndPseudoNormals(const Vec3d & first_g, const Vec3d & second_g, const Vec3d & third_g, const Vec3d pseudoNormals[7]);

  inline const Vec3d & pseudoNormal(int closestFeature) const { return pseudoNormal_[closestFeature]; }
  inline Vec3d interpolatedVertexPseudoNormal(double alpha, double beta, double gamma) const { return norm(alpha * pseudoNormal_[0] + beta * pseudoNormal_[1] + gamma * pseudoNormal_[2]); }

  // for vertices, returns the vertex itself
  // for edges, returns the midpoint of the edge (for the purposes of distance sign test, this could be any point on the edge)
  // for faces, returns the face centroid
  inline const Vec3d & pseudoClosestPosition(int closestFeature) const { return pseudoClosestPosition_[closestFeature]; }

protected:
  Vec3d pseudoNormal_[7];
  Vec3d pseudoClosestPosition_[7];

};

#endif

