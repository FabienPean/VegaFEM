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

#ifndef __BOUNDINGBOX_H__
#define __BOUNDINGBOX_H__

//  Bounding Box
//  Author: Jernej Barbic, CMU

#include <vector>
#include "minivector.h"
#include <cfloat>
#include <cmath>
#include <iostream>

class TriangleBasic;
class TriangleWithCollisionInfo;
class TriangleWithCollisionInfoAndPseudoNormals;

class BoundingBox
{
public:
  BoundingBox() : bmin_(0.0), bmax_(1.0), center_(0.5), halfSides_(0.5) {}
  BoundingBox(Vec3d bmin_g, Vec3d bmax_g): bmin_(bmin_g), bmax_(bmax_g) { updateData();}

  BoundingBox(int numVertices, const Vec3d * vertices);
  template <class Vec3dContainer> explicit BoundingBox(const Vec3dContainer & vertices);
  template <class IntContainer> explicit BoundingBox(const Vec3d * allVertices, const IntContainer & vertexIDs);
  template <class IntContainer> explicit BoundingBox(const Vec3d * allVertices, const Vec3i * allTriangles, const IntContainer & triIDs);

  explicit BoundingBox(const std::vector<BoundingBox> & bbs);
  explicit BoundingBox(const std::vector<TriangleBasic> & bbs);
  explicit BoundingBox(const std::vector<TriangleWithCollisionInfo> & bbs);
  explicit BoundingBox(const std::vector<TriangleWithCollisionInfoAndPseudoNormals> & bbs);

  // accessors
  const Vec3d & bmin() const { return bmin_;}
  const Vec3d & bmax() const { return bmax_;}

  const Vec3d & center() const { return center_;}
  const Vec3d & halfSides() const { return halfSides_;}

  double diameter() const { return 2.0 * len(halfSides_); }
  double volume() const { return halfSides_[0] * halfSides_[1] * halfSides_[2] * 8; }

  // mutators
  void setbmin(const Vec3d & bmin_g) { bmin_=bmin_g; updateData();}
  void setbmin(double x, double y, double z) { bmin_=Vec3d(x,y,z); updateData();}
  void setbmax(const Vec3d & bmax_g) { bmax_=bmax_g; updateData();}
  void setbmax(double x, double y, double z) { bmax_=Vec3d(x,y,z); updateData();}
  void setbminmax(const Vec3d & bmin, const Vec3d & bmax) { bmin_=bmin; bmax_=bmax; updateData();}

  void render() const;

  double distanceToPoint(const Vec3d & point) const { return sqrt(distanceToPoint2(point)); }
  double furthestDistanceToPoint(const Vec3d & point) const { return sqrt(furthestDistanceToPoint2(point)); }

  // get squared distance
  // computing squared distance is much faster than computing distance because the later requires sqrt()
  double distanceToPoint2(const Vec3d & point) const;
  double furthestDistanceToPoint2(const Vec3d & point) const;

  // return true if point is inside or touching the bounding box
  bool checkInside(const Vec3d & point) const;
  // return true if bb is completely inside or touching from inside the bounding box
  bool checkInside(const BoundingBox & bb) const;

  // sanity check bmin <= bmax
  bool verifyBox() const;

  // expands from the center
  // factor of 1.0 indicates no expansion
  void expand(double expansionFactor);
  // expand the bounding box to include this point
  void expand(const Vec3d & point);
  void regularize(); // converts the box into one with all sides equal

  bool lineSegmentIntersection(const Vec3d & segmentStart, const Vec3d & segmentEnd, Vec3d * intersection) const;
  bool intersect(const BoundingBox & bb) const;

  // return the index of the longest side (bmax[i] - bmin[i]) and its value
  std::pair<int, double> longestSide() const;

  friend inline std::ostream & operator << (std::ostream & o, const BoundingBox & bb) { return o << "[ " << bb.bmin_ << " " << bb.bmax_ << " ]"; }
  void print() const; // cout << (*this) << endl

  // The ordering of children is as follows:
  // child ID -> < b2_b1_b0 > (binary notation), each bi is 0 or 1
  // bi = 0: in dimension i, this child is on the negative side
  // bi - 1: in dimension i, this child is on the positive side
  void createChildBoundingBoxes(BoundingBox children[8]) const;

  // get intersection of two bounding boxes
  // will return an invalid bounding box if this and bb does not intersect
  BoundingBox getIntersection(const BoundingBox & bb) const;

protected:
  template<class Triangle> void buildFromTriangles(const std::vector<Triangle> & tripool);

  void updateOnePoint(const Vec3d & p);
  void updateData(); // updates center and half-sides
  Vec3d bmin_,bmax_;
  Vec3d center_, halfSides_;
};

template <class Vec3dContainer>
BoundingBox::BoundingBox(const Vec3dContainer & vertices)
{
  // set bmin_, bmax_
  bmin_ = Vec3d(+DBL_MAX, +DBL_MAX, +DBL_MAX);
  bmax_ = Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for(const Vec3d & p : vertices)
  {
    updateOnePoint(p);
  }
  updateData();
}

template <class IntContainer>
BoundingBox::BoundingBox(const Vec3d * allVertices, const IntContainer & vertexIDs)
{
  bmin_ = Vec3d(+DBL_MAX, +DBL_MAX, +DBL_MAX);
  bmax_ = Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for(int ID : vertexIDs)
  {
    updateOnePoint(allVertices[ID]);
  }
  updateData();
}

template <class IntContainer>
BoundingBox::BoundingBox(const Vec3d * allVertices, const Vec3i * allTriangles, const IntContainer & triIDs)
{
  bmin_ = Vec3d(+DBL_MAX, +DBL_MAX, +DBL_MAX);
  bmax_ = Vec3d(-DBL_MAX, -DBL_MAX, -DBL_MAX);
  for(int triID : triIDs)
    for(int vtxID : allTriangles[triID])
    {
      updateOnePoint(allVertices[vtxID]);
    }
  updateData();
}
#endif
