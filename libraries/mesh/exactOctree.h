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

#ifndef EXACT_OCTREE_H
#define EXACT_OCTREE_H

// An octree implementation using exact predicates for precise geometry queries.
// The base class is ExactOctreeBase. It has several derived classes:
// ExactVertexOctree, ExactTriMeshOctree and ExactTetMeshOctree.

#include "minivector.h"
#include "boundingBox.h"
#include "simpleSphere.h"
#include "halfSpace.h"
#include "triMeshGeo.h"
#include "tetMeshGeo.h"
#include <cstring>
#include <vector>
#include <array>
#include <cassert>
#include <functional>
#include <cfloat>

// Base class for exact octrees
// defines node struct inside the tree
// defines range query and nearest query algorithm

class ExactOctreeBase
{
public:
  int getDepth() const { return depth; } // get depth of the tree, which equals to the maximum Node::depth in all nodes
  int getNumNodes() const { return nodes.size(); }

  // return true if the queried object intersects bb
  using BBFilter = std::function<bool(const BoundingBox & bb)>;
  // return true if the queried object intersects the element
  using ElementFilter = std::function<bool(int eleID)>;
  // process each element, usually is defined to check if the queried object intersects the element
  using ElementProcess = std::function<void(int eleID)>;

  // Query the octree:
  // Starting at the root, if the bounding box query toBB returns true,
  // recursively visit each children node, until reaching a leaf node,
  // where the element query toEle is called to each elementIDs stored in this leaf node.
  void rangeQuery(BBFilter toBB, ElementProcess toEle) const;

  // Query the octree and return hit elementID list:
  // Starting at the root, if the bounding box query toBB returns true,
  // recursively visit each children node, until reaching a leaf node,
  // where if the element query toEle returns true on one elementID stored in this leaf node,
  // the elementID is then pushed into elementIDList.
  // Note: elementIDList is not sorted or deduplicated.
  inline void rangeQuery(BBFilter toBB, ElementFilter toEle, std::vector<int> & elementIDList) const;

  // return the near/far distance of the bounding box bb towards the queried object
  using DistanceToBB = std::function<std::pair<double, double>(const BoundingBox & bb)>;
  // return the distance of the element with elementID towards the queried object
  using DistanceToElement = std::function<double(int elementID, double minDistanceSoFar)>;

  // Nearest object query:
  // Starting at the root, call toBB to check the near/far distance of the bounding box towards the queried object.
  // If the near distance is farther than the closest distance found so far, then skip.
  // Otherwise, recursively visit each children node, until reaching a leaf node,
  // where if calling toElement to compute the distance between the element and the queried object gives
  // a distance smaller than the closest distance found so far. In this case, the closest distance will be updated.
  // distanceHi: the value to initialze the "closest distance found so far" value in the algorithm.
  // Ususally it is set to be DBL_MAX. But it can be other values to cull away elements that are too far away as a priori.
  void nearestQuery(DistanceToBB toBB, DistanceToElement toElement, double & minDistance, double distanceHi = DBL_MAX) const;

protected:
  ExactOctreeBase() {}
  ~ExactOctreeBase() {}

  using Subdivide = std::function<void(const std::vector<int> & elements, const BoundingBox & bb,
      std::vector<int> childIDs[8], BoundingBox subBBs[8])>;

  // derived class calls this function to build the octree tree
  void buildFromRoot(Subdivide divideNode, int maxDepth, int maxNumElementsPerNode);

  // whether a bounding box intersects an element
  using ElementBBIntersect = std::function<bool(int eleID, const BoundingBox & bb)>;
  // given element ID list, create a bounding box that covers them
  using BBFromElementList = std::function<BoundingBox(const std::vector<int> & elementList)>;

  // an implementation to divide the elements within the bounding box bb into eight sub groups, store them
  // in childIDs and their bounding boxes in subBBs
  void boundingBoxPartition(const std::vector<int> & elements, const BoundingBox & bb,
      std::vector<int> childIDs[8], BoundingBox subBBs[8],
      const std::vector<BoundingBox> & elementBBs, BBFromElementList, ElementBBIntersect=nullptr);

  struct Node
  {
    int depth; // for the root node of the tree, depth = 0
    std::vector<int> indices;
    BoundingBox bb;
    std::array<int, 8> childrenIDs{ {-1,-1,-1,-1,-1,-1,-1,-1} };
    Node(int depth, std::vector<int> indices)  : depth(depth), indices(std::move(indices)) {}
    bool isLeaf() const { return indices.size() > 0; }
    // get the next available child index in childrenIDs starting at nextChildIndex
    // return 8 if not available
    // return nextChildIndex if childrenIDs[nextChildIndex] is available
    // nextChildIndex: [0, 8]
    int getNextChild(int nextChildIndex) const { for(; nextChildIndex < 8 && childrenIDs[nextChildIndex] < 0; nextChildIndex++) {} return nextChildIndex; }
  };

  int depth{0}; // how deep in the hierarchy is this octree

  std::vector<Node> nodes;
};

// octree to query vertices

class ExactVertexOctree : public ExactOctreeBase
{
public:
  // make empty octree
  ExactVertexOctree() {}
  virtual ~ExactVertexOctree() { clear(); }

  // maxDepth: maximum tree depth allowed
  // build the octree
  void build(const std::vector<Vec3d> & verticesList, int maxDepth, int maxNumVerticesPerNode);

  // query vertices inside simpleSphere/halfSpace with double-precision
  // vertexIDList are not sorted or deduplicated
  void rangeQuery(const SimpleSphere & simpleSphere, std::vector<int> & vertexIDList) const;
  void rangeQuery(const HalfSpace & halfSpace, std::vector<int> & vertexIDList) const;

protected:
  std::vector<Vec3d> vertices;

  void clear();
};

// octree to query a triangle mesh

class ExactTriMeshOctree : public ExactOctreeBase
{
public:
  ExactTriMeshOctree() : numVertices(0), numTriangles(0) {}
  virtual ~ExactTriMeshOctree() { clear(); }

  void build(const TriMeshRef triMesh, int maxDepth, int maxNumTrianglesPerNode);

  const std::vector<BoundingBox> & triangleBoundingBoxes() const { return triBBs; }

  int numMeshVertices() const { return numVertices; }
  int numMeshTriangles() const { return numTriangles; }

  // void rangeQuery(const SimpleSphere & simpleSphere, std::vector<int> & vertexIDList);

  // do exact self intersection
  // input triMesh should be the same mesh as the one used in build()
  // output triangleID pairs are sorted and deduplicated
  void selfIntersectionExact(const TriMeshRef triMesh, std::vector<std::pair<int, int>> & triangleIDList) const;

  // do exact intersection with another triangle mesh
  // output triangleID pairs are sorted and deduplicated
  void intersectionExact(const TriMeshRef triMesh, const ExactTriMeshOctree & otherOctree, const TriMeshRef & otherMesh,
      std::vector<std::pair<int,int>> & triangleIDList);
  void intersectionExact(const TriMeshRef triMesh, const TriMeshRef & otherMesh,
        std::vector<std::pair<int,int>> & triangleIDList);

  // do exact intersection with line segment / triangle
  void lineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> & triangleIDList) const;
  void triangleIntersectionExact(const TriMeshRef triMesh, Vec3d t0, Vec3d t1, Vec3d t2, std::vector<int> & triangleIDList) const;

  // detect whether the line segment hits a triangle exactly.
  // if it does, find the segment parameter segWeight[2] to locate the first intersection point from segStart to segEnd by:
  // intersection_point = segStart * segWeight[0] + segEnd * segWeight[1]
  bool lineSegmentFirstIntersectionPoint(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, double segWeight[2]) const;

  // get the closest triangle to queryPosition, double-precision
  int getClosestTriangle(const TriMeshRef triMesh, const Vec3d & queryPosition, int & feature) const;

  // check whether bounding boxes are fine
  bool sanityCheck(const TriMeshRef triMesh) const;

protected:
  int numVertices, numTriangles;
  std::vector<BoundingBox> triBBs;

  void clear();
  std::vector<int> trianglesUnderNode(int nodeID) const;
};

// octree to query a tet mesh

class ExactTetMeshOctree : public ExactOctreeBase
{
public:
  ExactTetMeshOctree() : numVertices(0), numTets(0) {}

  void build(const TetMeshRef & tetMesh, int maxDepth, int maxNumTetsPerNode);

//  int getClosestTet(const Vec3d & queryPosition, double * minDistance = nullptr, double distanceHi = DBL_MAX) const;

protected:
  int numVertices, numTets;
  std::vector<BoundingBox> tetBBs;

  void clear();
};


// =====================================================
// ============= BELOW ARE IMPLEMENTATION ==============
// =====================================================

inline void ExactOctreeBase::rangeQuery(BBFilter toBB, ElementFilter toEle, std::vector<int> & elementIDList) const
{
  auto elementProcess = [&](int eleID)
   {
    if (toEle(eleID))
      elementIDList.push_back(eleID);
  };
  rangeQuery(toBB, elementProcess);
}


#endif

