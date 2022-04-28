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

#include <iostream>
#include "openGL-headers.h"
#include "predicates.h"
#include "triple.h"
#include <vector>
#include <stack>
#include <cstring>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <queue>
#include <fstream>
#include "basicAlgorithms.h"
#include "exactOctree.h"
#include "geometryQuery.h"
#include "containerHelper.h"
#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif

using namespace std;

/////////////////////////////////////////////////////////////
//                   ExactOctreeBase
/////////////////////////////////////////////////////////////


void ExactOctreeBase::rangeQuery(BBFilter toBB, ElementProcess toEle) const
{
  assert(nodes.empty() == false);
  stack<int> nodeStack;
  nodeStack.push(0);
  while(!nodeStack.empty())
  {
    int nodeID = nodeStack.top();
    nodeStack.pop();
    const auto & node = nodes[nodeID];
    if (toBB(node.bb) == false) continue;

    const auto & indices = node.indices;
    if(indices.size() > 0)
    { // leaf node
      for(int eleID : indices) { toEle(eleID); }
    }
    else
    {
      for(int childID : node.childrenIDs)
        if(childID >= 0) { nodeStack.push(childID); }
    }
  }
}

void ExactOctreeBase::nearestQuery(DistanceToBB toBB, DistanceToElement toElement,
    double & minDistance, double distanceHi) const
{
  assert(nodes.empty() == false);
  // distanceHi stores the current closest distances
  // traverse the hierarchy
  minDistance = distanceHi;

  // each entry in the stack is for one node: <higher/lower bound of distanceToNode, node ID>
  stack<tuple<double, double, int>> nodeStack;
  auto nodeBound = toBB(nodes[0].bb); // get the lower and higher bound of this BV node
  nodeStack.push(make_tuple(nodeBound.first, nodeBound.second, 0));

  while(nodeStack.size() > 0)
  {
    const auto & top = nodeStack.top();
    int nodeID = get<2>(top);
    double nodeLowerBound = get<0>(top), nodeHigherBound = get<1>(top);
    nodeStack.pop();
    const Node & node = nodes[nodeID];
    if (nodeLowerBound <= minDistance) // if this BV node's distance lower bound is smaller than distanceHi
    {
      // this node needs to be inspected further
      if (nodeHigherBound < minDistance)
        minDistance = nodeHigherBound;

      if (node.isLeaf())
      {
        // leaf node, traverse all children
        for(int eleID : node.indices)
        {
          double dist = toElement(eleID, minDistance);

          if (dist <= minDistance)
            minDistance = dist;
        }
      }
      else
      {
        // non-leaf node
        // sort the children nodes according to their distance
        vector<tuple<double, double, int> > childrenSortBuffer; // <lowerBound, higherBound, nodeID>
        for(int childID : node.childrenIDs)
        {
          if (childID < 0) continue;
          auto distBound = toBB(nodes[childID].bb);
          childrenSortBuffer.emplace_back(distBound.first, distBound.second, childID);
        }

        sort(childrenSortBuffer.begin(), childrenSortBuffer.end());
        for(auto rit = childrenSortBuffer.rbegin(); rit != childrenSortBuffer.rend(); rit++)
        {
          nodeStack.push(*rit);
        }
      }
    }
  }
}

void ExactOctreeBase::buildFromRoot(Subdivide divideNode, int maxDepth, int maxNumElementsPerNode)
{
  depth = 0;

  stack<int> nodeStack;
  nodeStack.push(0);
  while(nodeStack.empty() == false)
  {
    int nodeID = nodeStack.top();
    nodeStack.pop();

    int nodeDepth = nodes[nodeID].depth;
    if (nodeDepth > depth)
      depth = nodeDepth;


    int numElements = (int)nodes[nodeID].indices.size();

    // if there are fewer elements than the threshold value, or if max depth has been reached,
    // this node becomes a leaf which stores the element indices
    // max depth checking is necessary, otherwise the box might get split forever
    if ((numElements <= maxNumElementsPerNode) || (nodeDepth >= maxDepth))
      continue;

    vector<int> indices = move(nodes[nodeID].indices);
    assert(nodes[nodeID].indices.empty()); // non-leaf nodes should not store vertex indices

    std::vector<int> childIDs[8];
    BoundingBox subBBs[8];

    // assert(nodes[nodeID].bb.verifyBox());
    divideNode(indices, nodes[nodeID].bb, childIDs, subBBs);

    for(int i = 0; i < 8; i++)
    {
      if (childIDs[i].size() > 0)
      {
        int childrenID = int(nodes.size());
        nodes[nodeID].childrenIDs[i] = childrenID;
        nodes.emplace_back(nodeDepth+1, move(childIDs[i]));
        auto & childNode = nodes.back();
        childNode.bb = subBBs[i];
        // assert(subBBs[i].verifyBox());
        nodeStack.push(childrenID);
      }
    }
  }
}

void ExactOctreeBase::boundingBoxPartition(const std::vector<int> & elementList, const BoundingBox & bb,
    std::vector<int> childIDs[8], BoundingBox subBBs[8], const std::vector<BoundingBox> & elementBBs,
    BBFromElementList createBBFromElements, ElementBBIntersect intersectElementBB)
{
  Vec3d center = bb.center();

  BoundingBox spatialChildrenBBs[8];
  // partition node.bb into 8 children spatially
  bb.createChildBoundingBoxes(spatialChildrenBBs);

  // for each child, generate a list of vertices that intersect the child bounding box
  for(int eleID : elementList)
  {
    const auto & triBB = elementBBs[eleID];

    unsigned char lowerChildrenBBCode[3], higherChildrenBBCode[3]; // stores [0, 1]
    // lowerChildrenBBCCode stores which side bmin of triBB is in the childrenBBs on all dimenstions
    for(int dim = 0; dim < 3; dim++)
    {
      lowerChildrenBBCode[dim] = (triBB.bmin()[dim] > center[dim]);
      higherChildrenBBCode[dim] = (triBB.bmax()[dim] > center[dim]);
    }

    for(unsigned char ii = lowerChildrenBBCode[0]; ii <= higherChildrenBBCode[0]; ii++)
      for(unsigned char jj = lowerChildrenBBCode[1]; jj <= higherChildrenBBCode[1]; jj++)
        for(unsigned char kk = lowerChildrenBBCode[2]; kk <= higherChildrenBBCode[2]; kk++)
        {
          unsigned char code = ii + (jj << 1) + (kk << 2);
          // cout << spatialChildrenBBs[code] << endl;
          if (!intersectElementBB || intersectElementBB(eleID, spatialChildrenBBs[code]))
            childIDs[code].push_back(eleID);
//          assert(spatialChildrenBBs[code].intersect(triBB));
        }
  }

  for(int i = 0; i < 8; i++)
  {
    if (childIDs[i].size() > 0)
    {
      BoundingBox bbFromElements = createBBFromElements(childIDs[i]);

      // we can use spatialChildrenBBs[i] as childNode.bb
      // but it might be larger than what triangles occupy
      // we can shrink this bb by doing intersection with allTriVtxBB
      subBBs[i] = spatialChildrenBBs[i].getIntersection(bbFromElements);
//      assert(subBBs[i].verifyBox());
    }
  }
}


/////////////////////////////////////////////////////////////
//                   ExactVertexOctree
/////////////////////////////////////////////////////////////

void  ExactVertexOctree::clear()
{
  vertices.clear();
  nodes.clear();
}

void ExactVertexOctree::build(const std::vector<Vec3d> & verticesList, int maxDepth, int maxNumVerticesPerNode)
{
  nodes.clear();
  vertices = verticesList;
  vector<int> rootIndices(vertices.size());
  iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, move(rootIndices));
  nodes[0].bb = BoundingBox(verticesList);

  auto subdivide = [&](const vector<int> & elementList, const BoundingBox & bb,
      vector<int> childIDs[8], BoundingBox subBBs[8]) -> void
  {
    Vec3d center = bb.center();
    int numVertices = elementList.size();

    // for each child, generate a list of vertices that intersect the child bounding box
    for(int i = 0; i < numVertices; i++)
    {
      const Vec3d & v = verticesList[elementList[i]];
      unsigned char code = 0;
      code += (v[0] > center[0]);
      code += ((v[1] > center[1]) << 1);
      code += ((v[2] > center[2]) << 2);
      assert(code < 8);
      childIDs[code].push_back(elementList[i]);
    }

    for(int i = 0; i < 8; i++)
    {
      if (childIDs[i].size() > 0)
      {
        subBBs[i] = BoundingBox(verticesList.data(), childIDs[i]);
      }
    }
  };
  buildFromRoot(subdivide, maxDepth, maxNumVerticesPerNode);
}

void ExactVertexOctree::rangeQuery(const SimpleSphere & simpleSphere, std::vector<int> & vertexIDList) const
{
  auto toBB = [&](const BoundingBox & bb) -> bool
  {
    return simpleSphere.doesBoundingBoxIntersect(bb);
  };

  auto toEle = [&](int idx) -> bool
  {
    return simpleSphere.inside(vertices[idx]);
  };

  ExactOctreeBase::rangeQuery(toBB, toEle, vertexIDList);
}

void ExactVertexOctree::rangeQuery(const HalfSpace & halfSpace, std::vector<int> & vertexIDList) const
{
  auto toBB = [&](const BoundingBox & bb) -> bool
  {
    return halfSpace.intersect(bb);
  };

  auto toEle = [&](int idx) -> bool
  {
    return halfSpace.outside(vertices[idx]) <= 0;
  };

  ExactOctreeBase::rangeQuery(toBB, toEle, vertexIDList);
}

/////////////////////////////////////////////////////////////
//                   ExactTriMeshOctree
/////////////////////////////////////////////////////////////

void ExactTriMeshOctree::clear()
{
  nodes.clear();
  triBBs.clear();
}

void ExactTriMeshOctree::build(const TriMeshRef triMesh, int maxDepth, int maxNumTrianglesPerNode)
{
  nodes.clear();
  numVertices = triMesh.numVertices();
  numTriangles = triMesh.numTriangles();
  assert(numTriangles > 0);
  triBBs = triMesh.getTriangleBoundingBoxes();
  vector<int> rootIndices(triMesh.numTriangles());
  iota(rootIndices.begin(), rootIndices.end(), 0);
  BoundingBox rootBB(triMesh.positions(), triMesh.triangles(), rootIndices);
  nodes.emplace_back(0, move(rootIndices));
  // We don't build a bounding box out of all triMesh's positions because some vertices might not
  // be used by triangles
  nodes[0].bb = rootBB;

  // cout << streamRange(rootIndices) << endl;
  // cout << nodes[0].bb << endl;
  // assert(nodes[0].bb.verifyBox());

  auto bbFromTris = [&](const vector<int> & triList) -> BoundingBox
  {
    // get the bounding box of all triangles inside this childNode
    auto bb = BoundingBox(triMesh.positions(), triMesh.triangles(), triList);
    // assert(bb.verifyBox());
    return bb;
  };
  auto triBBIntersect = [&](int triID, const BoundingBox & bb) -> bool
  {
    return intersectTriAABB(triMesh.pos(triID, 0), triMesh.pos(triID, 1), triMesh.pos(triID, 2), bb.bmin(), bb.bmax());
  };

  auto subdivide = [&](const vector<int> & elementList, const BoundingBox & bb,
        vector<int> childIDs[8], BoundingBox subBBs[8]) -> void
  {
    // assert(bb.verifyBox());
    boundingBoxPartition(elementList, bb, childIDs, subBBs, triBBs, bbFromTris, triBBIntersect);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTrianglesPerNode);
}

void ExactTriMeshOctree::selfIntersectionExact(const TriMeshRef triMesh, std::vector<std::pair<int, int>> & triangleIDList) const
{
  // cout << "selfIntersectionExact" << endl;
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());
  assert(nodes.empty() == false);

  unordered_set<UEdgeKey> candidateSet;

  // BVTT is such that it always compares BVs at the same level,
  //   or at immediately adjacent levels, with recursing first into the first tree node

  // build a stack that stores <BVNodeA, BVNodeB, childIndex>
  // we will compare intersection between BVNodeA and BVNodeB
  // we make sure BVNodeA <= BVNodeB to avoid redundant comparison

  auto traverseTree = [&]()
  {
    stack<triple<int, int, int> > traversalStack;
    traversalStack.push(triple<int,int,int>(0, 0, 0));
    while (true)
    {
      triple<int,int,int> stackEntry = traversalStack.top();

      //printf("Top off stack: %d vs %d. Child index: %x . Stack size: %d\n", stackEntry.first, stackEntry.second, stackEntry.third, (int)traversalStack.size());

      int nodeAID = stackEntry.first;
      int nodeBID = stackEntry.second;

  //    cout << "stack top: " << nodeAID << " " << nodeBID << " " << stackEntry.third << endl;

      assert(nodeAID >= 0 && nodeAID < (int)nodes.size());
      assert(nodeBID >= 0 && nodeBID < (int)nodes.size());
      const auto * nodeA = &nodes[nodeAID];
      const auto * nodeB = &nodes[nodeBID];

      bool neighboringNodes = (nodeAID == nodeBID);

      //printf("Neighboring nodes: %s\n", neighboringNodes ? "TRUE" : "FALSE");

      bool potentialCollision = true;

      // for non-neighboring nodes, we can rule out collisions by checking the BVs
      // for neighboring nodes, this test would be useless as BVs always intersect
      if (!neighboringNodes)
      {
        potentialCollision = nodeA->bb.intersect(nodeB->bb);
      }

      // potentialCollision flag is now set (in either case; neighboring nodes or not)

      //printf("BVs %d %d potential collision is subtree: %s\n", nodeA, nodeB, potentialCollision ? "TRUE" : "FALSE");

      // are the nodes on deepest level
      // negative sign by convention indicates leaf
      bool deepestLevelA = nodeA->isLeaf();
      bool deepestLevelB = nodeB->isLeaf();

      bool BVVT_leaf = deepestLevelA && deepestLevelB; // if both nodes are leaves

      //printf("BVVT leaf: %s\n", BVVT_leaf ? "TRUE" : "FALSE");

      if (BVVT_leaf && potentialCollision)
      {
        //printf("Checking all triangle pairs...\n");
        // check all triangle pairs; note: nodeAp and nodeBp could be the same
        auto processTriPair = [&](int triIDA, int triIDB)
        {
          if (triMesh.tri(triIDA).intersect(triMesh.tri(triIDB))) return; // if two triangles are neighbors
          UEdgeKey ue(triIDA, triIDB);
          if (usetFind(candidateSet, ue)) return;
          // do tri vs tri intersection
          if (triBBs[triIDA].intersect(triBBs[triIDB]) == false) return;

          candidateSet.insert(ue);
        };

        if (nodeA == nodeB)
        {
          const auto & indices = nodeA->indices;
          for(size_t i = 0; i < indices.size(); i++)
            for(size_t j = i+1; j < indices.size(); j++)
            {
              processTriPair(indices[i], indices[j]);
            }
        }
        else
        {
          for(int triIDA : nodeA->indices)
            for(int triIDB : nodeB->indices)
            {
              if (triIDA == triIDB) continue; // same triangle, continue
              processTriPair(triIDA, triIDB);
            }
        }
      }

      int childIndex = stackEntry.third;

      // for intra-node tests, we store (i,j) indices of current child pair in upper and lower parts of childIndex;
      // we only traverse children with i <= j; this avoids checking pairs of nodes twice
      short int * childIndexB = (short int*)&childIndex;
      short int * childIndexA = childIndexB + 1;
      bool recurseOnFirstNode = true;
      int parentChildIndex;

      // go deeper from <nodeA, nodeB> to visit one more pair
      // if nodeA is nodeB, then the next would be their children pair at location childIndexA and childIndexB
      // else, either nodeA or nodeB will go deeper to its child at location childIndex
      auto drillDown = [&]()
      {
        if (nodeAID == nodeBID)
        {
          nodeAID = nodeA->childrenIDs[*childIndexA];
          nodeBID = nodeB->childrenIDs[*childIndexB];
        }
        else
        {
          if (recurseOnFirstNode)
            nodeAID = nodeA->childrenIDs[childIndex];
          else
            nodeBID = nodeB->childrenIDs[childIndex];
        }
        assert(nodeAID >= 0 && nodeAID < (int)nodes.size());
        assert(nodeBID >= 0 && nodeBID < (int)nodes.size());
        traversalStack.push(triple<int,int,int>(nodeAID, nodeBID, childIndex));
      };

      auto goToNextChild = [&]()
      {
        if (nodeA == nodeB) (*childIndexB)++;
        else childIndex++;
      };

      // modify childIndex to get the next available child index
      // return true if we can find one
      auto goToValidChild = [&]()
      {
        bool hasChild = true;
        if (nodeA == nodeB)
        {
          // find an available childIndexA
          for(; *childIndexA < 8 && nodeA->childrenIDs[*childIndexA] < 0; (*childIndexB) = ++(*childIndexA)) {}

          // find an available childIndexB
          for(; *childIndexB < 8 && nodeB->childrenIDs[*childIndexB] < 0; (*childIndexB)++) {}

          if (*childIndexB >= 8)
          {
            // find next available childIndexA
            for((*childIndexA)++; *childIndexA < 8 && nodeA->childrenIDs[*childIndexA] < 0; (*childIndexA)++) {}
            (*childIndexB) = (*childIndexA);
            hasChild = (*childIndexA < 8);
          }
          else
          {
            assert(nodeA->childrenIDs[*childIndexA] >= 0);
          }
        }
        else
        {
          recurseOnFirstNode = (nodeB->isLeaf() || ((nodeA->depth <= nodeB->depth) && (nodeA->isLeaf() == false) ));
          const auto recurseNode = (recurseOnFirstNode ? nodeA : nodeB);

          for(; childIndex < 8 && recurseNode->childrenIDs[childIndex] < 0; childIndex++) {}

          hasChild = (childIndex < 8);
        }
        return hasChild;
      };

      if (!potentialCollision || BVVT_leaf)
      {
        // we are traversing the octree
        // Since there's no collision or it's two-leaf case where collision has been checked
        // we stop this search direction and pop from the stack

        parentChildIndex = childIndex;
        bool hasChild = true;
        do
        {
          childIndex = parentChildIndex;

          traversalStack.pop();
          if (traversalStack.empty())
          {
            return;
          }
          triple<int, int, int> parentEntry = traversalStack.top();
          parentChildIndex = parentEntry.third;
          nodeAID = parentEntry.first;
          nodeBID = parentEntry.second;
          nodeA = &nodes[nodeAID];
          nodeB = &nodes[nodeBID];
          goToNextChild();
          hasChild = goToValidChild();
        }
        while (hasChild == false);
        drillDown();
      }
      else
      {
        // Here, nodeA/B cannot be both leaves
        // recursively go down one level
        childIndex = 0;
        goToValidChild();
        drillDown();
      }
    }
  };
  traverseTree();

  // cout << "Finish finding candidate pairs: " << candidatePairs.size() << endl;

  //  sortAndDeduplicate(candidatePairs);
  vector<UEdgeKey> candidatePairs(candidateSet.begin(), candidateSet.end()); // store all those candidate triangle pairs for parallel evaluations
  vector<char> intersected(candidatePairs.size(), 0);

#ifdef USE_TBB
  tbb::parallel_for(tbb::blocked_range<int>(0, candidatePairs.size()), [&](const tbb::blocked_range<int> & rng)
  {
    for (int i = rng.begin(); i != rng.end(); ++i)
    {
#else
    for(size_t i = 0; i < candidatePairs.size(); ++i)
    {
#endif
      int triIDA = candidatePairs[i][0], triIDB = candidatePairs[i][1];
      if (intersectTriTri(triMesh.pos(triIDA, 0), triMesh.pos(triIDA, 1), triMesh.pos(triIDA,2),
          triMesh.pos(triIDB, 0), triMesh.pos(triIDB, 1), triMesh.pos(triIDB,2)))
      {
        intersected[i] = 1;
      }
    }
#ifdef USE_TBB
  }, tbb::auto_partitioner()); //end for locations
#endif

  for(size_t i = 0; i < candidatePairs.size(); i++)
  {
    if (intersected[i])
    {
      triangleIDList.push_back(make_pair(candidatePairs[i][0], candidatePairs[i][1]));
    }
  }

  // cout << "selfIntersectionExact Done" << endl;
}

void ExactTriMeshOctree::intersectionExact(const TriMeshRef triMesh, const TriMeshRef & otherMesh,
      std::vector<pair<int,int>> & triangleIDList)
{
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());
  assert(nodes.empty() == false);

  vector<pair<int,int>> allPairList;

#ifdef USE_TBB
  struct ThreadLocalData
  {
    vector<int> IDlist;
    vector<pair<int,int>> pairList;
  };
  tbb::enumerable_thread_specific<ThreadLocalData> threadLocalData;
  tbb::parallel_for(tbb::blocked_range<int>(0, otherMesh.numTriangles()), [&](const tbb::blocked_range<int> & rng)
  {
//  for(int oID = 0; oID < otherMesh.numTriangles(); oID++)
    auto & local = threadLocalData.local();
    auto & IDlist = local.IDlist;
    auto & pairList = local.pairList;
    for(int oID = rng.begin(); oID != rng.end(); oID++)
#else
    vector<int> IDlist;
    auto & pairList = allPairList;
    for(int oID = 0; oID < otherMesh.numTriangles(); oID++)
#endif
    {
      IDlist.clear();
      triangleIntersectionExact(triMesh, otherMesh.pos(oID,0), otherMesh.pos(oID,1), otherMesh.pos(oID,2), IDlist);
      sortAndDeduplicate(IDlist);
      for(int triID : IDlist)
        pairList.emplace_back(triID, oID);
    }
#ifdef USE_TBB
  });
  for(const auto & local : threadLocalData)
    vectorInsertRangeBack(allPairList, local.pairList);
#endif

  sortAndDeduplicate(allPairList);
  vectorInsertRangeBack(triangleIDList, allPairList);
}

void ExactTriMeshOctree::intersectionExact(const TriMeshRef triMesh, const ExactTriMeshOctree & otherOctree, const TriMeshRef & otherMesh,
    std::vector<std::pair<int,int>> & triangleIDList)
{
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());
  assert(nodes.empty() == false);
  assert(otherOctree.numVertices == otherMesh.numVertices());
  assert(otherOctree.numTriangles == otherMesh.numTriangles());
  assert(otherOctree.nodes.empty() == false);

//  ofstream fout("bbs.txt");
  vector<pair<int,int>> candidatePairs;

  // BVTT is such that it always compares spheres at the same level,
  //   or at immediately adjacent levels, with recursing first into the first tree node

  stack<triple<int, int, int> > traversalStack;

  traversalStack.push(triple<int,int,int>(0, 0, 0));

//  int iter = 0;
  while (true)
  {
    triple<int,int,int> stackEntry = traversalStack.top();

//    cout << "====== " << iter++ << " ======" << endl;
//    printf("Top off stack: %d vs %d. Child index: %x . Stack size: %d\n", stackEntry.first, stackEntry.second, stackEntry.third, (int)traversalStack.size());

    int nodeAID = stackEntry.first;
    int nodeBID = stackEntry.second;
    assert(nodeAID >= 0 && nodeAID < sizei(nodes));
    assert(nodeBID >= 0 && nodeBID < sizei(otherOctree.nodes));
    const auto * nodeA = &nodes[nodeAID];
    const auto * nodeB = &otherOctree.nodes[nodeBID];

    bool potentialCollision = nodeA->bb.intersect(nodeB->bb);
//    cout << "children " << streamRange(nodeA->childrenIDs) << " " << streamRange(nodeB->childrenIDs) << endl;
//    fout << nodeA->bb << " " << nodeB->bb << endl;
//    cout << "bb collision: " << potentialCollision << endl;

    // are the nodes on the deepest level
    // negative sign by convention indicates leaf
    bool deepestLevelA = nodeA->isLeaf();
    bool deepestLevelB = nodeB->isLeaf();

    bool BVVT_leaf = deepestLevelA && deepestLevelB;

    if (BVVT_leaf && potentialCollision)
    {
      for(int triA : nodeA->indices)
      {
        const BoundingBox & bbA = triBBs[triA];
        for(int triB : nodeB->indices)
        {
          const BoundingBox & bbB = otherOctree.triBBs[triB];

          // do tri vs tri intersection
          if (bbA.intersect(bbB) == false) continue;
//          cout << "add candidate pair " << triA << " " << triB << endl;
          candidatePairs.push_back(make_pair(triA, triB));
        }
      }
    }

    if (!potentialCollision || BVVT_leaf) // find next node
    {
      int childIndex = 0; // index of this stack entry's node in its parent's childIDs array
      bool endTraversal = false;
      bool recurseOnFirstNode = false;
      do
      {
        childIndex = traversalStack.top().third;
        traversalStack.pop();
        if (traversalStack.empty())
        {
          endTraversal = true;
          break;
        }
        triple<int, int, int> parentEntry = traversalStack.top();
        nodeAID = parentEntry.first; // now it stores the parentNodeAID
        nodeBID = parentEntry.second;// now it stores the parentNodeBID
//        cout << "  back to parent " << nodeAID << " " << nodeBID << endl;
        assert(nodeAID >= 0 && nodeAID < sizei(nodes));
        assert(nodeBID >= 0 && nodeBID < sizei(otherOctree.nodes));
        nodeA = &nodes[nodeAID];             // now it refers to parent nodeA
        nodeB = &otherOctree.nodes[nodeBID]; // now it refers to parent nodeB

        childIndex += 1; // get next child
        recurseOnFirstNode = (nodeB->isLeaf() || ((nodeA->depth <= nodeB->depth) && nodeA->isLeaf() == false ));
        if (recurseOnFirstNode) { childIndex = nodeA->getNextChild(childIndex); }
        else {childIndex = nodeB->getNextChild(childIndex); }
      }
      while (childIndex >= 8);
      if (endTraversal) break;

//      cout << "  at parent " << nodeAID << " " << nodeBID << ", ";
      // push the new node
      // get the child nodeID from parent node ID
      if (recurseOnFirstNode)
        nodeAID = nodeA->childrenIDs[childIndex];
      else
        nodeBID = nodeB->childrenIDs[childIndex];
      // now nodeAID and nodeBID stores the node IDs for the next stack entry
//      cout << "found child " << nodeAID << " " << nodeBID << " " << endl;

      traversalStack.push(triple<int,int,int>(nodeAID, nodeBID, childIndex));
    }
    else
    {
      bool recurseOnFirstNode = (nodeB->isLeaf() || ((nodeA->depth <= nodeB->depth) && nodeA->isLeaf() == false ));
      int childIndex = 0;
      // push the new node
      if (recurseOnFirstNode)
      {
        childIndex = nodeA->getNextChild(0);
        nodeAID = nodeA->childrenIDs[childIndex];
      }
      else
      {
        childIndex = nodeB->getNextChild(0);
        nodeBID = nodeB->childrenIDs[childIndex];
      }
      assert(childIndex < 8);
//      cout << "at a new entry " << nodeAID << " " << nodeBID << endl;
      traversalStack.push(triple<int,int,int>(nodeAID, nodeBID, childIndex));
    }
  } // end traversal loop

  // cout << "Finish finding candidate pairs: " << candidatePairs.size() << endl;
  sortAndDeduplicate(candidatePairs);

  vector<char> intersected(candidatePairs.size(), 0);

#ifdef USE_TBB
  tbb::parallel_for(tbb::blocked_range<int>(0, candidatePairs.size()), [&](const tbb::blocked_range<int> & rng)
  {
    for (int i = rng.begin(); i != rng.end(); ++i)
    {
#else
    for(size_t i = 0; i < candidatePairs.size(); ++i)
    {
#endif
      int triIDA = candidatePairs[i].first, triIDB = candidatePairs[i].second;
      assert(triIDA >= 0 && triIDA < numTriangles);
      assert(triIDB >= 0 && triIDB < otherMesh.numTriangles());
      if (intersectTriTri(triMesh.pos(triIDA, 0), triMesh.pos(triIDA, 1), triMesh.pos(triIDA,2),
          otherMesh.pos(triIDB, 0), otherMesh.pos(triIDB, 1), otherMesh.pos(triIDB,2)))
      {
        intersected[i] = 1;
      }
    }
#ifdef USE_TBB
  }, tbb::auto_partitioner()); //end for locations
#endif

  for(size_t i = 0; i < candidatePairs.size(); i++)
  {
    if (intersected[i])
      triangleIDList.push_back(candidatePairs[i]);
  }
}

void ExactTriMeshOctree::lineSegmentIntersectionExact(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, std::vector<int> & triangleIDList) const
{
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox & bb) -> bool
  {
    return intersectSegAABB(&segStart[0], &segEnd[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  auto toEle = [&](int idx) -> bool
  {
    return intersectSegTri(&segStart[0], &segEnd[0], &triMesh.pos(idx, 0)[0], &triMesh.pos(idx, 1)[0], &triMesh.pos(idx, 2)[0]);
  };

  ExactOctreeBase::rangeQuery(toBB, toEle, triangleIDList);
}

// note: no tbb threading should be used in triangleIntersectionExact
// because I use tbb on intersectionExact with another TriMesh
// and this function calls triangleIntersectionExact in an parallelized loop
// nested multi-threading can be a problem when I use thread-local data
void ExactTriMeshOctree::triangleIntersectionExact(const TriMeshRef triMesh, Vec3d t0, Vec3d t1, Vec3d t2, std::vector<int> & triangleIDList) const
{
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox & bb) -> bool
  {
    return intersectTriAABB(&t0[0], &t1[0], &t2[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  auto toEle = [&](int idx) -> bool
  {
    return intersectTriTri(&t0[0], &t1[0], &t2[0], &triMesh.pos(idx, 0)[0], &triMesh.pos(idx, 1)[0], &triMesh.pos(idx, 2)[0]);
  };

  ExactOctreeBase::rangeQuery(toBB, toEle, triangleIDList);
}

bool ExactTriMeshOctree::lineSegmentFirstIntersectionPoint(const TriMeshRef triMesh, Vec3d segStart, Vec3d segEnd, double segWeight[2]) const
{
  assert(numVertices == triMesh.numVertices());
  assert(numTriangles == triMesh.numTriangles());

  auto toBB = [&](const BoundingBox & bb) -> bool
  {
    return intersectSegAABB(&segStart[0], &segEnd[0], &bb.bmin()[0], &bb.bmax()[0]);
  };

  segWeight[0] = segWeight[1] = DBL_MAX;
  auto toEle = [&](int idx) -> void
  {
    double alpha[5];
    if (intersectSegTri(segStart, segEnd, triMesh.pos(idx,0), triMesh.pos(idx,1), triMesh.pos(idx,2), alpha, alpha+2))
    {
      if (alpha[0] < segWeight[0])
      {
        memcpy(segWeight, alpha, sizeof(double) * 2);
      }
    }
  };

  ExactOctreeBase::rangeQuery(toBB, toEle);
  return segWeight[0] < DBL_MAX;
}

vector<int> ExactTriMeshOctree::trianglesUnderNode(int nodeID) const
{
  vector<int> ret;
  vector<int> candidates = { nodeID };
  int candBegin = 0, candEnd = 1;
  while(candBegin < candEnd)
  {
    for(int i = candBegin; i < candEnd; i++)
    {
      int ID = candidates[i];
      if (nodes[ID].isLeaf())
      {
        ret.insert(ret.end(), nodes[ID].indices.begin(), nodes[ID].indices.end());
      }
      else
      {
        for(int childID : nodes[ID].childrenIDs)
        {
          if (childID < 0)
            continue;
          candidates.push_back(childID);
        }
      }
    }
    candBegin = candEnd;
    candEnd = int(candidates.size());
  }

  sortAndDeduplicate(ret);
  return ret;
}


int ExactTriMeshOctree::getClosestTriangle(const TriMeshRef triMesh, const Vec3d & queryPosition, int & feature) const
{
  assert(triMesh.numVertices() == numVertices && triMesh.numTriangles() == numTriangles);

  int closestTriID = -1;
  int closestFeature = -1;

  auto toBB = [&](const BoundingBox & bb) -> pair<double, double>
  {
    double dist2 = bb.distanceToPoint2(queryPosition);
    double furDist2 = bb.furthestDistanceToPoint2(queryPosition);
    return { dist2, furDist2 };
  };


  auto toElement = [&](int triID, double minDistance) -> double
  {
    const Vec3d & v0 = triMesh.pos(triID, 0);
    const Vec3d & v1 = triMesh.pos(triID, 1);
    const Vec3d & v2 = triMesh.pos(triID, 2);

    int feature = -1;
    double dist2 = getSquaredDistanceToTriangle(queryPosition, v0, v1, v2, feature);

    if (dist2 <= minDistance)
    {
      closestTriID = triID;
      closestFeature = feature;
    }
    return dist2;
  };

  double minDist = DBL_MAX;
  nearestQuery(toBB, toElement, minDist);

  assert(closestTriID >= 0);
  assert(closestFeature >= 0);
  feature = closestFeature;
  return closestTriID;
}

bool ExactTriMeshOctree::sanityCheck(const TriMeshRef triMesh) const
{
  assert(triMesh.numVertices() == numVertices && triMesh.numTriangles() == numTriangles);

  for(size_t i = 0; i < nodes.size(); i++)
  {
    const BoundingBox & bb = nodes[i].bb;
    for(size_t j = 0; j < nodes[i].childrenIDs.size(); j++)
    {
      int childID = nodes[i].childrenIDs[j];
      if (childID < 0) continue;
      const BoundingBox & childBB = nodes[childID].bb;
      if (bb.checkInside(childBB) == false)
      {
        cout << "Error: boundingBox ID " << childID << " is not inside its parent: ID " << i << endl;
        bb.print();
        childBB.print();
        return false;
      }
    }
    // we don't check whether triangles are inside the box because our construction
    // allow triangles to be partially outside as long as the outside part are covered by other boxes
  }
  return true;
}


/////////////////////////////////////////////////////////////
//                   ExactTetMeshOctree
/////////////////////////////////////////////////////////////

void ExactTetMeshOctree::clear()
{
  nodes.clear();
  tetBBs.clear();
}

void ExactTetMeshOctree::build(const TetMeshRef & tetMesh, int maxDepth, int maxNumTetsPerNode)
{
  nodes.clear();
  numVertices = tetMesh.numVertices();
  numTets = tetMesh.numTets();
  tetBBs = tetMesh.getTetBoundingBoxes();
  vector<int> rootIndices(numTets);
  iota(rootIndices.begin(), rootIndices.end(), 0);

  nodes.emplace_back(0, move(rootIndices));
  nodes[0].bb = BoundingBox(tetMesh.numVertices(), tetMesh.positions());

  auto bbFromTris = [&](const vector<int> & tetList) -> BoundingBox
  {
    // get the bounding box of all triangles inside this childNode
    vector<int> allVtxIDsInThisNode;
    tetMesh.getVerticesInTets(tetList, allVtxIDsInThisNode);
    return BoundingBox (tetMesh.positions(), allVtxIDsInThisNode);
  };
  auto subdivide = [&](const vector<int> & elementList, const BoundingBox & bb,
        vector<int> childIDs[8], BoundingBox subBBs[8]) -> void
  {
    boundingBoxPartition(elementList, bb, childIDs, subBBs, tetBBs, bbFromTris);
  };

  buildFromRoot(subdivide, maxDepth, maxNumTetsPerNode);
}

//int ExactTetMeshOctree::getClosestTet(const Vec3d & queryPosition, double * minDistance, double distanceHi) const
//{
//  auto toBB = [&](const BoundingBox & bb) -> pair<double, double>
//  {
//    double dist = bb.distanceToPoint(queryPosition);
//    return { dist, dist + bb.diameter() };
//  };
//
//  auto toElement = [&](int childIndex, double minDistance) -> double
//  {
//    pointInTet(&queryPosition[0], )
//    Vec3d triangleVertices[3];
//    for(int vtx=0; vtx<3; vtx++)
//    {
//      int vertex = triangles[3*childIndex+vtx];
//      for(int dim=0; dim<3; dim++)
//        triangleVertices[vtx][dim] = vertices[3*vertex+dim];
//    }
//
//    TriangleWithCollisionInfo triangle(triangleVertices[0], triangleVertices[1], triangleVertices[2]);
//
//    int feature;
//    double alpha, beta, gamma;
//    Vec3d queryPosV(queryPos);
//    double distToTriangle2 = triangle.distanceToPoint2(queryPosV, &feature, &alpha, &beta, &gamma);
//
//    //printf("child %d: exactd=%G\n", ch, distToBV);
//    if (distToTriangle2 <= minDistance * minDistance)
//    {
//      *nearestTriangle = childIndex;
//      *nearestFeature = feature;
//      *nearestAlpha = alpha;
//      *nearestBeta = beta;
//      *nearestGamma = gamma;
//      return sqrt(distToTriangle2);
//    }
//    return DBL_MAX;
//  };
//
//  double minDist = 0.0;
//  nearestQuery(toBB, toElement, minDist, distanceHi);
//  if (minDistance) *minDistance = minDist;
//}
