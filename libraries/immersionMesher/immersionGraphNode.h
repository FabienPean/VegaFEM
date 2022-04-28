/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "immersionMesher" library , Copyright (C) 2018 USC                    *
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

#ifndef IMMERSIONGRAPHNODE_H
#define IMMERSIONGRAPHNODE_H

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <tuple>

/*
  Implements a node in our immersion graph.
*/

// The ownership of a node to a B-patch.
enum ImmOwnership
{
  OWNED,
  DECLINED,
  UNDECIDED
};

std::ostream & operator << (std::ostream & o, ImmOwnership os);

// whether the owerships of the two geometric neighboring B-patches on a node is compatible by Rule 5
inline bool neighboringOwernshipValid(ImmOwnership o1, ImmOwnership o2)
{
  if (o1 == OWNED && o2 == OWNED) return false;
  return true;
}

class ImmersionGraphNode
{
public:
  // cellPatches is the var cellPatches stored in class ImmersionMeshing,
  // it defines a mapping: cellID -> patchID -> whether this patch points outward for this cell
  // the constructor initializes nbrs to be nbrs[patchID] = -1, where -1 means no neighbors
  // patchOwnership is initialized to: patchOwnership[patchID] = DECLINED, if the orientation of the B-patch and the node don't agree
  //                                                             UNDECIDED, otherwise
  ImmersionGraphNode(int cellID, int nodeID, const std::vector<std::map<int, bool>> & cellPatches);

  int getNodeID() const { return nodeID; }
  int getCellID() const { return cellID; }
  const std::map<int, int> & getNbrs() const { return nbrs; }
  const std::map<int, ImmOwnership> & getOwnerships() const { return patchOwnership; }

  // whether the node has connected to another node across patchID
  bool hasNbrAtPatch(int patchID) const;

  // get the node ID of the neighboring node across patchID
  // return -1 if no neighboring node connected
  int getNbrIDAtPatch(int patchID) const;

  // set a neighboring node with nbrID to this node across patchID
  void setNbrIDAtPatch(int patchID, int nbrID);

  // return true if and only if there is a neighboring node across patchID but this neighbor is not nodeID
  bool hasNbrAtPatchAndNotNode(int patchID, int nodeID) const;

  // get the ownership at patchID
  ImmOwnership getPatchOwnership(int patchID) const;

  // set the ownership at patchID to o
  void setPatchOwnership(int patchID, ImmOwnership o);

  // integrity check using Rule 5 on the node's B-patches
  // input cellPatchBouNbrs is a member var of class ImmersionMeshing,
  // it defines a mapping: cellID -> patchID -> bouID -> <nbr bouID, nbr patchID >
  bool checkPatchValid(const std::vector<std::map<int, std::map<int , std::pair<int, int>>>> & cellPatchBouNbrs, bool verbose = false) const;

  // print patch ownership and neighboring nodes
  void print() const;

  bool operator == (const ImmersionGraphNode & node2) const;
  bool operator != (const ImmersionGraphNode & node2) const { return ! (*this == node2); }

protected:
  int nodeID = -1;
  int cellID = -1; // the cell this node belongs to
  std::map<int, int> nbrs;                    // patchID -> nbr nodeID, in the constructor, we initialize: nbrs[patchID] = -1
  std::map<int, ImmOwnership> patchOwnership; // patchID -> ownership
};

// print a vector of nodes for debugging
void printNodes(const std::vector<ImmersionGraphNode> & nodes);

// When ambiguity shows up during the search for an immersion,
// we store relevant data onto a stack and pick a potential direction to continue searching.
// This structure is used to hold the data of our immersion algorithm.
struct ImmStackEntry
{
  std::vector<ImmersionGraphNode> nodes;               // the nodes vector, storing all the nodes
  std::vector<int> patchOwnerID;            // patchID -> nodeID which owns this patch, -1 if no node owns it
  std::vector<std::set<int>> usedNodes;     // cellID -> nodeIDs, stores the current nodes that are visited by the search algorithm
  std::vector<std::set<int>> unusedNodes;   // cellID -> nodeIDs, stores those nodes not visited yet
  std::vector<int> nodeOrder;               // record the order of nodes added into the graph, used only for illustrating the process of the algorithm
  std::set<std::tuple<int,int,int>> triedDirections; // the search direction already tried at the state when this data is pushed to the stack
};

#endif

