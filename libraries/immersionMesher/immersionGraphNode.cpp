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

#include "immersionGraphNode.h"

#include <cassert>
using namespace std;

ostream & operator << (ostream & o, ImmOwnership os)
{
  switch (os)
  {
  case OWNED:
    cout << "OWN"; break;
  case DECLINED:
    cout << "DEC"; break;
  case UNDECIDED:
    cout << "UND"; break;
  }
  return o;
}

ImmersionGraphNode::ImmersionGraphNode(int cellID, int nodeID, const vector<map<int, bool>> & cellPatches) : nodeID(nodeID), cellID(cellID)
{
  for(auto p : cellPatches[cellID])
  {
    int patchID = p.first;
    nbrs[patchID] = -1;
    patchOwnership[patchID] = (p.second ? UNDECIDED : DECLINED);
  }
}

bool ImmersionGraphNode::hasNbrAtPatch(int patchID) const
{
  auto it = nbrs.find(patchID);
  assert(it != nbrs.end());
  return (it->second >= 0);
}

int ImmersionGraphNode::getNbrIDAtPatch(int patchID) const
{
  auto it = nbrs.find(patchID);
  assert(it != nbrs.end());
  return it->second;
}

void ImmersionGraphNode::setNbrIDAtPatch(int patchID, int nbrID)
{
  assert(nbrs.find(patchID) != nbrs.end());
  assert(nbrs[patchID] == -1);
  nbrs[patchID] = nbrID;
}

bool ImmersionGraphNode::hasNbrAtPatchAndNotNode(int patchID, int nodeID) const
{
  auto it = nbrs.find(patchID);
  assert(it != nbrs.end());
  return (it->second >= 0 && it->second != nodeID);
}

ImmOwnership ImmersionGraphNode::getPatchOwnership(int patchID) const
{
  assert(nbrs.find(patchID) != nbrs.end());
  auto it = patchOwnership.find(patchID);
  assert(it != patchOwnership.end());
  return patchOwnership.find(patchID)->second;
}

void ImmersionGraphNode::setPatchOwnership(int patchID, ImmOwnership o)
{
  assert(o != UNDECIDED);
  assert(nbrs.find(patchID) != nbrs.end());
  assert(patchOwnership[patchID] == UNDECIDED || patchOwnership[patchID] == o);
  patchOwnership[patchID] = o;
}

// cellPatchBouNbrs: cellID -> patchID -> bouID -> <nbr bouID, nbr patchID >
bool ImmersionGraphNode::checkPatchValid(const vector<map<int, map<int , pair<int, int>>>> & cellPatchBouNbrs, bool verbose) const
{
  for(auto p : patchOwnership) // for each patch
  {
    int patchID = p.first;
    ImmOwnership o = p.second;

    auto patchNbrIt = cellPatchBouNbrs[cellID].find(patchID);
    if (patchNbrIt == cellPatchBouNbrs[cellID].end()) continue; // if this cell has only one B-patch, then this B-patch has no geo nbrs
    const auto & nbrPatchMap = patchNbrIt->second; // bouID -> <nbr bouID, nbr patchID >
    for(auto p2: nbrPatchMap)
    {
      int nbrPatchID = p2.second.second;
      if (nbrPatchID == patchID) continue; // if patch self-intersect
      if (neighboringOwernshipValid(o, getPatchOwnership(nbrPatchID)) == false)
      {
        if (verbose)
        {
          cout << "checkPatchValid fails on patch " << patchID << " at a node " << nodeID << " at cell " << cellID << endl;
          cout << "nbr patch " << nbrPatchID << endl;
          print();
        }
        return false;
      }
    }
  }
  if (verbose)
  {
    cout << "checkPatchValid finishes" << endl;
  }
  return true;
}

void ImmersionGraphNode::print() const
{
  cout << "node " << nodeID << " at cell " << cellID << " ";
  for(auto p : nbrs)
  {
    cout << p.first << "->(" << getPatchOwnership(p.first) << " " << p.second << ") ";
  }
  cout << endl;
}

void printNodes(const vector<ImmersionGraphNode> & nodes)
{
  for(auto & node : nodes)
  {
    node.print();
  }
}

bool ImmersionGraphNode::operator == (const ImmersionGraphNode & node2) const
{
  if (nodeID != node2.nodeID) return false;
  if (cellID != node2.cellID) return false;
  if (nbrs != node2.nbrs) return false;
  if (patchOwnership != node2.patchOwnership) return false;
  return true;
}

