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

#include "containerHelper.h"
#include "basicAlgorithms.h"
#include "tetKey.h"
#include "vec4d.h"
#include "listIO.h"
#include <cassert>
#include "immersionMesher.h"
#include "immersionGraphNode.h"
using namespace std;

// update patch owership and connection of nodeA based on some prior changes, according to several rules
// optional input: patchIDToAdd: >=0 perform an action of letting nodeA own this patch first
//                               -1  no such action
bool ImmersionMesher::tryUpdate(ImmStackEntry * entry, int nodeA, bool nodeHasChanged, int patchIDToOwn, int patchIDToDecline)
{
  vector<ImmersionGraphNode> & nodes = entry->nodes;
  vector<int> & patchOwnerID = entry->patchOwnerID;
  vector<set<int>> & usedNodes = entry->usedNodes;
  const auto & cellIDsAtPatch = selfCutMesh.cellIDsAtPatch;

  // newlyOwnedPatchIDs stores the patchIDs that will be owned in this function or from input parameter patchIDToOwn
  // newlyDeclinedPatchIDs stores the patchIDs that will be declined in this function or from input parameter patchIDToDecline
  // the usage of these two vars is to propogate update to other nodes who has the same patches
  vector<int> newlyOwnedPatchIDs, newlyDeclinedPatchIDs;
  if (debugImmersion) cout << "updating node " << nodeA << endl;

  // use Rule 5 to check whether B-patches are compatible on nodeA
  if (nodes[nodeA].checkPatchValid(cellPatchBouNbrs) == false)
  {
    cout << "Error: node is broken before calling tryUpdate" << endl;
    return false;
  }

  if (patchIDToOwn >= 0)  // peform the action of nodeA owning patchID
  {
    if (nodes[nodeA].getPatchOwnership(patchIDToOwn) != UNDECIDED) return false;
    assert(patchOwnerID[patchIDToOwn] < 0);
    nodeHasChanged = true;
    nodes[nodeA].setPatchOwnership(patchIDToOwn, OWNED);
    patchOwnerID[patchIDToOwn] = nodeA;
    newlyOwnedPatchIDs.push_back(patchIDToOwn);
  }

  if (patchIDToDecline >= 0)
  {
    if (nodes[nodeA].getPatchOwnership(patchIDToDecline) == OWNED) return false;
    nodes[nodeA].setPatchOwnership(patchIDToDecline, DECLINED);
    nodeHasChanged = true;
    newlyDeclinedPatchIDs.push_back(patchIDToDecline);
  }

  int cellA = nodes[nodeA].getCellID();

  // apply Rule 4 to update patch ownerships on node A
  for(auto p : nodes[nodeA].getNbrs()) // for all nbring node
  {
    int patchID = p.first;
    int nodeB = p.second;
    if (nodeB < 0) continue; // nodeB < 0 means there is no neighboring node across patchID
    int cellB = nodes[nodeB].getCellID();
    const set<int> & arcsAroundPatchID = patchNbrArcs[patchID];

    for(auto p2 : nodes[nodeB].getOwnerships())
    {
      int patchB = p2.first;
      if (patchB == patchID) continue; // nodeA is connected to B across patchID, nothing can be done to patchID now
      ImmOwnership ownershipB = p2.second;
      if (ownershipB == UNDECIDED) continue;

      // Rule 4 only applies on patches sharing the arc with patchID
      // this condition helps on the case of self-touching cells
      if (intersect(patchNbrArcs[patchB], arcsAroundPatchID) == false) continue;

      for(auto q : patchBouNbrs[patchB]) // for each topological neighbor of patchB
      {
        int bouID = q.first;
        int arcID = bou2Arcs[bouID];
        int topoNbrPatchA = q.second; // get one topological neighbor of patchB
        // skip if topoNbrPatchA is not a B-patch of cellA
        if (cellPatches[cellA].find(topoNbrPatchA) == cellPatches[cellA].end()) continue;

        // Rule 4 only applies on patches sharing the arc with patchID
        // this condition helps on the case of self-touching cells
        if (arcsAroundPatchID.find(arcID) == arcsAroundPatchID.end()) continue;

        if (topoNbrPatchA == patchID) continue; // nodeA is connected to B across patchID, nothing can be done to patchID now

        ImmOwnership existingA = nodes[nodeA].getPatchOwnership(topoNbrPatchA);

        if (existingA == ownershipB) continue;

        // now the two owernships are different on node A and B
        if (existingA != UNDECIDED) // Rule 4 does not hold here!
        {
          if (debugImmersion)
          {
            cout << "update fail because nodeA has patch not compatible with a nbr node " << nodeB << " from cell " << cellB << endl;
            cout << "patch on A " << topoNbrPatchA << " " << existingA << endl;
            cout << "patch on B " << patchB << " " << ownershipB << endl;
          }
          return false;
        }
        // then existingA must be UNDECIDE, apply Rule 4 to update patch ownership on A

        if (debugImmersion)
        {
          cout << "set patch " << topoNbrPatchA << " on node " << nodeA << " to " << ownershipB <<
              " because nbring node " << nodeB << " has a topo nbr patch " << patchB << endl;
        }

        nodeHasChanged = true;
        nodes[nodeA].setPatchOwnership(topoNbrPatchA, ownershipB);
        if (ownershipB == OWNED)
        {
          patchOwnerID[topoNbrPatchA] = nodeA;
          newlyOwnedPatchIDs.push_back(topoNbrPatchA);
        }
        else
        {
          assert(ownershipB == DECLINED);
          newlyDeclinedPatchIDs.push_back(topoNbrPatchA);
        }
      }
    }
  } // end apply Rule 4 to update patch ownerships on node A

  if (debugImmersion) cout << "after patch pass, checking whether patches are OK" << endl;

  // check whether patches are OK on nodeA
  if (nodes[nodeA].checkPatchValid(cellPatchBouNbrs, debugImmersion) == false) return false;

  if (debugImmersion)  cout << "checking patches passed" << endl;

  // Apply Rule 5, go through all UNDECIDED patches, decline them if necessary.
  for(auto p : nodes[nodeA].getOwnerships())
  {
    int patchID = p.first;
    ImmOwnership o = p.second;
    if (o != UNDECIDED) continue; // we visit all UNDECIDED patch in this loop

    // if this UNDECIDED patch has been owned by others, then we decline it
    if (patchOwnerID[patchID] >= 0)
    {
      assert(patchOwnerID[patchID] != nodeA);
      if (debugImmersion)
      {
        cout << "set patch " << patchID << " on node " << nodeA << " to DECLINED "
            "because this patch has been owned by others" << endl;
      }
      nodes[nodeA].setPatchOwnership(patchID, DECLINED);
      nodeHasChanged = true;
      // here we don't store this patchID into newlyDeclinedPatchIDs because
      // we want to use newlyDeclinedPatchIDs to propogate update to other nodes
      // by checking whether all nodes have declined this patchID except one node
      // in this case, this node will own this patch
      // since here in this if branch the ower of patchID is already known,
      // we cannot do any propogation on it, so we don't store this patch to newlyDeclinedPatchIDs
      continue;
    }

    // if this UNDECIDED patch has a geometric patch nbr on the node that is owned
    // then we decline the patch (Rule 5)
    bool hasNbrOwned = false;
    for(auto p2: cellPatchBouNbrs[cellA][patchID])
    {
      int nbrPatchID = p2.second.second;
      if (nbrPatchID == patchID) continue; // if patch self-intersect
      if (nodes[nodeA].getPatchOwnership(nbrPatchID) == OWNED)
      {
        hasNbrOwned = true;
        break;
      }
    }
    if (hasNbrOwned)
    {
      nodeHasChanged = true;
      nodes[nodeA].setPatchOwnership(patchID, DECLINED);
      newlyDeclinedPatchIDs.push_back(patchID);
    }
  } // end apply Rule 5, go through all UNDECIDED patches, decline them if necessar

  if (debugImmersion)  cout << "declining patches finished" << endl;

  // check whether patches are OK
  if (nodes[nodeA].checkPatchValid(cellPatchBouNbrs) == false) return false;

  vector<int> otherChangedNodes; // nodes other than nodeA that are changed in the following update
  // apply Rule 7 to decline patches on other nodes
  for(int ownedPatchID : newlyOwnedPatchIDs)
  {
    int cellID = cellIDsAtPatch[ownedPatchID].second; // get the cellID whose orientation agrees with the patch
    assert(cellID == cellA);
    for(int node = 0; node < sizei(nodes); node++)
    {
      if (node == nodeA) continue;
      if (nodes[node].getCellID() != cellID) continue;
      auto o = nodes[node].getPatchOwnership(ownedPatchID);
      if (o == OWNED) return false; // this patch is owned by both node A and another node, which by Rule 7, is not an immersion
      if (o == DECLINED) continue;
      // now o can only be UNDECIDED, we will assign this to be DECLINED according to Rule 7
      nodes[node].setPatchOwnership(ownedPatchID, DECLINED);

      if (debugImmersion)
      {
        cout << "Since patch " << ownedPatchID << " is owned by " << nodeA << ", this node " << node <<
        " has to DECLINE it" << endl;
      }

      if (setFind(usedNodes[cellID], node)) // we leave unusedNodes because they are not connected yet
        otherChangedNodes.push_back(node);
    }
  }

  for(int declinedPatchID : newlyDeclinedPatchIDs)
  {
    if (patchOwnerID[declinedPatchID] >= 0) continue; // already have an owner node, skip
    if (cellIDsAtPatch[declinedPatchID].second != cellA) continue;
    int numNodesUndecided = 0;
    int targetNode = -1;
    for(int nodeID : cellID2NodeIDs[cellA])
    {
      if (nodes[nodeID].getPatchOwnership(declinedPatchID) == UNDECIDED)
      {
        numNodesUndecided++;
        targetNode = nodeID;
      }
    }
//
//    cout << "checking declined patchID " << numNodesUndecided << endl;
//    cout << "Cell " << cellA << endl;
//    cout << streamRange(cellID2NodeIDs[cellA]) << endl;
    if (numNodesUndecided == 0)
    {
//      cout << "target patch  " << declinedPatchID << endl;
//      printNodes(nodes);
      return false; // all nodes decline this patch, not an immersion
    }
    if (numNodesUndecided == 1)
    {
      if (setNotFind(usedNodes[cellA],targetNode)) // if this target node is in unused node
      {
        usedNodes[cellA].insert(targetNode);
        entry->unusedNodes[cellA].erase(targetNode);
        entry->nodeOrder.push_back(targetNode);
      }
      bool success = tryUpdate(entry, targetNode, true, declinedPatchID);
      if (success == false) return false;
    }
  }

  if (nodeHasChanged) // since patches have been changed on A, we should propogate update to nbring nodes, according to Rule 4
  {
    for(auto p : nodes[nodeA].getNbrs())
    {
      int nbrNodeID = p.second;
      if (nbrNodeID < 0) continue; // if not connected across p, skip
      bool success = tryUpdate(entry, nbrNodeID, false);
      if (success == false) return false;
    }
  }

  for(int node : otherChangedNodes)
  {
    bool success = tryUpdate(entry, node, true); // propogate update()
    if (success == false) return false;
  }

  if (debugImmersion) { cout << "updating node " << nodeA << " end" << endl; }
  return true;
}

// Use Rule 4 is to connect nodes if their patches are topological neighbors
// Use Rule 6 to connect nodes around an arc
bool ImmersionMesher::tryGlobalUpdate(ImmStackEntry * entry)
{
  vector<ImmersionGraphNode> & nodes = entry->nodes;
  vector<int> & patchOwnerID = entry->patchOwnerID;
  const vector<set<int>> & usedNodes = entry->usedNodes;
  const auto & cellIDsAtPatch = selfCutMesh.cellIDsAtPatch;
  // try to find any unconnected node pairs that can be connect via Rule 4 and connect them
  for(int patchID = 0; patchID < numPatches; patchID++) // for each owned patch
  {
    int nodeID = patchOwnerID[patchID];
    if (nodeID < 0) continue; // if no node owns patchID, continue
    assert(nodeID >= 0 && nodeID < sizei(nodes));
    int cellID = nodes[nodeID].getCellID();
    for(auto p : patchBouNbrs[patchID]) // for each topological nbr of patchID
    {
      int bouID = p.first;
      int topoNbrPatchID = p.second;
      if (topoNbrPatchID < patchID) continue; // we only visit each unordered pair once
      int nbrNodeID = patchOwnerID[topoNbrPatchID];
      if (nbrNodeID < 0) continue; // skip if no node owns topoNbrNodeID
      int nbrCellID = nodes[nbrNodeID].getCellID();
      if (cellID == nbrCellID) return false; // two nodes from one cell should not own a topo-nbr patch pair
      assert(cellID != nbrCellID);

      // There should be a third patch between the two nodes which own patchID and topoNbrPatchID.
      // By Rule 4, the two nodes connect to each other when the third patch shares the same arc with the other two patches
      // Here we are checking this: geoNbrPatchID is the third patch.
      UEdgeKey cellIDKey(cellID, nbrCellID);
      // cellPatchBouNbrs: cellID -> patchID -> bouID -> <nbr bouID, nbr patchID>
      auto p2 = cellPatchBouNbrs[cellID][patchID][bouID];
      int geoNbrPatchID = p2.second;
      int geoNbrBou = p2.first;
      assert(bou2Arcs[bouID] == bou2Arcs[geoNbrBou]); // assert bouID and geoNbrBou are on the same arc
      assert(mapFind(cellPatchBouNbrs[nbrCellID][topoNbrPatchID], bouID));
      auto p3 = cellPatchBouNbrs[nbrCellID][topoNbrPatchID][bouID];
      if (p3.second != geoNbrPatchID) return false;
      if (p3.first != geoNbrBou) return false;

      auto cellIDs = cellIDsAtPatch[geoNbrPatchID];
      UEdgeKey foundCellIDKey(cellIDs.first, cellIDs.second);
      if (cellIDKey != foundCellIDKey) continue;
      int thisNbrID = nodes[nodeID].getNbrIDAtPatch(geoNbrPatchID); // we get the node connected to nodeID across geoNbrPatchID
      if (thisNbrID == nbrNodeID) continue; // the two already connected
      if (thisNbrID >= 0) return false; // another node connected to nodeID, Rule 4 violated, return false

      if (debugImmersion)
      { cout << "patch pair " << patchID << " " << topoNbrPatchID <<
          " found a patch leading to connection: " << nodeID << " " <<
          nbrNodeID << " patch " << geoNbrPatchID << endl;
      }
      bool success = tryConnect(entry, nodeID, nbrNodeID, geoNbrPatchID);
      if (debugImmersion) { cout << "connection " << (success ? "finished" : "failed") << endl; }
      if (success == false) return false;
    }
  }

  // try to apply Rule 6:
  //    d    q   e
  //         |
  // ---p--------r-----
  //         |              p, q, r, s: patches
  //    d    s   f          c, d, e, f: cells
  for(int arcID = 0; arcID < numArcs; arcID++)
  {
    assert(arc2Bous[arcID].size() == 2);
    vector<int> bous(arc2Bous[arcID].begin(), arc2Bous[arcID].end());

    // get patchIDs p, q, r, s, which all neighbors arcID
    int p = bouPatchIDs[bous[0]][0];
    int q = bouPatchIDs[bous[1]][0];
    int r = bouPatchIDs[bous[0]][1];
    int s = bouPatchIDs[bous[1]][1];

    // get the cellIDs c, d, e, f around the arcs
    int c = cellIDsAtPatch[p].first;
    int d = cellIDsAtPatch[p].second;
    int f = cellIDsAtPatch[r].first;
    int e = cellIDsAtPatch[r].second;
    if (c == 0 || d == 0 || f == 0 || e == 0) continue;
    if (cellWindingNumbers[c] == 0 || cellWindingNumbers[d] == 0 || cellWindingNumbers[e] == 0 || cellWindingNumbers[f] == 0)
      continue;  // there is an air cell among them, skip

    // cellIDsAtPatch: patchID -> <cellID on front side of patch, cellID on back side of patch>
    // adjust c,d,e,f so that they match the picture above
    assert(c != f);
    if ((cellIDsAtPatch[s].first != c && cellIDsAtPatch[s].second != c) ||
        (cellIDsAtPatch[s].first != f && cellIDsAtPatch[s].second != f))
    {
      swap(c, d);
      swap(e, f);
    }

    for(int nodeC : usedNodes[c]) // for each node C on cell c,
    {
      int nodeD = nodes[nodeC].getNbrIDAtPatch(p);
      int nodeF = nodes[nodeC].getNbrIDAtPatch(s);
      // first, let's check whether Rule 6 can be used to connect E and F across r, or connect D and E across q
      // connect on r or q
      if (nodeD >= 0 && nodeF >= 0)
      {
        int nodeE = nodes[nodeD].getNbrIDAtPatch(q);
        if (nodeE >= 0)
        { // check if we can connect E and F across r
          if(nodes[nodeE].getPatchOwnership(r) == DECLINED && nodes[nodeF].getPatchOwnership(r) == DECLINED
              && nodes[nodeE].getNbrIDAtPatch(r) != nodeF)
          {
            // if node E/F already connected to another node not F/E, then Rule 6 is broken
            if (nodes[nodeE].hasNbrAtPatchAndNotNode(r, nodeF)) return false;
            if (nodes[nodeF].hasNbrAtPatchAndNotNode(r, nodeE)) return false;
            if (tryConnect(entry, nodeE, nodeF, r) == false) return false;
            // we return true here because tryConnect will call this tryGlobalUpdate() again recursively
            // so when tryConnect returns, there should be no more arcs to apply Rule 6 to, so we can safely return true
            return true;
          }
        }
        else
        { // check if we can connect D and E across q
          nodeE = nodes[nodeF].getNbrIDAtPatch(r);
          if (nodeE >= 0 && nodes[nodeE].getPatchOwnership(q) == DECLINED && nodes[nodeD].getPatchOwnership(q) == DECLINED)
            // we don't need to check nodes[nodeE].getNbrIDAtPatch(q) != nodeD here because if this condition is true, we
            // will enter the branch above instead
          {
            if (nodes[nodeE].hasNbrAtPatchAndNotNode(q, nodeD)) return false;
            if (nodes[nodeD].hasNbrAtPatchAndNotNode(q, nodeE)) return false;
            if (tryConnect(entry, nodeE, nodeD, q) == false) return false;
            return true;
          }
        }
      }
      else if (nodeF < 0 && nodeD >= 0)
      { // check if we can connect C and F across s
        int nodeE = nodes[nodeD].getNbrIDAtPatch(q);
        if (nodeE >= 0)
        {
          nodeF = nodes[nodeE].getNbrIDAtPatch(r);
          if (nodeF >= 0 && nodes[nodeF].getPatchOwnership(s) == DECLINED && nodes[nodeC].getPatchOwnership(s) == DECLINED)
            // we don't need to check nodes[nodeF].getNbrIDAtPatch(s) != nodeC here because if this condition is true, we
            // will enter the branch: if (nodeD >= 0 && nodeF >= 0)
          {
            if (nodes[nodeC].hasNbrAtPatchAndNotNode(s, nodeF)) return false;
            if (nodes[nodeF].hasNbrAtPatchAndNotNode(s, nodeC)) return false;
            if (tryConnect(entry, nodeC, nodeF, s) == false) return false;
            return true;
          }
        }
      }
      else if (nodeD < 0 && nodeF >= 0)
      { // check if we can connect C and D across p
        int nodeE = nodes[nodeF].getNbrIDAtPatch(r);
        if (nodeE >=0)
        {
          nodeD = nodes[nodeE].getNbrIDAtPatch(q);
          if (nodeD >= 0 && nodes[nodeD].getPatchOwnership(p) == DECLINED && nodes[nodeC].getPatchOwnership(p) == DECLINED)
          {
            if (nodes[nodeC].hasNbrAtPatchAndNotNode(p, nodeD)) return false;
            if (nodes[nodeD].hasNbrAtPatchAndNotNode(p, nodeC)) return false;
            if (tryConnect(entry, nodeC, nodeD, p) == false) return false;
            return true;
          }
        }
      }
    } // end for nodeC
  }

  return true;
}

// connect node A to nodeB across patchID
// return whether it succeeds
bool ImmersionMesher::tryConnect(ImmStackEntry * entry, int nodeA, int nodeB, int patchID)
{
  vector<ImmersionGraphNode> & nodes = entry->nodes;

  if (debugImmersion)
  {
    cout << "node A and B: " << nodeA << " " << nodeB << " connecting"<< endl;
    nodes[nodeA].print();
    nodes[nodeB].print();
    cout << "patch " << patchID << endl;
  }

  nodes[nodeA].setNbrIDAtPatch(patchID, nodeB);
  nodes[nodeB].setNbrIDAtPatch(patchID, nodeA);

  // check ownership of patchID for node A and B
  // by rule 2, they should both decline this patch in order to connect to each other
  auto oA = nodes[nodeA].getPatchOwnership(patchID);
  if(oA == OWNED) return false;
  nodes[nodeA].setPatchOwnership(patchID, DECLINED);
  auto oB = nodes[nodeB].getPatchOwnership(patchID);
  if(oB == OWNED) return false;
  nodes[nodeB].setPatchOwnership(patchID, DECLINED);

  if (debugImmersion) { cout << "updating nodeA..." << endl; }

  if (tryUpdate(entry, nodeA, true, -1, patchID) == false) return false;

  if (debugImmersion) { cout << "updating nodeB..." << endl;}

  if (tryUpdate(entry, nodeB, true, -1, patchID) == false) return false;

  if (debugImmersion) { cout << "finish two nodes update" << endl; }

  if (tryGlobalUpdate(entry) == false) return false;

  return true;
}

// find a suitable seed node from unused nodes to begin the search
// return -1 if no nodes found, otherwise return seed nodeID
int ImmersionMesher::addSeed(ImmStackEntry * entry)
{
  vector<ImmersionGraphNode> & nodes = entry->nodes;
  vector<int> & patchOwnerID = entry->patchOwnerID;
  vector<set<int>> & usedNodes = entry->usedNodes;
  vector<set<int>> & unusedNodes = entry->unusedNodes;

  int numNodes = nodes.size();
  // find a seed
  for(int nodeID = 0; nodeID < numNodes; nodeID++) // for each unused node
  {
    int cellID = nodes[nodeID].getCellID();
    if (setFind(usedNodes[cellID], nodeID)) continue; // if the node is used, then skip
    for(auto p : cellPatches[cellID])
    {
      int patchID = p.first;
      // to find a useful seed, we should find a node which can own a B-patch
      // since new nodes are initialized to have patch owernships either DECLINED or UNDECIDED
      // if we meet a node with a UNDECIDED patch, we will set this patch to OWNED to
      // make this node a seed
      if (nodes[nodeID].getPatchOwnership(patchID) == UNDECIDED) // find a node with a UNDECIDED patch
      {
        assert(patchOwnerID[patchID] == -1);
        usedNodes[cellID].insert(nodeID);                // add this node to usedNodes
        unusedNodes[cellID].erase(nodeID);               // remove this node from unusedNodes

        // update patch owernships to reflect owning patchID
        bool success = tryUpdate(entry, nodeID, true, patchID);
        if (success == false)
        {
          cout << "Error, update fail on the seed node!" << endl;
          return -1;
        }
        entry->nodeOrder.push_back(nodeID);
        return nodeID;
      }
    }
  }
  return -1;
}

// find an incomplete patch
// incomplete patch is defined as a declined patch of a node with no neighboring node connected across the patch
// if there is an incomplete patch, we can try to connect a node across this patch to advance our search for immersion
// return: a tuple of node A, node B and patch p, where A can connect to B across p
//         numIncompletePatchesVisited: #incomplete patches visited in this function
//         hasAmbiguityInSearch: whether all incomplete patches have ambiguity on how to connect another node
// if triedDir == nullptr, then the code uses a heuristic to only advance at the incomplete patch where
// there is no ambigutiy on witch node to connect to.
// if triedDir != nullptr, it stores the search direction we already tried,
// therefore, we shall avoid those search directions, then we just return the first valid search direction
// that has not appeared in triedDir. Note in this case the heuristic does not apply
// Here search direction is defined as a potential action, either as a tuple of <A, B, p> which meanings node A can
// connect to B across patch p
tuple<int,int,int> ImmersionMesher::getIncompletePatchHeuristically(const ImmStackEntry * entry,
    int & numIncompletePatchesVisited, bool & hasAmbiguityInSearch, set<tuple<int,int,int>> * triedDir)
{
  const vector<ImmersionGraphNode> & nodes = entry->nodes;
  const vector<set<int>> & usedNodes = entry->usedNodes;
  const vector<set<int>> & unusedNodes = entry->unusedNodes;
  const int numNodes = nodes.size();

  bool debugOuterLoop = false;
  numIncompletePatchesVisited = 0;
  hasAmbiguityInSearch = false;

  for(int nodeID = 0; nodeID < numNodes; nodeID++) // for each used node
  {
    int cellID = nodes[nodeID].getCellID();
    if (setNotFind(usedNodes[cellID], nodeID)) continue;
    for(auto p : cellPatches[cellID]) // check each B-patch of the node
    {
//        cout << "check open patch at node " << nodeID << " cell " << cellID << " patch " << p.first << endl;
      int patchID = p.first;
      // we look for incomplete patch
      if (nodes[nodeID].getPatchOwnership(patchID) != DECLINED || nodes[nodeID].hasNbrAtPatch(patchID)) continue;

      numIncompletePatchesVisited++;
      // found an incomplete patch:
      // check whether we can connect a new node to it
      assert(cellNeighborsAtPatch[cellID].find(patchID) != cellNeighborsAtPatch[cellID].end());
      int nbrCellID = cellNeighborsAtPatch[cellID][patchID]; // get the neighboring cell ID across this patch
      if (cellWindingNumbers[nbrCellID] == 0) continue; // if nbrCellID is the outer cell (empty space), skip

      bool ambiguityFound = false;
      int targetNbrNodeID = -1; // the target nodeID we should use to connect at the incomplete patch
      // first, we go through all used nodes on this nbrCellID to check whether it can be connected to nodeID
      // without violating any rules
      for(int nbrNodeID : usedNodes[nbrCellID])
      {
        assert(nbrCellID == nodes[nbrNodeID].getCellID());
        // if nbrNodeID already connects to sth. acorss patchID, then skip (Rule 1)
        if (nodes[nbrNodeID].getNbrIDAtPatch(patchID) >= 0) continue;

        auto tmpEntry = *entry;

        if (debugOuterLoop)
        {
          cout << "check open patch at node " << nodeID << " cell " << cellID << " patch " << patchID << " to node " << nbrNodeID <<
              " cell " << nbrCellID << endl;
        }
        // connect nodeID and nbrNodeID across patchID
        bool success = tryConnect(&tmpEntry, nodeID, nbrNodeID, patchID);
        if (debugOuterLoop)
        {
          cout << "result is " << success << endl;
        }
        if (success)
        {
          if (triedDir)  // we should avoid the tried direction
          {
            tuple<int,int,int> dir(nodeID, nbrNodeID, patchID);
            if (setFind(*triedDir, dir)) continue; // this connection was in tried direction
            return dir;
          }
          if (targetNbrNodeID < 0)
          {
            targetNbrNodeID = nbrNodeID; // store the nbrNodeID
          }
          else // targetNbrNodeID >= 0, we already found one valid nbrNodeID at this incomplete patch, ambiguity!
          {
            ambiguityFound = true;
            break; // break the search on nbrNodeID across the patchID
          }
        }
      } // end searching for nbrNodeID across the patchID

      if (triedDir) // we tried all the used nodes at this nbrCellID, but we cannot connect them
      {             // either because the search direction is in triedDir, or the connection is not valid (i.e., breaks rules)
        if(unusedNodes[nbrCellID].size() > 0) // we can still try connecting to an unused node
        {
          int nbrNodeID = *unusedNodes[nbrCellID].begin();
          tuple<int,int,int> dir(nodeID, nbrNodeID, patchID);
          if (setFind(*triedDir, dir)) continue;
          return dir;
        }
        continue;
      }

      // if there is ambiguity on which used node to connect, or ambiguity on whether to connect to one used node or
      // to an unused node, then we apply the heuristic and try the next incomplete patch
      if (ambiguityFound || (targetNbrNodeID >= 0 && unusedNodes[nbrCellID].size() > 0))
      {
        if (verbose)
          cout << "---------------------------- AMB ---------------------------" << endl;
        hasAmbiguityInSearch = true;
        continue; // let 's try the next imcomplete Patch
      }

      if (targetNbrNodeID >= 0) // we find one used no to connect to
      {
        return tuple<int,int,int>(nodeID, targetNbrNodeID, patchID);
      }
      else // targetNbrNodeID < 0 but unusedNodes[nbrCellID].size() > 0,
      {    // we can connect to one unused node
        int nbrNodeID = *unusedNodes[nbrCellID].begin();
        return tuple<int,int,int>(nodeID, nbrNodeID, patchID);
      }
    } // for all patch on node
  } // for all node

  return tuple<int,int,int>(-1,-1,-1);
}

// find one B-patch that is UNDECIDED, skip search directions stored in triedDir
// the format of tuples in the triedDir is: <node A, node B, patch p>
// if node B >= 0, this search direction represents a connection from A to B across p
// else, this search direction represents A owning p
// return the tuple representing a valid search direction
// no heuristic applied here, unlike the function getIncompletePatchHeuristically
tuple<int,int,int> ImmersionMesher::getAvailableUndecidedPatch(const ImmStackEntry * entry, set<tuple<int,int,int>> * triedDir)
{
  const vector<ImmersionGraphNode> & nodes = entry->nodes;
  const vector<set<int>> & usedNodes = entry->usedNodes;

  bool debugOuterLoop = false;
  const int numNodes = nodes.size();

  for(int nodeID = 0; nodeID < numNodes; nodeID++) // for all used nodes
  {
    int cellID = nodes[nodeID].getCellID();
    if (setNotFind(usedNodes[cellID], nodeID)) continue;
    for(auto p : cellPatches[cellID])
    {
//        cout << "check undecided patch at node " << nodeID << " cell " << cellID << " patch " << p.first << endl;
      int patchID = p.first;
      // we are seraching for a UNDECIDED B-patch
      if (nodes[nodeID].getPatchOwnership(patchID) != UNDECIDED) continue;

      // found one!
      // check whether this node can own this patch
      {
        tuple<int,int,int> dir(nodeID, -1, patchID);
        if (setNotFind(*triedDir, dir)) // this direction not tried yet
        {
          auto tmpEntry = *entry;

          assert(tmpEntry.patchOwnerID[patchID] < 0);
          bool success = tryUpdate(&tmpEntry, nodeID, true, patchID);
          success = success && tryGlobalUpdate(&tmpEntry);
          // if we can successfully own this patch, then we have found a good search direction
          if (success) { return dir; }
        }
      }

      // check whether we can connect another node to this patch
      assert(cellNeighborsAtPatch[cellID].find(patchID) != cellNeighborsAtPatch[cellID].end());
      int nbrCellID = cellNeighborsAtPatch[cellID][patchID];
      if (cellWindingNumbers[nbrCellID] == 0) continue; // if nbr cell is the outer space, skip

      for(int nbrNodeID : usedNodes[nbrCellID])
      {
        assert(nbrCellID == nodes[nbrNodeID].getCellID());
        if (nodes[nbrNodeID].getNbrIDAtPatch(patchID) >= 0) continue; // if nbrNodeID has connected to sth. already across patchID
        auto tmpEntry = *entry;

        if (debugOuterLoop)
        {
          cout << "check open patch at node " << nodeID << " cell " << cellID << " patch " << patchID << " to node " << nbrNodeID <<
              " cell " << nbrCellID << endl;
        }
        tuple<int,int,int> dir(nodeID, nbrNodeID, patchID);
        if (setFind(*triedDir, dir)) continue;

        bool success = tryConnect(&tmpEntry, nodeID, nbrNodeID, patchID);
        if (debugOuterLoop) { cout << "result is " << success << endl; }
        if (success) { return dir; }
      }
    } // for all patch on node
  } // for all node

  return tuple<int,int,int>(-1,-1,-1);
}

void ImmersionMesher::printGraphInfo(const ImmStackEntry * curEntry) const
{
  cout << "==================================" << endl;
  printNodes(curEntry->nodes);
  cout << "used nodes: ";
  for(int cellID = 1; cellID < sizei(curEntry->usedNodes); cellID++)
  {
    cout << cellID << "->" << streamRange(curEntry->usedNodes[cellID]) << " ";
  }
  cout << endl;
  cout << "unused nodes: ";
  int numUnused = 0;
  for(int cellID = 1; cellID < sizei(curEntry->unusedNodes); cellID++) { numUnused += curEntry->usedNodes[cellID].size(); }
  cout << numUnused << ", ";
  for(int cellID = 1; cellID < sizei(curEntry->unusedNodes); cellID++)
  {
    cout << cellID << "->" << streamRange(curEntry->usedNodes[cellID]) << " ";
  }
  cout << endl;
  cout << "==================================" << endl;
}

void ImmersionMesher::runNodeSearchMethod(vector<vector<ImmersionGraphNode>> & graphs)
{
//  verbose = true;
  cout << "=======================================================" << endl;
  cout << "          Start to building immersion graphs" << endl;
  cout << "=======================================================" << endl;
  vector<ImmStackEntry> searchStack; // define a stack
  // When ambiguity shows up during the search for an immersion,
  // we store relevant data onto a stack and pick a potential direction to continue searching
  searchStack.emplace_back(); // first, push a stack entry representing the currect search
  ImmStackEntry * curEntry = &searchStack[0];

  curEntry->patchOwnerID.resize(numPatches, -1);

  // initialize usedNodes, unusedNodes, nodes
  cout << "Winding number for each cell: ";

  for(int cellID = 1; cellID < numCells; cellID++) // we skip cellID == 0, because cellID 0 is the outer cell, representing the empty space
  {
    cout << "(cell " << cellID << ": " << cellWindingNumbers[cellID] << ") ";
    if (cellWindingNumbers[cellID] < 0) // winding number < 0 on this cell, there is inversion, not valid input
    {
      cout << "Error: winding number is negative! Immersion not possible!" << endl;
      return;
    }

    // create some nodes at this cell. the number of the nodes is the cell's winding number
    for(int i = 0; i < cellWindingNumbers[cellID]; i++)
    {
      int nodeID = curEntry->nodes.size();
      curEntry->nodes.emplace_back(cellID, nodeID, cellPatches);
    }
  }
  cout << endl;

  curEntry->usedNodes.resize(numCells);
  curEntry->unusedNodes.resize(numCells);
  for(int nodeID = 0; nodeID < sizei(curEntry->nodes); nodeID++)
  {
    int cellID = curEntry->nodes[nodeID].getCellID();
    curEntry->unusedNodes[cellID].insert(nodeID);
  }
  cellID2NodeIDs = curEntry->unusedNodes; // initialize cellID2NodeIDs

  const int numNodes = curEntry->nodes.size();
  cout << "We have " << numNodes << " nodes to connect in the graph." << endl;

  // pick a seed to start
  int seedNode = addSeed(curEntry);
  if (seedNode == -1)
  {
    cout << "Error: cannot find a seed to start! Input not immersion!" << endl;
    return;
  }

  // pop the search stack, return true if success
  auto processRetreat = [&]() -> bool
  {
    if (searchStack.size() > 1)
    {
      cout << "Retreat to the last fork" << endl;
      searchStack.pop_back();
      curEntry = &(searchStack.back());

      cout << "Now the current graph is: " << endl;
      printGraphInfo(curEntry);
      return true;
    }
    else
    {
      cout << "=======================================================" << endl;
      cout << "         Nothing left in stack, search over" << endl;
      cout << "=======================================================" << endl;
      return false;
    }
  };

//  bool debugOuterLoop = false;
  bool hasAmbiguityInSearch = false;
  for(int iter = 0; ; iter++)
  {
    if (verbose)
      cout << "****************iter: " << iter << endl;
//    char tmp = 0;
//    cin >> tmp;
//    cout << "input  tmp is " << tmp << endl;

    bool hasFinished = true; // whether the search for immersion is finished
    hasAmbiguityInSearch = false;

    // if there is a patch not owned by any node, the search is not done yet
    if(find(curEntry->patchOwnerID.begin(), curEntry->patchOwnerID.end(), -1) != curEntry->patchOwnerID.end())
      hasFinished = false;

    int numIncompletePatchesVisited = 0;
    // try to find one incomplete patch and a node to connect across the patch
    // incomplete patch is defined as a declined patch of a node with no neighboring node connected across the patch
    // if there is an incomplete patch, we can try to connect a node across this patch to advance our search for immersion
    // the function getIncompletePatchHeuristically retruns a tuple of node A, node B and patch p, where A can connect to B across p
    //         numIncompletePatchesVisited: #incomplete patches visited in this function
    //         hasAmbiguityInSearch: whether all incomplete patches have ambiguity on how to connect another node
    // Here triedDir == nullptr, the code uses a heuristic to only advance at the incomplete patch where
    // there is no ambigutiy on witch node to connect to.
    // If all incomplete patches have ambiguity, the funtion returns <-1,-1,-1>, representing failure in using the heuristic
    tuple<int,int,int> incompletePatchResult = getIncompletePatchHeuristically(curEntry, numIncompletePatchesVisited, hasAmbiguityInSearch);
    if (numIncompletePatchesVisited > 0) { hasFinished = false; }

    if (hasFinished)
    {
      cout << "=======================================================" << endl;
      cout << "            SUCCESS in finding an immersion" << endl;
      cout << "=======================================================" << endl;
      int graphID = graphs.size();
      graphs.push_back(curEntry->nodes);
      if (verbose)
        ListIO::save(("nodeOrder" + to_string(graphID) + ".txt").c_str(), curEntry->nodeOrder, 0);

      if (processRetreat()) { continue; }
      else break;
    }

    if (get<0>(incompletePatchResult) >= 0) // we have found an incomplete B-patch and it is unambiguous to connect
    {
      int nodeID = get<0>(incompletePatchResult);
      int nbrNodeID = get<1>(incompletePatchResult);
      int nbrCellID = curEntry->nodes[nbrNodeID].getCellID();
      int patchID = get<2>(incompletePatchResult);
      if (verbose)
      {
        cout << "connect node " << nodeID << " cell " << curEntry->nodes[nodeID].getCellID() << " patch " << patchID <<
            " to node " << nbrNodeID << " cell " << curEntry->nodes[nbrNodeID].getCellID() << endl;
      }

      // if the node to connect to, nbrNodeID is in unusedNodes, then we move it from unusedNodes to usedNodes
      if (setFind(curEntry->unusedNodes[nbrCellID], nbrNodeID))
      {
        // grab an unused node
        curEntry->unusedNodes[nbrCellID].erase(nbrNodeID);
        curEntry->usedNodes[nbrCellID].insert(nbrNodeID);
        curEntry->nodeOrder.push_back(nbrNodeID);
        if (verbose)
          cout << "grab an unused node" << endl;
      }

      bool success = tryConnect(curEntry, nodeID, nbrNodeID, patchID);
      if (success)
      {
        // successfully connected to nbrNodeID
        if (verbose)
        {
          cout << "connect node " << nodeID << " cell " << curEntry->nodes[nodeID].getCellID() << " patch " << patchID <<
              " to node " << nbrNodeID << " cell " << curEntry->nodes[nbrNodeID].getCellID() << endl;
          printGraphInfo(curEntry);
        }
        continue; // go the next iteration
      }

      // fail to connect
      cout << "Failed to make this connection!" << endl;
      // If nbrNodeID is in usedNodes in the last iteration, then getIncompletePatchHeuristically() already called
      // tryConnect() to determine that connecting nbrNodeID is fine.
      // So we enter this "else" branch only when nbrNodeID has just been added into usedNodes in this iteration.
      // This means the existing graph failed to add an unused node.
      // Then, the existing graph must have been broken, we should stop this search direction and go the previous one
      // on the stack
      if (processRetreat()) { continue; } // if we succeed in retreating, then go to the next iteration
      else break;                         // else, nothing left in the stack, break the search
    } // end if (get<0>(incompletePatchResult) >= 0)

    // else, no available unabmbiguous, incomplete patch to connect, and the algorithm has not finished

    // if no ambiguity on incomplete B-patches, and we have found at least one incomplete B-patch
    // this then means that on that incomplete B-patch, no nodes can connect to it!
    // then clearly this current search direction is broken, we have to retreat
    if (hasAmbiguityInSearch == false && numIncompletePatchesVisited > 0)
    {
      if (processRetreat()) { continue; } // if we succeed in retreating, then go to the next iteration
      else break;                         // else, nothing left in the stack, break the search
    }

    // then it means we either:
    // a) have another connected component to discover, which we will find by adding a new seed on the component
    // or
    // b) have ambiguity on which node to own/connect to an undecided B-patch, which we will determine by branching and searching
    // or
    // c) have ambiguity on which node to connect to an incomplete B-patch, which we will determine by branching and searching
    if (numIncompletePatchesVisited == 0)
    {
      cout << "finished one connected component, now goes to the next" << endl;
      // first we try case (a), searching for a new connected component in the graph
      ImmStackEntry tmpEntry = *curEntry;
      int newSeed = addSeed(&tmpEntry);
      if (newSeed >= 0) // found one
      {
        *curEntry = tmpEntry;
        continue;
      }
    }

    // if we have ambiguity during the search,
    // or if no incomplete patches found and no new connected componet can be found
    // then we save current state in the stack and branch the search direction
    cout << "=======================================================" << endl;
    cout << " Save current state on stack and begin a new direction" << endl;
    cout << "=======================================================" << endl;
    if (hasAmbiguityInSearch)
      cout << "Because ambiguity found in the search" << endl;
    else
    {
      cout << "Because num incomplete patches visited is 0" << endl;
    }

    { // push the entry to the stack
      // actually in our implementation, the current active state is always in the var searchStack,
      // so technically it's not "pushing" but "copying"...  but you got the idea
      ImmStackEntry newEntry = *curEntry; // copy the current state
      searchStack.emplace_back(move(newEntry));
      curEntry = &(searchStack.back());
    }

    // tried direction saves which branching direction we have tried at this point
    set<tuple<int,int,int>> & triedDir = searchStack[searchStack.size()-2].triedDirections;
    cout << "#Tried directions: " << triedDir.size() << endl;
    tuple<int,int,int> result(-1,-1,-1);
    if (numIncompletePatchesVisited == 0)
    {
      // we are in case (b)
      cout << "All left are undecided patches, try one" << endl;
      result = getAvailableUndecidedPatch(curEntry, &triedDir);
    }
    else
    {
      // case (c)
      result = getIncompletePatchHeuristically(curEntry, numIncompletePatchesVisited, hasAmbiguityInSearch, &triedDir);
    }

    int nodeID = get<0>(result);
    int nbrNodeID = get<1>(result);
    int nbrCellID = curEntry->nodes[nbrNodeID].getCellID();
    int patchID = get<2>(result);
    if (nodeID >= 0) // we find a possible search direction
    {
      triedDir.insert(result); // save this direction to triedDir so that we won't enter again

      if (nbrNodeID >= 0) // the action of the direction is to connect to nbrNodeID
      {
        if (setFind(curEntry->unusedNodes[nbrCellID], nbrNodeID)) // is nbrNodeID is in unusedNodes, move it to usedNodes first
        {
          // grab an unused node
          curEntry->unusedNodes[nbrCellID].erase(nbrNodeID);
          curEntry->usedNodes[nbrCellID].insert(nbrNodeID);
          curEntry->nodeOrder.push_back(nbrNodeID);
          if (verbose)
            cout << "grab an unused node" << endl;
        }

        bool success = tryConnect(curEntry, nodeID, nbrNodeID, patchID);
        if (success == false)
        {
          cout << "Fail to connect to a possible search direction" << endl;
          if (processRetreat()) continue;
          else break;
        }

        if (verbose)
        {
          cout << "connect node " << nodeID << " cell " << curEntry->nodes[nodeID].getCellID() << " patch " << patchID <<
              " to node " << nbrNodeID << " cell " << curEntry->nodes[nbrNodeID].getCellID() << endl;
          printGraphInfo(curEntry);
        }
      }
      else // action is to own a patch
      {
        bool success = tryUpdate(curEntry, nodeID, true, patchID);
        success = success && tryGlobalUpdate(curEntry);
        assert(success); // this action has been tried in getAvailableUndecidedPatch(), so it must be success
        if (verbose)
        {
          cout << "assign node " << nodeID << " cell " << curEntry->nodes[nodeID].getCellID() << " an owned patch " << patchID << endl;
          printGraphInfo(curEntry);
        }
      }
      continue; // finish the action, go to the next iteration
    } // end if (nodeID >= 0) // we find a possible search direction
    else // no possible direction can be found on this state
    {
      cout << "Finish searching available directions at this fork, have to retreat twice" << endl;
      if (processRetreat() && processRetreat()) continue; // we retreat twice because no more possible direction
      else break;
    }
  } // end algorithm loop

  cout << "Finished building immersion graph." << endl;
  if (graphs.size() == 0) { cout << "Found no immersion!" << endl; }
  else { cout << "Found " << graphs.size() << " possible immersion graph" << (graphs.size() == 1 ? "" : "s") << "." << endl; }
}
