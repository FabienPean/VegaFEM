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

#ifndef IMMERSIONMESHER_H
#define IMMERSIONMESHER_H

#include "triKey.h"
#include "tetMeshGeo.h"
#include "triMeshGeo.h"
#include "triMeshNeighbor.h"
#include "barycentricCoordinates.h"
#include "immersionGraphNode.h"
#include "iglRemeshSelfIntersection.h"
#include "profiler.h"
#include <iostream>
#include <vector>
#include <map>
#include <array>

/*
  This class implements:

  Yijing Li, Jernej Barbic: Immersion of Self-Intersecting Solids and Surfaces,
  ACM Transactions on Graphics 37(4) (SIGGRAPH 2018), Vancouver, Canada, 2018.

  It can compute a tet mesh with properly duplicated tets to embed
  a self-intersecting triangle mesh.
*/

class ImmersionMesher
{
public:

  ImmersionMesher();

  // Set verbosity.
  void setVerbose(bool verbose) { this->verbose = verbose; }

  // Use CSG to generate virtual tets, as opposed to the virtual tets algorithm of our paper (default: false, i.e., use our paper).
  void useCSGForVirtualTets(bool useCSGInVT) { this->useCSGInVT = useCSGInVT; }

  // The main routine. It computes an ``unglued'' tet mesh by properly duplicating the tets of the input tet mesh, and embeds
  // the input triangle mesh into it (please see paper).
  // Input:
  //   triangle mesh (manifold and potentially self-intersecting)
  //   tetMesh covering the space occupied by the triangle mesh
  // Output:
  //   Note: Each output is an std::vector because there may be multiple tet mesh outputs, due to ambiguities (please see paper).
  //   outputTetMeshes: output tet meshes
  //   embeddingWeights: barycentric coordinates of input triangles embedded into each output tet mesh
  //   cellMeshes (optional; for visualization purposes): The triangle mesh of each cell. Each std::vector element is the triangle mesh of one cell.
  //   cellMeshWeights (optional; for visualization purposes): Barycentric coordinates to embed cellMeshes into the output tet mesh. Each std::vector element gives the barycentric coordinates of the embedding of the triangle mesh of one cell into the output tet mesh.
  void run(const TriMeshGeo & triMesh,
           const TetMeshGeo & tetMesh,
           std::vector<TetMeshGeo> & outputTetMeshes,
           std::vector<BarycentricCoordinates> & embeddingWeights,
           std::vector<std::vector<TriMeshGeo>> * cellMeshes = nullptr,
           std::vector<std::vector<BarycentricCoordinates>> * cellMeshWeights = nullptr);

  // Self-intersection divides the space into several R^3 cells (see paper).
  int getNumCells() const { return numCells; }
  // Patches are subsets of the input triangle surface that form the cells' boundary (see paper).
  int getNumPatches() const { return numPatches; }
  // Arcs are common boundaries between patches (see paper).
  int getNumArcs() const { return numArcs; }

  // To obtain profiling data after calling "run".
  Profiler & getProfiler() { return profiler; }

  // The mesh obtained after the input mesh self-cuts.
  const TriMeshRef getCutMesh() const { return selfCutMesh.cutMesh; }

  // patchID -> triangles in cutVtxIDs of this patch
  const std::vector<std::vector<Vec3i>> & getPatchTriangles() const { return patchTris; }

  // patchID -> nbring bouID -> topological nbring patchID sharing the bou
  const std::vector<std::map<int, int>> & getPatchNbrPatches() const { return patchBouNbrs; }

  // arcID -> nbring patchIDs
  const std::vector<std::set<int>> & getArcNbrPatches() const { return arcNbrPatches; }

  // cellID -> triangles in cutVtxIDs of this cell,
  // triangle orientations are adjusted so that their normals are pointing outward from the cell
  const std::vector<std::vector<Vec3i>> & getCellTriangles() const { return cellTris; }

  // arcID -> stitchID UEdgeKeys on the arc
  const std::vector<std::set<UEdgeKey>> & getArcs() const { return arcs; }

  // stitchID -> position in world space of this stitchID
  const std::vector<Vec3d> & getStitchPositions() const { return stitchPositions; }

protected:

  void buildBasicData();
  void prepareDataForSolve();
  void buildCellSurfaceMeshes();
  void runNodeSearchMethod(std::vector<std::vector<ImmersionGraphNode>> &);
  void buildCellTetMeshes();
  bool tryUpdate(ImmStackEntry * entry, int nodeA, bool nodeAlreadyChagned, int patchIDToOwn = -1, int patchIDToDecline = -1);
  bool tryGlobalUpdate(ImmStackEntry * entry);
  bool tryConnect(ImmStackEntry * entry, int nodeA, int nodeB, int patchID);
  int addSeed(ImmStackEntry * entry);
  std::tuple<int,int,int> getIncompletePatchHeuristically(const ImmStackEntry * entry,
      int & numOpenPatchChecked, bool & hasAmbiguityInSearch, std::set<std::tuple<int,int,int>> * triedDir = nullptr);
  std::tuple<int,int,int> getAvailableUndecidedPatch(const ImmStackEntry * entry, std::set<std::tuple<int,int,int>> * triedDir = nullptr);
  void printGraphInfo(const ImmStackEntry * curEntry) const;
  void produceFinalTetMesh(std::vector<ImmersionGraphNode> & cellNodes, int graphID,
      TetMeshGeo & outputTetMesh, BarycentricCoordinates & outputInterpWeight,
      std::vector<TriMeshGeo> * allCellMesh, std::vector<BarycentricCoordinates> * allCellInerpWeight);

  // The following three function implement the three stages of our work. They must be called one after another.
  // Call this function first.
  void buildTopologyData(const TriMeshGeo & triMesh, const TetMeshGeo & tetMesh);
  // Then, call this function.
  void runImmersionAlgorithm();
  // Finally, call this function, which gives you the results.
  void generateImmersedTetMesh(std::vector<TetMeshGeo> & tetMeshes, std::vector<BarycentricCoordinates> & embeddingWeights,
      std::vector<std::vector<TriMeshGeo>> * allCellMeshes = nullptr, std::vector<std::vector<BarycentricCoordinates>> * allCellMeshWeights = nullptr);

  // If two triangles intersect in a normal way (no degeneracy),
  // in world space where the triangles intersect, there is one line segement l_w in R^3 where the intersection happens.
  // In abstract manifold space, there are two line segements l_t{0,1} that correspond to the intersection on the two triangles, respectively.
  // When mapping the texture space to the world space, the two line segments l_t{0,1} are mapped to the same l_w.
  // On the entire mesh, those l_w form curves where self-intersection happens.
  // These curves are the boundaries of the patches. We call them arcs.
  // "Bous" are the arcs represented in the abstract manifold space (see paper).
  // In most cases two bous are mapped to one arc, when mapping from abstract manifold space to world space.
  int getNumBous() const { return numBous; }
  // Note that after cutting by self-intersection, the input triangle mesh becomes a new mesh called 'cutMesh'.
  // The vertices (named by cutVtxID) of the cutMesh are built in such a way that all cutVtx triangles do not intersect with
  // each other except on the edges.
  // Also, if two bous l_t0 and l_t1 correspond to the same arc, the cutVtxIDs on l_t0 and l_t1 are DIFFERENT.
  // In this way, we are able to differentiate between l_t0 and l_t1.
  // To visualize arcs, we have another vertex ID called 'stitchID', where two cutVtxIDs that correspond to the same point on an arc
  // have the same stitchID. In this way, we can count and visualize arcs.

  // bouID -> cutVtxID UEdgeKeys on the bou
  const std::vector<std::set<UEdgeKey>> & getBous() const { return bous; }

  Profiler profiler;

  TetMeshGeo inputTetMeshGeo;
  TriMeshGeo inputTriMeshGeo;

  bool cutMeshSaved = false;
  bool useCSGInVT = false;
  bool verbose = false;

  // output
  iglInterface::SelfCutMeshData selfCutMesh;

  int numCells = 0;
  int numPatches = 0;
  int numArcs = 0;
  int numBous = 0;

  std::vector<Vec3d> stitchPositions; // stitchID -> pos
  TriMeshNeighbor cutMeshNbr;
  std::vector<std::vector<Vec3i>> patchTris; // patchID -> cut triangles in cutVtxIDs of this patch
  std::vector<std::vector<int>> patchTriIDs; // patchID -> cutTriIDs of the patch
  std::map<UTriKey, int> cutMeshTriKeyIDs;   // cut UTriKey -> cutTriID of the patch

  std::vector<std::vector<Vec3i>> cellTris; // cellID -> ourward cut triangles in cutVtxIDs of this cell
  std::vector<std::map<int, bool>> cellPatches; // cellID -> patchID -> whether this patch points outward for this cell
  std::vector<int> cellWindingNumbers; // cellID -> windingNumber
  std::vector<std::map<int, int>> cellNeighborsAtPatch; // cellID -> patchID -> nbr cellID

  std::vector<std::map<int, int>> patchBouNbrs; // patchID -> bouID -> topological nbring patchID sharing the bou
  std::vector<std::set<UEdgeKey>> bous; // bouID -> cutVtxID UEdgeKeys on the bou
  std::vector<std::set<int>> bouVertices;    // bouID -> vtx on the bou
  std::vector<std::array<int,2>> bouPatchIDs; // bouID -> two patchIDs

  std::vector<std::set<UEdgeKey>> arcs; // arcID -> stitchID UEdgeKeys on the arc
  std::vector<std::set<int>> arc2Bous; // arcID -> bouID
  std::vector<int> bou2Arcs; // bouID -> arcID

  std::vector<std::set<int>> patchNbrArcs; // patchID -> nbring arcIDs
  std::vector<std::set<int>> arcNbrPatches; // arcID -> nbring patchIDs;
  std::vector<std::map<int, std::set<int>>> patchNbrArcBous; // patchID -> arcID as part of the patch boundary -> bouID as the patch boundary

  // On a cell, boundary patches (B-patches) form a closed "manifold" shape.
  // we wish to find which B-patch neighbors to which B-patch on a cell, so we use the following two vars:
  // cellID -> patchID -> bouID -> <nbr bouID, nbr patchID >
  std::vector<std::map<int, std::map<int , std::pair<int, int>>>> cellPatchBouNbrs;
  // When two B-patches neighbors on a cell, they share at least one arc.
  // on this arc, the bou vertices from each patch corresponds to each other:
  // cellID -> <unordered bouID pair> -> vector of vtxIDs on both bous that are correspondent
  std::vector<std::map<UEdgeKey, std::vector<std::pair<int, int>>>> cellBouCorrespondence;

  std::vector<std::map<int, int>> cellPatchStart; // cellID -> patchID -> the start of the triangles from this patch in this cell mesh (cellTriID)
  std::vector<TriMeshGeo> manifoldCellMesh; // cellID -> the manifold cell surface mesh
  std::vector<std::vector<std::set<int>>> manifoldCellMeshOriVtxIDs; // cellID -> mnCellVtxID -> cutVtxID
  std::vector<std::map<int,std::map<int, int>>> manifoldCellMeshPatchOri2NewVtxIDMap; // cellID -> patchID -> cutVtxID -> mnCellVtxID
  std::vector<std::vector<int>> manifoldCellOriTriIDs; // cellID -> cellTriID -> cutTriID
  std::vector<std::vector<Vec3i>> manifoldCellOriTris; // cellID -> cellTriID -> cut triangles in cutVtxIDs (this could be removed)
  std::vector<std::map<int, std::pair<int, bool>>> manifoldCellOri2NewTriIDMap; // cellID -> cutTriID -> <cellTriID, orientation>

  std::vector<TetMeshGeo> cellTetMeshes;
  std::vector<BarycentricCoordinates> cellTetMeshInterps; // embed: cellTetMeshes -> manifoldCellMesh
  std::vector<std::vector<std::vector<int>>> cellTetMeshTri2TetIDs; // cellID -> cellTriID -> cellTetIDs intersecting the triangle
  std::vector<std::vector<int>> cellTetMeshOriTetIDs;               // cellID -> cellTetID -> inputTetID
  // When creating a manifold cell mesh, we need to know which part in the final cell mesh is from which patch.

  std::vector<std::set<int>> cellID2NodeIDs; // cellID -> nodeIDs of this cellID
  std::vector<std::vector<ImmersionGraphNode>> computedCellNodes;
  bool debugImmersion = false;
};

#endif

