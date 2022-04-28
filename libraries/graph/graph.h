/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "graph" library , Copyright (C) 2018 USC                              *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
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

#ifndef _GRAPH_H_
#define _GRAPH_H_

#include <vector>
#include <set>
#include <map>
#include "sparseMatrix.h"

/*
  A class to store an undirected graph (nodes connected with edges).
*/

class Graph
{
public:

  Graph(); 
  // if sortEdgeVertices=1, each edge (v0, v1) in the file will be sorted to ensure v0 < v1 
  // before added into internal data
  Graph(const char * filename, int sortEdgeVertices=1); // load graph from file
  // each edge is given by two integers; length of "edges" should be 2xnumEdges
  // if sortEdgeVertices=1, each edge (v0, v1) in edges will be sorted to ensure v0 < v1 
  // before added into internal data
  Graph(int numVertices, int numEdges, const int * edges, int sortEdgeVertices=1);

  // convert a matrix into a graph; numVertices = matrix->GetNumRows()
  // two vtx (v0, v1) share an edge if they have sparse entries at matrix(v0,v1) or matrix(v1,v0)
  Graph(const SparseMatrix * matrix);

  Graph(const Graph & graph);
  Graph & operator=(const Graph & graph);
  virtual ~Graph();

  void Save(const char * filename) const; // save graph to file

  int GetNumVertices() const;
  int GetNumEdges() const;

  int GetNumNeighbors(int vertex) const;
  int GetNeighbor(int vertex, int i) const;
  // return 0 if vtx1 and vtx2 are not neighbors
  // return neighbor index + 1 (range: [1, #neighbor]) otherwise
  int IsNeighbor(int vtx1, int vtx2) const;
  // neighborhoodSize >= 0
  // find all vertices whose distance to vertex <= neighborhoodSize
  // return their vertex ID and distance
  std::map<int,int> GetNeighborhoodWithDistance(int vertex, int neighborhoodSize);
  std::map<int,int> GetNeighborhoodWithDistance(const std::set<int> & seedVertices, int neighborhoodSize);

  int GetMinDegree() const;
  int GetMaxDegree() const;
  double GetAvgDegree() const;
  double GetStdevDegree() const;

  void ExpandNeighbors(); // connects every node to all the neighbors of every neighbor
  void PrintInfo() const;

  // if scaleRows == 1, each row will be scaled to sum to one
  void GetLaplacian(SparseMatrix ** L, int scaleRows=0) const;

  // returns the Cartesian graph product of "this" and graph2
  Graph * CartesianProduct(Graph & graph2) const;
  // return the index of vertex (vertex1, vertex2) in the cartesian product of "this" with another graph (which is not needed explicitly)
  int GetCartesianProductVertexIndex(int vertex1, int vertex2) const;
  // converts in the opposite direction
  void GetCartesianProductVertexIndexComponents(int productVertex, int * vertex1, int * vertex2) const;

  // clusters given vertices into connected components
  // previous data stored in clusters will be cleared
  void Cluster(const std::set<int> & vertices, std::vector<std::set<int> > & clusters) const;

  // computes the shortest distance from the given seed vertices (distance of zero) to all the graph vertices
  // input: seed vertices
  // output: distance to the set of seed vertices, for each mesh vertex
  void ShortestDistance(const std::set<int> & seedVertices, std::vector<int> & distances) const;

  // computes one of the shortest path from seedVertices to one of the destinationVertices
  // return whether the path exists.
  // if path != NULL, previous data stored in path will be cleared. On return, it stores the path from
  // one of the seed vertices to one of the destination vertices, both ends included
  bool FindShortestPath(const std::set<int> & seedVertices, const std::set<int> & destinationVertices, std::vector<int> * path = NULL) const;

  // find shortest path only in the local vertices set
  bool FindLocalShortestPath(const std::set<int> & localVertices, const std::set<int> & seedVertices, const std::set<int> & destinationVertices, std::vector<int> * path = NULL) const;

  // find the connected component of vtx
  void GetConnectedComponent(int vtx, std::set<int> & connectedVertices) const;

protected:
  int numVertices, numEdges; // num vertices, num edges
  std::set< std::pair<int, int> > edges;
  // for each vtx, mapping: vtx index of neighbor [0, numVertices) -> neighbor index [0, #neighbor)
  std::vector< std::map<int, int> > vertexNeighbors;
  std::vector< std::vector<int> > vertexNeighborsVector;

  void BuildVertexNeighbors();
  void BuildVertexNeighborsVector();
};

inline int Graph::GetNumVertices() const
{
  return numVertices;
}

inline int Graph::GetNumEdges() const
{
  return numEdges;
}

inline int Graph::GetNumNeighbors(int vertex) const
{
  return (int) vertexNeighborsVector[vertex].size();
}

inline int Graph::GetNeighbor(int vertex, int i) const
{
  return vertexNeighborsVector[vertex][i];
}

#endif

