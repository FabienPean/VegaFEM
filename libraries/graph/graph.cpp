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

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include "matrixIO.h"
#include "graph.h"
using namespace std;

Graph::Graph() 
{ 
  numVertices = 0; 
  numEdges = 0; 
}

Graph::Graph(const Graph & graph) 
{
  (*this) = graph;
}

Graph & Graph::operator=(const Graph & graph) 
{
  numVertices = graph.numVertices;
  numEdges = graph.numEdges;
  edges = graph.edges;
  vertexNeighbors = graph.vertexNeighbors;
  vertexNeighborsVector = graph.vertexNeighborsVector;
  return *this;
}

Graph::~Graph()
{
}

Graph::Graph(int numVertices_, int numEdges_, const int * edges_, int sortEdgeVertices): numVertices(numVertices_), numEdges(numEdges_)
{
  //printf("num vertices: %d\n", numVertices);
  for(int i=0; i<numEdges; i++)
  {
    //printf("Edge: %d %d\n", edges_[2*i+0], edges_[2*i+1]);
    const int * edge = edges_ + 2*i;
    if (sortEdgeVertices && (edge[0] > edge[1])) // keep the two indices in each edge sorted
      edges.insert(make_pair(edge[1], edge[0])); 
    else
      edges.insert(make_pair(edge[0], edge[1])); 
  }

  BuildVertexNeighbors();
}

Graph::Graph(const char * filename, int sortEdgeVertices)
{
  FILE * fin;
  OpenFile_(filename, &fin, "r");
  int code = fscanf(fin, "%d %d\n", &numVertices, &numEdges);
  if (code != 2)
    throw 1;

  //printf("num vertices: %d\n", numVertices);
  for(int i=0; i<numEdges; i++)
  {
    int vtxA, vtxB;
    int code = fscanf(fin, "%d %d\n", &vtxA, &vtxB);
    if (code != 2)
      throw 2;

    //printf("Edge: %d %d\n", vtxA, vtxB);
    // keep the two indices in each edge sorted
    if (sortEdgeVertices && (vtxA > vtxB))
      edges.insert(make_pair(vtxB, vtxA));
    else
      edges.insert(make_pair(vtxA, vtxB));
  }

  fclose(fin);

  BuildVertexNeighbors();
}

Graph::Graph(const SparseMatrix * matrix)
{
  numVertices = matrix->GetNumRows();
  for(int row = 0; row < matrix->GetNumRows(); row++)
  {
    int rowLen = matrix->GetRowLength(row);
    for(int j = 0; j < rowLen; j++)
    {
      int col = matrix->GetColumnIndex(row, j);
      if (row < col)
        edges.insert(make_pair(row, col));
      else
        edges.insert(make_pair(col, row));
    }
  }
  numEdges = edges.size();

  BuildVertexNeighbors();
}

void Graph::Save(const char * filename) const
{
  FILE * fout;
  OpenFile_(filename, &fout, "w");
  fprintf(fout, "%d %d\n", numVertices, numEdges);

  for(set<pair<int, int> > :: iterator iter = edges.begin(); iter != edges.end(); iter++)
  {
    int vtxA = iter->first;
    int vtxB = iter->second;
    fprintf(fout, "%d %d\n", vtxA, vtxB);
  }

  fclose(fout);
}

void Graph::BuildVertexNeighbors()
{
  vertexNeighbors.assign(numVertices, map<int, int>());

  for(set<pair<int,int> > :: iterator iter = edges.begin(); iter != edges.end(); iter++)
  {
    vertexNeighbors[iter->first].insert(make_pair(iter->second,0)); 
    vertexNeighbors[iter->second].insert(make_pair(iter->first,0)); 
  }

  // number the neighbors
  for(int i=0; i<numVertices; i++)
  {
    int count = 0;
    for(map<int,int> :: iterator iter = vertexNeighbors[i].begin(); iter != vertexNeighbors[i].end(); iter++)
    {
      iter->second = count;
      count++;
    }
  }

  BuildVertexNeighborsVector();
}

void Graph::BuildVertexNeighborsVector()
{
  vertexNeighborsVector.assign(numVertices, vector<int>());

  // create a copy of the data (in a vector), so that can access ith element fast
  for(int i=0; i<numVertices; i++)
    for(map<int,int> :: iterator iter = vertexNeighbors[i].begin(); iter != vertexNeighbors[i].end(); iter++)
      vertexNeighborsVector[i].push_back(iter->first);
}

int Graph::GetMaxDegree() const
{
  int maxDegree = 0;
  for(int vtx=0; vtx<numVertices; vtx++)
    if ((int)vertexNeighbors[vtx].size() > maxDegree)
      maxDegree = vertexNeighbors[vtx].size();
  return maxDegree;
}

int Graph::GetMinDegree() const
{
  int minDegree = INT_MAX;
  for(int vtx=0; vtx<numVertices; vtx++)
    if ((int)vertexNeighbors[vtx].size() < minDegree)
      minDegree = vertexNeighbors[vtx].size();
  return minDegree;
}

double Graph::GetAvgDegree() const
{
  double avgDegree = 0;
  for(int vtx=0; vtx<numVertices; vtx++)
    avgDegree += vertexNeighbors[vtx].size();
  return avgDegree / numVertices;
}

double Graph::GetStdevDegree() const
{
  double avgDegree_ = GetAvgDegree();
  double std = 0;
  for(int vtx=0; vtx<numVertices; vtx++)
    std += (vertexNeighbors[vtx].size() - avgDegree_) * (vertexNeighbors[vtx].size() - avgDegree_);
  return sqrt(std / numVertices);
}

int Graph::IsNeighbor(int vtx1, int vtx2) const
{
  map<int,int> :: const_iterator iter = vertexNeighbors[vtx1].find(vtx2);
  if (iter == vertexNeighbors[vtx1].end())
    return 0;
  else
    return iter->second + 1;
}

std::map<int,int> Graph::GetNeighborhoodWithDistance(int vertex, int neighborhoodSize)
{
  set<int> seed = {vertex};
  return GetNeighborhoodWithDistance(seed, neighborhoodSize);
}

std::map<int,int> Graph::GetNeighborhoodWithDistance(const std::set<int> & seedVertices, int neighborhoodSize)
{
  map<int, int> foundVertices;
  vector<int> lastLayerVertices;;
  for(int vtx : seedVertices)
  {
    foundVertices.emplace(vtx, 0);
    lastLayerVertices.push_back(vtx);
  }

  vector<int> newAffectedVertices;
  for(int i=1; i<=neighborhoodSize; i++)
  {
    for(size_t j = 0; j < lastLayerVertices.size(); j++)
    {
      // traverse all neighbors and check if they were already previously inserted
      int vtx = lastLayerVertices[j];
      int deg = GetNumNeighbors(vtx);
      for(int k=0; k<deg; k++)
      {
        int vtxNeighbor = GetNeighbor(vtx, k);
        if (foundVertices.find(vtxNeighbor) == foundVertices.end())
        {
          // discovered new vertex
          newAffectedVertices.push_back(vtxNeighbor);
          foundVertices.emplace(vtxNeighbor, i);
        }
      }
    }

    lastLayerVertices.swap(newAffectedVertices);
    newAffectedVertices.clear();
  }
  return foundVertices;
}

void Graph::ExpandNeighbors()
{
  // over all edges:
  // insert neigbors of every vtx into the edges

  set<pair<int, int> > expandedEdges = edges;

  for(set<pair<int, int> > :: iterator iter = edges.begin(); iter != edges.end(); iter++)
  {
    int vtxA = iter->first; 
    int vtxB = iter->second; 

    // connect all neighbors of A to B
    for(map<int,int> :: iterator mapIter = vertexNeighbors[vtxA].begin(); mapIter != vertexNeighbors[vtxA].end(); mapIter++)
    {
      if (vtxB < mapIter->first)
        expandedEdges.insert(make_pair(vtxB, mapIter->first)); 
    }

    // connect all neigbhors of B to A
    for(map<int,int> :: iterator mapIter = vertexNeighbors[vtxB].begin(); mapIter != vertexNeighbors[vtxB].end(); mapIter++)
      if (vtxA < mapIter->first)
        expandedEdges.insert(make_pair(vtxA, mapIter->first)); 
  }
 
  edges = expandedEdges;
  numEdges = edges.size();
  BuildVertexNeighbors();
}

void Graph::PrintInfo() const
{
  printf("Graph vertices: %d\n", numVertices);
  printf("Graph edges: %d\n", numEdges);
  printf("Graph min degree: %d\n", GetMinDegree());
  printf("Graph max degree: %d\n", GetMaxDegree());
  printf("Graph avg degree: %G\n", GetAvgDegree());
  printf("Graph degree stdev: %G\n", GetStdevDegree());
}

void Graph::GetLaplacian(SparseMatrix ** L, int scaleRows) const
{
  SparseMatrixOutline outline(3*numVertices);
  for(int i=0; i<numVertices; i++)
  {
    int numNeighbors = (int)vertexNeighborsVector[i].size();
    if (numNeighbors == 0)
      continue;

    for(int k=0; k<3; k++)
      outline.AddEntry(3 * i + k, 3 * i + k, (scaleRows != 0) ? 1.0 : numNeighbors);

    double weight;
    if (scaleRows != 0)
      weight = -1.0 / numNeighbors;
    else
      weight = -1.0;

    for(int j=0; j<numNeighbors; j++)
      for(int k=0; k<3; k++)
        outline.AddEntry(3 * i + k, 3 * vertexNeighborsVector[i][j] + k, weight);
  }

  *L = new SparseMatrix(&outline);
}

Graph * Graph::CartesianProduct(Graph & graph2) const
{
  int numProductVertices = numVertices * graph2.numVertices;
  int numProductEdges = numEdges * graph2.numVertices + numVertices * graph2.numEdges;
  int * productEdges = (int*) malloc (sizeof(int) * 2 * numProductEdges);
 
  printf("Num space-time graph vertices: %d\n", numProductVertices);
  printf("Num space-time graph edges: %d\n", numProductEdges);

  int edge = 0;
  for(int j=0; j<graph2.numVertices; j++)
  {
    for(int i=0; i<numVertices; i++)
    {
      // connect every vertex of graph1 to its neighbors
      //std::vector< std::vector<int> > vertexNeighborsVector;
      for(int k=0; k<(int)vertexNeighborsVector[i].size(); k++)
      {  
        if (i > vertexNeighborsVector[i][k])
        {
          productEdges[2*edge+0] = GetCartesianProductVertexIndex(i, j);
          productEdges[2*edge+1] = GetCartesianProductVertexIndex(vertexNeighborsVector[i][k], j);
          edge++;
        }
      }
      // connect every vertex of graph2 to its neighbors
      for(int k=0; k<(int)(graph2.vertexNeighborsVector[j].size()); k++)
      {  
        if (j > graph2.vertexNeighborsVector[j][k])
        {
          productEdges[2*edge+0] = GetCartesianProductVertexIndex(i, j);
          productEdges[2*edge+1] = GetCartesianProductVertexIndex(i, graph2.vertexNeighborsVector[j][k]);
          edge++;
        }
      }
    }
  }

  Graph * graph = new Graph(numProductVertices, numProductEdges, productEdges);
  free(productEdges);
  return graph;
}

// cluster given vertices into connected components
void Graph::Cluster(const set<int> & vertices, vector<set<int> > & clusters) const
{
  clusters.clear();

  set<int> remainingVertices = vertices;

  while(remainingVertices.size() > 0)
  {
    int seed = *remainingVertices.begin();
    clusters.resize(clusters.size() + 1);
    set<int> & curCluster = clusters.back();

    vector<int> oldFront, front;

    curCluster.insert(seed);
    remainingVertices.erase(seed);
    oldFront.push_back(seed);

    while (oldFront.size() > 0)
    {
      // create the front
      front.clear();
      for(size_t i = 0; i < oldFront.size(); i++)
      {
        int node = oldFront[i];
        for(size_t j = 0; j < vertexNeighborsVector[node].size(); j++)
        {
          int neighbor = vertexNeighborsVector[node][j];
          if (remainingVertices.find(neighbor) != remainingVertices.end())
          {
            front.push_back(neighbor);
            remainingVertices.erase(neighbor);
            curCluster.insert(neighbor);
          }
        }
      }
      oldFront = front;
    }
  }
}

int Graph::GetCartesianProductVertexIndex(int vertex1, int vertex2) const
{
  return vertex2 * numVertices + vertex1;
}

void Graph::GetCartesianProductVertexIndexComponents(int productVertex, int * vertex1, int * vertex2) const
{
  *vertex2 = productVertex / numVertices;
  *vertex1 = productVertex % numVertices;
}

void Graph::ShortestDistance(const std::set<int> & seedVertices, std::vector<int> & distances) const
{
  distances.assign(numVertices, INT_MAX);
  for(set<int> :: const_iterator iter = seedVertices.begin(); iter != seedVertices.end(); iter++)
    distances[*iter] = 0;

  int distance = 0;
  vector<int> oldFront, front;
  oldFront.reserve(numVertices);
  front.reserve(numVertices);
  oldFront.insert(oldFront.end(), seedVertices.begin(), seedVertices.end());
  while (oldFront.size() > 0)
  {
    distance++;

    // create the front
    front.clear();
    for(size_t i = 0; i < oldFront.size(); i++)
    {
      int node = oldFront[i];
      for(size_t j = 0; j < vertexNeighborsVector[node].size(); j++)
      {
        int neighbor = vertexNeighborsVector[node][j];
        if (distances[neighbor] == INT_MAX)
        {
          front.push_back(neighbor);
          distances[neighbor] = distance;
        }
      }
    }

    oldFront = front; 
  }
}

bool Graph::FindShortestPath(const std::set<int> & seedVertices, const std::set<int> & destinationVertices, std::vector<int> * path) const
{
  if (path)
    path->clear();
  int distance = 0;
  vector<int> oldFront, front;
  oldFront.insert(oldFront.end(), seedVertices.begin(), seedVertices.end());
  map<int, int> parentNode;
  for(set<int>::const_iterator it = seedVertices.begin(); it != seedVertices.end(); it++)
  {
    parentNode[*it] = -1;
    if (destinationVertices.find(*it) != destinationVertices.end()) // one of the start locations is the destination
    {
      if (path)
        path->push_back(*it);
      return true;
    }
  }

  while (oldFront.size() > 0)
  {
    distance++;

    // create the front
    front.clear();
    for(size_t i = 0; i < oldFront.size(); i++)
    {
      int node = oldFront[i];
      for(size_t j = 0; j < vertexNeighborsVector[node].size(); j++)
      {
        int neighbor = vertexNeighborsVector[node][j];
        if (parentNode.find(neighbor) == parentNode.end())
        {
          front.push_back(neighbor);
          parentNode.insert(pair<int, int>(neighbor, node));
          if (destinationVertices.find(neighbor) != destinationVertices.end()) // we reach the destination
          {
            if (path)
            {
              int parent = neighbor, curNode = 0;
              do
              {
                curNode = parent;
                path->push_back(curNode); // push nodes into path in a reverse order
                parent = parentNode[curNode];
              } while(parent >= 0);
              reverse(path->begin(), path->end()); // recover the correct order
            }
            return true;
          }
        }
      }
    }

    oldFront = front;
  }

  return false;
}

bool Graph::FindLocalShortestPath(const std::set<int> & localVertices, const std::set<int> & seedVertices, const std::set<int> & destinationVertices, std::vector<int> * path) const
{
  if (path)
    path->clear();
  int distance = 0;
  vector<int> oldFront, front;
  map<int, int> parentNode;

  for(set<int>::const_iterator it = seedVertices.begin(); it != seedVertices.end(); it++)
  {
    if (localVertices.find(*it) == localVertices.end())
      continue; // the seed vertex is not inside localVertices
    oldFront.push_back(*it);
    parentNode[*it] = -1;
    if (destinationVertices.find(*it) != destinationVertices.end()) // one of the start locations is the destination
    {
      if (path)
        path->push_back(*it);
      return true;
    }
  }

  while (oldFront.size() > 0)
  {
    distance++;

    // create the front
    front.clear();
    for(size_t i = 0; i < oldFront.size(); i++)
    {
      int node = oldFront[i];
      for(size_t j = 0; j < vertexNeighborsVector[node].size(); j++)
      {
        int neighbor = vertexNeighborsVector[node][j];
        if (localVertices.find(neighbor) == localVertices.end())
          continue;
        if (parentNode.find(neighbor) == parentNode.end())
        {
          front.push_back(neighbor);
          parentNode.insert(pair<int, int>(neighbor, node));
          if (destinationVertices.find(neighbor) != destinationVertices.end()) // we reach the destination
          {
            if (path)
            {
              int parent = neighbor, curNode = 0;
              do
              {
                curNode = parent;
                path->push_back(curNode); // push nodes into path in a reverse order
                parent = parentNode[curNode];
              } while(parent >= 0);
              reverse(path->begin(), path->end()); // recover the correct order
            }
            return true;
          }
        }
      }
    }

    oldFront = front;
  }

  return false;
}


void Graph::GetConnectedComponent(int vtx, std::set<int> &connectedVertices) const
{
  std::vector<int> Q;
  std::set<int> visited;
  size_t start = 0;

  Q.push_back(vtx);
  visited.insert(vtx);

  while (Q.size() - start > 0) {
    int v = Q[start];
    start++;

    for (int i = 0; i < (int)vertexNeighborsVector[v].size(); i++) {
      int c = vertexNeighborsVector[v][i];
      if (visited.find(c) == visited.end()) {
        Q.push_back(c);
        visited.insert(c);
      }
    }
  }

  connectedVertices = visited;
}
