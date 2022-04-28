/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "interpolationCoordinates" library , Copyright (C) 2018 USC           *
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

#include "barycentricCoordinates.h"
#include "generateInterpolationMatrix.h"
#ifdef CONSTRUCTOR_WITH_SEED
  #include "objMeshGraph.h"
  #include "tetMeshManifold.h"
#endif
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
#include <climits>
using namespace std;

BarycentricCoordinates::BarycentricCoordinates(int numLocations_, int numElementVertices_, const int* indices_, const double* weights_,
  const int * elementIndices_)
{
  numLocations = numLocations_;
  numElementVertices = numElementVertices_;
  indices.resize(numElementVertices * numLocations);
  weights.resize(numElementVertices * numLocations);
  memcpy(&indices[0], indices_, sizeof(int) * indices.size());
  memcpy(&weights[0], weights_, sizeof(double) * weights.size());
  elements.resize(numLocations, -1);
  if (elementIndices_)
    elements.assign(elementIndices_, elementIndices_ + numLocations);
}

BarycentricCoordinates::BarycentricCoordinates(const std::vector<Vec4i> & tetVtxIndices, const std::vector<Vec4d> & tetWeights,
  const int * elementIndices) : BarycentricCoordinates(min(tetVtxIndices.size(), tetWeights.size()), 4,
  (const int*)tetVtxIndices.data(), (const double*) tetWeights.data(), elementIndices)
{
  assert(tetVtxIndices.size() == tetWeights.size());
}

#ifdef CONSTRUCTOR_WITH_SEED

BarycentricCoordinates::BarycentricCoordinates(const ObjMesh * objMesh,  const TetMesh * volumetricMesh, int numThreads, int seed)
{
  numLocations = objMesh->getNumVertices();
  numElementVertices = volumetricMesh->getNumElementVertices();
  indices.resize(numElementVertices * numLocations, 0);
  weights.resize(numElementVertices * numLocations, 0);
  elements.resize(numLocations, -1);

  vector<double> locations;
  exportObjMeshGeometry(objMesh, locations);
  TetMeshManifold manifold;

  map<UTetKey, int> index;
  for(int i = 0; i < volumetricMesh->getNumElements(); i++)
  {
    assert(manifold.add(volumetricMesh->getElementVertexIndices(i)));
    UTetKey key(volumetricMesh->getElementVertexIndices(i));
    index.insert(pair<UTetKey,int>(key, i));
  }
  const TetMeshManifold::TetMap & tetMap = manifold.getTetMap();


  vector<Vec3d> vs;
  for(int i = 0; i < volumetricMesh->getNumVertices(); i++)
    vs.push_back(*volumetricMesh->getVertex(i));
  VerticesQuery query(vs.size(), vs.data());

  assert((int)seed < numLocations);
  FastInterpolationWeights interpolation(volumetricMesh);
  interpolation.generateInterpolationWeights(&locations[3*seed], &indices[numElementVertices * seed], &weights[numElementVertices * seed], -1, &elements[seed]);

  // Creates a graph where the nodes are obj mesh vertices. Two nodes are connected if they are adjacent in the mesh.
  Graph * graph = ObjMeshGraph::GenerateVertexGraph(objMesh, 0);
  set<int> visited;
  queue<int> candidates;
  candidates.push(seed);
  visited.insert(seed);

  while(candidates.size() > 0)
  {
    int vtx = candidates.front();
    int element = elements[vtx];
    candidates.pop();
    if (element < 0)
      continue;

    int nn = graph->GetNumNeighbors(vtx);
    for(int i = 0; i < nn; i++)
    {
      int nbr = graph->GetNeighbor(vtx, i);
      if (visited.insert(nbr).second == false) // this nbr is already visited before
        continue;

      Vec3d p = objMesh->getPosition(nbr);

      // find elements containing nbr
      int curEle = element;
      while (true)
      {
        UTetKey key(volumetricMesh->getElementVertexIndices(curEle));
        TetMeshManifold::TetCIter it = tetMap.find(key);
        assert(it != tetMap.end());
        const TetMeshManifold::Tetrahedron * t = it->second;
        const TetMeshManifold::Tetrahedron * next = NULL;

        int j = 0;
        for(j = 0; j < 4; j++)
        {
          int f[3];
          t->getOppositeFace(j, f);

          if (query.toPlane(p, f[0], f[1], f[2]) > 0)
          {
            //            assert(weight[j] < 0);
            const TetMeshManifold::Tetrahedron * nbr = t->getNeighbor(j);
            if (nbr)
            {
              next = nbr;
              break;
            }
          }
        }
        if (j == 4)  // the point is inside this tet
          break;
        else if (next)
        {
          int vtx[4];
          for(int j = 0; j < 4; j++)
            vtx[j] = next->getVtx(j);
          UTetKey nbkey(vtx[0], vtx[1], vtx[2], vtx[3]);
          assert(nbkey < key || key < nbkey);
          curEle = index[nbkey];
        }
        else // the point is outside of a surface face on the tet
        {
          curEle = -1;
          break;
        }
      } // end while (true)

      if (curEle >= 0)
      {
        candidates.push(nbr);
        elements[nbr] = curEle;
        double * weight = &weights[numElementVertices * nbr];
        int * vtxIndex = &indices[numElementVertices * nbr];
        volumetricMesh->computeBarycentricWeights(curEle, p, weight);
        for(int j = 0; j < 4; j++)
          vtxIndex[j] = volumetricMesh->getVertexIndex(curEle,j);
      }
    }
  }

  for(int i = 0; i < numLocations; i++)
  {
    if (elements[i] < 0)
    {
      // use fastInterpolationWeights to compute interp data
      interpolation.generateInterpolationWeights(&locations[3*i], &indices[numElementVertices * i], &weights[numElementVertices * i], -1, &elements[i]);
    }
  }
}

#endif

void BarycentricCoordinates::deform(const double * verticesDisp, double* locationDisp) const
{
  VolumetricMesh::interpolate(verticesDisp, locationDisp, numLocations, numElementVertices, &indices[0], &weights[0]);
}

int BarycentricCoordinates::saveInterpolationWeights(const string & filename) const
{
  //saveInterpolationWeights(const char * filename, int numTargetLocations, int numElementVertices, int * vertices, double * weights);
  int ret = VolumetricMesh::saveInterpolationWeights(filename.c_str(), numLocations, numElementVertices, &indices[0], &weights[0]);
  cout << (ret == 0 ? "Saved" : "Failed to save");
  cout << " interpolation weights (numLocations: " << numLocations << ", numElementVertices: " << numElementVertices << ") to " << filename << "." << endl;
  return ret;
}

SparseMatrix * BarycentricCoordinates::generateInterpolationMatrix() const
{
  SparseMatrix * A = NULL;
  GenerateInterpolationMatrix::generate(numLocations, numElementVertices, &indices[0], &weights[0], &A);
  return A;
}

#define READ_ONE_PAIR \
  do { \
    int index = 0; \
    double w = 0; \
    ss >> index; \
    ss >> w; \
    if (ss.fail()) { \
      cerr << index << " " << w << endl; \
      cerr << "Error: incorrect interp file format at \"" << buffer << "\" in interp file " << filename << "." << endl; \
      throw 1; \
    } \
    if (index < 0) { \
      cerr << "Error: invalid index at \"" << buffer << "\" in interp file " << filename << "." << endl; \
      throw 1; \
    } \
    indices.push_back(index); \
    weights.push_back(w); \
  } while(0)

BarycentricCoordinates::BarycentricCoordinates(const std::string & filename)
{
  ifstream fin(filename.c_str(), ios::binary);
  if (!fin)
  {
    cerr << "Error: cannot open interp file " << filename << "." << endl;
    throw 1;
  }
  fin >> ws;

  int count = 0;
  string buffer;
  stringstream ss;
  numElementVertices = INT_MAX;
  while(!fin.eof())
  {
    getline(fin, buffer);
    fin >> ws;
    if (buffer.size() == 0)
      continue;

//    istringstream ss(buffer);
    ss.clear();
    ss.str(buffer);
    ss.seekg(0);
    int c = 0;
    ss >> c;
//    cout << "At line: " << buffer << endl;
    if (c != count)
    {
      cerr << "Warning: wrong line index at \"" << buffer << "\" in interp file " << filename << "." << endl;
    }
    if (numElementVertices == INT_MAX)
    {
      // determine #element vertices
      ss >> ws;
//      cout << "First line " << endl;
      numElementVertices = 0;
      while(!ss.eof())
      {
//        cout << "RRR" << endl;
        READ_ONE_PAIR;
        numElementVertices++;
//        cout << "numElementVertices now: " << numElementVertices << endl;
        ss >> ws;
      }
    }
    else
    {
      for(int i = 0; i < numElementVertices; i++)
      {
        READ_ONE_PAIR;
      }
    }
    count++;
  }
  numLocations = count;

  elements.resize(numLocations, -1); // no information for elements read from file

  cout << "Loaded " << numLocations << " locations from file " << filename << " with numElementVertices = " << numElementVertices << "." << endl;
  assert((int)indices.size() == numLocations * numElementVertices);
  assert((int)weights.size() == numLocations * numElementVertices);
}
