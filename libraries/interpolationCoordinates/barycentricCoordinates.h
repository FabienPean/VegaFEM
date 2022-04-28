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

#ifndef BARYCENTRICCOORDINATES_H
#define BARYCENTRICCOORDINATES_H

#include "volumetricMesh.h"
#include "interpolationCoordinates.h"
#include "sparseMatrix.h"
#include "vec4d.h"
#include "vec4i.h"
#include "tetMesh.h"
#include "objMesh.h"
#include <vector>
#include <string>

// Compute and store barycentric coordinates.
// Optionally, compute the indices of the elements that the embedded vertices belong to.

class BarycentricCoordinates : public InterpolationCoordinates
{
public:

  BarycentricCoordinates() {}

  // Load from file. No element indices are available.
  BarycentricCoordinates(const std::string & filename);

  // Initialize with only element vertex indices and weights. No element indices are available.
  BarycentricCoordinates(int numLocations, int numElementVertices, const int * elementVertexIndices, const double * weights,
    const int * elementIndices = nullptr);

  BarycentricCoordinates(const std::vector<Vec4i> & tetVtxIndices, const std::vector<Vec4d> & tetWeights, const int * elementIndices = nullptr);

  virtual ~BarycentricCoordinates() {}

  virtual void deform(const double * verticesDisp, double * locationDisp) const override;

  int getNumLocations() const { return numLocations; }
  int getNumElementVertices() const { return numElementVertices; }

  const int * getEmbeddingVertexIndices(int embeddedVtx) const { return &indices[numElementVertices * embeddedVtx]; }
  const double * getEmbeddingWeights(int embeddedVtx) const { return &weights[numElementVertices * embeddedVtx]; }

  // Get internal data.
  std::vector<int> & getEmbeddingVertexIndices() { return indices; }
  std::vector<double> & getEmbeddingWeights() { return weights; }
  std::vector<int> & getElements() { return elements; }

  // Return embedding element index if it's available; -1 otherwise.
  int getEmbeddingElement(int embeddedVtx) const { return elements[embeddedVtx]; }

  // Save weights to file.
  int saveInterpolationWeights(const std::string & filename) const;

  SparseMatrix * generateInterpolationMatrix() const;

protected:
  int numLocations = 0;
  int numElementVertices = 0;
  std::vector<int> indices;
  std::vector<double> weights;
  std::vector<int> elements;
};

#endif /* BARYCENTRICCOORDINATES_H */

