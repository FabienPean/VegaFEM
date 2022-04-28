/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Stencil Force Model" library , Copyright (C) 2018 USC                *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Bohan Wang, Jernej Barbic                               *
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

#ifndef _STENCILFORCEMODEL_H_
#define _STENCILFORCEMODEL_H_

#include <vector>

/*
  To simulate a deformable object, we usually discretize it into elements.
  For example, a 3D solid can be represented by a set of tetrahedra,
  and a cloth can be represented by a set of triangles. When we simulate 
  an object, we need to evaluate the elastic energy, internal forces
  and sometimes the tangent stiffness matrix. These quantities are usually
  evaluated by computing them for each individual element
  and then assembled into a single global per-object vector.
  For example, when we compute the elastic energy of a tetrahedral mesh,
  we can compute the energy of each tetrahedron first and then add them together.
  In summary, the entire evaluation process can be considered as
  an assembling of quantities at the individual units. We call the units "stencils".
  In cloth simulation, a stencil is a single triangle for in-plane stretch and shear.
  Another stencil type is an edge between two triangles to model the bending energy.
  In FEM solid simulation, each element can be treated as a stencil.
  In mass-spring systems, each spring can be treated as a stencil.
  In this file, we define a class called StencilForceModel that provides an 
  abstract interface to stencils.
*/

class StencilForceModel
{
public:
  StencilForceModel() {}
  virtual ~StencilForceModel() {}

  // Get the number of degrees of freedom (DOF) of the object (= # vertices x 3).
  int Getn3() const { return n3; }
  // Return the number of different types of stencils within a single object.
  int GetNumStencilTypes() const { return (int)numStencilsInDifferentTypes.size(); }
  // Return the number of stencils for each type.
  int GetNumStencils(int stencilType) const { return numStencilsInDifferentTypes[stencilType]; }
  // Return the number of vertices involved in a single stencil.
  int GetNumStencilVertices(int stencilType) const { return numStencilVerticesInDifferentTypes[stencilType]; }
  // Return the number of DOFs of the internal forces in a single stencil.
  int GetStencilInternalForceSize(int stencilType) const { return numStencilVerticesInDifferentTypes[stencilType] * 3; }
  // Return the size of the stiffness matrix of a single stencil.
  int GetStencilStiffnessMatrixSize(int stencilType) const { return GetStencilInternalForceSize(stencilType) * GetStencilInternalForceSize(stencilType); }
  // Compute the energy, internal forces, tangent stiffness matrix of stencil stencilId in type stencilType.
  // Parameter u is the displacement vector of object vertices in R^n3.
  // Energy points to a double value.
  // internalForces points to a dense vector.
  // tangentStiffnessMatrix points to a dense column-major matrix
  // The pointers energy, internalForces and tangentStiffnessMatrix can be nullptr, in which case the corresponding quantity will not be computed.
  virtual void GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix) = 0;
  
  // Return an array of vertex indices that a stencil 'stencilId' in type 'stencilType' involves.
  // Typically, a tetrahedron involves 4 vertices. The return pointer will point to an array with 4 integers.
  virtual const int *GetStencilVertexIndices(int stencilType, int stencilId) const = 0;

  // Vertex gravity.
  virtual void GetVertexGravityForce(int vertexId, double gravity[3]) { gravity[0] = gravity[1] = gravity[2] = 0.0; }

protected:
  std::vector<int> numStencilsInDifferentTypes;
  std::vector<int> numStencilVerticesInDifferentTypes;
  int n3;
};

#endif

