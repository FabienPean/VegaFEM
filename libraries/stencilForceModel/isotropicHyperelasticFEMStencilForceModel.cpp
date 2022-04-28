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

#include "isotropicHyperelasticFEMStencilForceModel.h"
#include <cassert>

using namespace std;

IsotropicHyperelasticFEMStencilForceModel::IsotropicHyperelasticFEMStencilForceModel(IsotropicHyperelasticFEM * isotropicHyperelasticFEM_): isotropicHyperelasticFEM(isotropicHyperelasticFEM_)
{
  n3 = 3 * isotropicHyperelasticFEM->GetTetMesh()->getNumVertices();
  numStencilsInDifferentTypes.push_back(isotropicHyperelasticFEM->GetTetMesh()->getNumElements());
  numStencilVerticesInDifferentTypes.push_back(isotropicHyperelasticFEM->GetTetMesh()->getNumElementVertices());
}

IsotropicHyperelasticFEMStencilForceModel::~IsotropicHyperelasticFEMStencilForceModel() {}

void IsotropicHyperelasticFEMStencilForceModel::GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix)
{
  assert(stencilType == 0);

  isotropicHyperelasticFEM->GetElementLocalEnergyAndForceAndMatrix(stencilId, u, energy, internalForces, tangentStiffnessMatrix);
}

const int *IsotropicHyperelasticFEMStencilForceModel::GetStencilVertexIndices(int stencilType, int stencilId) const
{
  assert(stencilType == 0);

  return isotropicHyperelasticFEM->GetTetMesh()->getVertexIndices(stencilId);
}

