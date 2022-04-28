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

#include "StVKStencilForceModel.h"
#include <cassert>
using namespace std;

StVKStencilForceModel::StVKStencilForceModel(StVKFEM * fem): stvkFEM(fem)
{
  n3 = stvkFEM->GetVolumetricMesh()->getNumVertices() * 3;
  numStencilsInDifferentTypes.push_back(stvkFEM->GetVolumetricMesh()->getNumElements());
  numStencilVerticesInDifferentTypes.push_back(stvkFEM->GetVolumetricMesh()->getNumElementVertices());
}

StVKStencilForceModel::~StVKStencilForceModel() {}

void StVKStencilForceModel::GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix)
{
  assert(stencilType == 0);
  stvkFEM->ComputeElementLocalEnergyAndInternalForcesAndStiffnessMatrix(u, stencilId, energy, internalForces, tangentStiffnessMatrix);
}

const int *StVKStencilForceModel::GetStencilVertexIndices(int stencilType, int stencilId) const
{
  assert(stencilType == 0);
  return stvkFEM->GetVolumetricMesh()->getVertexIndices(stencilId);
}
