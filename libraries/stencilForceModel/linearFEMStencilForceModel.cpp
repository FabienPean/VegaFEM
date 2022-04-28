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

#include "linearFEMStencilForceModel.h"

#ifdef USE_TBB
#include <tbb/tbb.h>
#endif

#include <cassert>
#include <cstring>
#include <iostream>

using namespace std;

LinearFEMStencilForceModel::LinearFEMStencilForceModel(StencilForceModel * fem): stencilForceModel(fem)
{
  assert(dynamic_cast<LinearFEMStencilForceModel*>(stencilForceModel) == nullptr);

  n3 = fem->Getn3();
  elementK.resize(fem->GetNumStencilTypes());

  std::vector<double> u(n3, 0.0);
  for (int eltype = 0; eltype < fem->GetNumStencilTypes(); eltype++) 
  {
    numStencilsInDifferentTypes.push_back(stencilForceModel->GetNumStencils(eltype));
    numStencilVerticesInDifferentTypes.push_back(stencilForceModel->GetNumStencilVertices(eltype));

    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    int nele = stencilForceModel->GetNumStencils(eltype);

    elementK[eltype].resize(nele * nelev * nelev * 9);

    cout << "Computing element stiffness matrices.." << endl;

#ifdef USE_TBB
    tbb::parallel_for(0, nele, [&] (int el) 
    {
      stencilForceModel->GetStencilLocalEnergyAndForceAndMatrix(eltype, el, u.data(), nullptr, nullptr, elementK[eltype].data() + el * nelev * nelev * 9);
    }, tbb::static_partitioner());
#else
    for (int el = 0; el < nele; el++)
    {
      stencilForceModel->GetStencilLocalEnergyAndForceAndMatrix(eltype, el, u.data(), nullptr, nullptr, elementK[eltype].data() + el * nelev * nelev * 9);
    }
#endif

    cout << "Done." << endl;
  }
}

LinearFEMStencilForceModel::~LinearFEMStencilForceModel() {}

void LinearFEMStencilForceModel::GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix)
{
  int nelev = numStencilVerticesInDifferentTypes[stencilType];
  int dof = nelev * 3;
  int dof2 = dof * dof;
  const double *K0 = elementK[stencilType].data() + stencilId * dof2;
  const int *vtxIdx = stencilForceModel->GetStencilVertexIndices(stencilType, stencilId);

  if (tangentStiffnessMatrix)
    memcpy(tangentStiffnessMatrix, K0, sizeof(double) * dof2);

  if (internalForces) 
  {
    memset(internalForces, 0, sizeof(double) * dof);
    for (int i = 0; i < dof; i++)
    {
      for (int j = 0; j < dof; j++)
      {
        internalForces[j] += K0[i * dof + j] * u[vtxIdx[i / 3] * 3 + (i % 3)];
      }
    }
  }

  if (energy)
  {
    *energy = 0.0;
    for (int i = 0; i < dof; ++i)
      for (int j = 0; j < dof; ++j)
        *energy += u[vtxIdx[i / 3] * 3 + (i % 3)] * u[vtxIdx[j / 3] * 3 + (j % 3)] * K0[i * dof + j] * 0.5;
  }

}

const int *LinearFEMStencilForceModel::GetStencilVertexIndices(int stencilType, int stencilId) const
{
  return stencilForceModel->GetStencilVertexIndices(stencilType, stencilId);
}
