/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "forceModel" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC     *
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

#include "corotationalLinearFEMForceModel.h"

CorotationalLinearFEMForceModel::CorotationalLinearFEMForceModel(CorotationalLinearFEM * corotationalLinearFEM_, int warp_): corotationalLinearFEM(corotationalLinearFEM_), warp(warp_)
{
  r = 3 * corotationalLinearFEM->GetVolumetricMesh()->getNumVertices();
}

CorotationalLinearFEMForceModel::~CorotationalLinearFEMForceModel() {}

double CorotationalLinearFEMForceModel::GetElasticEnergy(const double * u)
{
  double energy = 0.0;
  corotationalLinearFEM->ComputeEnergyAndForceAndStiffnessMatrix(u, &energy, NULL, NULL, warp);
  return energy;
}

void CorotationalLinearFEMForceModel::GetInternalForce(const double * u, double * internalForces)
{
  corotationalLinearFEM->ComputeEnergyAndForceAndStiffnessMatrix(u, NULL, internalForces, NULL, warp);
}

void CorotationalLinearFEMForceModel::GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
  corotationalLinearFEM->GetStiffnessMatrixTopology(tangentStiffnessMatrix);
}

void CorotationalLinearFEMForceModel::GetTangentStiffnessMatrix(const double * u, SparseMatrix * tangentStiffnessMatrix)
{
  corotationalLinearFEM->ComputeEnergyAndForceAndStiffnessMatrix(u, NULL, NULL, tangentStiffnessMatrix, warp);
} 

void CorotationalLinearFEMForceModel::GetForceAndMatrix(const double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
  corotationalLinearFEM->ComputeEnergyAndForceAndStiffnessMatrix(u, NULL, internalForces, tangentStiffnessMatrix, warp);
}

