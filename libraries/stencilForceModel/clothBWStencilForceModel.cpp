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

#include "clothBWStencilForceModel.h"
#include <iostream>
#include <cstdlib>

ClothBWStencilForceModel::ClothBWStencilForceModel(ClothBW * clothBW_): clothBW(clothBW_)
{ 
  numStencilsInDifferentTypes.push_back(clothBW->GetNumTriangles());
  numStencilsInDifferentTypes.push_back(clothBW->GetNumQuads());

  numStencilVerticesInDifferentTypes.push_back(3);
  numStencilVerticesInDifferentTypes.push_back(4);

  n3 = clothBW->GetNumVertices() * 3;
}

ClothBWStencilForceModel::~ClothBWStencilForceModel()
{
}

void ClothBWStencilForceModel::GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix)
{
  if (stencilType == 0)
    clothBW->ComputeTriangleElement(stencilId, u, energy, internalForces, tangentStiffnessMatrix);
  else if (stencilType == 1)
    clothBW->ComputeQuadElement(stencilId, u, energy, internalForces, tangentStiffnessMatrix);
  else {
    std::cerr << "Wrong type of element!" << std::endl;
    abort();
  }
}

const int *ClothBWStencilForceModel::GetStencilVertexIndices(int stencilType, int stencilId) const
{
  if (stencilType == 0)
    return clothBW->GetTriangleVertexIndices(stencilId);
  else if (stencilType == 1)
    return clothBW->GetQuadVertexIndices(stencilId);
  else {
    std::cerr << "Wrong type of element!" << std::endl;
    abort();
  }

  return nullptr;
}

void ClothBWStencilForceModel::GetVertexGravityForce(int vid, double gravity[3])
{
  clothBW->ComputeVertexGravity(vid, gravity);
}

