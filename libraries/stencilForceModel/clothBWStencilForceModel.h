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


#ifndef _CLOTHBW_STENCIL_FORCEMODEL_H_
#define _CLOTHBW_STENCIL_FORCEMODEL_H_

#include "clothBW.h"
#include "stencilForceModel.h"

// Stencils for cloth simulation. There are two different stencils here.
// One is a triangle (# vertices = 3) (for in-place stretch and shear);
// and the other is an edge joining two triangles (# vertices = 4) (for bending).
// See comments in the parent class.
class ClothBWStencilForceModel : public StencilForceModel
{
public:
  ClothBWStencilForceModel(ClothBW * clothBW);
  virtual ~ClothBWStencilForceModel();

  virtual void GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix) override;
  virtual const int *GetStencilVertexIndices(int stencilType, int stencilId) const override;

  virtual void GetVertexGravityForce(int vertexId, double gravity[3]) override;

  ClothBW * GetForceModelHandle() { return clothBW; }
protected:
  ClothBW * clothBW;
};

#endif

