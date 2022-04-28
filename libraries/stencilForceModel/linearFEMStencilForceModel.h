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


#ifndef _LINEARFEM_STENCIL_FORCEMODEL_H_
#define _LINEARFEM_STENCIL_FORCEMODEL_H_

#include "stencilForceModel.h"

#include <vector>

// Stencils for linear FEM.
// A stencil is one FEM element.
// See comments in the parent class.
class LinearFEMStencilForceModel : public StencilForceModel
{
public:
  // this class is able to linearize any non-linear force model
  // therefore, any StencilForceModel (except itself) can be used to initialize this class.
  LinearFEMStencilForceModel(StencilForceModel * stencilForceModel);
  virtual ~LinearFEMStencilForceModel();

  virtual const int *GetStencilVertexIndices(int stencilType, int stencilId) const override;
  virtual void GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix) override;

  StencilForceModel * GetForceModelHandle() { return stencilForceModel; }

protected:
  StencilForceModel * stencilForceModel;

  std::vector<std::vector<double>> elementK;
};

#endif

