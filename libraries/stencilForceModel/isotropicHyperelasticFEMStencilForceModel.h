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

#ifndef _ISOTROPICHYPERELASTICFEM_STENCIL_FORCEMODEL_H_
#define _ISOTROPICHYPERELASTICFEM_STENCIL_FORCEMODEL_H_

#include "isotropicHyperelasticFEM.h"
#include "stencilForceModel.h"

// Stencils for isotropic hyperelastic FEM.
// A stencil here refers to one FEM element.
// See comments in the parent class.
class IsotropicHyperelasticFEMStencilForceModel : public StencilForceModel
{
public:
  IsotropicHyperelasticFEMStencilForceModel(IsotropicHyperelasticFEM * isotropicHyperelasticFEM);
  virtual ~IsotropicHyperelasticFEMStencilForceModel();

  virtual const int *GetStencilVertexIndices(int stencilType, int stencilId) const override;
  virtual void GetStencilLocalEnergyAndForceAndMatrix(int stencilType, int stencilId, const double * u, double * energy, double * internalForces, double * tangentStiffnessMatrix) override;

  IsotropicHyperelasticFEM * GetForceModelHandle() { return isotropicHyperelasticFEM; }

protected:
  IsotropicHyperelasticFEM * isotropicHyperelasticFEM;
};

#endif

