/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "isotropic hyperelastic FEM" library , Copyright (C) 2018 USC         *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Fun Shing Sin                            *
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

#ifndef VEGAFEM_NEOHOOKEANISOTROPICMATERIAL_H
#define VEGAFEM_NEOHOOKEANISOTROPICMATERIAL_H

#include "tetMesh.h"
#include "isotropicMaterialWithCompressionResistance.h"

namespace vegafem
{

/*
   Neo-Hookean material. Material properties are read from the tet mesh, and can be heterogeneous.

   The implemented neo-Hookean material is described in:
   BONET J., WOOD R. D.: Nonlinear Continuum Mechanics
   for Finite Element Analysis, 2nd Ed. Cambridge University
   Press, 2008, page 162
*/

class NeoHookeanIsotropicMaterial : public IsotropicMaterialWithCompressionResistance
{
public:
  NeoHookeanIsotropicMaterial(TetMesh * tetMesh, int enableCompressionResistance=0, double compressionResistance=0.0);
  virtual ~NeoHookeanIsotropicMaterial();

  virtual double ComputeEnergy(int elementIndex, double * invariants);
  virtual void ComputeEnergyGradient(int elementIndex, double * invariants, double * gradient); // invariants and gradient are 3-vectors
  virtual void ComputeEnergyHessian(int elementIndex, double * invariants, double * hessian); // invariants is a 3-vector, hessian is a 3x3 symmetric matrix, unrolled into a 6-vector, in the following order: (11, 12, 13, 22, 23, 33).

protected:
  double * lambdaLame; 
  double * muLame; 

  double compressionResistance;
  double * EdivNuFactor;
  virtual double GetCompressionResistanceFactor(int elementIndex);
};


}//namespace vegafem

#endif


