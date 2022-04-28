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

/*
  Force model for the Corotational linear FEM material.
*/

#ifndef _COROTATIONALLINEARFEMFORCEMODEL_H_
#define _COROTATIONALLINEARFEMFORCEMODEL_H_

#include "corotationalLinearFEM.h"
#include "forceModel.h"

class CorotationalLinearFEMForceModel : public ForceModel
{
public:
  CorotationalLinearFEMForceModel(CorotationalLinearFEM * corotationalLinearFEM, int warp=1);
  virtual ~CorotationalLinearFEMForceModel(); 

  virtual double GetElasticEnergy(const double * u) override;
  virtual void GetInternalForce(const double * u, double * internalForces) override;
  virtual void GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix) override;
  virtual void GetTangentStiffnessMatrix(const double * u, SparseMatrix * tangentStiffnessMatrix) override; 

  virtual void GetForceAndMatrix(const double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix) override;

  inline void SetWarp(int warp) { this->warp = warp; }
  inline int GetWarp() { return warp; }
  static inline int GetMaxWarp() { return 2; }

protected:
  CorotationalLinearFEM * corotationalLinearFEM;
  int warp;
};

#endif

