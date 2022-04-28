/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "StVK" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Bohan Wang, Jernej Barbic                                *
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

#ifndef _STVK_FEM_H_
#define _STVK_FEM_H_

#include "volumetricMesh.h"
#include "StVKElementABCD.h"

#include <vector>


// This class implements the StVK force model in a element style.
// It specifically computes the simulation quantities in a scale of a FEM element,
// using 'ComputeElementLocalEnergyAndInternalForcesAndStiffnessMatrix'.
// If you want to compute the quantities in object scale,
// please look into stencilForceModel lib,
// or use stVKInternalForces.* and stVKStiffnessMatrix.*.
// The later files are the same as the previous version only for backward compatibility.
// The same functionality can be achieved by using stencilForceModel/forceModelAssembler.*.

class StVKFEM
{
public:
  StVKFEM(VolumetricMesh * volumetricMesh, StVKElementABCD * precomputedABCDIntegrals, bool addGravity=false, double g=9.81);
  virtual ~StVKFEM();  

  // enables or disables the gravity (note: you can also set this in the constructor; 
  // use this routine to turn the gravity on/off during the simulation)
  // if addGravity is enabled, ComputeForces will subtract the gravity force from the internal forces 
  // (note: subtraction, not addition, is used because the internal forces are returned with the sign as described in the f_int(x) comment above)
  void SetGravity(bool addGravity) { this->addGravity = addGravity; InitGravity(); } 

  // get simulation mesh
  const VolumetricMesh * GetVolumetricMesh() const { return volumetricMesh; }

  // get integral class object
  StVKElementABCD * GetPrecomputedIntegrals() const { return precomputedIntegrals; }

  // compute the energy E, internal forces f and stiffness matrix K of a single FEM element.
  // u is the input displacement with dimention (# all vtx x 3)
  // E is a pointer to a double
  // f is a pointer to a dense vector
  // K is a pointer to a dense matrix in column major
  // E, f, K can be nullptr. If it is, the function will not compute it.
  void ComputeElementLocalEnergyAndInternalForcesAndStiffnessMatrix(const double *u, int el, double * energy, double * fint, double * K);

protected:
  VolumetricMesh * volumetricMesh;
  StVKElementABCD * precomputedIntegrals;

  std::vector<double> gravityForce;
  double g;
  bool addGravity;  
  void InitGravity(); // aux function

  int numElementVertices;
  void ResetVector(double * vec); // aux function

  std::vector<double> lambdaLame;
  std::vector<double> muLame;
  std::vector<void*> internalElementData;

  int dof, dof2;

  inline void AddMatrix3x3Block(int c, int a, const Mat3d & matrix, double * K)
  {
    for (int k = 0; k < 3; k++) 
    {
      for (int l = 0; l < 3; l++)
      {
        int row = c * 3 + k;
        int col = a * 3 + l;

        K[col * dof + row] += matrix[k][l];
      }
    }
  }

  inline void AddMatrix3x3Block(int c, int a, const double matrix[9], double * K)
  {
    for (int k = 0; k < 3; k++) 
    {
      for (int l = 0; l < 3; l++)
      {
        int row = c * 3 + k;
        int col = a * 3 + l;

        K[col * dof + row] += matrix[k * 3 + l];
      }
    }
  }
};

#endif
