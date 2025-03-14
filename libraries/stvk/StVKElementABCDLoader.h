/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "StVK" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
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
  A loader class for StVK cubic (unreduced) coefficients.
  See also StVKReducedInternalForces.h .
*/

#ifndef VEGAFEM_STVKELEMENTABCDLOADER_H
#define VEGAFEM_STVKELEMENTABCDLOADER_H

#include "volumetricMesh.h"
#include "StVKElementABCD.h"

namespace vegafem
{

class StVKElementABCDLoader
{
public:
  // creates the StVKElementABCD coefficients, given a volumetric mesh
  // loadingFlag: 
  //   0 : use the low-memory version (default)
  //   1 : use the high-memory version (only applies with tet meshes); with this setting, computation speeds will be higher, at the expense of more memory (however, difference is typically not large and speeds might even decrease with large meshes when running out of memory)
  static StVKElementABCD * load(VolumetricMesh * volumetricMesh, unsigned int loadingFlag=0); 
};


}//namespace vegafem

#endif


