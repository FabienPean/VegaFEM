/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "massSpringSystem" library, Copyright (C) 2007 CMU, 2009 MIT,         *
 *                                           2018 USC                    *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic, Daniel Schroeder                          *
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

#ifndef VEGAFEM_MASSSPRINGSYSTEMFROMCUBEMESHCONFIGFILE_H
#define VEGAFEM_MASSSPRINGSYSTEMFROMCUBEMESHCONFIGFILE_H

#include "massSpringSystem.h"

namespace vegafem
{

class MassSpringSystemCubicMeshConfiguration
{
public:
  MassSpringSystemCubicMeshConfiguration() : cubicMeshFilename(NULL), surfaceMeshFilename(NULL) {}
  ~MassSpringSystemCubicMeshConfiguration() { free(cubicMeshFilename); free(surfaceMeshFilename); }

  char * cubicMeshFilename;
  char * surfaceMeshFilename;
  double density;
  double tensileStiffness;
  double damping;
  int addGravity;
};

class MassSpringSystemFromCubicMeshConfigFile
{
public:

  // generates a mass spring system from the given configuration file
  // if massSpringSystemConfiguration is not NULL, it also returns the mass spring parameters read from the file (massSpringSystemConfiguration is output parameter only)
  int GenerateMassSpringSystem(const char * configFilename, MassSpringSystem ** massSpringSystem, MassSpringSystemCubicMeshConfiguration * massSpringSystemCubicMeshConfiguration = NULL);

protected:
  char * DuplicateString(const char*);
};


}//namespace vegafem

#endif


