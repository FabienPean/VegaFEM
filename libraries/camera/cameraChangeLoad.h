/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "camera" library , Copyright (C) 2018 USC                             *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Jernej Barbic                                *
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

#ifndef CAMERACHANGELOAD_H
#define CAMERACHANGELOAD_H

#include "camera.h"
#include <string>
#include <vector>

// load a script to control camera from disk

// file format
// each line is:
// <frameID> <dr> <dphi> <dtheta> (<duration>)
// comment start with '#'
// dr: change of spherical camera radius
// dphi: change of spherical camera phi angle, the longitude angle
// dtheta: change of spherical camera theta angle, the lattitude angle
// duration: the length of the change, default is 1 frame

class CameraChangeLoad
{
public:
  CameraChangeLoad() {}
  
  CameraChangeLoad(const std::string & filename);

  int load(const std::string & filename);

  void controlCamera(int frameID, SphericalCamera & camera) const;

protected:
  struct Entry
  {
    Entry() {}
    int frameID = 0;
    double dr = 0.0, dphi = 0.0, dtheta = 0.0;
    int duration = 1;
  };
  std::vector<Entry> entries;
};

#endif
