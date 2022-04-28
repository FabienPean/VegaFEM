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

#include "cameraChangeLoad.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;

#define CHECK_FAILURE(cond) \
  if (cond) \
  { \
    cerr << "Error CameraChangeLoad file format at line " << lineCount << endl; \
    return 2; \
  }

CameraChangeLoad::CameraChangeLoad(const std::string & filename)
{
  if (load(filename) != 0) throw 1;
}

int CameraChangeLoad::load(const std::string & filename)
{
  ifstream fin(filename);
  if (!fin)
  {
    cerr << "Cannot open HandleScript file " << filename << endl;
    return 1;
  }

  string line;
  fin >> ws;

  int lastTime = 0;
  int lineCount = 0;
  while(!fin.eof())
  {
    lineCount++;
    getline(fin, line);
    if (line.size() == 0)
      continue;
    if (line[0] == '#')
      continue;

    Entry c;
    istringstream is(line);
    is >> c.frameID;
    CHECK_FAILURE(is.fail());
    if (c.frameID < lastTime)
    {
      cerr << "Error CameraChangeLoad file format at line " << lineCount << ", time becomes smaller" << endl;
      return 2;
    }
    is >> c.dr >> c.dphi >> c.dtheta >> ws;
    CHECK_FAILURE(is.fail());
    if (!is.eof()) is >> c.duration;
    CHECK_FAILURE(is.fail() || c.duration <= 0);

    lastTime = c.frameID+c.duration;
    c.dphi *= M_PI / 180.0;
    c.dtheta *= M_PI / 180.0;

    entries.push_back(c);

    fin >> ws;
  }
  return 0;
}

void CameraChangeLoad::controlCamera(int timestepCount, SphericalCamera & camera) const
{
  int curIndex = 0;
  for(; curIndex < (int)entries.size(); curIndex++) {
    const auto & e = entries[curIndex];
    if (e.frameID <= timestepCount && timestepCount < e.frameID+e.duration) {
      if (e.dphi != 0.0) camera.MoveRight(e.dphi);
      if (e.dr != 0.0) camera.ZoomIn(-e.dr);
      if (e.dtheta != 0.0) camera.MoveUp(e.dtheta);
    }
  }
}

