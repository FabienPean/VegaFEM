/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "animationHelper" library , Copyright (C) 2018 USC                    *
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

#include "simulationRecorder.h"
#include <iostream>
#include <cstdio>
#include "matrixIO.h"
using namespace std;

FrameRecorder::FrameRecorder(int size, const std::string & fn, int mf, int numSkippedIntervalFrames) :
  r(size), numFrames(0), filename(fn), maxFrames(mf), saved(false), numSkipIter(numSkippedIntervalFrames+1), skipCounter(numSkipIter) {}

void FrameRecorder::addFrame(const double * frame)
{
  if (skipCounter >= numSkipIter)
  {
    if (numFrames < maxFrames)
    {
      //memcpy(u_record + timestepCounter * 3 * n, integratorBase->Getq(), sizeof(double) * 3 * n);
      for(int i = 0; i < r; i++)
        u.push_back(frame[i]);
      numFrames++;
    }

    if (!saved && numFrames >= maxFrames)
      save();

    skipCounter = 0;
  }
  skipCounter++;
}

bool FrameRecorder::save()
{
  printf("Writing %d frames into %s... ", numFrames, filename.c_str());
  int ret = WriteMatrixToDisk(filename.c_str(), r, numFrames, &u[0]);
  if(ret == 0) {
    cout << "Done." << endl;
    saved = true;
    return true;
  }
  cout << "Failed." << endl;
  return false;
}

void FrameRecorder::clear()
{
  saved = false;
  u.clear();
  numFrames = 0;
  skipCounter = numSkipIter;
}

////////////////////////////////////////////////

SimulationRecorder::SimulationRecorder(size_t numVtx, std::string & fileBaseName, size_t maxFrames) :
uFrames(3*numVtx, fileBaseName + "d.u", maxFrames),
vFrames(3*numVtx, fileBaseName + "v.u", maxFrames),
fFrames(3*numVtx, fileBaseName + "f.u", maxFrames) {}

SimulationRecorder::~SimulationRecorder() {}

void SimulationRecorder::addFrame(const double * u, const double * v, const double * f)
{
  uFrames.addFrame(u);
  vFrames.addFrame(v);
  fFrames.addFrame(f);
}

bool SimulationRecorder::save()
{
  bool uret = uFrames.save();
  bool vret = vFrames.save();
  bool fret = fFrames.save();
  return uret && vret && fret;
}

void SimulationRecorder::clear()
{
  uFrames.clear();
  vFrames.clear();
  fFrames.clear();
}
