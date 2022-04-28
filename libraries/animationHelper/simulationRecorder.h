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

#ifndef SIMULATIONRECORDER_H
#define SIMULATIONRECORDER_H
#include <climits>
#include <string>
#include <vector>

// record frames of data and save to disk
class FrameRecorder
{
public:
  // r: frame size, 
  // If numSkippedIntervalFrames == 0, every frame passed to addFrame() is stored.
  // otherwise, the recorder skips numSkippedIntervalFrames before accepts a frame in addFrame()
  FrameRecorder(int r, const std::string & filename, int maxFrames_ = INT_MAX, int numSkippedIntervalFrames = 0);

  // add one frame and save all the frames to disk if maxFrames is reached
  // it will not accept more frames after it has saved the frames to disk
  // u may not be stored depending on numSkippedIntervalFrames and maxFrames
  void addFrame(const double * u);

  // save the current frames to disk
  // return false if it cannot be saved
  bool save();

  void setFilename(const std::string & newName) { filename = newName; }

  void setMaxFrames(int f) { maxFrames = f; }

  // return true if it has been saved 
  bool hasSaved() const { return saved;}

  // restart the recorder. Clear all the stored frames and reset saved flag
  void clear();

  // get # frames stored so far
  int getNumFrames() const { return numFrames; }

protected:
  int r; // size of each frame
  int numFrames;
  std::string filename;
  int maxFrames;
  std::vector<double> u;
  bool saved;
  int numSkipIter;
  int skipCounter;
};

// save displacement (u), velocity (v) and force (f) to disk and name them as <basename>d.u, <basename>v.u and <basename>f.u

class SimulationRecorder
{
public:
  SimulationRecorder(size_t numVtx, std::string & basename, size_t maxFrames = INT_MAX);
  virtual ~SimulationRecorder();

  void addFrame(const double * u, const double * v, const double * f);
  size_t getNumFrames() const { return uFrames.getNumFrames(); }

  bool save();
  bool hasSaved() const { return uFrames.hasSaved(); }

  void clear();

protected:

  FrameRecorder uFrames, vFrames, fFrames;
};

#endif
