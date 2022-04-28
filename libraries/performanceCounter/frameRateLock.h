/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "performanceCounter" library , Copyright (C) 2018 USC                 *
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

#include "performanceCounter.h"
#include <cassert>

// used to lock frame-rate
class FrameRateLock
{
public:
  FrameRateLock(double targetFrameRate = 30.0) : frameRate(targetFrameRate), invFrameRate(1.0 / frameRate) { assert(frameRate > 0.0); }
  virtual ~FrameRateLock() {}

  // start internal counter
  void start() { counter.StartCounter(); }

  // called in openGL idle function, pause the code until target frame rate is reached
  inline void lock();

protected:
  double frameRate, invFrameRate;
  PerformanceCounter counter;
};


inline void FrameRateLock::lock()
{
  double elapsedTime = 0.0;
  do
  {
    counter.StopCounter();
    elapsedTime = counter.GetElapsedTime();
  }
  while (elapsedTime <= invFrameRate);

  counter.StartCounter();
}
