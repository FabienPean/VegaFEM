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

#ifndef STOPWATCH_H
#define STOPWATCH_H

/*
  === A stopwatch to measure code execution time. ===
  Modified from performanceCounter.h

  Designed for real-time system (e.g., real-time computer animation, haptics),
  but useful in general. You can time arbitrary segments of your code.
  Same interface under Windows, Linux and Mac OS X.

  Under Linux/MAC OS X, accuracy is that of gettimeofday function,
  which gives time in seconds and microseconds.
  In practice, it has been accurate down to microsecond range.

  Under Windows, the counter uses the QueryPerformanceCounter Windows API call.
  Again, accuracy has been in the microsecond range in practice.
*/

#if defined(_WIN32) || defined(WIN32)
#include <windows.h>
#include <cstdlib>
#include <cstdio>
#else
#include <cstdlib>
#include <cstdio>
#include "sys/time.h"
#endif

class StopWatch
{
public:
  // initialize a stopped watch 
  inline StopWatch();

  // start counting
  // have no effect if the watch has already started
  inline void start();

  // stop counting and store the elapsed time
  // can be resumed by calling start() again
  // have no effect if the watch has already stopped
  inline void stop();

  // stop the watch and reset internal timer
  inline void reset();

  // If the watch has stopped, return the stored time
  // otherwise, return the stored time + the time from the recent start() to this getElapsedTime() call
  inline double getElapsedTime();

  inline bool hasStarted() const { return stopped == 0; }
  inline bool hasStopped() const { return stopped; }

protected:
  inline double getTime();

#if defined(_WIN32) || defined(WIN32)
  LARGE_INTEGER timerFrequency;
  LARGE_INTEGER startCount,stopCount;
#else
  long startCountSec,stopCountSec,startCountMicroSec,stopCountMicroSec;
#endif
  double savedTime;
  int stopped;
};

// helper class for measuring the time of a scope
class FunctionTimer
{
public:
  FunctionTimer(StopWatch & w) : watch(&w) { watch->start(); }
  ~FunctionTimer() { watch->stop(); }
protected:
  StopWatch * watch;
};

// helper class to print time of a function
class ScopeTimeReporter
{
public:
  ScopeTimeReporter(const char * scopeName) : scopeName(scopeName) { watch.start(); }
  ~ScopeTimeReporter() { watch.stop(); printf("%s: %lfs\n", scopeName, watch.getElapsedTime());}
protected:
  const char * scopeName;
  StopWatch watch;
};

#define REPORT_FUNCTION_TIME ScopeTimeReporter __function_timer(__func__)


//////////////////////////////////////////////////////////
//                   IMPLEMENTATION
//////////////////////////////////////////////////////////

inline StopWatch::StopWatch()
{
#if defined(_WIN32) || defined(WIN32)
  QueryPerformanceFrequency(&timerFrequency);
  startCount.QuadPart = stopCount.QuadPart = 0;
#else
  startCountSec = stopCountSec = startCountMicroSec = stopCountMicroSec = 0.0;
#endif
  savedTime = 0.0;
  stopped = 1;
}

inline void StopWatch::start()
{
  if (stopped == 0)
    return;
  stopped = 0;

#if defined(_WIN32) || defined(WIN32)
  LARGE_INTEGER count;
  QueryPerformanceCounter(&count);
#else
  struct timeval tv;
  gettimeofday(&tv,NULL);
#endif
#if defined(_WIN32) || defined(WIN32)
  startCount = count;
#else
  startCountSec = tv.tv_sec;
  startCountMicroSec = tv.tv_usec;
#endif
}

inline double StopWatch::getElapsedTime() 
{
  if (stopped)
    return savedTime;

  double currentTime = getTime();
  return savedTime + currentTime;
}

inline void StopWatch::stop() 
{
  double currentTime = getTime();
  if (stopped)
    return;

  stopped = 1;
  savedTime += currentTime;
}

inline void StopWatch::reset()
{
  savedTime = 0.0;
  stopped = 1;
}

inline double StopWatch::getTime() 
{
#if defined(_WIN32) || defined(WIN32)
  QueryPerformanceCounter(&stopCount);
  return ((double)(stopCount.QuadPart - startCount.QuadPart))
        / ((double)timerFrequency.QuadPart);
#else
  struct timeval tv;
  gettimeofday(&tv,NULL);
  stopCountSec = tv.tv_sec;
  stopCountMicroSec = tv.tv_usec;
  return 1.0 * (stopCountSec-startCountSec) + 1E-6 * (stopCountMicroSec - startCountMicroSec);
#endif
}


#endif /* STOPWATCH_H */
