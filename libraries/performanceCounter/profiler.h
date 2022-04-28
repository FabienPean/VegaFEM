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

#ifndef PROFILER_H
#define PROFILER_H

#include "stopWatch.h"
#include "averagingBuffer.h"
#include <map>
#include <string>

// a profiling tool that can record timings in multiple parts of the code
class Profiler
{
public:
  Profiler();
  virtual ~Profiler() {}

  // start the timer on one section
  // it does nothing if the timer has already started
  void startTimer(const std::string & sectionName);

  // stop the timer on one section
  // it does nothing if there is no such section or the timer has already stopped
  // when you stop the timer and start it again, new time is added onto the previous time measured on the last stopTimer call()
  void stopTimer(const std::string & sectionName);
  // stop the timer that started last
  void stopLastTimer();

  // return the elapsed time on one section
  // return 0.0 if that section does not exist
  double getTime(const std::string & sectionName);

  void startExtraTimer(const std::string & sectionName);
  void stopExtraTimer(const std::string & sectionName);

  // clear all stored timers
  void clear();

  // return a string that writes the timing
  // the format is:
  // =============== PROFILING ================
  // <SECTION A>: 10.00 s                 5.00%
  // <SECTION B>: 20.00 s                10.00%
  // <SECTION C>: 170.00 s               85.00%
  // Overall: 200.00 s
  // ==========================================
  std::string toString();

protected:
  std::map<std::string, StopWatch> timer;
  std::map<std::string, StopWatch>::iterator lastTimer;
  std::map<std::string, StopWatch> extraTimer;
};

// a profiling tool that records timing on multiple sections in a iterative program
class IterationProfiler : public Profiler
{
public:
  IterationProfiler(int averageBufferSize = 20, bool keepHistory = false);
  virtual ~IterationProfiler() {}

  // tell the profiler that one iteration is finished
  // calling this function will reset all timers and their time is stored in 
  // corresponding averaging buffer and history buffer (if keepHistory is on) and maxTime is updated
  virtual void finishOneIteration();
  // set the number of recent iterations used to calculate average timing; default is 20
  void setAverageBufferSize(int bufferSize);

  // clear all stored data
  void clear();

  // return a string that writes the timing
  // the format is:
  // =============== PROFILING ================
  // <SECTION A>: 10.00 ms                 5.00%
  // <SECTION B>: 20.00 ms                10.00%
  // <SECTION C>: 170.00 ms               85.00%
  // Overall: 200.00 ms
  // ==========================================
  std::string toString();

  bool saveHistory(const std::string & filename) const;

protected:
  int bufferSize;
  struct ProfilingData
  {
    ProfilingData(int bufferSize) : avgBuffer(bufferSize), maxTime(0) {}
    AveragingBuffer avgBuffer;
    double maxTime;
    std::vector<double> history;
  };
  std::map<std::string, ProfilingData> profilingData;

  bool keepHistory;
  int numIterInHistory;
  int numIterRecorded;
};

// helper profiler measure timing on a section
// by using destructor function, it stops the timer whenever it goes out of the code block

class ProfilerSection
{
public:
  ProfilerSection(Profiler * profiler, const std::string & sectionName) : p(profiler) { if (p) { name = sectionName; } }
  virtual ~ProfilerSection() { stop(); }

  void start() { if (p) p->startTimer(name); }
  void stop() { if (p) p->stopTimer(name); }

protected:
  Profiler * p;
  std::string name;
};

class ProfilerExtraSection
{
public:
  ProfilerExtraSection(Profiler * profiler, const std::string & sectionName) : p(profiler) { if (p) { name = sectionName; } }
  virtual ~ProfilerExtraSection() { stop(); }

  void start() { if (p) p->startExtraTimer(name); }
  void stop() { if (p) p->stopExtraTimer(name); }

protected:
  Profiler * p;
  std::string name;
};
#endif
