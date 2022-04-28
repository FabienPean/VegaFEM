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

#include "profiler.h"
#include "matrixIO.h"
#include <iomanip>
#include <sstream>
#include <cassert>
#include <cstring>
#include <algorithm>
using namespace std;


Profiler::Profiler()
{
  lastTimer = timer.end();
}

void Profiler::startTimer(const std::string & name)
{
  timer[name].start();
  lastTimer = timer.find(name);
}

void Profiler::stopTimer(const std::string & name)
{
  if (timer.find(name) != timer.end())
    timer[name].stop();
}

void Profiler::stopLastTimer()
{
  if (lastTimer != timer.end())
    lastTimer->second.stop();
}

void Profiler::startExtraTimer(const std::string & name)
{
  extraTimer[name].start();
}

void Profiler::stopExtraTimer(const std::string & name)
{
  if (extraTimer.find(name) != extraTimer.end())
    extraTimer[name].stop();
}

double Profiler::getTime(const std::string & name)
{
  if (timer.find(name) != timer.end())
    return timer[name].getElapsedTime();
  else
    return 0.0;
}

std::string Profiler::toString()
{
  const string title = " PROFILING ";
  const size_t titleLineLen = 30; // number of "=" in one side of ======= PROFILING =======
  const size_t maxSize = titleLineLen * 2 + title.size();

  map<string, double> times, extraTimes; // first record the timing on each timer
  double overallTime = 0.0;
  for(auto & p : timer)
  {
    double t = p.second.getElapsedTime();
    times[p.first] = t;
    overallTime += t;
  }
  for(auto & p: extraTimer) {
    extraTimes[p.first] = p.second.getElapsedTime();
  }

  ostringstream os;
  // print ======= PROFILING =======
  for(size_t i = 0; i < titleLineLen; i++)
    os << '=';
  os << title;
  for(size_t i = 0; i < titleLineLen; i++)
    os << '=';
  os << endl;

  // print each timer
  for(map<string, double>::iterator it = times.begin(); it != times.end(); it++)
  {
    ostringstream left, right;
    // left                                 right
    // <SECTION A>: 10.00 s                 5.00%
    left << it->first << ": " << fixed << setprecision(2) << it->second << " s";
    double per = it->second * 100 / overallTime;
    right << fixed << setprecision(2) << per << "%";
    string leftStr = left.str(), rightStr = right.str();

    int numSpace = 0; // number of space character to connect left and right parts
    if (leftStr.size() + rightStr.size() < maxSize)
      numSpace = maxSize - leftStr.size() - rightStr.size();

    os << leftStr;
    for(int i = 0; i < numSpace; i++)
      os << ' ';
    os << rightStr << endl;
  }
  if (extraTimer.size() > 0) {
    os << "extra ---\n";
    for(auto & p : extraTimes) {
      os << "  " << p.first << ": " << fixed << setprecision(2) << p.second << " s\n";
    }
  }
  os << "Overall: " << overallTime << " s" << endl;

  // print the ending =========================
  for(size_t i = 0; i < maxSize; i++) os << '=';
  os << endl;
  return os.str();
}

void Profiler::clear()
{
  timer.clear();
  lastTimer = timer.end();
}

////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////

IterationProfiler::IterationProfiler(int size, bool keep)
{
  assert(size > 0);
  bufferSize = size;
  numIterRecorded = 0;
  keepHistory = keep;
  numIterInHistory = 0;
}

void IterationProfiler::finishOneIteration()
{
  // add each stopWatch's time into averagingBuffer
  // notice that we need to ensure that each stopWatch has a one-to-one relationship with the averagingBuffers
  // because currently Profiler and IterationProfiler has no function to remove individual timers, 
  //   we are sure that each averagingBuffer corresponds to one stopWatch
  // but each stopWatch may not corresponds to an averagingBuffer
  for(map<string, StopWatch>::iterator it = timer.begin(); it != timer.end(); it++)
  {
    double t = it->second.getElapsedTime();
    it->second.reset();
    map<string, ProfilingData>::iterator foundBuffer = profilingData.find(it->first);
    if (foundBuffer == profilingData.end())
    {
      // because it's a new stopWatch, we create a new struct of ProfilingData.
      foundBuffer = profilingData.insert(pair<string, ProfilingData>(it->first, ProfilingData(bufferSize))).first;
      // since this profilingData is new and stores nothing on previous iterations, we should fill in previous timing of 0.0
      for(int i = 0; i < numIterRecorded; i++)
        foundBuffer->second.avgBuffer.addValue(0.0); // we only fill in numIterRecorded times
      if (keepHistory)
        foundBuffer->second.history.resize(numIterInHistory, 0.0); // if we keep history, we put all previous timing of 0.0 into history buffer
    }
    ProfilingData & data = foundBuffer->second;
    data.avgBuffer.addValue(t);
    if (data.maxTime < t)
      data.maxTime = t;
    if (keepHistory) // if we keep history, we add the new time
      data.history.push_back(t);
  }
  if (numIterRecorded < bufferSize)
    numIterRecorded++;

  if (keepHistory)
    numIterInHistory++;
}

void IterationProfiler::setAverageBufferSize(int newSize)
{
  assert(newSize > 0);
  if (bufferSize == newSize)
    return;

  bufferSize = newSize;
  for(map<string, ProfilingData>::iterator it = profilingData.begin(); it != profilingData.end(); it++)
    it->second.avgBuffer.setBufferSize(bufferSize);
  numIterRecorded = min(numIterRecorded, bufferSize);
}

void IterationProfiler::clear()
{
  profilingData.clear();
  numIterRecorded = 0;
  numIterInHistory = 0;
  Profiler::clear();
}

std::string IterationProfiler::toString()
{
  const string title = " PROFILING ";
  const size_t titleLineLen = 30;
  const size_t maxSize = titleLineLen * 2 + title.size();

  map<string, double> times; // first record the timing on each timer
  double overallTime = 0.0;
  for(map<string, ProfilingData>::iterator it = profilingData.begin(); it != profilingData.end(); it++)
  {
    double t = it->second.avgBuffer.getAverage() * 1000.0;  // unit to ms
    if (t < 0.0) t = 0.0;
    times[it->first] = t;
    overallTime += t;
  }


  ostringstream os;

  for(size_t i = 0; i < titleLineLen; i++)
    os << '=';
  os << title;
  for(size_t i = 0; i < titleLineLen; i++)
    os << '=';
  os << endl;


  
  for(map<string, double>::iterator it = times.begin(); it != times.end(); it++)
  {
    ostringstream left, right;
    double av = it->second;
    left << it->first << ": " << av << " ms, max: " << (profilingData.find(it->first)->second.maxTime*1000.0) << " ms";
    double per = av * 100 / overallTime;
    right << fixed << setprecision(2) << per << "%";
    string leftStr = left.str(), rightStr = right.str();

    int numSpace = 0;
    if (leftStr.size() + rightStr.size() < maxSize)
      numSpace = maxSize - leftStr.size() - rightStr.size();
    os << leftStr;
    for(int i = 0; i < numSpace; i++)
      os << ' ';
    os << rightStr << endl;
  }
  os << "Overall: " << overallTime << " ms" << endl;

  for(size_t i = 0; i < maxSize; i++)
    os << '=';
  os << endl;
  return os.str();
}

bool IterationProfiler::saveHistory(const std::string & filename) const
{
  if (keepHistory == false)
    return true;

  map<string, ProfilingData>::const_iterator it = profilingData.begin();
  if (it == profilingData.end())
    return true;
  int m = it->second.history.size();
  int n = profilingData.size();

  vector<double> matrix(m*n);
  int index = 0;
  for(map<string, ProfilingData>::const_iterator it = profilingData.begin(); it != profilingData.end(); it++)
  {
    memcpy(&matrix[m*index], it->second.history.data(), sizeof(double) * m);
    index++;
  }
  int ret = WriteMatrixToDisk(filename.c_str(), m, n, &matrix[0]);
  return ret == 0;
}
