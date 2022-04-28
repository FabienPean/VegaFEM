/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "basicAlgorithms" library , Copyright (C) 2018 USC                    *
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

#include "fileIO.h"
#include "stringHelper.h"
#include <cstring>
#include <iostream>
using namespace std;

int readEachLine(const std::string & filename, std::function<int(std::string & line)> processLine, const char * comment)
{
  fstream fin(filename.c_str());
  if (!fin) return 1;
  return readEachLine(fin, processLine, comment);
}

int readEachLine(std::istream & fin, std::function<int(std::string & line)> processLine, const char * comment)
{
  fin >> std::ws;
  std::string line;
  int commentStrLen = (comment ? strlen(comment) : 0);
  while(fin.eof() == false)
  {
    line.clear();
    std::getline(fin, line);
    fin >> std::ws;
    line = strip(line);
    if (line.size() == 0) continue;
    if (comment && strncmp(line.c_str(), comment, commentStrLen) == 0) continue;
    int ret =processLine(line);
    if (ret != 0) return ret;
  }
  return 0;
}
