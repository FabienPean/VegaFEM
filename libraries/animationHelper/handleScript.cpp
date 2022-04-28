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

#include "handleScript.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
using namespace std;

const string HandleScript::commandName[HandleScript::NUM_COMMANDS] = {
  "VOID",
  "ADD",
  "REMOVE",
  "MOVE",
  "SPECIAL",
  "END"
};


const HandleScript::Command HandleScript::endCommand(END);

#define CHECK_FAILURE(cond) \
  if (cond) \
  { \
    cerr << "Error HandleScript file format at line " << lineCount << endl; \
    throw 2; \
  }

HandleScript::HandleScript() : nextCommand(0)
{
}

HandleScript::HandleScript(const char * filename) : nextCommand(0)
{
  ifstream fin(filename);
  if (!fin)
  {
    cerr << "Cannot open HandleScript file " << filename << endl;
    throw 1;
  }

  string line;
  fin >> ws;

  string typeName;
  int lineCount = 0;
  while(!fin.eof())
  {
    lineCount++;
    getline(fin, line);
    if (line.size() == 0)
      continue;
    if (line[0] == '#')
      continue;

    Command c;
    istringstream is(line);
    typeName.clear();
    is >> typeName;
    CHECK_FAILURE(is.fail());

    if (typeName == "VOID")
    {
      // do nothing
    }
    else if (typeName == "END")
    {
      c.type = END;
    }
    else if (typeName == "SPECIAL")
    {
      c.type = SPECIAL;
      is >> c.vtx;
      if (is.fail()) 
        c.vtx = 0;
    }
    else if (typeName == "ADD")
    {
      c.type = ADD;
      is >> c.vtx;
      CHECK_FAILURE(is.fail());
      CHECK_FAILURE(c.vtx < 0);
    }
    else if (typeName == "REMOVE")
    {
      c.type = REMOVE;
      is >> c.vtx;
      CHECK_FAILURE(is.fail());
      CHECK_FAILURE(c.vtx < 0);
    }
    else if (typeName == "MOVE")
    {
      c.type = MOVE;
      is >> c.vtx;
      CHECK_FAILURE(is.fail());
      CHECK_FAILURE(c.vtx < 0);
      is >> ws;
      for(int i = 0; i < 3; i++)
      {
        char s = is.peek();
        if (s == ',')
          is >> s;
        is >> c.vec[i];
        CHECK_FAILURE(is.fail());
      }
    }

    commands.push_back(c);

    fin >> ws;
  }
}

const HandleScript::Command & HandleScript::getNextCommand()
{
  if (nextCommand >= commands.size())
    return endCommand;

  return commands[nextCommand++];
}

void HandleScript::addCommand(const Command & c)
{
  commands.push_back(c);
}

bool HandleScript::save(const char * filename)
{
  ofstream fout(filename);
  if(!fout)
    return false;

  for(size_t i = 0; i < commands.size(); i++)
  {
    Command & c = commands[i];
    CommandType type = c.type;
    assert(type < NUM_COMMANDS);
    fout << commandName[(int)type]; //output command name
    if (type == ADD || type == REMOVE || type == SPECIAL)
    {
      fout << " " << c.vtx;
    }
    else if (type == MOVE)
    {
      fout << " " << c.vtx << " " << c.vec[0] << ", " << c.vec[1] << ", " << c.vec[2];
    }
    fout << endl;
    if (fout.fail())
      return false;
  }

  fout.close();
  return true;
}
