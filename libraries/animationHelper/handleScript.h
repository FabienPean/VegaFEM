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

#ifndef HANDLESCRIPT_H
#define HANDLESCRIPT_H

#include "vec3d.h"
#include <vector>
#include <string>

// a simple class to load/save a script of handle creation/deletion/movement

class HandleScript
{
public:

  enum CommandType
  {
    VOID, ADD, REMOVE, MOVE, SPECIAL, END, NUM_COMMANDS
  };

  struct Command
  {
    CommandType type;
    int vtx;
    Vec3d vec;
    Command(CommandType t=VOID) : type(t), vtx(0), vec(0.0) {}
  };

  HandleScript();
  HandleScript(const char * filename);
  virtual ~HandleScript() {}

  const Command & getNextCommand();

  void addCommand(const Command & c);
  bool save(const char * filename);

  static const std::string commandName[NUM_COMMANDS];

protected:
  std::vector<Command> commands;
  size_t nextCommand;
  const static Command endCommand;
};


#endif
