/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "getopts" library , Copyright (C) 2018 USC                            *
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

#ifndef COMMANDLINEPARSER_H
#define COMMANDLINEPARSER_H

#include "vegalong.h"
#include <string>
#include <vector>

// parser for command line arguments
// same argument syntax as in getopts.h:
//
// switches come with either "-" or "/", optionally followed by "=". e.g. -a  /h  -x=
// the parser accepts:
//
// string:
//   -s abc
//   -sabc
// whitespace must follow a string argument, no -sstring-a is allowed
// number
//   -a123.4
//   -a 123.4
//   -v-57
// there must be no whitespace between the switch and a negative number, no -v -57 is allowed
// boolean
//   -b+
//   -b
//   -b-
// the parser accepts a "false" boolean value only when a '-' after the switch without whitespace, like -b-
// all the other situations are considered as "true"

class CommandLineParser
{
public:
  CommandLineParser();
  virtual ~CommandLineParser() {}

  void addOption(const std::string & name, int & value);
  void addOption(const std::string & name, char * value);
  void addOption(const std::string & name, bool & value);
  void addOption(const std::string & name, vegalong & value);
  void addOption(const std::string & name, float & value);
  void addOption(const std::string & name, double & value);
  void addOption(const std::string & name, std::string & value);

  // parse args from argv. numSkipArg tells how many entries in argv are skipped before args are reached
  // on success, return argc
  // on failure, return the 0-index of the arg that caused failure
  int parse(int argc, char ** argv, int numSkipArg = 1);
  int parse(int argc, const char ** argv, int numSkipArg = 1);

protected:

  enum Type
  {
    INT = 1,
    CSTR = 2,
    BOOL = 3,
    LONG = 4,
    FLOAT = 5,
    DOUBLE = 6,
    STRING = 7
  };

  struct Entry
  {
    std::string name;
    Type type;
    void * value;
    Entry(const std::string & name, Type type, void * value);
  };
  std::vector<Entry> entries;
};




#endif
