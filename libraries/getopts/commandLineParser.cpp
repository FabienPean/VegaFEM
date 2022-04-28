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

#include "commandLineParser.h"
#include <iostream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
using namespace std;

CommandLineParser::CommandLineParser() {}

void CommandLineParser::addOption(const string & name, int & value)
{
  entries.push_back(Entry(name, INT, (void*)&value));
}

void CommandLineParser::addOption(const string & name, char * value)
{
  entries.push_back(Entry(name, CSTR, (void*)value));
}

void CommandLineParser::addOption(const string & name, bool & value)
{
  entries.push_back(Entry(name, BOOL, (void*)&value));
}

void CommandLineParser::addOption(const std::string & name, vegalong & value)
{
  entries.push_back(Entry(name, LONG, (void*)&value));
}

void CommandLineParser::addOption(const std::string & name, float & value)
{
  entries.push_back(Entry(name, FLOAT, (void*)&value));
}

void CommandLineParser::addOption(const std::string & name, double & value)
{
  entries.push_back(Entry(name, DOUBLE, (void*)&value));
}

void CommandLineParser::addOption(const std::string & name, std::string & value)
{
  entries.push_back(Entry(name, STRING, (void*)&value));
}

int CommandLineParser::parse(int argc, char ** argv, int numSkipArg)
{
  return parse(argc, (const char**)argv, numSkipArg);
}

int CommandLineParser::parse(int argc, const char ** argv, int numSkipArg)
{
  for (int i = numSkipArg; i < argc; i++)
  {
    const char * arg = argv[i];
    // option begins with '-' or '/'
    if ((*arg != '-') && (*arg != '/'))
      return (i);
    arg++; // skip '-' or '/'

    for (size_t j = 0; j < entries.size(); j++) //find options for this arg
    {
      //options are case sensitive!!
      size_t l = entries[j].name.size();

      //printf("%s\n", opttable[j].sw);
      if (strncmp(arg, entries[j].name.data(), l) == 0) //found this option
      {
        const char * var = nullptr; //pointers to the input data for this option
//        bool standAloneData = false; // whether the data for this arg is in a separate string, like -v 123

        // There's ambiguity when a negative number arg is next to a switch, like -f -48.9
        // For now we only intepret it as two successive switches
        // A smarter implementation could be to forbid switch name to start with numbers, and check whether the second "switch"
        // starts with "-\d"

        // if a data string follows this arg string, like: -v 123
        if ( (strlen(arg) == l) && ( (i+1 < argc) && (argv[i+1] != nullptr) && ( argv[i+1][0] != '-' && argv[i+1][0] != '/' ) ) )
        {
          // go to the next data string
          i++;
          var = argv[i];
//          standAloneData = true;
        }
        else
        {
          //copy the rest in this arg string, like: -v123
          if ( ( *(arg+l)=='=' ) || (*(arg+l)==' ') ) // accept '=' or ' ' after the switch
            var = arg + l + 1;
          else
            var = arg + l;
        }

        // if (*var == '\0')
        //   continue;
        //printf("%d %s\n",i, var);
        switch (entries[j].type)
        {
        case INT :
          if (*var != '\0')
          {
            *((int *)entries[j].value) = (int)strtol(var,nullptr,10);
            if (errno == ERANGE)
              return (i);
          }
          else // the data is absent. We consider it an error
            return i;
          break;

        case LONG :
          if (*var != '\0')
          {
            *((vegalong *)entries[j].value) = strtol(var,nullptr,10);
            if (errno == ERANGE)
              return (i);
          }
          else
            return i;
          break;

        case FLOAT :
          if (*var != '\0')
          {
            *((float *)entries[j].value) = (float)strtod(var,nullptr);
            if (errno == ERANGE)
              return (i);
          }
          else
            return i;
          break;

        case DOUBLE :
          if (*var != '\0')
          {
            *((double *)entries[j].value) = strtod(var,nullptr);
            if (errno == ERANGE)
              return (i);
          }
          else
            return i;
          break;

          // case OPTSUBOPT :
          //   //option with more than 1 char as suboptions:
          //   // /oail or /op
          //   strcpy((char *)opttable[j].var, arg+l);
          //   arg += strlen(arg)-1;
          //   break;

        case CSTR :
          if (*var != '\0')
            strcpy((char *)entries[j].value, var);
          else
            return i;
          break;

        case STRING :
          if (*var != '\0')
            *((string*)(entries[j].value)) = var;
          else
            return i;
          break;

        case BOOL :
          //check a + or - after the sw: "/a- /b+"
          {
            bool boolValue = true;
            if ( *var == '-' )
              boolValue = false;
            else
              boolValue = true;

            *((bool *)entries[j].value) = boolValue;
            break;
          }
        }
        break;      //break the for
      } // end if (strncmp(arg, entries[j].name.data(), l) == 0) //found this option
    } // end for (size_t j = 0; j < entries.size(); j++)
  }
  return argc;
}

CommandLineParser::Entry::Entry(const std::string & n, Type t, void * v) : name(n), type(t), value(v) {}
