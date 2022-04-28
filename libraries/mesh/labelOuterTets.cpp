/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "mesh" library , Copyright (C) 2018 USC                               *
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

#include "labelOuterTets.h"
#include <cassert>
using namespace std;

vector<bool> labelOuterTets(const TetMeshRef & tetMesh, const TetNeighbor & tetNeighbor,
    function<bool(int tetID)> isTetBoundary, function<bool(int tetID)> istetOuter)
{
  enum class TetLabel : char
  {
    IN,
    OUT,
    UNKNOWN
  };

  auto tetBoundaries = tetNeighbor.findTetBoundaries(tetMesh.numTets(), tetMesh.tets());
  vector<int> airTetIDs;

  vector<TetLabel> tetLabel(tetMesh.numTets(), TetLabel::UNKNOWN);

  //  airProfiler.startTimer("initialAirTetLabel");
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    if (isTetBoundary(tetID))
    {
      tetLabel[tetID] = TetLabel::IN;
    }
  }
  for(auto p : tetBoundaries)
  {
    int tetID = p.first;
    assert(tetID >= 0 && tetID < tetMesh.numTets());
    if (tetLabel[tetID] == TetLabel::UNKNOWN)
    {
      airTetIDs.push_back(tetID);
      tetLabel[tetID] = TetLabel::OUT;
    }
  }
  //  airProfiler.stopLastTimer();

//  cout << "Try finding air tets..." << endl;
  size_t candidateBegin = 0, candidateEnd = airTetIDs.size();

  auto floodFillAir = [&]()
  {
    while(candidateBegin != candidateEnd)
    {
      for(size_t i = candidateBegin; i < candidateEnd; i++)
      {
        int tetID = airTetIDs[i];
        for(int nbr : tetNeighbor.getTetNeighbors(tetID))
        {
          if (nbr < 0) continue;
          assert(nbr < tetMesh.numTets());
          if (tetLabel[nbr] != TetLabel::UNKNOWN) continue;
          tetLabel[nbr] = TetLabel::OUT;
          airTetIDs.push_back(nbr);
        }
      }
      candidateBegin = candidateEnd;
      candidateEnd = airTetIDs.size();
    }
  };

  auto floodFillInterior = [&](int seedTetID)
  {
    vector<int> buffer = { seedTetID }, nextBuffer;
    tetLabel[seedTetID] = TetLabel::IN;
    while(buffer.size() > 0)
    {
      for(int tetID : buffer)
      {
        for(int nbr : tetNeighbor.getTetNeighbors(tetID))
        {
          if (nbr < 0) continue;
          assert(nbr < tetMesh.numTets());
          if (tetLabel[nbr] != TetLabel::UNKNOWN) continue;
          tetLabel[nbr] = TetLabel::IN;
          nextBuffer.push_back(nbr);
        }
      }
      buffer.swap(nextBuffer);
      nextBuffer.clear();
    }
  };
  //  airProfiler.startTimer("initialFloodFill");
  floodFillAir();

  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
  {
    if (tetLabel[tetID] != TetLabel::UNKNOWN) continue;
    if (istetOuter(tetID)) // outside
    {
      candidateBegin = airTetIDs.size();
      airTetIDs.push_back(tetID);
      tetLabel[tetID] = TetLabel::OUT;
      candidateEnd = airTetIDs.size();
      floodFillAir();
    }
    else // inside
    {
      floodFillInterior(tetID);
    }
  }

  vector<bool> ret(tetMesh.numTets());
  for(int tetID = 0; tetID < tetMesh.numTets(); tetID++)
    ret[tetID] = (tetLabel[tetID] == TetLabel::OUT);

  return ret;
}
