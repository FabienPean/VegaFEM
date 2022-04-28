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

#include "meshIntersection.h"
#include "basicAlgorithms.h"
#include "exactOctree.h"
#include "predicates.h"
#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif
using namespace std;

vector<vector<int>> computeTrianglesIntersectingEachTetExact(const TetMeshRef tetMesh, const TriMeshRef triMesh)
{
  ExactTriMeshOctree octree;
  octree.build(triMesh, 5, 10);

  vector<vector<int>> tetEmbedTri(tetMesh.numTets());

#ifdef USE_TBB
  tbb::parallel_for(tbb::blocked_range<int>(0, tetMesh.numTets()), [&](const tbb::blocked_range<int> & rng)
  {
    for (int tetID = rng.begin(); tetID != rng.end(); ++tetID)
    {
#else
    for(int tetID = 0; tetID < tetMesh.numTets(); ++tetID)
    {
#endif
      array<Vec3d, 4> tet;
      for(int j = 0; j < 4; j++) tet[j] = tetMesh.pos(tetID, j);
      BoundingBox tetbb(tet);

      auto toBB = [&](const BoundingBox & bb)
      {
        return (tetbb.intersect(bb));
      };
      auto toTri = [&](int tri)
      {
        return intersectTriTet(triMesh.pos(tri, 0), triMesh.pos(tri, 1), triMesh.pos(tri, 2),
            tet[0], tet[1], tet[2], tet[3]);
      };
      octree.rangeQuery(toBB, toTri, tetEmbedTri[tetID]);
      sortAndDeduplicate(tetEmbedTri[tetID]);
    }
#ifdef USE_TBB
  }, tbb::auto_partitioner()); //end for locations
#endif
  return tetEmbedTri;
}
