/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "exactArithmetic" library , Copyright (C) 2018 USC                    *
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

#include "vec3ER.h"

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

/////////////////////////////////////////////////////////////////
//                             V3
/////////////////////////////////////////////////////////////////

const Vec3ER VEC3ER_NULL(0.0, 0.0, 0.0);

#ifdef VEGA_USE_CGAL_HEADER
template <class CGAL_Kernel_Vector_3>
Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal)
{
  Vec3ER ret;
  ret[0] = cgal[0];
  ret[1] = cgal[1];
  ret[2] = cgal[2];
  return ret;
}

template Vec3ER assignCGALToVec3ER<K::Vector_3>(const K::Vector_3 & cgal);


bool intersectTriTri(const Vec3ER & pa, const Vec3ER & pb, const Vec3ER & pc, const Vec3ER & qa, const Vec3ER & qb, const Vec3ER & qc)
{
  K::Point_3 p0(pa[0], pa[1], pa[2]);
  K::Point_3 p1(pb[0], pb[1], pb[2]);
  K::Point_3 p2(pc[0], pc[1], pc[2]);

  K::Point_3 q0(qa[0], qa[1], qa[2]);
  K::Point_3 q1(qb[0], qb[1], qb[2]);
  K::Point_3 q2(qc[0], qc[1], qc[2]);

  K::Triangle_3 p(p0, p1, p2), q(q0, q1, q2);
  auto result =  CGAL::intersection(p, q);
  if(result) return true;
  return false;
}

#else


namespace
{
//  const K::RT & RT(const void * p) {
//    return *(const K::RT*)p;
//  }
  K::RT & RT(void * p) {
    return *(K::RT*)p;
  }
}


template <class CGAL_Kernel_Vector_3>
Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal)
{
  ER::UninitializedLabel l;
  Vec3ER ret(l);
  ret.elt[0].p = new K::RT(cgal[0]);
  ret.elt[1].p = new K::RT(cgal[1]);
  ret.elt[2].p = new K::RT(cgal[2]);
  return ret;
}

template Vec3ER assignCGALToVec3ER<K::Vector_3>(const K::Vector_3 & cgal);

bool intersectTriTri(const Vec3ER & pa, const Vec3ER & pb, const Vec3ER & pc, const Vec3ER & qa, const Vec3ER & qb, const Vec3ER & qc)
{
  K::Point_3 p0(RT(pa[0].p), RT(pa[1].p), RT(pa[2].p));
  K::Point_3 p1(RT(pb[0].p), RT(pb[1].p), RT(pb[2].p));
  K::Point_3 p2(RT(pc[0].p), RT(pc[1].p), RT(pc[2].p));

  K::Point_3 q0(RT(qa[0].p), RT(qa[1].p), RT(qa[2].p));
  K::Point_3 q1(RT(qb[0].p), RT(qb[1].p), RT(qb[2].p));
  K::Point_3 q2(RT(qc[0].p), RT(qc[1].p), RT(qc[2].p));

  K::Triangle_3 p(p0, p1, p2), q(q0, q1, q2);
  auto result =  CGAL::intersection(p, q);
  if(result) return true;
  return false;
}
#endif

bool intersectTriAABB(const Vec3ER & ta, const Vec3ER tb, const Vec3ER tc, const Vec3ER bmin, const Vec3ER bmax)
{
  assert(bmin[0] <= bmax[0] && bmin[1] <= bmax[1] && bmin[2] <= bmax[2]);
  if (bmin[0] <= ta[0] && ta[0] <= bmax[0] && bmin[1] <= ta[1] && ta[1] <= bmax[1] && bmin[2] <= ta[2] && ta[2] <= bmax[2]) return true;
  if (bmin[0] <= tb[0] && tb[0] <= bmax[0] && bmin[1] <= tb[1] && tb[1] <= bmax[1] && bmin[2] <= tb[2] && tb[2] <= bmax[2]) return true;
  if (bmin[0] <= tc[0] && tc[0] <= bmax[0] && bmin[1] <= tc[1] && tc[1] <= bmax[1] && bmin[2] <= tc[2] && tc[2] <= bmax[2]) return true;
  // create a box mesh for AABB
  const Vec3ER p[8] = { {bmin[0], bmin[1], bmin[2]}, { bmax[0], bmin[1], bmin[2] }, { bmax[0], bmax[1], bmin[2] }, { bmin[0], bmax[1], bmin[2] },
                  { bmin[0], bmin[1], bmax[2] }, { bmax[0], bmin[1], bmax[2] }, {bmax[0], bmax[1], bmax[2]}, { bmin[0], bmax[1], bmax[2] } };
  const int t[12][3] = { {0,3,2}, {0,2,1}, {4,5,6}, {4,6,7}, {0,1,5}, {0,5,4},
                                 {3,7,6}, {3,6,2}, {1,2,6}, {1,6,5}, {0,4,7}, {0,7,3} };
  for(int i = 0; i < 12; i++) {
    if (intersectTriTri(ta, tb, tc, p[t[i][0]], p[t[i][1]], p[t[i][2]])) return true;
  }
  return false;
}

