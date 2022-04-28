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

#include "predicates.h"
#include <cassert>
#include <functional>
#include <cstring>
#include <algorithm>
using namespace std;

// use predicates from Shewchuk's predicates
extern "C" void exactinit();
extern "C" double orient2d(double* pa, double *pb, double *pc);
extern "C" double orient3d(double* pa, double *pb, double *pc, double *pd);
extern "C" double incircle(double* pa, double *pb, double *pc, double *pd);
extern "C" double insphere(double* pa, double *pb, double *pc, double *pd, double *pe);

extern "C" double orient2dfast(double* pa, double *pb, double *pc);
extern "C" double orient3dfast(double* pa, double *pb, double *pc, double *pd);
extern "C" double incirclefast(double* pa, double *pb, double *pc, double *pd);
extern "C" double inspherefast(double* pa, double *pb, double *pc, double *pd, double *pe);


// The following routines are based on Robert Bridson's Tunicate library

namespace
{
bool same_sign(double a, double b)
{
   return (a<=0 && b<=0) || (a>=0 && b>=0);
}
}

int simplex_intersection2d(int k, const double* x0, const double* x1, const double* x2, const double* x3,
    double* alpha0, double* alpha1, double* alpha2, double* alpha3);

int simplex_intersection3d(int k, const double* x0, const double* x1, const double* x2, const double* x3, const double* x4, double* alpha0,
    double* alpha1, double* alpha2, double* alpha3, double* alpha4);

// degenerate test in 3d - assumes four points lie on the same plane
int simplex_intersection3d(int k, const double* x0, const double* x1, const double* x2, const double* x3,
    double* alpha0, double* alpha1, double* alpha2, double* alpha3);

// -------------------------------------------------------------

void initPredicates()
{
  exactinit();
}

double orient2d(const double pa[2], const double pb[2], const double pc[2])
{
  return orient2d((double*)pa, (double*)pb, (double*)pc);
}

double orient3d(const double pa[3], const double pb[3], const double pc[3], const double pd[3])
{
  return orient3d((double*)pa, (double*)pb, (double*)pc, (double*)pd);
}

double incircle(const double pa[2], const double pb[2], const double pc[2], const double pd[2])
{
  return incircle((double*)pa, (double*)pb, (double*)pc, (double*)pd);
}

double insphere(const double pa[3], const double pb[3], const double pc[3], const double pd[3], const double pe[3])
{
  return insphere((double*)pa, (double*)pb, (double*)pc, (double*)pd, (double*)pe);
}

double orient2dfast(const double pa[2], const double pb[2], const double pc[2])
{
  return orient2dfast((double*)pa, (double*)pb, (double*)pc);
}

double orient3dfast(const double pa[3], const double pb[3], const double pc[3], const double pd[3])
{
  return orient3dfast((double*)pa, (double*)pb, (double*)pc, (double*)pd);
}

double incirclefast(const double pa[2], const double pb[2], const double pc[2], const double pd[2])
{
  return incirclefast((double*)pa, (double*)pb, (double*)pc, (double*)pd);
}

double inspherefast(const double pa[3], const double pb[3], const double pc[3], const double pd[3], const double pe[3])
{
  return inspherefast((double*)pa, (double*)pb, (double*)pc, (double*)pd, (double*)pe);
}

// osega = orient3d(tria, trib, tric, sega)
// osegb = orient3d(tria, trib, tric, segb)
bool intersectSegTriWithOrientOnSeg(const double sega[3], const double segb[3], const double trip[3], const double triq[3], const double trir[3],
    double osega, double osegb)
{
  if ((osega > 0 && osegb > 0) || (osega < 0 && osegb < 0)) 
    return false;
  double alpha2 = orient3d(sega, segb, trip, triq);
  double alpha3 = orient3d(sega, segb, triq, trir);
  if (!same_sign(alpha2, alpha3)) 
    return false;
  double alpha4 = orient3d(sega, segb, trir, trip);
  if (!same_sign(alpha2, alpha4)) 
    return false;
  if (!same_sign(alpha3, alpha4)) 
    return false;

  double sum2 =alpha2+alpha3+alpha4;

  if (osega != osegb && sum2 != 0.0) 
    return true;

  // degenerate: segment and triangle in same plane
  double a[4];
  if (simplex_intersection3d(1, segb, trip, triq, trir, a, a+1, a+2, a+3)) 
    return true;
  if (simplex_intersection3d(1, sega, trip, triq, trir, a, a+1, a+2, a+3)) 
    return true;
  if (simplex_intersection3d(2, sega, segb, triq, trir, a, a+1, a+2, a+3)) 
    return true;
  if (simplex_intersection3d(2, sega, segb, trip, trir, a, a+1, a+2, a+3)) 
    return true;
  if (simplex_intersection3d(2, sega, segb, trip, triq, a, a+1, a+2, a+3)) 
    return true;
  return false;
}

bool pointInTet(const double point[3], const double teta[3], const double tetb[3], const double tetc[3], const double tetd[3])
{
  double alpha[5];
  return simplex_intersection3d(1, point, teta, tetb, tetc, tetd, alpha, alpha+1, alpha+2, alpha+3, alpha+4) != 0;
}

bool intersectSegTri(const double sega[3], const double segb[3], const double tria[3], const double trib[3], const double tric[3])
{
  double osega = orient3d(tria, trib, tric, sega);
  double osegb = orient3d(tria, trib, tric, segb);
  return intersectSegTriWithOrientOnSeg(sega, segb, tria, trib, tric, osega, osegb);
}

bool intersectSegTri(const double sega[3], const double segb[3], const double tria[3], const double trib[3], const double tric[3],
    double segWeight[2], double triangleWeight[3])
{
  return simplex_intersection3d(2, sega, segb, tria, trib, tric, segWeight, segWeight+1, triangleWeight, triangleWeight+1, triangleWeight+2) != 0;
}

// check each edge against the other triangle
bool intersectTriTri(const double pa[3], const double pb[3], const double pc[3], const double qx[3], const double qy[3], const double qz[3])
{
  double oqa = orient3d(pa, pb, pc, qx);
  double oqb = orient3d(pa, pb, pc, qy);
  double oqc = orient3d(pa, pb, pc, qz);

  // if qa, qb, qc are all on the same side of triangle p, then disjoint
  if ((oqa > 0 && oqb > 0 && oqc > 0) || (oqa < 0 && oqb < 0 && oqc < 0)) 
    return false;

  double opa = orient3d(qx, qy, qz, pa);
  double opb = orient3d(qx, qy, qz, pb);
  double opc = orient3d(qx, qy, qz, pc);

  // if pa, pb, pc are all on the same side of triangle q, then disjoint
  if ((opa > 0 && opb > 0 && opc > 0) || (opa < 0 && opb < 0 && opc < 0)) return false;

  return intersectSegTriWithOrientOnSeg(qx, qy, pa, pb, pc, oqa, oqb) ||
         intersectSegTriWithOrientOnSeg(qy, qz, pa, pb, pc, oqb, oqc) ||
         intersectSegTriWithOrientOnSeg(qz, qx, pa, pb, pc, oqc, oqa) ||
         intersectSegTriWithOrientOnSeg(pa, pb, qx, qy, qz, opa, opb) ||
         intersectSegTriWithOrientOnSeg(pb, pc, qx, qy, qz, opb, opc) ||
         intersectSegTriWithOrientOnSeg(pc, pa, qx, qy, qz, opc, opa);

  ////////
//  double otherOrient[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//  bool otherOrientComputed[9];
//  memset(otherOrientComputed, 0, sizeof(otherOrientComputed));
//
//  auto intersectSegTriWithOrient = [&](const double sega[3], const double segb[3], const double trip[3], const double triq[3], const double trir[3],
//      double osega, double osegb, int oabpqID, int oabqrID, int oabrpID)->bool {
//    if ((osega > 0 && osegb > 0) || (osega < 0 && osegb < 0)) return false;
//
//    auto computeOrient = [&otherOrientComputed, &otherOrient](const double a[3], const double b[3], const double c[3], const double d[3], int ID) -> double {
//      double alpha = 0.0;
//      if (otherOrientComputed[ID] == false) {
//        alpha = orient3d(a, b, c, d);
//        otherOrient[ID] = alpha;
//        otherOrientComputed[ID] = true;
//      }
//      else { alpha = otherOrient[ID]; }
//      return alpha;
//    };
//    double alpha2 = computeOrient(sega, segb, trip, triq, oabpqID);
//    double alpha3 = computeOrient(sega, segb, triq, trir, oabqrID);
//
//    if(!same_sign(alpha2, alpha3)) return 0;
//
//    double alpha4 = computeOrient(sega, segb, trir, trip, oabrpID);
//
//    if(!same_sign(alpha2, alpha4)) return 0;
//    if(!same_sign(alpha3, alpha4)) return 0;
//
//    double sum2=alpha2+alpha3+alpha4;
//
//    if(osega != osegb && sum2) { return true; }
//
//    // degenerate: segment and triangle in same plane
//    double a[4];
//    if(simplex_intersection3d(1, segb, trip, triq, trir,
//        a, a+1, a+2, a+3)) { return true; }
//    if(simplex_intersection3d(1, sega, trip, triq, trir,
//        a, a+1, a+2, a+3)) { return true; }
//    if(simplex_intersection3d(2, sega, segb, triq, trir,
//        a, a+1, a+2, a+3)) { return true; }
//    if(simplex_intersection3d(2, sega, segb, trip, trir,
//        a, a+1, a+2, a+3)){ return true; }
//    if(simplex_intersection3d(2, sega, segb, trip, triq,
//        a, a+1, a+2, a+3)){ return true; }
//    return false;
//  };
//
//  // order of otherOrient:
////           0     1     2     3     4     5     6     7     8
////  double xyab, xybc, xyca, yzab, yzbc, yzca, zxab, zxbc, zxca;
//  //         0     3     6     1     4     7     2     5     8
//  //       abxy, abyz, abzx, bcxy, bcyz, bczx, caxy, cayz, cazx
//
//  return intersectSegTriWithOrient(qx, qy, pa, pb, pc, oqa, oqb, 0, 1, 2) ||
//         intersectSegTriWithOrient(qy, qz, pa, pb, pc, oqb, oqc, 3, 4, 5) ||
//         intersectSegTriWithOrient(qz, qx, pa, pb, pc, oqc, oqa, 6 ,7, 8) ||
//         intersectSegTriWithOrient(pa, pb, qx, qy, qz, opa, opb, 0, 3, 6) ||
//         intersectSegTriWithOrient(pb, pc, qx, qy, qz, opb, opc, 1, 4, 7) ||
//         intersectSegTriWithOrient(pc, pa, qx, qy, qz, opc, opa, 2, 5, 8);
}

bool intersectTriTet(const double tria[3], const double trib[3], const double tric[3],
    const double teta[3], const double tetb[3], const double tetc[3], const double tetd[3])
{
  return pointInTet(tria, teta, tetb, tetc, tetd) ||
         // triangle edges against tet faces
         intersectSegTri(tria, trib, teta, tetb, tetc) || intersectSegTri(tria, trib, teta, tetb, tetd) ||
         intersectSegTri(tria, trib, teta, tetc, tetd) || intersectSegTri(tria, trib, tetb, tetc, tetd) ||
         intersectSegTri(tria, tric, teta, tetb, tetc) || intersectSegTri(tria, tric, teta, tetb, tetd) ||
         intersectSegTri(tria, tric, teta, tetc, tetd) || intersectSegTri(tria, tric, tetb, tetc, tetd) ||
         intersectSegTri(trib, tric, teta, tetb, tetc) || intersectSegTri(trib, tric, teta, tetb, tetd) ||
         intersectSegTri(trib, tric, teta, tetc, tetd) || intersectSegTri(trib, tric, tetb, tetc, tetd) ||
         // tet edges against triangle
         intersectSegTri(teta, tetb, tria, trib, tric) || intersectSegTri(teta, tetc, tria, trib, tric) ||
         intersectSegTri(teta, tetd, tria, trib, tric) || intersectSegTri(tetb, tetc, tria, trib, tric) ||
         intersectSegTri(tetb, tetd, tria, trib, tric) || intersectSegTri(tetc, tetd, tria, trib, tric);
}

// d: [0,3), direction (X,Y,Z) which this face is perpendicular to
bool intersectSegAABBFace(int d, const double sa[3], const double sb[3],
    const double ra[3], const double rb[3], const double rc[3], const double rd[3])
{
  double x = ra[d];
  if (sa[d] < x && sb[d] < x) return 0;
  if (sa[d] > x && sb[d] > x) return 0;
  double a2 = orient3d(sa, sb, rc, rd);
  double a3 = orient3d(sa, sb, rd, ra);
  if(!same_sign(a2, a3)) return 0;
  double a4 = orient3d(sa, sb, ra, rb);
  if(!same_sign(a2, a4)) return 0;
  if(!same_sign(a3, a4)) return 0;
  double a5 = orient3d(sa, sb, rb, rc);
  if(!same_sign(a2, a5)) return 0;
  if(!same_sign(a3, a5)) return 0;
  if(!same_sign(a4, a5)) return 0;

  double sum1 = a2 + a3 + a4 + a5;
  if ((sa[d] != x || sb[d] != x) && sum1 != 0) // non-degenerate
    return 1;
  // process degeneracy

  if(simplex_intersection3d(1, sa, ra, rb, rc, &a2, &a3, &a4, &a5)) // whether sa inside triangle ra,rb,rc
  {
     return 1;
  }
  if(simplex_intersection3d(1, sa, ra, rc, rd, &a2, &a3, &a4, &a5)) // whether sa inside triangle ra,rc,rd
  {
     return 1;
  }
  if(simplex_intersection3d(1, sb, ra, rb, rc, &a2, &a3, &a4, &a5)) // whether sb inside triangle ra,rb,rc
  {
     return 1;
  }
  if(simplex_intersection3d(1, sb, ra, rc, rd, &a2, &a3, &a4, &a5)) // whether sb inside triangle ra,rc,rd
  {
     return 1;
  }
  if(simplex_intersection3d(2, sa, sb, ra, rb, &a2, &a3, &a4, &a5)) // whether segment sa, sb intersect segment ra, rb
  {
     return 1;
  }
  if(simplex_intersection3d(2, sa, sb, rb, rc, &a2, &a3, &a4, &a5)) // whether segment sa, sb intersect segment rb, rc
  {
     return 1;
  }
  if(simplex_intersection3d(2, sa, sb, rc, rd, &a2, &a3, &a4, &a5)) // whether segment sa, sb intersect segment rc, rd
  {
     return 1;
  }
  if(simplex_intersection3d(2, sa, sb, rd, ra, &a2, &a3, &a4, &a5)) // whether segment sa, sb intersect segment rd, ra
  {
     return 1;
  }
  return 0;
}

// d: [0,3), direction (X,Y,Z) which this face is perpendicular to
bool intersectTriAABBFace(int d, const double ta[3], const double tb[3], const double tc[3], // t for triangle
    const double ra[3], const double rb[3], const double rc[3], const double rd[3])          // r for rectangle
{
//  return intersectSegTri(pa, pb, qa, qb, qc) || intersectSegTri(pa, pc, qa, qb, qc) || intersectSegTri(pb, pc, qa, qb, qc) ||
//      intersectSegTri(qa, qb, pa, pb, pc) || intersectSegTri(qa, qc, pa, pb, pc) || intersectSegTri(qb, qc, pa, pb, pc);
  return intersectSegAABBFace(d, ta, tb, ra, rb, rc, rd) ||
         intersectSegAABBFace(d, ta, tc, ra, rb, rc, rd) ||
         intersectSegAABBFace(d, tb, tc, ra, rb, rc, rd) ||
         intersectSegTri(ra, rb, ta, tb, tc) ||
         intersectSegTri(rb, rc, ta, tb, tc) ||
         intersectSegTri(rc, rd, ta, tb, tc) ||
         intersectSegTri(rd, ra, ta, tb, tc);
}

// local variable const double p[8][3] in the following functions represent the AABB vertic4es
// order of p is:
//     3 - - - 2
//    /|      /|
//   7 - - - 6 |       y
//   | |     | |       |
//   | 0 - - | 1       |_ _ _x
//   |/      |/       /
//   4 - - - 5       z

bool intersectTriAABB(const double ta[3], const double tb[3], const double tc[3], const double bmin[3], const double bmax[3])
{
  assert(bmin[0] <= bmax[0] && bmin[1] <= bmax[1] && bmin[2] <= bmax[2]);
  if (bmin[0] <= ta[0] && ta[0] <= bmax[0] && bmin[1] <= ta[1] && ta[1] <= bmax[1] && bmin[2] <= ta[2] && ta[2] <= bmax[2]) return true;
  if (bmin[0] <= tb[0] && tb[0] <= bmax[0] && bmin[1] <= tb[1] && tb[1] <= bmax[1] && bmin[2] <= tb[2] && tb[2] <= bmax[2]) return true;
  if (bmin[0] <= tc[0] && tc[0] <= bmax[0] && bmin[1] <= tc[1] && tc[1] <= bmax[1] && bmin[2] <= tc[2] && tc[2] <= bmax[2]) return true;

  double minx = min(min(ta[0], tb[0]), tc[0]);
  double miny = min(min(ta[1], tb[1]), tc[1]);
  double minz = min(min(ta[2], tb[2]), tc[2]);
  double maxx = max(max(ta[0], tb[0]), tc[0]);
  double maxy = max(max(ta[1], tb[1]), tc[1]);
  double maxz = max(max(ta[2], tb[2]), tc[2]);

  if ( ( bmax[0] < minx ) || ( maxx < bmin[0] ) || ( bmax[1] < miny ) || ( maxy < bmin[1] ) || ( bmax[2] < minz ) || ( maxz < bmin[2] ) )
    return false; //no overlap

  // create a box mesh for AABB
  const double p[8][3] = { {bmin[0], bmin[1], bmin[2]}, { bmax[0], bmin[1], bmin[2] }, { bmax[0], bmax[1], bmin[2] }, { bmin[0], bmax[1], bmin[2] },
                            { bmin[0], bmin[1], bmax[2] }, { bmax[0], bmin[1], bmax[2] }, {bmax[0], bmax[1], bmax[2]}, { bmin[0], bmax[1], bmax[2] } };

  if (intersectSegAABBFace(0, ta, tb, p[0], p[4], p[7], p[3]) ||
      intersectSegAABBFace(0, ta, tc, p[0], p[4], p[7], p[3]) ||
      intersectSegAABBFace(0, tb, tc, p[0], p[4], p[7], p[3])) return true;
  if (intersectSegAABBFace(0, ta, tb, p[1], p[2], p[6], p[5]) ||
      intersectSegAABBFace(0, ta, tc, p[1], p[2], p[6], p[5]) ||
      intersectSegAABBFace(0, tb, tc, p[1], p[2], p[6], p[5])) return true;

  if (intersectSegAABBFace(1, ta, tb, p[0], p[1], p[5], p[4]) ||
      intersectSegAABBFace(1, ta, tc, p[0], p[1], p[5], p[4]) ||
      intersectSegAABBFace(1, tb, tc, p[0], p[1], p[5], p[4])) return true;
  if (intersectSegAABBFace(1, ta, tb, p[2], p[3], p[7], p[6]) ||
      intersectSegAABBFace(1, ta, tc, p[2], p[3], p[7], p[6]) ||
      intersectSegAABBFace(1, tb, tc, p[2], p[3], p[7], p[6])) return true;

  if (intersectSegAABBFace(2, ta, tb, p[0], p[3], p[2], p[1]) ||
      intersectSegAABBFace(2, ta, tc, p[0], p[3], p[2], p[1]) ||
      intersectSegAABBFace(2, tb, tc, p[0], p[3], p[2], p[1])) return true;
  if (intersectSegAABBFace(2, ta, tb, p[4], p[5], p[6], p[7]) ||
      intersectSegAABBFace(2, ta, tc, p[4], p[5], p[6], p[7]) ||
      intersectSegAABBFace(2, tb, tc, p[4], p[5], p[6], p[7])) return true;

  if (intersectSegTri(p[0], p[4], ta, tb, tc) || intersectSegTri(p[4], p[7], ta, tb, tc) ||
      intersectSegTri(p[7], p[3], ta, tb, tc) || intersectSegTri(p[3], p[0], ta, tb, tc) ||
      intersectSegTri(p[1], p[2], ta, tb, tc) || intersectSegTri(p[2], p[6], ta, tb, tc) ||
      intersectSegTri(p[6], p[5], ta, tb, tc) || intersectSegTri(p[5], p[1], ta, tb, tc) ||
      intersectSegTri(p[0], p[1], ta, tb, tc) || intersectSegTri(p[4], p[5], ta, tb, tc) ||
      intersectSegTri(p[6], p[7], ta, tb, tc) || intersectSegTri(p[2], p[3], ta, tb, tc)) return true;
//
//  if (intersectTriAABBFace(0, ta, tb, tc, p[0], p[4], p[7], p[3])) return true; // check intersection with face 0-4-7-3
//  if (intersectTriAABBFace(0, ta, tb, tc, p[1], p[2], p[6], p[5])) return true;
//
//  if (intersectTriAABBFace(1, ta, tb, tc, p[0], p[1], p[5], p[4])) return true;
//  if (intersectTriAABBFace(1, ta, tb, tc, p[2], p[3], p[7], p[6])) return true;
//
//  if (intersectTriAABBFace(2, ta, tb, tc, p[0], p[3], p[2], p[1])) return true;
//  if (intersectTriAABBFace(2, ta, tb, tc, p[4], p[5], p[6], p[7])) return true;
  return false;
}


bool intersectSegAABB(const double sa[3], const double sb[3], const double bmin[3], const double bmax[3])
{
  assert(bmin[0] <= bmax[0] && bmin[1] <= bmax[1] && bmin[2] <= bmax[2]);
  auto getOutcode = [](const double pos[3], const double bmin[3], const double bmax[3])
  {
    char code = 0;
    code += (pos[0] < bmin[0]);
    code += (pos[0] > bmax[0]) << 1;
    code += (pos[1] < bmin[1]) << 2;
    code += (pos[1] > bmax[1]) << 3;
    code += (pos[2] < bmin[2]) << 4;
    code += (pos[2] > bmax[2]) << 5;
    return code;
  };
  char codea = getOutcode(sa, bmin, bmax), codeb = getOutcode(sb, bmin, bmax);
  if (codea == 0 || codeb == 0) return true;
  if ((codea & codeb) != 0) return false;

  const double p[8][3] = { {bmin[0], bmin[1], bmin[2]}, { bmax[0], bmin[1], bmin[2] }, { bmax[0], bmax[1], bmin[2] }, { bmin[0], bmax[1], bmin[2] },
      { bmin[0], bmin[1], bmax[2] }, { bmax[0], bmin[1], bmax[2] }, {bmax[0], bmax[1], bmax[2]}, { bmin[0], bmax[1], bmax[2] } };
//  const int t[12][3] = { {0,3,2}, {0,2,1}, {4,5,6}, {4,6,7}, {0,1,5}, {0,5,4},
//      {3,7,6}, {3,6,2}, {1,2,6}, {1,6,5}, {0,4,7}, {0,7,3} };
//  for(int i = 0; i < 12; i++) {
//    if (intersectSegTri(sa, sb, p[t[i][0]], p[t[i][1]], p[t[i][2]])) return true;
//  }
//  return false;

  if (intersectSegAABBFace(0, sa, sb, p[0], p[4], p[7], p[3])) return true; // check intersection with face 0-4-7-3
  if (intersectSegAABBFace(0, sa, sb, p[1], p[2], p[6], p[5])) return true;
  if (intersectSegAABBFace(1, sa, sb, p[0], p[1], p[5], p[4])) return true;
  if (intersectSegAABBFace(1, sa, sb, p[2], p[3], p[7], p[6])) return true;
  if (intersectSegAABBFace(2, sa, sb, p[0], p[3], p[2], p[1])) return true;
  if (intersectSegAABBFace(2, sa, sb, p[4], p[5], p[6], p[7])) return true;
  return false;
}

namespace
{

void getXY(const double v[3], double o[2]) { o[0] = v[0]; o[1] = v[1]; }
void getYZ(const double v[3], double o[2]) { o[0] = v[1]; o[1] = v[2]; }
void getZX(const double v[3], double o[2]) { o[0] = v[2]; o[1] = v[0]; }
std::function<void(const double v[3], double o[2])> getComp[3] = { getXY, getYZ, getZX };

}

bool isTriangleDegenerate(const double ta[3], const double tb[3], const double tc[3])
{
  for(int i = 0; i < 3; i++)
  {
    auto f = getComp[i];
    double a[2], b[2], c[2];
    f(ta, a);
    f(tb, b);
    f(tc, c);

    if (orient2d(a,b,c) != 0.0) return false;
  }
  return true;
}

int inCircumsphereOnPlane(const double ta[3], const double tb[3], const double tc[3], const double td[3])
{
  for(int i = 0; i < 3; i++)
  {
    auto f = getComp[i];
    double a[2], b[2], c[2], d[2];
    f(ta, a);
    f(tb, b);
    f(tc, c);
    f(td, d);
    double ret = incircle(a, b, c, d);
    if (ret < 0) return -1;
    if (ret > 0) return 1;
  }
  return 0;
}

bool intersectSegSeg2d(const double sa[2], const double sb[2], const double ta[2], const double tb[2])
{
  double alpha[4];
  return simplex_intersection2d(2, sa, sb, ta, tb, alpha, alpha+1, alpha+2, alpha+3);
}

bool intersectSegSeg2d(const double sa[2], const double sb[2], const double ta[2], const double tb[2], double sw[2], double tw[2])
{
  return simplex_intersection2d(2, sa, sb, ta, tb, sw, sw+1, tw, tw+1);
}

