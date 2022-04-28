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

#include "er.h"

#ifndef VEGA_USE_CGAL_HEADER

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/intersections.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;

namespace
{
//  const K::RT & RT(const void * p) {
//    return *(const K::RT*)p;
//  }
  K::RT & RT(void * p) {
    return *(K::RT*)p;
  }
}

ER::ER(double s)
{
  p = new K::RT(s);
}

ER::ER(const ER & s)
{
  p = new K::RT(RT(s.p));
}

ER::ER(ER && s)
{
  p = s.p;
  s.p = new K::RT(0.0);
}

ER::~ER() { delete (K::RT*)p; }

ER & ER::operator = (const ER & s)
{
  RT(p) = RT(s.p);
  return *this;
}

ER & ER::operator = (ER && s)
{
  p = s.p;
  s.p = new K::RT(0.0);
  return *this;
}

ER ER::operator + (const ER & s) const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(RT(p) + RT(s.p));
  return ret;
}

ER ER::operator - (const ER & s) const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(RT(p) - RT(s.p));
  return ret;
}

ER ER::operator * (const ER & s) const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(RT(p) * RT(s.p));
  return ret;
}

ER ER::operator / (const ER & s) const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(RT(p) / RT(s.p));
  return ret;
}

ER ER::operator - () const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(-RT(p));
  return ret;
}

bool ER::operator == (const ER & s) const
{
  return RT(p) == RT(s.p);
}

bool ER::operator != (const ER & s) const
{
  return RT(p) != RT(s.p);
}

bool ER::operator <= (const ER & s) const
{
  return RT(p) <= RT(s.p);
}

bool ER::operator < (const ER & s) const
{
  return RT(p) < RT(s.p);
}

bool ER::operator >= (const ER & s) const
{
  return RT(p) >= RT(s.p);
}

bool ER::operator > (const ER & s) const
{
  return RT(p) > RT(s.p);
}

ER & ER::operator += (const ER & s)
{
  RT(p) += RT(s.p);
  return *this;
}

ER & ER::operator -= (const ER & s)
{
  RT(p) -= RT(s.p);
  return *this;
}

ER & ER::operator *= (const ER & s)
{
  RT(p) *= RT(s.p);
  return *this;
}

ER & ER::operator /= (const ER & s)
{
  RT(p) /= RT(s.p);
  return *this;
}


ER ER::abs() const
{
  UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(CGAL::abs(RT(p)));
  return ret;
}

char ER::sign() const
{
  auto s = CGAL::sign(RT(p));
  static_assert((char)CGAL::POSITIVE == 1, "Error: CGAL::POSITIVE is not 1");
  static_assert((char)CGAL::NEGATIVE == -1, "Error: CGAL::NEGATIVE is not -1");
  static_assert((char)CGAL::ZERO == 0, "Error: CGAL::ZERO is not 0");
  return (char)s;
}

double ER::toDouble() const
{
  // following code from libigl
  // FORCE evaluation of the exact type otherwise interval might be huge.
  const auto cgal = RT(p).exact();
  const auto interval = CGAL::to_interval(cgal);
  double d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
  return d;
}

template <class CGAL_Kernel_RT>
ER assignCGALToER(const CGAL_Kernel_RT & cgal)
{
  ER::UninitializedLabel l;
  ER ret(l);
  ret.p = new K::RT(cgal);
  return ret;
}

template ER assignCGALToER<K::RT>(const K::RT & cgal);

#endif


