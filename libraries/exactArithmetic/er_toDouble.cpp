// This file is modified from libigl, a simple c++ geometry processing library.
// libigl is under the Mozilla Public License v. 2.0.
// According to the license restriction, we release this code including the changes
// we made to libigl code under the same license.
// 
// Copyright of our changes: (C) 2018 USC
// Code authors of our changes: Yijing Li, Jernej Barbic 
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.

#include "er.h"

#ifdef VEGA_USE_CGAL_HEADER

double ER_toDouble(const ER & er)
{
  // following code from libigl
  // FORCE evaluation of the exact type otherwise interval might be huge.
  const auto cgal = er.exact();
  const auto interval = CGAL::to_interval(cgal);
  double d = interval.first;
  do
  {
    const double next = nextafter(d, interval.second);
    if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next))
      break;
    d = next;
  } while (d < interval.second);
  return d;
}

#endif
