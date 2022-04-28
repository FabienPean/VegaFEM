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

#ifndef CGALKERNELEXACT_H_
#define CGALKERNELEXACT_H_

#include <iostream>

#define VEGA_USE_CGAL_HEADER
//#define VEGA_ER_INHERIT_FROM_CGAL

// A exact arithmetic real number
// implemented using CGAL's exact kernel without sqrt
// the benefit is that the code that relies on it won't be affected by the sloooooow CGAL header compilation

#ifdef VEGA_USE_CGAL_HEADER
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>


#ifdef VEGA_ER_INHERIT_FROM_CGAL

class Vec3ER;
class ER : public CGAL::Exact_predicates_exact_constructions_kernel::RT
{
public:
  typedef CGAL::Exact_predicates_exact_constructions_kernel K;
  typedef K::RT S;

  ER(double s = 0.0) : S(s) {}
  ER(const ER & o) : S(o) {}
  ER(ER && o) : S(move(o)) {}
  ER(const S & s) : S(s) {}
  ~ER() {}

  ER & operator = (const ER & o) { S::operator =(o); return *this; }
  ER & operator = (ER && o) { S::operator =(move(o)); return *this; }

  // arithmetic operations
  ER operator + (const ER & o) const { return *(const S *)this + *(const S *)&o; }
  ER operator - (const ER & o) const { return *(const S *)this - *(const S *)&o; }
  ER operator * (const ER & o) const { return *(const S *)this * *(const S *)&o; }
  ER operator / (const ER & o) const { return *(const S *)this / *(const S *)&o; }
  ER operator - () const { return -*(const S *)this; }

  ER & operator += (const ER & o) { S::operator+=(o); return *this; }
  ER & operator -= (const ER & o) { S::operator-=(o); return *this; }
  ER & operator *= (const ER & o) { S::operator*=(o); return *this; }
  ER & operator /= (const ER & o) { S::operator/=(o); return *this; }

  // absolute
  ER abs() const { return CGAL::abs(*(const S *)this); }
  // return {-1, 0, +1}
  char sign() const { return char(CGAL::sign(*(const S*)this)); }

  friend double ER_toDouble(const ER & er);
  double toDouble() const { return ER_toDouble(*this); }

protected:
  friend class Vec3ER;
  template <class CGAL_Kernel_Vector_3>
  friend Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal);
  friend bool intersectTriTri(const Vec3ER & pa, const Vec3ER & pb, const Vec3ER & pc, const Vec3ER & qa, const Vec3ER & qb, const Vec3ER & qc);
};

inline ER ER_abs(const ER & er)
{
  return er.abs();
}

inline char ER_sign(const ER & er)
{
  return er.sign();
}


#else // VEGA_ER_INHERIT_FROM_CGAL

typedef CGAL::Exact_predicates_exact_constructions_kernel::RT ER;
inline ER ER_abs(const ER & er)
{
  return CGAL::abs(er);
}

inline char ER_sign(const ER & er)
{
  return char(CGAL::sign(er));;
}

#endif

double ER_toDouble(const ER & er);

template <class CGAL_Kernel_RT>
ER assignCGALToER(const CGAL_Kernel_RT & cgal) { return cgal; }


#else
class ER
{
public:
  ER(double s = 0.0);
  ER(const ER &);
  ER(ER && );
  ~ER();

  ER & operator = (const ER &);
  ER & operator = (ER &&);

  // arithmetic operations
  ER operator + (const ER &) const;
  ER operator - (const ER &) const;
  ER operator * (const ER &) const;
  ER operator / (const ER &) const;
  ER operator - () const;

  // comparison
  bool operator == (const ER &) const;
  bool operator != (const ER &) const;
  bool operator <= (const ER &) const;
  bool operator <  (const ER &) const;
  bool operator >= (const ER &) const;
  bool operator >  (const ER &) const;

  // arithmetic updates
  ER & operator += (const ER &);
  ER & operator -= (const ER &);
  ER & operator *= (const ER &);
  ER & operator /= (const ER &);

  // absolute
  ER abs() const;
  char sign() const; // return {-1, 0, +1}

  // conversion
  double toDouble() const;

  template <class CGAL_Kernel_RT>
  friend ER assignCGALToER(const CGAL_Kernel_RT & cgal);

  struct UninitializedLabel {}; // used for constructing a S with uninitialized data for performance
  ER(UninitializedLabel) {}

protected:
  void * p = nullptr;

  friend class Vec3ER;
  template <class CGAL_Kernel_Vector_3>
  friend Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal);
  friend bool intersectTriTri(const Vec3ER & pa, const Vec3ER & pb, const Vec3ER & pc, const Vec3ER & qa, const Vec3ER & qb, const Vec3ER & qc);
};



// used to convert CGAL Exact kernel scalar to ES without the need to include CGAL headers
template <class CGAL_Kernel_RT>
ER assignCGALToER(const CGAL_Kernel_RT & cgal);

#endif

#endif
