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

#ifndef VEC3ER_H
#define VEC3ER_H

#include "er.h"

// vector of 3 ES
class Vec3ER
{
public:
  inline Vec3ER() {}
  inline Vec3ER(ER x, ER y, ER z) : elt { x, y, z } {}
  inline Vec3ER(double x, double y, double z) : elt { x, y, z } {}
  // create a vector with all entries "entry" (can create zero vector for entry=0.0)
  inline explicit Vec3ER(ER entry) : elt { entry, entry, entry } {}
  inline explicit Vec3ER(double entry) : elt { entry, entry, entry } {}
  // create a vector from the array of three Ss pointed to by "vec"
  inline Vec3ER(const ER vec[3]) : elt { vec[0], vec[1], vec[2] } {}
  inline Vec3ER(const double vec[3]) : elt { vec[0], vec[1], vec[2] } {}
  inline Vec3ER(const Vec3ER & vec);

  inline void set(ER x0, ER x1, ER x2); // assign vector [x0, x1, x2]
  inline void set(ER value); // set all elements to value

  inline Vec3ER & operator=(const Vec3ER & source);
  inline bool operator==(const Vec3ER & vec2) const;
  inline bool operator!=(const Vec3ER & vec2) const;

  inline const Vec3ER operator+ (const Vec3ER & vec2) const;
  inline Vec3ER & operator+= (const Vec3ER & vec2);

  inline const Vec3ER operator- (const Vec3ER & vec2) const;
  inline Vec3ER & operator-= (const Vec3ER & vec2);

  inline const Vec3ER operator* (ER scalar) const;
  inline Vec3ER & operator*= (ER scalar);

  inline Vec3ER operator/ (ER scalar) const;
  inline Vec3ER & operator/= (ER scalar);

  // operator for V3 to be used as a key in std::set, std::map, etc.
  inline bool operator < (const Vec3ER & vec2) const;

  friend inline Vec3ER operator* (ER scalar, const Vec3ER & vec2);
  friend inline Vec3ER operator/ (ER scalar, const Vec3ER & vec2);
  friend inline Vec3ER operator- (const Vec3ER & vec1);

  friend inline ER dot(const Vec3ER & vec1, const Vec3ER & vec2); // dot product

  friend inline Vec3ER cross(const Vec3ER & vec1, const Vec3ER & vec2); // cross product

  friend inline std::ostream & operator << (std::ostream &s, const Vec3ER &v);
  void print() const;

  inline ER & operator[] (int index); // v[i] returns i-th entry of v
  inline const ER & operator[] (int index) const;

  // copies the vector into an array of length 3
  inline void convertToArray(ER vecArray[3]) const;
  // adds the vector into an array of length 3
  inline void addToArray(ER vecArray[3]) const;

  template <class CGAL_Kernel_Vector_3>
  friend Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal);

protected:
  ER elt[3];
#ifndef VEGA_USE_CGAL_HEADER
  inline Vec3ER(ER::UninitializedLabel);
#endif
};

extern const Vec3ER VEC3ER_NULL;

template <class CGAL_Kernel_Vector_3>
Vec3ER assignCGALToVec3ER(const CGAL_Kernel_Vector_3 & cgal);

// triangle-triangle intersection
bool intersectTriTri(const Vec3ER & pa, const Vec3ER & pb, const Vec3ER & pc, const Vec3ER & qa, const Vec3ER & qb, const Vec3ER & qc);
bool intersectTriAABB(const Vec3ER & tria, const Vec3ER trib, const Vec3ER tric, const Vec3ER bmin, const Vec3ER bmax);

// === below is the implementation ===

inline Vec3ER::Vec3ER(const Vec3ER & vec)
{
  elt[0] = vec.elt[0];
  elt[1] = vec.elt[1];
  elt[2] = vec.elt[2];
}

inline Vec3ER & Vec3ER::operator=(const Vec3ER & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];

  return *this;
}

inline bool Vec3ER::operator==(const Vec3ER & vec2) const
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]) &&
          (elt[2] == vec2[2]));
}

inline bool Vec3ER::operator!=(const Vec3ER & vec2) const
{
  return ((elt[0] != vec2[0]) ||
          (elt[1] != vec2[1]) ||
          (elt[2] != vec2[2]));
}

inline bool Vec3ER::operator<(const Vec3ER & vec2) const
{
  if(elt[0] < vec2[0])
    return true;
  if(elt[0] > vec2[0])
    return false;
  if(elt[1] < vec2[1])
    return true;
  if(elt[1] > vec2[1])
    return false;
  return elt[2] < vec2[2];
}

inline Vec3ER operator* (ER scalar, const Vec3ER & vec2)
{
  Vec3ER result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;

  return result;
}

inline Vec3ER operator/ (ER scalar, const Vec3ER & vec2)
{
  Vec3ER result = vec2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;
  result.elt[2] /= scalar;

  return result;
}

inline Vec3ER operator- (const Vec3ER & vec1)
{
  Vec3ER result = vec1;
  result.elt[0] *= -1.0;
  result.elt[1] *= -1.0;
  result.elt[2] *= -1.0;

  return result;
}

inline const Vec3ER Vec3ER::operator+ (const Vec3ER & vec2) const
{
  Vec3ER sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];
  sum.elt[2] += vec2.elt[2];

  return sum;
}

inline Vec3ER & Vec3ER::operator+= (const Vec3ER & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];
  elt[2] += vec2.elt[2];

  return *this;
}

inline const Vec3ER Vec3ER::operator- (const Vec3ER & vec2) const
{
  Vec3ER sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];
  sum.elt[2] -= vec2.elt[2];

  return sum;
}

inline Vec3ER & Vec3ER::operator-= (const Vec3ER & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];
  elt[2] -= vec2.elt[2];

  return *this;
}

inline ER & Vec3ER::operator[] (int index)
{
  return elt[index];
}

inline const ER & Vec3ER::operator[] (int index) const
{
  return elt[index];
}

inline Vec3ER & Vec3ER::operator*= (ER scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;
  return *this;
}

inline const Vec3ER Vec3ER::operator* (ER scalar) const
{
  return (Vec3ER(elt[0] * scalar, elt[1] * scalar, elt[2] * scalar));
}

inline Vec3ER Vec3ER::operator/ (ER scalar) const
{
  return (Vec3ER(elt[0] / scalar, elt[1] / scalar, elt[2] / scalar));
}

inline Vec3ER & Vec3ER::operator/= (ER scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;
  return *this;
}

inline std::ostream &operator << (std::ostream & s, const Vec3ER & v)
{
  return(s << '[' << ER_toDouble(v[0]) << ' ' << ER_toDouble(v[1]) << ' ' << ER_toDouble(v[2]) << ']');
}

inline void Vec3ER::convertToArray(ER vecArray[3]) const
{
  vecArray[0] = elt[0];
  vecArray[1] = elt[1];
  vecArray[2] = elt[2];
}

inline void Vec3ER::addToArray(ER vecArray[3]) const
{
  vecArray[0] += elt[0];
  vecArray[1] += elt[1];
  vecArray[2] += elt[2];
}

inline void Vec3ER::print() const
{
  printf("[%G %G %G]\n", ER_toDouble(elt[0]), ER_toDouble(elt[1]), ER_toDouble(elt[2]));
}

inline void Vec3ER::set(ER x0, ER x1, ER x2) // assign vector [x0, x1, x2]
{
  elt[0] = x0;
  elt[1] = x1;
  elt[2] = x2;
}

inline void Vec3ER::set(ER value) // set all elements to value
{
  elt[0] = value;
  elt[1] = value;
  elt[2] = value;
}

inline ER dot(const Vec3ER & vec1, const Vec3ER & vec2)
{
  return (vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]);
}

inline Vec3ER cross(const Vec3ER & vec1, const Vec3ER & vec2)
{
  Vec3ER result(vec1[1] * vec2[2] - vec2[1] * vec1[2],
              -vec1[0] * vec2[2] + vec2[0] * vec1[2],
               vec1[0] * vec2[1] - vec2[0] * vec1[1]);

  return result;
}

inline ER len2(const Vec3ER & vec1)
{
  return(dot(vec1,vec1));
}

#ifndef VEGA_USE_CGAL_HEADER
inline Vec3ER::Vec3ER(ER::UninitializedLabel) : elt {ER::UninitializedLabel(), ER::UninitializedLabel(), ER::UninitializedLabel()}
{
}
#endif



#endif
