/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "minivector" library , Copyright (C) 2018 USC                         *
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


#ifndef VEC3I_H
#define VEC3I_H

#include <stdio.h>
#include <math.h>
#include <ostream>
#include "vec3d.h"

class Vec3i
{
public:
  inline Vec3i() {}
  inline Vec3i(int x, int y, int z) {elt[0]=x; elt[1]=y; elt[2]=z;}
  inline Vec3i(int entry); // create a vector with all entries "entry" (can create zero vector for entry=0)
  inline Vec3i(const int vec[3]); // create a vector from the array of three ints pointed to by "vec"
  inline Vec3i(const Vec3i & vec);

  inline void set(int x0, int x1, int x2); // assign vector [x0, x1, x2]
  inline void set(int value); // set all elements to value

  inline Vec3i & operator=(const Vec3i & source);
  inline bool operator==(const Vec3i & vec2) const;
  inline bool operator!=(const Vec3i & vec2) const;

  inline const Vec3i operator+ (const Vec3i & vec2) const;
  inline Vec3i & operator+= (const Vec3i & vec2);

  inline const Vec3i operator- (const Vec3i & vec2) const;
  inline Vec3i & operator-= (const Vec3i & vec2);

  inline const Vec3i operator* (int scalar) const;
  inline Vec3i & operator*= (int scalar);

  inline const Vec3i operator/ (int scalar) const;
  inline Vec3i & operator/= (int scalar);

  inline const Vec3d operator* (double scalar) const;
  inline const Vec3d operator/ (double scalar) const;

  // operator for Vec3i to be used as a key in std::set, std::map, etc.
  inline bool operator < (const Vec3i & vec2) const;

  friend inline Vec3i operator* (int scalar, const Vec3i & vec2);
  friend inline Vec3i operator- (const Vec3i & vec1);

  friend inline int dot(const Vec3i & vec1, const Vec3i & vec2); // dot product

  friend inline Vec3i cross(const Vec3i & vec1, const Vec3i & vec2); // cross product

  friend inline std::ostream &operator << (std::ostream & s, const Vec3i & v);
  void print() const;

  inline int & operator[] (int index); // v[i] returns i-th entry of v
  inline const int & operator[] (int index) const;

  // copy the vector into an array of length 3
  inline void convertToArray(int vecArray[3]) const;
  // add the vector into an array of length 3
  inline void addToArray(int vecArray[3]) const;

  // find the first index in elt which equals to value; return -1 if not found
  inline int getInvertedIndex(int value) const;
 
  // rotate elt[0] to elt[3] so that elt[newStartIndex] is moved to elt[0]
  inline void rotate(int newStartIndex);

  // do set intersection; return whether this and vec2 share at least one element
  inline bool intersect(const Vec3i & vec2) const;

  const int * begin() const { return elt; }
  int * begin() { return elt; }
  const int * end() const { return elt + 3; }
  int * end() { return elt + 3; }

protected:
  int elt[3];
};

// === below is the implementation ===

inline Vec3i::Vec3i(int entry)
{
  elt[0] = entry;
  elt[1] = entry;
  elt[2] = entry;
}

inline Vec3i::Vec3i(const int vec[3])
{
  elt[0] = vec[0];
  elt[1] = vec[1];
  elt[2] = vec[2];
}

inline Vec3i::Vec3i(const Vec3i & vec)
{
  elt[0] = vec.elt[0];
  elt[1] = vec.elt[1];
  elt[2] = vec.elt[2];
}

inline Vec3i & Vec3i::operator=(const Vec3i & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];

  return *this;
}

inline bool Vec3i::operator==(const Vec3i & vec2) const
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]) &&
          (elt[2] == vec2[2]));
}

inline bool Vec3i::operator!=(const Vec3i & vec2) const
{
  return ((elt[0] != vec2[0]) ||
          (elt[1] != vec2[1]) ||
          (elt[2] != vec2[2]));
}

inline bool Vec3i::operator<(const Vec3i & vec2) const
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

inline Vec3i operator* (int scalar, const Vec3i & vec2)
{
  Vec3i result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;

  return result;
}

inline Vec3i operator- (const Vec3i & vec1)
{
  Vec3i result = vec1;
  result.elt[0] *= -1;
  result.elt[1] *= -1;
  result.elt[2] *= -1;

  return result;
}

inline const Vec3i Vec3i::operator+ (const Vec3i & vec2) const
{
  Vec3i sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];
  sum.elt[2] += vec2.elt[2];

  return sum;
}

inline Vec3i & Vec3i::operator+= (const Vec3i & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];
  elt[2] += vec2.elt[2];

  return *this;
}

inline const Vec3i Vec3i::operator- (const Vec3i & vec2) const
{
  Vec3i sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];
  sum.elt[2] -= vec2.elt[2];

  return sum;
}

inline Vec3i & Vec3i::operator-= (const Vec3i & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];
  elt[2] -= vec2.elt[2];

  return *this;
}

inline int & Vec3i::operator[] (int index)
{
  return elt[index];
}

inline const int & Vec3i::operator[] (int index) const
{
  return elt[index];
}

inline int dot(const Vec3i & vec1, const Vec3i & vec2)
{
  return (vec1.elt[0] * vec2.elt[0] + vec1.elt[1] * vec2.elt[1] + vec1.elt[2] * vec2.elt[2]);
}

inline Vec3i cross(const Vec3i & vec1, const Vec3i & vec2)
{
  Vec3i result(vec1.elt[1] * vec2.elt[2] - vec2.elt[1] * vec1.elt[2],
              -vec1.elt[0] * vec2.elt[2] + vec2.elt[0] * vec1.elt[2],
               vec1.elt[0] * vec2.elt[1] - vec2.elt[0] * vec1.elt[1]);

  return result;
}

inline Vec3i & Vec3i::operator*= (int scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;
  return *this;
}

inline const Vec3i Vec3i::operator* (int scalar) const
{
  return (Vec3i(elt[0] * scalar, elt[1] * scalar, elt[2] * scalar));
}

inline const Vec3i Vec3i::operator/ (int scalar) const
{
  return (Vec3i(elt[0] / scalar, elt[1] / scalar, elt[2] / scalar));
}

inline const Vec3d Vec3i::operator* (double scalar) const
{
  return (Vec3d(elt[0] * scalar, elt[1] * scalar, elt[2] * scalar));
}

inline const Vec3d Vec3i::operator/ (double scalar) const
{
  return (Vec3d(elt[0] / scalar, elt[1] / scalar, elt[2] / scalar));
}


inline Vec3i & Vec3i::operator/= (int scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;
  return *this;
}

inline std::ostream &operator << (std::ostream &s, const Vec3i &v)
{
  return(s << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ']');
}

inline void Vec3i::convertToArray(int vecArray[3]) const
{
  vecArray[0] = elt[0];
  vecArray[1] = elt[1];
  vecArray[2] = elt[2];
}

inline void Vec3i::addToArray(int vecArray[3]) const
{
  vecArray[0] += elt[0];
  vecArray[1] += elt[1];
  vecArray[2] += elt[2];
}

inline void Vec3i::print() const
{
  printf("[%d %d %d]\n", elt[0], elt[1], elt[2]);
}

inline void Vec3i::set(int x0, int x1, int x2) // assign vector [x0, x1, x2]
{
  elt[0] = x0;
  elt[1] = x1;
  elt[2] = x2;
}

inline void Vec3i::set(int value) // set all elements to value
{
  elt[0] = value;
  elt[1] = value;
  elt[2] = value;
}

inline int Vec3i::getInvertedIndex(int value) const
{
  if (value == elt[0]) return 0;
  if (value == elt[1]) return 1;
  if (value == elt[2]) return 2;
  return -1;
}

inline void Vec3i::rotate(int newStartIndex)
{
  if (newStartIndex == 1) // rotate left one time
  {
    int tmp = elt[0];
    elt[0] = elt[1];
    elt[1] = elt[2];
    elt[2] = tmp;
  }
  else if (newStartIndex == 2) // rotate right one time
  {
    int tmp = elt[2];
    elt[2] = elt[1];
    elt[1] = elt[0];
    elt[0] = tmp;
  }
}

inline bool Vec3i::intersect(const Vec3i & vec2) const
{
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      if (elt[i] == vec2.elt[j])
        return true;
  return false;
}

#endif
