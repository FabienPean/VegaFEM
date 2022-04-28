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


#ifndef VEC4I_H
#define VEC4I_H


#include <stdio.h>
#include <math.h>
#include <ostream>
#include <algorithm>

class Vec4i
{
public:
  inline Vec4i() {}
  inline Vec4i(int x, int y, int z, int w) {elt[0]=x; elt[1]=y; elt[2]=z; elt[3]=w; }
  inline Vec4i(int entry); // create a vector with all entries "entry" (can create zero vector for entry=0)
  inline Vec4i(const int vec[4]); // create a vector from the array of four ints pointed to by "vec"
  inline Vec4i(const Vec4i & vec);

  inline void set(int x0, int x1, int x2, int x3); // assign vector [x0, x1, x2, x3]
  inline void set(int value); // set all elements to value

  inline Vec4i & operator=(const Vec4i & source);
  inline bool operator==(const Vec4i & vec2) const;
  inline bool operator!=(const Vec4i & vec2) const;

  inline const Vec4i operator+ (const Vec4i & vec2) const;
  inline Vec4i & operator+= (const Vec4i & vec2);

  inline const Vec4i operator- (const Vec4i & vec2) const;
  inline Vec4i & operator-= (const Vec4i & vec2);

  inline const Vec4i operator* (int scalar) const;
  inline Vec4i & operator*= (int scalar);

  inline Vec4i operator/ (int scalar) const;
  inline Vec4i & operator/= (int scalar);

  // operator for Vec4i to be used as a key in std::set, std::map, etc.
  inline bool operator < (const Vec4i & vec2) const;

  friend inline Vec4i operator* (int scalar, const Vec4i & vec2);
  friend inline Vec4i operator/ (int scalar, const Vec4i & vec2);
  friend inline Vec4i operator- (const Vec4i & vec1);

  friend inline int dot(const Vec4i & vec1, const Vec4i & vec2); // dot product

  friend inline std::ostream &operator << (std::ostream & s, const Vec4i & v);
  void print() const;

  inline int & operator[] (int index); // v[i] returns i-th entry of v
  inline const int & operator[] (int index) const;

  const int * data() const { return &elt[0]; }
  int * data() { return &elt[0]; }

  // copy the vector into an array of length 4
  inline void convertToArray(int vecArray[4]) const;
  // add the vector into an array of length 4
  inline void addToArray(int vecArray[4]) const;

  // find the first index in elt which equals to value; return -1 if not found
  inline int getInvertedIndex(int value) const;
  // rotate elt[0] to elt[4] so that elt[newStartIndex] is moved to elt[0]
  inline void rotate(int newStartIndex);

  // do set intersection; return whether this and vec2 share at least one element
  inline bool intersect(const Vec4i & vec2) const;

  const int * begin() const { return elt; }
  int * begin() { return elt; }
  const int * end() const { return elt + 4; }
  int * end() { return elt + 4; }

protected:
  int elt[4];
};

// === below is the implementation ===

inline Vec4i::Vec4i(int entry)
{
  elt[0] = entry;
  elt[1] = entry;
  elt[2] = entry;
  elt[3] = entry;
}

inline Vec4i::Vec4i(const int vec[4])
{
  elt[0] = vec[0];
  elt[1] = vec[1];
  elt[2] = vec[2];
  elt[3] = vec[3];
}

inline Vec4i::Vec4i(const Vec4i & vec)
{
  elt[0] = vec.elt[0];
  elt[1] = vec.elt[1];
  elt[2] = vec.elt[2];
  elt[3] = vec.elt[3];
}

inline Vec4i & Vec4i::operator=(const Vec4i & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];
  elt[2] = source.elt[2];
  elt[3] = source.elt[3];

  return *this;
}

inline bool Vec4i::operator==(const Vec4i & vec2) const
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]) &&
          (elt[2] == vec2[2]) &&
          (elt[3] == vec2[3]));
}

inline bool Vec4i::operator!=(const Vec4i & vec2) const
{
  return ((elt[0] != vec2[0]) ||
          (elt[1] != vec2[1]) ||
          (elt[2] != vec2[2]) ||
          (elt[3] != vec2[3]));
}

inline bool Vec4i::operator<(const Vec4i & vec2) const
{
  if(elt[0] < vec2[0])
    return true;
  if(elt[0] > vec2[0])
    return false;
  if(elt[1] < vec2[1])
    return true;
  if(elt[1] > vec2[1])
    return false;
  if(elt[2] < vec2[2])
      return true;
  if(elt[2] > vec2[2])
      return false;
  return elt[3] < vec2[3];
}

inline Vec4i operator* (int scalar, const Vec4i & vec2)
{
  Vec4i result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;
  result.elt[2] *= scalar;
  result.elt[3] *= scalar;

  return result;
}

inline Vec4i operator/ (int scalar, const Vec4i & vec2)
{
  Vec4i result = vec2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;
  result.elt[2] /= scalar;
  result.elt[3] /= scalar;

  return result;
}

inline Vec4i operator- (const Vec4i & vec1)
{
  Vec4i result = vec1;
  result.elt[0] *= -1;
  result.elt[1] *= -1;
  result.elt[2] *= -1;
  result.elt[3] *= -1;

  return result;
}

inline const Vec4i Vec4i::operator+ (const Vec4i & vec2) const
{
  Vec4i sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];
  sum.elt[2] += vec2.elt[2];
  sum.elt[3] += vec2.elt[3];

  return sum;
}

inline Vec4i & Vec4i::operator+= (const Vec4i & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];
  elt[2] += vec2.elt[2];
  elt[3] += vec2.elt[3];

  return *this;
}

inline const Vec4i Vec4i::operator- (const Vec4i & vec2) const
{
  Vec4i sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];
  sum.elt[2] -= vec2.elt[2];
  sum.elt[3] -= vec2.elt[3];

  return sum;
}

inline Vec4i & Vec4i::operator-= (const Vec4i & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];
  elt[2] -= vec2.elt[2];
  elt[3] -= vec2.elt[3];

  return *this;
}

inline int & Vec4i::operator[] (int index)
{
  return elt[index];
}

inline const int & Vec4i::operator[] (int index) const
{
  return elt[index];
}

inline int dot(const Vec4i & vec1, const Vec4i & vec2)
{
  return (vec1.elt[0] * vec2.elt[0] + vec1.elt[1] * vec2.elt[1] + vec1.elt[2] * vec2.elt[2] + vec1.elt[3] * vec2.elt[3]);
}

inline Vec4i & Vec4i::operator*= (int scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  elt[2] *= scalar;
  elt[3] *= scalar;
  return *this;
}

inline const Vec4i Vec4i::operator* (int scalar) const
{
  return (Vec4i(elt[0] * scalar, elt[1] * scalar, elt[2] * scalar, elt[3] * scalar));
}

inline Vec4i Vec4i::operator/ (int scalar) const
{
  return (Vec4i(elt[0] / scalar, elt[1] / scalar, elt[2] / scalar, elt[3] / scalar));
}

inline Vec4i & Vec4i::operator/= (int scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  elt[2] /= scalar;
  elt[3] /= scalar;
  return *this;
}

inline std::ostream &operator << (std::ostream &s, const Vec4i &v)
{
  return(s << '[' << v[0] << ' ' << v[1] << ' ' << v[2] << ' ' << v[3] << ']');
}

inline void Vec4i::convertToArray(int vecArray[4]) const
{
  vecArray[0] = elt[0];
  vecArray[1] = elt[1];
  vecArray[2] = elt[2];
  vecArray[3] = elt[3];
}

inline void Vec4i::addToArray(int vecArray[4]) const
{
  vecArray[0] += elt[0];
  vecArray[1] += elt[1];
  vecArray[2] += elt[2];
  vecArray[3] += elt[3];
}

inline void Vec4i::print() const
{
  int a = elt[0];
  int b = elt[1];
  int c = elt[2];
  int d = elt[3];

  printf("[%d %d %d %d]\n", a, b, c, d);
}

inline void Vec4i::set(int x0, int x1, int x2, int x3) // assign vector [x0, x1, x2]
{
  elt[0] = x0;
  elt[1] = x1;
  elt[2] = x2;
  elt[3] = x3;
}

inline void Vec4i::set(int value) // set all elements to value
{
  elt[0] = value;
  elt[1] = value;
  elt[2] = value;
  elt[3] = value;
}

inline int Vec4i::getInvertedIndex(int value) const
{
  if (value == elt[0]) return 0;
  if (value == elt[1]) return 1;
  if (value == elt[2]) return 2;
  if (value == elt[3]) return 3;
  return -1;
}

inline void Vec4i::rotate(int newStartIndex)
{
  if (newStartIndex == 1) // rotate left one time, (1, 2, 3, 0)
  {
    int tmp = elt[0];
    elt[0] = elt[1];
    elt[1] = elt[2];
    elt[2] = elt[3];
    elt[3] = tmp;
  }
  else if (newStartIndex == 2) // rotate left two times (2, 3, 0, 1)
  {
    std::swap(elt[0], elt[2]);
    std::swap(elt[1], elt[3]);
  }
  else if (newStartIndex == 3) // rotate right one time (3, 0, 1, 2)
  {
    int tmp = elt[3];
    elt[3] = elt[2];
    elt[2] = elt[1];
    elt[1] = elt[0];
    elt[0] = tmp;
  }
}

inline bool Vec4i::intersect(const Vec4i & vec2) const
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      if (elt[i] == vec2.elt[j])
        return true;
  return false;
}

#endif
