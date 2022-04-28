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
 * Funding: National Science Foundation                                  *
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


#ifndef VEC2I_H
#define VEC2I_H

#include <stdio.h>
#include <math.h>
#include <ostream>

class Vec2i
{
public:
  inline Vec2i() {}
  inline Vec2i(int x, int y) {elt[0]=x; elt[1]=y;}
  inline Vec2i(int entry); // create a vector with all entries "entry" (can create zero vector for entry=0)
  inline Vec2i(const int vec[2]); // create a vector from the array of three ints pointed to by "vec"
  inline Vec2i(const Vec2i & vec);

  inline void set(int x0, int x1); // assign vector [x0, x1]
  inline void set(int value); // set all elements to value

  inline Vec2i & operator=(const Vec2i & source);
  inline bool operator==(const Vec2i & vec2) const;
  inline bool operator!=(const Vec2i & vec2) const;

  inline const Vec2i operator+ (const Vec2i & vec2) const;
  inline Vec2i & operator+= (const Vec2i & vec2);

  inline const Vec2i operator- (const Vec2i & vec2) const;
  inline Vec2i & operator-= (const Vec2i & vec2);

  inline const Vec2i operator* (int scalar) const;
  inline Vec2i & operator*= (int scalar);

  inline Vec2i operator/ (int scalar) const;
  inline Vec2i & operator/= (int scalar);

  // operator for Vec2i to be used as a key in std::set, std::map, etc.
  inline bool operator < (const Vec2i & vec2) const;

  friend inline Vec2i operator* (int scalar, const Vec2i & vec2);
  friend inline Vec2i operator/ (int scalar, const Vec2i & vec2);
  friend inline Vec2i operator- (const Vec2i & vec1);

  friend inline int dot(const Vec2i & vec1, const Vec2i & vec2); // dot product

  friend inline Vec2i cross(const Vec2i & vec1, const Vec2i & vec2); // cross product

  friend inline std::ostream &operator << (std::ostream & s, const Vec2i & v);
  void print() const;

  inline int & operator[] (int index); // v[i] returns i-th entry of v
  inline const int & operator[] (int index) const;

  // copy the vector into an array of length 2
  inline void convertToArray(int vecArray[2]) const;
  // add the vector into an array of length 2
  inline void addToArray(int vecArray[2]) const;

  // find the first index in elt which equals to value; return -1 if not found
  inline int getInvertedIndex(int value) const;
 
  // rotate elt[0] to elt[1] so that elt[newStartIndex] is moved to elt[0]
  inline void rotate(int newStartIndex);

  // do set intersection; return whether this and vec2 share at least one element
  inline bool intersect(const Vec2i & vec2) const;

  const int * begin() const { return elt; }
  int * begin() { return elt; }
  const int * end() const { return elt + 2; }
  int * end() { return elt + 2; }

protected:
  int elt[2];
};

// === below is the implementation ===

inline Vec2i::Vec2i(int entry)
{
  elt[0] = entry;
  elt[1] = entry;
}

inline Vec2i::Vec2i(const int vec[2])
{
  elt[0] = vec[0];
  elt[1] = vec[1];
}

inline Vec2i::Vec2i(const Vec2i & vec)
{
  elt[0] = vec.elt[0];
  elt[1] = vec.elt[1];
}

inline Vec2i & Vec2i::operator=(const Vec2i & source)
{
  elt[0] = source.elt[0];
  elt[1] = source.elt[1];

  return *this;
}

inline bool Vec2i::operator==(const Vec2i & vec2) const
{
  return ((elt[0] == vec2[0]) &&
          (elt[1] == vec2[1]));
}

inline bool Vec2i::operator!=(const Vec2i & vec2) const
{
  return ((elt[0] != vec2[0]) ||
          (elt[1] != vec2[1]));
}

inline bool Vec2i::operator<(const Vec2i & vec2) const
{
  if(elt[0] < vec2[0])
    return true;
  if(elt[0] > vec2[0])
    return false;
  return elt[1] < vec2[1];
}

inline Vec2i operator* (int scalar, const Vec2i & vec2)
{
  Vec2i result = vec2;
  result.elt[0] *= scalar;
  result.elt[1] *= scalar;

  return result;
}

inline Vec2i operator/ (int scalar, const Vec2i & vec2)
{
  Vec2i result = vec2;
  result.elt[0] /= scalar;
  result.elt[1] /= scalar;

  return result;
}

inline Vec2i operator- (const Vec2i & vec1)
{
  Vec2i result = vec1;
  result.elt[0] *= -1;
  result.elt[1] *= -1;

  return result;
}

inline const Vec2i Vec2i::operator+ (const Vec2i & vec2) const
{
  Vec2i sum = *this;
  sum.elt[0] += vec2.elt[0];
  sum.elt[1] += vec2.elt[1];

  return sum;
}

inline Vec2i & Vec2i::operator+= (const Vec2i & vec2)
{
  elt[0] += vec2.elt[0];
  elt[1] += vec2.elt[1];

  return *this;
}

inline const Vec2i Vec2i::operator- (const Vec2i & vec2) const
{
  Vec2i sum = *this;
  sum.elt[0] -= vec2.elt[0];
  sum.elt[1] -= vec2.elt[1];

  return sum;
}

inline Vec2i & Vec2i::operator-= (const Vec2i & vec2)
{
  elt[0] -= vec2.elt[0];
  elt[1] -= vec2.elt[1];

  return *this;
}

inline int & Vec2i::operator[] (int index)
{
  return elt[index];
}

inline const int & Vec2i::operator[] (int index) const
{
  return elt[index];
}

inline int dot(const Vec2i & vec1, const Vec2i & vec2)
{
  return (vec1.elt[0] * vec2.elt[0] + vec1.elt[1] * vec2.elt[1]);
}

inline Vec2i & Vec2i::operator*= (int scalar)
{
  elt[0] *= scalar;
  elt[1] *= scalar;
  return *this;
}

inline const Vec2i Vec2i::operator* (int scalar) const
{
  return (Vec2i(elt[0] * scalar, elt[1] * scalar));
}

inline Vec2i Vec2i::operator/ (int scalar) const
{
  return (Vec2i(elt[0] / scalar, elt[1] / scalar));
}

inline Vec2i & Vec2i::operator/= (int scalar)
{
  elt[0] /= scalar;
  elt[1] /= scalar;
  return *this;
}

inline std::ostream &operator << (std::ostream &s, const Vec2i &v)
{
  return(s << '[' << v[0] << ' ' << v[1] << ']');
}

inline void Vec2i::convertToArray(int vecArray[2]) const
{
  vecArray[0] = elt[0];
  vecArray[1] = elt[1];
}

inline void Vec2i::addToArray(int vecArray[2]) const
{
  vecArray[0] += elt[0];
  vecArray[1] += elt[1];
}

inline void Vec2i::print() const
{
  printf("[%d %d]\n", elt[0], elt[1]);
}

inline void Vec2i::set(int x0, int x1) // assign vector [x0, x1, x2]
{
  elt[0] = x0;
  elt[1] = x1;
}

inline void Vec2i::set(int value) // set all elements to value
{
  elt[0] = value;
  elt[1] = value;
}

inline int Vec2i::getInvertedIndex(int value) const
{
  if (value == elt[0]) return 0;
  if (value == elt[1]) return 1;
  return -1;
}

inline void Vec2i::rotate(int newStartIndex)
{
  if (newStartIndex == 1) // rotate left one time
  {
    std::swap(elt[0], elt[1]);
  }
}

inline bool Vec2i::intersect(const Vec2i & vec2) const
{
  for(int i = 0; i < 2; i++)
    for(int j = 0; j < 2; j++)
      if (elt[i] == vec2.elt[j])
        return true;
  return false;
}

#endif
