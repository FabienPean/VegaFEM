/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "basicAlgorithms" library , Copyright (C) 2018 USC                    *
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

#ifndef UNIQUERECENTQUEUE_H_
#define UNIQUERECENTQUEUE_H_

#include "range.h"
#include <list>

template<class T>
class UniqueRecentQueue
{
public:
  UniqueRecentQueue(int maxSize) : maxSize(maxSize) {}

  void push(const T & v);
  void push(T && v);
  template<typename... _Args>
  void emplace(_Args &&... __args);

  typename std::list<T>::const_iterator begin() const { return l.begin(); }
  typename std::list<T>::const_iterator end() const { return l.end(); }

  std::size_t size() const { return l.size(); }
  Range<typename std::list<T>::const_iterator> last(int numElements) const;

protected:
  std::list<T> l;
  int maxSize = 0;
};

template<class T>
void UniqueRecentQueue<T>::push(const T & v)
{
  if (l.empty())
  {
    l.push_back(v);
  }
  else if (l.back() != v)
  {
    l.push_back(v);
    if (l.size() > size_t(maxSize))
      l.pop_front();
  }
}

template<class T>
void UniqueRecentQueue<T>::push(T && v)
{
  if (l.empty())
  {
    l.push_back(std::move(v));
  }
  else if (l.back() != v)
  {
    l.push_back(std::move(v));
    if (l.size() > maxSize)
      l.pop_front();
  }
}

template<class T> template<typename... _Args>
void UniqueRecentQueue<T>::emplace(_Args &&... __args)
{
  l.emplace_back(std::forward<_Args>(__args)...);
}

template<class T>
Range<typename std::list<T>::const_iterator> UniqueRecentQueue<T>::last(int numElements) const
{
  int n = l.size() - numElements;
  typename std::list<T>::const_iterator it = l.begin();
  for(int i = 0; i < n; i++, it++) {}
  return { it, l.end() };
}

#endif /* UNIQUERECENTQUEUE_H_ */
