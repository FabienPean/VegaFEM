/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "hashTable" library , Copyright (C) 2007 CMU                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
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

/*
  A simple hash table. 
  Keys are of type 'unsigned int', and datatype can be arbitrary (templated).

  The hash function is:
  key * (key + 3) % p,
  where appropriate p is chosen automatically by the constructor.
*/

#ifndef VEGAFEM_HASHTABLE_H
#define VEGAFEM_HASHTABLE_H

#include <vector>
#include <list>

namespace vegafem
{

template<class Data>
class HashTable
{
public:
  HashTable(unsigned int suggestedSize);

  // clears the hash table
  void clear();
  
  // insert (key, data)
  void insert(unsigned int key, Data & data); // makes an internal copy of data; does not check for key duplication; if you want to use data that you allocated externally, pass pointers to your data as 'data'

  // erase entry with the given key
  void erase(unsigned int key); // does not check for duplicate entries (only erases the first matching entry)

  Data * find(unsigned int key); // returns pointer to internal copy of data if found, or NULL otherwise

  bool contains(unsigned int key, Data & data); // does (key,Data) pair appear in the hash table

  void printInfo(); // prints info on hash table efficiency

protected:
  inline unsigned int hashFunction(unsigned int key);
  bool isPrime(unsigned int n);
  
  typedef std::pair<unsigned int, Data > listEntry;
  std::vector<std::list<listEntry> > buckets;
  typedef typename std::list<listEntry>::iterator listIterator;

  unsigned int p;
};

template<class Data>
inline unsigned int HashTable<Data>::hashFunction(unsigned int key)
{
  return (key*(key+3) % p);
}


}//namespace vegafem

#endif


