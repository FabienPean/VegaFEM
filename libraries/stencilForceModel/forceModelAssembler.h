/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Stencil Force Model" library , Copyright (C) 2018 USC                *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Bohan Wang, Jernej Barbic                               *
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

#ifndef _FORCEMODEL_ASSEMBLER_H_
#define _FORCEMODEL_ASSEMBLER_H_

#include "forceModel.h"
#include "stencilForceModel.h"

#ifdef USE_TBB
  #include <tbb/tbb.h>
#endif

/*
  For each stencil type, this class assembles values at individual stencils into global object quantities.
  E.g., form the global internal force vector from stencil force vectors, or form the global tangent stiffness matrix
  from stencil tangent stiffness matrices.
  Please see also stencilForceModel.h for more information.
  This class inherits from ForceModel and implements the necessary interfaces.
  If Intel TBB is provided, assembly will be performed in parallel.
  The number of threads can be controlled outside the class using the Intel TBB APIs.
  If Intel TBB is not provided, the computation will be single-threaded.
*/

class ForceModelAssembler : public ForceModel
{
public:
  ForceModelAssembler(StencilForceModel *stencilForceModel);
  virtual ~ForceModelAssembler();

  // See comments in the parent class for the following functions.
  virtual double GetElasticEnergy(const double * u) override;
  virtual void GetInternalForce(const double * u, double * internalForces) override;
  virtual void GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix) override;
  virtual void GetTangentStiffnessMatrix(const double * u, SparseMatrix * tangentStiffnessMatrix) override;
  virtual void GetForceAndMatrix(const double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix) override;

  // This function computes the energy, internal forces and tangent stiffness matrix of the 'object'.
  // u is the displacement vector of all vertices.
  // energy points to a double variable.
  // internalForces points to a double array which dimension is the same as u.
  // tangentStiffnessMatrix point to a sparse matrix object.
  // energy, internalForces, tangentStiffnessMatrix can be nullptr. If nullptr, the corresponding quantity will not be computed.
  virtual void GetEnergyAndForceAndMatrix(const double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix);

protected:
  StencilForceModel * stencilForceModel = nullptr;
  SparseMatrix * Ktemplate = nullptr;
  std::vector<std::vector<int>> inverseIndices;

  // data structures for parallelism
#ifdef USE_TBB
  typedef tbb::cache_aligned_allocator<double> BufferAllocator;
  typedef std::vector<double, BufferAllocator> Buffer;
  std::vector<tbb::enumerable_thread_specific<Buffer>*> localBuffers;
  tbb::affinity_partitioner * partitioners = nullptr;
  tbb::enumerable_thread_specific<double> energyLocalBuffer;
  tbb::spin_mutex * internalForceVertexLocks, *stiffnessMatrixVertexRowLocks = nullptr;
#else
  // data structures for single-threaded computation
  typedef std::allocator<double> BufferAllocator;
  typedef std::vector<double, BufferAllocator> Buffer;
#endif

  std::vector<Buffer> bufferExamplars;
};

#endif

