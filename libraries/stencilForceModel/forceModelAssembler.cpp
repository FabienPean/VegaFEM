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

#include "forceModelAssembler.h"
#include <cassert>

using namespace std;

ForceModelAssembler::ForceModelAssembler(StencilForceModel *eleFM) : stencilForceModel(eleFM)
{
  r = stencilForceModel->Getn3();
  SparseMatrixOutline *smo = new SparseMatrixOutline(r);
  SparseMatrixOutline *smo1 = new SparseMatrixOutline(r / 3);

  for (int eltype = 0; eltype < stencilForceModel->GetNumStencilTypes(); eltype++) 
  {
    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    for (int ele = 0; ele < stencilForceModel->GetNumStencils(eltype); ele++) 
    {
      const int *vertexIndices = stencilForceModel->GetStencilVertexIndices(eltype, ele);

      for (int vi = 0; vi < nelev; vi++) 
      {
        for (int vj = 0; vj < nelev; vj++) 
        {
          smo->AddBlock3x3Entry(vertexIndices[vi], vertexIndices[vj]);
          smo1->AddEntry(vertexIndices[vi], vertexIndices[vj]);
        }
      }
    }
  }

  // compute stiffness matrix topology
  Ktemplate = new SparseMatrix(smo);
  delete smo;

  SparseMatrix *vertexK = new SparseMatrix(smo1);
  delete smo1;

  inverseIndices.resize(stencilForceModel->GetNumStencilTypes());
  for (int eltype = 0; eltype < stencilForceModel->GetNumStencilTypes(); eltype++) 
  {
    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    std::vector<int> &indices = inverseIndices[eltype];
    indices.resize(nelev * nelev * stencilForceModel->GetNumStencils(eltype));

    for (int ele = 0; ele < stencilForceModel->GetNumStencils(eltype); ele++) 
    {
      const int *vertexIndices = stencilForceModel->GetStencilVertexIndices(eltype, ele);
      int *inverseVtxIdx = indices.data() + ele * nelev * nelev;

      for (int vi = 0; vi < nelev; vi++) 
      {
        for (int vj = 0; vj < nelev; vj++) 
        {
          int vtxColIdx = vertexK->GetInverseIndex(vertexIndices[vi], vertexIndices[vj]);
          assert(vtxColIdx >= 0);
          inverseVtxIdx[vj * nelev + vi] = vtxColIdx;
        } // end vi
      } // end vj
    } // end ele
  } // end ele type

  delete vertexK;

  // initialize all necessary buffers
  bufferExamplars.resize(stencilForceModel->GetNumStencilTypes());
  for (int eltype = 0; eltype < stencilForceModel->GetNumStencilTypes(); eltype++) 
  {
    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    bufferExamplars[eltype].resize(nelev * 3 + nelev * nelev * 9);
  }

#ifdef USE_TBB
  localBuffers.resize(stencilForceModel->GetNumStencilTypes());
  for (int eltype = 0; eltype < stencilForceModel->GetNumStencilTypes(); eltype++) 
  {
    localBuffers[eltype] = new tbb::enumerable_thread_specific<Buffer>(bufferExamplars[eltype]);
  }

  partitioners = new tbb::affinity_partitioner[stencilForceModel->GetNumStencilTypes()];

  internalForceVertexLocks = new tbb::spin_mutex[r / 3];
  stiffnessMatrixVertexRowLocks = new tbb::spin_mutex[r / 3];
#endif
}

ForceModelAssembler::~ForceModelAssembler()
{
#ifdef USE_TBB
  delete[] internalForceVertexLocks;
  delete[] stiffnessMatrixVertexRowLocks;
  delete[] partitioners;
#endif
}

void ForceModelAssembler::GetTangentStiffnessMatrixTopology(SparseMatrix ** tangentStiffnessMatrix)
{
  *tangentStiffnessMatrix = new SparseMatrix(*Ktemplate);
  std::vector<double> zeros(stencilForceModel->Getn3(), 0.0);
  GetEnergyAndForceAndMatrix(zeros.data(), nullptr, nullptr, *tangentStiffnessMatrix);
}

void ForceModelAssembler::GetEnergyAndForceAndMatrix(const double * u, double * energy, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
  // reset to zero
  if (internalForces)
    memset(internalForces, 0, sizeof(double) * r);

  if (tangentStiffnessMatrix)
    tangentStiffnessMatrix->ResetToZero();

#ifdef USE_TBB
  for (auto itt = energyLocalBuffer.begin(); itt != energyLocalBuffer.end(); ++itt)
    *itt = 0.0;

  tbb::parallel_for(0, stencilForceModel->GetNumStencilTypes(), 1, [&] (int eltype) 
  {
    tbb::enumerable_thread_specific<Buffer> &tls = *localBuffers[eltype];
    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    int nele = stencilForceModel->GetNumStencils(eltype);

    tbb::parallel_for(0, nele, 1, [&] (int ele) 
    {
      Buffer &localBuffer = tls.local();
      double *fEle = localBuffer.data();
      double *KEle = localBuffer.data() + nelev * 3;

      double energyEle = 0;

      stencilForceModel->GetStencilLocalEnergyAndForceAndMatrix(eltype, ele, u,
        (energy ? &energyEle : nullptr),
        (internalForces ? fEle : nullptr),
        (tangentStiffnessMatrix ? KEle : nullptr)
      );

      const int *vIndices = stencilForceModel->GetStencilVertexIndices(eltype, ele);

      if (internalForces) 
      {
        for (int v = 0; v < nelev; v++) 
        {
          internalForceVertexLocks[vIndices[v]].lock();

          internalForces[vIndices[v] * 3] += fEle[v * 3];
          internalForces[vIndices[v] * 3 + 1] += fEle[v * 3 + 1];
          internalForces[vIndices[v] * 3 + 2] += fEle[v * 3 + 2];

          internalForceVertexLocks[vIndices[v]].unlock();
        }
      }

      if (tangentStiffnessMatrix) 
      {
        const int *vtxColIndices = inverseIndices[eltype].data() + ele * nelev * nelev;

        // write matrices in place
        for (int va = 0; va < nelev; va++) 
        {
          int vIdxA = vIndices[va];
          stiffnessMatrixVertexRowLocks[vIdxA].lock();

          for (int vb = 0; vb < nelev; vb++) 
          {
            int columnIndexCompressed = vtxColIndices[vb * nelev + va];

            for (int i = 0; i < 3; i++) 
            {
              for (int j = 0; j < 3; j++) 
              {
                int row = 3 * vIdxA + i;
                int columnIndex = 3 * columnIndexCompressed + j;

                int local_row = 3 * va + i;
                int local_col = 3 * vb + j;

                tangentStiffnessMatrix->AddEntry(row, columnIndex, KEle[local_col * nelev * 3 + local_row]);
              } // i
            } // j
          } // vb
          stiffnessMatrixVertexRowLocks[vIdxA].unlock();
        } // va
      }

      if (energy) 
      {
        energyLocalBuffer.local() += energyEle;
      }
    }, partitioners[eltype]);
  });

  if (energy) 
  {
    *energy = 0;
    for (auto itt = energyLocalBuffer.begin(); itt != energyLocalBuffer.end(); ++itt) 
    {
      *energy += *itt;
    }
  }
#else
  for (int eltype = 0; eltype < stencilForceModel->GetNumStencilTypes(); eltype++) 
  {
    int nelev = stencilForceModel->GetNumStencilVertices(eltype);
    int nele = stencilForceModel->GetNumStencils(eltype);

    for (int ele = 0; ele < nele; ele++) 
    {
      double *fEle = bufferExamplars[eltype].data();
      double *KEle = bufferExamplars[eltype].data() + nelev * 3;

      double energyEle = 0;

      stencilForceModel->GetStencilLocalEnergyAndForceAndMatrix(eltype, ele, u,
        (energy ? &energyEle : nullptr),
        (internalForces ? fEle : nullptr),
        (tangentStiffnessMatrix ? KEle : nullptr)
      );

      const int *vIndices = stencilForceModel->GetStencilVertexIndices(eltype, ele);

      if (internalForces) 
      {
        for (int v = 0; v < nelev; v++) 
        {
          internalForces[vIndices[v] * 3] += fEle[v * 3];
          internalForces[vIndices[v] * 3 + 1] += fEle[v * 3 + 1];
          internalForces[vIndices[v] * 3 + 2] += fEle[v * 3 + 2];
        }
      }

      if (tangentStiffnessMatrix) 
      {
        const int *vtxColIndices = inverseIndices[eltype].data() + ele * nelev * nelev;

        // write matrices in place
        for (int va = 0; va < nelev; va++) 
        {
          int vIdxA = vIndices[va];
          for (int vb = 0; vb < nelev; vb++) 
          {
            int columnIndexCompressed = vtxColIndices[vb * nelev + va];

            for (int i = 0; i < 3; i++) 
            {
              for (int j = 0; j < 3; j++) 
              {
                int row = 3 * vIdxA + i;
                int columnIndex = 3 * columnIndexCompressed + j;

                int local_row = 3 * va + i;
                int local_col = 3 * vb + j;

                tangentStiffnessMatrix->AddEntry(row, columnIndex, KEle[local_col * nelev * 3 + local_row]);
              } // i
            } // j
          } // vb
        } // va
      }

      if (energy) 
      {
        *energy += energyEle;
      }
    }
  }
#endif

  if (internalForces)
  {
    for (int vi = 0; vi < stencilForceModel->Getn3() / 3; vi++) {
      double g[3];
      stencilForceModel->GetVertexGravityForce(vi, g);
      internalForces[vi * 3 + 1] += g[1];
    }
  }
}

double ForceModelAssembler::GetElasticEnergy(const double *u)
{
  double E = 0;
  GetEnergyAndForceAndMatrix(u, &E, nullptr, nullptr);

  return E;
}

void ForceModelAssembler::GetInternalForce(const double * u, double * internalForces)
{
  GetEnergyAndForceAndMatrix(u, nullptr, internalForces, nullptr);
}

void ForceModelAssembler::GetTangentStiffnessMatrix(const double * u, SparseMatrix *tangentStiffnessMatrix)
{
  GetEnergyAndForceAndMatrix(u, nullptr, nullptr, tangentStiffnessMatrix);
}

void ForceModelAssembler::GetForceAndMatrix(const double * u, double * internalForces, SparseMatrix * tangentStiffnessMatrix)
{
  GetEnergyAndForceAndMatrix(u, nullptr, internalForces, tangentStiffnessMatrix);
}
