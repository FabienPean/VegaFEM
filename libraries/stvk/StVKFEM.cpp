/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "StVK" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Bohan Wang, Jernej Barbic                                *
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

#include "StVKFEM.h"
#include "volumetricMeshENuMaterial.h"

#include <cassert>
#include <cstdio>

StVKFEM::StVKFEM(VolumetricMesh * volumetricMesh_, StVKElementABCD * precomputedABCDIntegrals_, bool addGravity_, double g_):
  volumetricMesh(volumetricMesh_), precomputedIntegrals(precomputedABCDIntegrals_), g(g_), addGravity(addGravity_)
{
  int numElements = volumetricMesh->getNumElements();

  lambdaLame.resize(numElements, 0.0);
  muLame.resize(numElements, 0.0);
  internalElementData.resize(numElements, nullptr);

  for(int el=0; el<numElements; el++)
  {
    VolumetricMesh::Material * material = volumetricMesh->getElementMaterial(el);
    VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
    if (eNuMaterial == NULL)
    {
      printf("Error: StVKInternalForces: mesh does not consist of E, nu materials.\n");
      throw 1;
    }

    lambdaLame[el] = eNuMaterial->getLambda();
    muLame[el] = eNuMaterial->getMu();

    void * elIter;
    precomputedIntegrals->AllocateElementIterator(&elIter);
    precomputedIntegrals->PrepareElement(el, elIter);
    internalElementData[el] = elIter;
  }

  numElementVertices = volumetricMesh->getNumElementVertices();

  dof = numElementVertices * 3;
  dof2 = dof * dof;

  InitGravity();
}

StVKFEM::~StVKFEM()
{
  for (int el=0; el < volumetricMesh->getNumElements(); el++) {
   precomputedIntegrals->ReleaseElementIterator(internalElementData[el]);
  }
}

void StVKFEM::InitGravity()
{
  if (addGravity)
  {
    gravityForce.resize(3 * volumetricMesh->getNumVertices(), 0.0);
    volumetricMesh->computeGravity(gravityForce.data(), g);
  }
}

void StVKFEM::ComputeElementLocalEnergyAndInternalForcesAndStiffnessMatrix(const double *u, int el, double * energy, double * fint, double * K)
{
  void *elIter = internalElementData[el];
  const int *vertices = volumetricMesh->getVertexIndices(el);
  double lambda = lambdaLame[el]; 
  double mu = muLame[el];

  if (energy)
    *energy = 0;

  if (fint)
  {
    memset(fint, 0, dof * sizeof(double));

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
    {
      Vec3d qc(u[3*vertices[c]+0], u[3*vertices[c]+1], u[3*vertices[c]+2]);

      // linear terms
      for (int a=0; a<numElementVertices; a++) // over all vertices
      {
        Vec3d qa(u[3*vertices[a]+0], u[3*vertices[a]+1], u[3*vertices[a]+2]);

        Vec3d force = lambda * (precomputedIntegrals->A(elIter,c,a) * qa) +
                      (mu * precomputedIntegrals->B(elIter,a,c)) * qa +
                      mu * (precomputedIntegrals->A(elIter,a,c) * qa);

        fint[3*c+0] += force[0];
        fint[3*c+1] += force[1];
        fint[3*c+2] += force[2];

        if (energy)
          *energy += 0.5 * dot(qc, force);
      }
    }

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
    {
      Vec3d qc(u[3*vertices[c]+0], u[3*vertices[c]+1], u[3*vertices[c]+2]);

      // quadratic terms
      for (int a=0; a<numElementVertices; a++) // over all vertices
      {
        for(int b=0; b<numElementVertices; b++)
        {
/*
          Vec3d force(0,0,0);
          Vec3d qa(vertexDisplacements[3*vertices[a]+0],
                   vertexDisplacements[3*vertices[a]+1],
                   vertexDisplacements[3*vertices[a]+2]);

          Vec3d qb(vertexDisplacements[3*vertices[b]+0],
                   vertexDisplacements[3*vertices[b]+1],
                   vertexDisplacements[3*vertices[b]+2]);

          double dotp = dot(qa,qb);

          force += 0.5 * lambda * dotp * precomputedIntegrals->C(el,c,a,b) +
                   mu * dotp * precomputedIntegrals->C(el,a,b,c);

          Vec3d C = lambda * precomputedIntegrals->C(el,a,b,c) +
                    mu * (precomputedIntegrals->C(el,c,a,b) + precomputedIntegrals->C(el,b,a,c)); 

          force += dot(C,qa) * qb;

          forces[3*vertices[c]+0] += force[0];
          forces[3*vertices[c]+1] += force[1];
          forces[3*vertices[c]+2] += force[2];
*/

          double qa[3] = { u[3*vertices[a]+0], u[3*vertices[a]+1], u[3*vertices[a]+2] };
          double qb[3] = { u[3*vertices[b]+0], u[3*vertices[b]+1], u[3*vertices[b]+2] };

          double dotp = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2];

          Vec3d forceTerm1 = 0.5 * lambda * dotp * precomputedIntegrals->C(elIter,c,a,b) +
                             mu * dotp * precomputedIntegrals->C(elIter,a,b,c);

          Vec3d C = lambda * precomputedIntegrals->C(elIter,a,b,c) +
                    mu * (precomputedIntegrals->C(elIter,c,a,b) + precomputedIntegrals->C(elIter,b,a,c)); 

          double dotCqa = C[0] * qa[0] + C[1] * qa[1] + C[2] * qa[2];

          double force[3] = 
          {
            forceTerm1[0] + dotCqa * qb[0],
            forceTerm1[1] + dotCqa * qb[1],
            forceTerm1[2] + dotCqa * qb[2],
          };

          fint[3*c+0] += force[0];
          fint[3*c+1] += force[1];
          fint[3*c+2] += force[2];

          if (energy)
            *energy += dot(qc, Vec3d(force)) / 3.0;
        }
      }
    }

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing force on vertex c
    {
      int vc = vertices[c];
      Vec3d qc(u[3*vc+0], u[3*vc+1], u[3*vc+2]);
      // cubic terms
      for(int a=0; a<numElementVertices; a++) // over all vertices
      {
        int va = vertices[a];
        for(int b=0; b<numElementVertices; b++)
        {
          int vb = vertices[b];
          for(int d=0; d<numElementVertices; d++)
          {
            int vd = vertices[d];
/*
            Vec3d qa(vertexDisplacements[3*va+0],
                     vertexDisplacements[3*va+1],
                     vertexDisplacements[3*va+2]);

            Vec3d qb(vertexDisplacements[3*vb+0],
                     vertexDisplacements[3*vb+1],
                     vertexDisplacements[3*vb+2]);

            Vec3d qd(vertexDisplacements[3*vd+0],
                     vertexDisplacements[3*vd+1],
                     vertexDisplacements[3*vd+2]);

            double dotp = dot(qa,qb);

            Vec3d force = 0.5 * lambda * dotp * precomputedIntegrals_->D(a,b,c,d) * qd +
                          mu * dotp * precomputedIntegrals_->D(a,c,b,d) * qd;

            forces[3*vertices[c]+0] += force[0];
            forces[3*vertices[c]+1] += force[1];
            forces[3*vertices[c]+2] += force[2];
*/
            const double * qa = &(u[3*va]);
            const double * qb = &(u[3*vb]);
            const double * qd = &(u[3*vd]);

            double dotp = qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2]; 
            double scalar = dotp * (0.5 * lambda * precomputedIntegrals->D(elIter,a,b,c,d) + mu * precomputedIntegrals->D(elIter,a,c,b,d) );

            double force[3] = { scalar * qd[0], scalar * qd[1], scalar * qd[2] };

            fint[3*c+0] += force[0];
            fint[3*c+1] += force[1];
            fint[3*c+2] += force[2];

            if (energy)
              *energy += dot(qc, Vec3d(force)) / 4.0;
          }
        }
      }
    }
  }

  if (K)
  {    
    memset(K, 0, dof2 * sizeof(double));

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing row of vertex c
    {
      // linear terms
      for (int a=0; a<numElementVertices; a++) // over all vertices
      {
        Mat3d matrix(1.0);
        matrix *= mu * precomputedIntegrals->B(elIter,a,c);
        matrix += lambda * precomputedIntegrals->A(elIter,c,a) +
                  mu * precomputedIntegrals->A(elIter,a,c);

        AddMatrix3x3Block(c, a, matrix, K);
      }
    }

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing row of vertex c
    {
      // quadratic terms
      for (int e=0; e<numElementVertices; e++) // compute contribution to block (c,e) of the stiffness matrix
      {
        double matrix[9];
        memset(matrix, 0, sizeof(double) * 9);

        for(int a=0; a<numElementVertices; a++)
        {
          double qa[3] = { u[3*vertices[a]+0], u[3*vertices[a]+1], u[3*vertices[a]+2] };

          Vec3d C0v = lambda * precomputedIntegrals->C(elIter,c,a,e) + mu * (precomputedIntegrals->C(elIter,e,a,c) + precomputedIntegrals->C(elIter,a,e,c));
          double C0[3] = {C0v[0], C0v[1], C0v[2]};

          // C0 tensor qa
          matrix[0] += C0[0] * qa[0]; matrix[1] += C0[0] * qa[1]; matrix[2] += C0[0] * qa[2];
          matrix[3] += C0[1] * qa[0]; matrix[4] += C0[1] * qa[1]; matrix[5] += C0[1] * qa[2];
          matrix[6] += C0[2] * qa[0]; matrix[7] += C0[2] * qa[1]; matrix[8] += C0[2] * qa[2];

          Vec3d C1v = lambda * precomputedIntegrals->C(elIter,e,a,c) + mu * (precomputedIntegrals->C(elIter,c,e,a) + precomputedIntegrals->C(elIter,a,e,c));
          double C1[3] = {C1v[0], C1v[1], C1v[2]};

          // qa tensor C1
          matrix[0] += qa[0] * C1[0]; matrix[1] += qa[0] * C1[1]; matrix[2] += qa[0] * C1[2];
          matrix[3] += qa[1] * C1[0]; matrix[4] += qa[1] * C1[1]; matrix[5] += qa[1] * C1[2];
          matrix[6] += qa[2] * C1[0]; matrix[7] += qa[2] * C1[1]; matrix[8] += qa[2] * C1[2];

          Vec3d C2v = lambda * precomputedIntegrals->C(elIter,a,e,c) + mu * (precomputedIntegrals->C(elIter,c,a,e) + precomputedIntegrals->C(elIter,e,a,c));
          double C2[3] = {C2v[0], C2v[1], C2v[2]};

          // qa dot C2
          double dotp = qa[0]*C2[0] + qa[1]*C2[1] + qa[2]*C2[2];
          matrix[0] += dotp; 
          matrix[4] += dotp; 
          matrix[8] += dotp;
        }

        AddMatrix3x3Block(c, e, matrix, K);
      }
    }

    for (int c=0; c<numElementVertices; c++) // over all vertices of the voxel, computing derivative on force on vertex c
    {
      // cubic terms
      for (int e=0; e<numElementVertices; e++) // compute contribution to block (c,e) of the stiffness matrix
      {
        double matrix[9];
        memset(matrix, 0, sizeof(double) * 9);
        for(int a=0; a<numElementVertices; a++)
        {
          int va = vertices[a];
          const double * qa = &(u[3*va]);
          for(int b=0; b<numElementVertices; b++)
          {
            int vb = vertices[b];

            const double * qb = &(u[3*vb]);

            double D0 = lambda * precomputedIntegrals->D(elIter,a,c,b,e) +
                        mu * ( precomputedIntegrals->D(elIter,a,e,b,c) + precomputedIntegrals->D(elIter,a,b,c,e) );

            matrix[0] += D0 * qa[0] * qb[0]; matrix[1] += D0 * qa[0] * qb[1]; matrix[2] += D0 * qa[0] * qb[2];
            matrix[3] += D0 * qa[1] * qb[0]; matrix[4] += D0 * qa[1] * qb[1]; matrix[5] += D0 * qa[1] * qb[2];
            matrix[6] += D0 * qa[2] * qb[0]; matrix[7] += D0 * qa[2] * qb[1]; matrix[8] += D0 * qa[2] * qb[2];

            double D1 = 0.5 * lambda * precomputedIntegrals->D(elIter,a,b,c,e) +
                        mu * precomputedIntegrals->D(elIter,a,c,b,e);

            double dotpD = D1 * (qa[0] * qb[0] + qa[1] * qb[1] + qa[2] * qb[2]);

            matrix[0] += dotpD; 
            matrix[4] += dotpD; 
            matrix[8] += dotpD; 
          }
        }

        AddMatrix3x3Block(c, e, matrix, K);
      }
    }
  }
}
