/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "corotational linear FEM" library , Copyright (C) 2018 USC            *
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

#include "corotationalLinearFEM.h"
#include "polarDecomposition.h"
#include "matrixMultiplyMacros.h"
#include "mat3d.h"
#include "volumetricMeshENuMaterial.h"
#include "volumetricMeshOrthotropicMaterial.h"
#include "cubicMesh.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
using namespace std;

// Corotational linear FEM for tet meshes is implemented by following:
// 
// M. Mueller, M. Gross: Interactive Virtual Materials.
// In Proc. of Graphics Interface 2004 (2004), pp. 239–246.
//
// Details can be found in the above paper.

// Corotational linear FEM for cubic meshes is implemented by following:
// Jesus Perez, Alvaro G. Perez and Miguel A. Otaduy:
// Simulation of Hyperelastic Materials Using Energy Constraints
// Proc. of Congreso Espanol de Informatica Grafica, 2013
//
// Details can be found in the above paper. Some additional notes here:
//
// We use natural coordinates s = (s_1, s_2, s_3)^T , where s_i belongs to [-1, 1].
// The natural coordinates for eight cubic nodes are: sn = ((-1)^i, (-1)^j, (-1)^k)^T, n = 4(i-1) + 2(j-1) + k, for i,j,k = {1,2}.
// So the shape functions are: N_n(s) = (1+sn_1 s_1) (1+sn_2 s_2) (1+sn_3 s_3) / 8 .
// The interpolation scheme based on shape functions is: y(s) = \sum_{n=1 to 8} yn N_n(s) .
// The element stiffness matrix at rest pose is:
// Ke = \integrate_X B'(X)^T E B'(X) dV_X , where X is the point in undeformed space, E is the 6x6 elasticity stiffness tensor
//    = \integratr_s B(s)^T E B(s) w dV_s, and s are the natural coordinates, and
// W = \partial{X} / \partial{s} is the scaling matrix from natural to material space.
// For cubic meshes, W can be just a scaled identity, W = w I.
// In the above, B(s) is a 6x24 matrix: B(s) = (diag(\partial{N_1}/\partial{s}), ..., diag(\partial{N_}/\partial{s}) ...).
// For a pair of corresponding X and s, B'(x) = B(s) * w^(-1).
// Due to the scaling from natural to material space, we have dV_X = w^3 dV_s .
// We approximate the integration of the elastic energy with second-order Gaussian quadrature:
// s_q = ((-1)^i, (-1)^j, (-1)^k)^T / \sqrt(3), q = 4(i-1) + 2(j-1) + k, for i,j,k = {1,2}
// weights w_q = 1
// So, Ke = \sum_q w_q B(s_q)^T E B(s_q) w .
//
// Shape functions are:
// n=1, i,j,k = 1,1,1, N_1 = (1-s1)(1-s2)(1-s3) / 8
// n=2, i,j,k = 1,1,2, N_2 = (1-s1)(1-s2)(1+s3) / 8
// n=3, i,j,k = 1,2,1, N_3 = (1-s1)(1+s2)(1-s3) / 8
// n=4, i,j,k = 1,2,2, N_4 = (1-s1)(1+s2)(1+s3) / 8
// n=5, i,j,k = 2,1,1, N_5 = (1+s1)(1-s2)(1-s3) / 8
// n=6, i,j,k = 2,1,2, N_6 = (1+s1)(1-s2)(1+s3) / 8
// n=7, i,j,k = 2,2,1, N_7 = (1+s1)(1+s2)(1-s3) / 8
// n=8, i,j,k = 2,2,2, N_8 = (1+s1)(1+s2)(1+s3) / 8
//
// B(s) = [-(1-s2)(1-s3) 0 0  ...  (1+s2)(1+s3) 0 0] / 8
//        [0 -(1-s1)(1-s3) 0  ... 0  (1+s1)(1+s3) 0]
//        [0 0 -(1-s1)(1-s2)  ... 0 0  (1+s1)(1+s2)]
//        [-(1-s1)(1-s3) -(1-s2)(1-s3)             ]
//        [...                                     ]
//        [...                                     ]
//
// Vertex order for cubic meshes (same as in the CubicMesh class) :
//
//     3 - - - 2
//    /|      /|
//   7 - - - 6 |       y
//   | |     | |       |
//   | 0 - - | 1       |_ _ _x
//   |/      |/       /
//   4 - - - 5       z
//
// We choose the tet (3,4,1,6) to compute the rotation component of the cube deformation:
static int cubeRotationTetIndex[4] = {3, 4, 1, 6};

// MInverse is pre-allocated
void CorotationalLinearFEM::buildMInverse(double * MInverse, VolumetricMesh * volumetricMesh) {

  int numElements = volumetricMesh->getNumElements();
  VolumetricMesh::elementType type = volumetricMesh->getElementType();
  int vtxIndex[4] = {0};
  for(int el = 0; el < numElements; el++)
  {
    // get the integer indices of the tet vertices
    if(type == VolumetricMesh::TET)
    {
      for(int vtx=0; vtx<4; vtx++)
        vtxIndex[vtx] = volumetricMesh->getVertexIndex(el, vtx);
    }
    else if(type == VolumetricMesh::CUBIC)
    {
      for(int i = 0; i < 4; i++)
        vtxIndex[i] = volumetricMesh->getVertexIndex(el, cubeRotationTetIndex[i]);
    }

    //  Form matrix:
    //  M = [ v0   v1   v2   v3 ]
    //      [  1    1    1    1 ]
    double M[16]; // row-major
    for(int vtx=0; vtx<4; vtx++)
      for(int dim=0; dim<3; dim++)
        M[4 * dim + vtx] = (volumetricMesh->getVertex(vtxIndex[vtx]))[dim]; //undeformedPositions[3 * vtxIndex[vtx] + dim];
    M[12] = M[13] = M[14] = M[15] = 1.0;

    inverse4x4(M, MInverse + el * 16);
  }
}

CorotationalLinearFEM::CorotationalLinearFEM(VolumetricMesh * volumetricMesh_) : volumetricMesh(volumetricMesh_)
{
  VolumetricMesh::elementType type = volumetricMesh->getElementType();
  if(type != VolumetricMesh::TET && type != VolumetricMesh::CUBIC)
  {
    printf("Error: CorotationalLinearFEM: unknown element type for the volumetric mesh.\n");
    throw 1;
  }

  numVertices = volumetricMesh->getNumVertices();
  int numElements = volumetricMesh->getNumElements();

  // store the undeformed positions
  undeformedPositions = (double*) malloc (sizeof(double) * 3 * numVertices);
  for(int i=0; i < numVertices; i++)
  {
    const Vec3d & v = volumetricMesh->getVertex(i);
    for(int j=0; j<3; j++)
      undeformedPositions[3*i+j] = v[j];
  }

  MInverse = (double**) malloc (sizeof(double*) * numElements);
  for(int el = 0; el < numElements; el++)
  {
    // get the integer indices of the tet vertices
    int vtxIndex[4];
    if(type == VolumetricMesh::TET)
    {
      for(int vtx=0; vtx<4; vtx++)
        vtxIndex[vtx] = volumetricMesh->getVertexIndex(el, vtx);
    }
    else if(type == VolumetricMesh::CUBIC)
    {
      for(int i = 0; i < 4; i++)
        vtxIndex[i] = volumetricMesh->getVertexIndex(el, cubeRotationTetIndex[i]);
    }

    /*
       Form matrix:
       M = [ v0   v1   v2   v3 ]
           [  1    1    1    1 ]
    */
    double M[16]; // row-major
    for(int vtx=0; vtx<4; vtx++)
      for(int dim=0; dim<3; dim++)
        M[4 * dim + vtx] = undeformedPositions[3 * vtxIndex[vtx] + dim];
    M[12] = M[13] = M[14] = M[15] = 1.0;

    // invert M and cache inverse (see [Mueller 2004])
    MInverse[el] = (double*) malloc (sizeof(double) * 16);
    inverse4x4(M, MInverse[el]);
  }

  // build acceleration indices for fast writing to the global stiffness matrix
  SparseMatrix * sparseMatrix = NULL;
  GetStiffnessMatrixTopology(&sparseMatrix);
  BuildRowColumnIndices(sparseMatrix);
  delete(sparseMatrix);

  // compute stiffness matrices for all the elements in the undeformed configuration
  KElementUndeformed = (double**) calloc (numElements, sizeof(double*));

  if(type == VolumetricMesh::TET)
  {
    for (int el = 0; el < numElements; el++)
    {
      double * MInv = MInverse[el];

      // Form stiffness matrix of the element in the undeformed configuration.
      // The procedure below is standard in FEM solid mechanics.
      // This code implements the equations given in Ahmed A. Shabana: Theory of Vibration, Volume II: Discrete and Continuous Systems, Springer--Verlag, New York, NY, 1990.
      // compute elasticity stiffness tensor
      double E[36]; // 6 x 6 matrix, stored row-major (symmetric, so row vs column-major storage makes no difference anyway)
      // theta = E * epsilon
      // epsilon = {e11, e22, e33, 2*e12, 2*e23, 2*e31}
      // theta = {t11, t22, t33, t12, t23, t31}
      if(computeElasticityStiffnessTensor(E, el) == false) // failure in constructing E
      {
        clear();  //clean data allocated so far
        throw 2;
      }

      double B[72] =
        { MInv[0], 0, 0,       MInv[4], 0, 0,       MInv[8], 0, 0,        MInv[12], 0, 0,
          0, MInv[1], 0,       0, MInv[5], 0,       0, MInv[9], 0,        0, MInv[13], 0,
          0, 0, MInv[2],       0, 0, MInv[6],       0, 0, MInv[10], 0,    0, MInv[14],
          MInv[1], MInv[0], 0, MInv[5], MInv[4], 0, MInv[9], MInv[8], 0,  MInv[13], MInv[12], 0,
          0, MInv[2], MInv[1], 0, MInv[6], MInv[5], 0, MInv[10], MInv[9], 0, MInv[14], MInv[13],
          MInv[2], 0, MInv[0], MInv[6], 0, MInv[4], MInv[10], 0, MInv[8], MInv[14], 0, MInv[12] };

      // EB = E * B
      double EB[72];
      memset(EB, 0, sizeof(double) * 72);
      for (int i=0; i<6; i++)
        for (int j=0; j<12; j++)
          for (int k=0; k<6; k++)
            EB[12 * i + j] += E[6 * i + k] * B[12 * k + j];

      // KElementUndeformed[el] = B^T * EB
      KElementUndeformed[el] = (double*) calloc (144, sizeof(double)); // element stiffness matrix
      for (int i=0; i<12; i++)
        for (int j=0; j<12; j++)
          for (int k=0; k<6; k++)
            KElementUndeformed[el][12 * i + j] += B[12 * k + i] * EB[12 * k + j];

      // KElementUndeformed[el] *= volume
      double volume = volumetricMesh->getElementVolume(el);

      for(int i=0; i<144; i++)
        KElementUndeformed[el][i] *= volume;
    }
  }
  else if(type == VolumetricMesh::CUBIC)
  {
    CubicMesh * cubicMesh = dynamic_cast<CubicMesh *>(volumetricMesh);
    double cubicSize = cubicMesh->getCubeSize();
    // w is the scale from nature to material coord.
    double w = cubicSize / 2;

    // create sampled B(s)
    double Bq[8][6*24];
    memset(Bq, 0, sizeof(Bq));
    assert(sizeof(Bq) == sizeof(double) * 8 * 6 * 24);

    const double sqrt3 = sqrt(3.);

    // the order of x[8],y[8],z[8] is the same as the vtx order in CUBIC:
    int x[8] = {-1, 1, 1,-1,-1, 1, 1,-1};
    int y[8] = {-1,-1, 1, 1,-1,-1, 1, 1};
    int z[8] = {-1,-1,-1,-1, 1, 1, 1, 1};
    for(int q = 0; q < 8; q++)
    {
      // s_q = ((-1)^i, (-1)^j, (-1)^k)^T / \sqrt(3), q = 4(i-1) + 2(j-1) + k, for i,j,k = {1,2}
      double s[3] = {x[q]/sqrt3, y[q]/sqrt3, z[q]/sqrt3};
      double s1 = s[0], s2 = s[1], s3 = s[2];
      double ms1 = 1-s1, ps1 = 1+s1, ms2 = 1-s2, ps2 = 1+s2, ms3 = 1-s3, ps3 = 1+s3;
      double s1c[2] = {ms1, ps1}, s2c[2] = {ms2, ps2}, s3c[2] = {ms3, ps3};
      // double N1x = -ms2*ms3, N1y = -ms1*ms3, N1z = -ms1*ms2;
      double dNds[8][3];
      for(int k = 0; k < 8; k++)
      {
        dNds[k][0] = x[k] * s2c[(y[k]+1)/2] * s3c[(z[k]+1)/2] / 8;
        dNds[k][1] = s1c[(x[k]+1)/2] * y[k] * s3c[(z[k]+1)/2] / 8;
        dNds[k][2] = s1c[(x[k]+1)/2] * s2c[(y[k]+1)/2] * z[k] / 8;
      }

      // n=1, i,j,k = 1,1,1, N_1 = (1-s1)(1-s2)(1-s3) / 8
      // n=2, i,j,k = 1,1,2, N_2 = (1-s1)(1-s2)(1+s3) / 8
      // ...

      // B(s) = [-(1-s2)(1-s3) 0 0  ...  (1+s2)(1+s3) 0 0] / 8
      //        [0 -(1-s1)(1-s3) 0  ... 0  (1+s1)(1+s3) 0]
      //        [0 0 -(1-s1)(1-s2)  ... 0 0  (1+s1)(1+s2)]
      //        [-(1-s1)(1-s3) -(1-s2)(1-s3)             ]
      //        [...                                     ]
      //        [...                                     ]

      // compute B(sq)
      for(int k = 0; k < 8; k++)
      {
        Bq[q][         k*3] = dNds[k][0];
        Bq[q][24   + 1+k*3] = dNds[k][1];
        Bq[q][24*2 + 2+k*3] = dNds[k][2];
        Bq[q][24*3 +   k*3] = dNds[k][1]; Bq[q][24*3 + 1+k*3] = dNds[k][0];
        Bq[q][24*4 + 1+k*3] = dNds[k][2]; Bq[q][24*4 + 2+k*3] = dNds[k][1];
        Bq[q][24*5 +   k*3] = dNds[k][2]; Bq[q][24*5 + 2+k*3] = dNds[k][0];
      }
    }

    for (int el = 0; el < numElements; el++)
    {
      // double * MInv = MInverse[el];
      double E[36]; // 6 x 6 matrix, stored row-major (symmetric, so row vs column-major storage makes no difference anyway)
      if(computeElasticityStiffnessTensor(E, el) == false) // failure in constructing E
      {
        clear();  //clean data allocated so far
        throw 1;
      }

      // Ke is of size 24 * 24
      KElementUndeformed[el] = (double*) calloc(576, sizeof(double));
      for(int q = 0; q < 8; q++)
      {
        double weight_q = 1.0;
        // Ke = \sum_q weight_q B(s_q)^T E B(s_q) w

        // EB = E * B
        double EB[6*24];
        memset(EB, 0, sizeof(double) * 6*24);
        for (int i = 0; i < 6; i++)
          for (int j = 0; j < 24; j++)
            for (int k = 0; k < 6; k++)
              EB[24 * i + j] += E[6 * i + k] * Bq[q][24 * k + j];

        // KElementUndeformed[el] += weight_q * B^T * EB
        for (int i = 0; i < 24; i++)
          for (int j = 0; j < 24; j++)
            for (int k = 0; k < 6; k++)
              KElementUndeformed[el][24 * i + j] += weight_q * Bq[q][24 * k + i] * EB[24 * k + j];
      }
      for(int i = 0; i < 576; i++)
        KElementUndeformed[el][i] *= w;

    }
  }
}

CorotationalLinearFEM::~CorotationalLinearFEM()
{
  clear();
}

void CorotationalLinearFEM::clear()
{
  free(undeformedPositions);
  undeformedPositions = NULL;
  for(int el=0; el < volumetricMesh->getNumElements(); el++)
  {
    free(KElementUndeformed[el]);
    free(MInverse[el]);
  }
  free(KElementUndeformed);
  KElementUndeformed = NULL;

  free(MInverse);
  MInverse = NULL;

  ClearRowColumnIndices();
}

// compute elasticity stiffness tensor
// E: 6 x 6 matrix, stored row-major (symmetric, so row vs column-major storage makes no difference anyway)
bool CorotationalLinearFEM::computeElasticityStiffnessTensor(double E[36], int el)
{
  VolumetricMesh::Material * material = volumetricMesh->getElementMaterial(el);

  // check if material is ENuMaterial (i.e., isotropic)
  VolumetricMesh::ENuMaterial * eNuMaterial = downcastENuMaterial(material);
  if (eNuMaterial != NULL)
  {
    // material is isotropic, specified by E, nu
    // compute Lame coefficients
    double lambda = eNuMaterial->getLambda();
    double mu = eNuMaterial->getMu();

    double Et[36] = { lambda + 2 * mu, lambda, lambda, 0, 0, 0,
        lambda, lambda + 2 * mu, lambda, 0, 0, 0,
        lambda, lambda, lambda + 2 * mu, 0, 0, 0,
        0, 0, 0, mu, 0, 0,
        0, 0, 0, 0, mu, 0,
        0, 0, 0, 0, 0, mu };

    memcpy(E, Et, sizeof(double) * 36);
  }
  else
  {
    // orthotropic material
    // we follow the following references:
    // Yijing Li and Jernej Barbic: Stable Orthotropic Materials, Symposium on Computer Animation 2014
    // http://en.wikipedia.org/wiki/Orthotropic_material
    // http://www.solidmechanics.org/text/Chapter3_2/Chapter3_2.htm

    // test if material is OrthotropicMaterial (i.e., orthotropic)
    VolumetricMesh::OrthotropicMaterial * orthotropicMaterial = downcastOrthotropicMaterial(material);
    if (orthotropicMaterial != NULL)
    {
      double E1 = orthotropicMaterial->getE1();
      double E2 = orthotropicMaterial->getE2();
      double E3 = orthotropicMaterial->getE3();
      double nu12 = orthotropicMaterial->getNu12();
      double nu23 = orthotropicMaterial->getNu23();
      double nu31 = orthotropicMaterial->getNu31();
      double G12 = orthotropicMaterial->getG12();
      double G23 = orthotropicMaterial->getG23();
      double G31 = orthotropicMaterial->getG31();

      double nu21 = nu12 * E2 / E1;
      double nu32 = nu23 * E3 / E2;
      double nu13 = nu31 * E1 / E3;

      double Y = 1.0 / (1.0 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 - 2.0 * nu21 * nu32 * nu13);

      double ELocal[36] = { E1 * (1.0 - nu23 * nu32) * Y, E1 * (nu21 + nu31 * nu23) * Y, E1 * (nu31 + nu21 * nu32) * Y, 0.0, 0.0, 0.0,
          E1 * (nu21 + nu31 * nu23) * Y, E2 * (1.0 - nu13 * nu31) * Y, E2 * (nu32 + nu12 * nu31) * Y, 0.0, 0.0, 0.0,
          E1 * (nu31 + nu21 * nu32) * Y, E2 * (nu32 + nu12 * nu31) * Y, E3 * (1.0 - nu12 * nu21) * Y, 0.0, 0.0, 0.0,
          0, 0, 0, G12, 0, 0,
          0, 0, 0, 0, G23, 0,
          0, 0, 0, 0, 0, G31 };

      //memcpy(E, ELocal, sizeof(double) * 36); // debug

      double R[9]; // row-major
      orthotropicMaterial->getR(R);

      // rotate Elocal into the basis given by the columns of R
      #define Relt(i,j) (R[3*(i)+(j)])
      double rotator[36];
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          rotator[6 * i + j] = Relt(i,j) * Relt(i,j);
          rotator[6 * i + 3 + j] = 2.0 * Relt(i, j) * Relt(i, (j+1) % 3);
          rotator[6 * (i + 3) + j] = Relt(i, j) * Relt((i+1) % 3, j);
          rotator[6 * (i + 3) + 3 + j] = Relt(i, j) * Relt((i+1) % 3, (j+1) % 3) + Relt(i, (j+1) % 3) * Relt((i+1) % 3, j);
        }
      #undef Relt

      // debug
      //memset(rotator, 0, sizeof(double) * 36);
      //for(int i=0; i<6; i++)
      //rotator[6*i+i] = 1.0;

      // E = rotator * ELocal * rotator^T
      double buffer[36];
      memset(buffer, 0, sizeof(double) * 36);
      // buffer = ELocal * rotator^T
      for(int i=0; i<6; i++)
        for(int j=0; j<6; j++)
          for(int k=0; k<6; k++)
            buffer[6 * i + j] += ELocal[6 * i + k] * rotator[6 * j + k];

      // E = rotator * buffer
      memset(E, 0, sizeof(double) * 36);
      for(int i=0; i<6; i++)
        for(int j=0; j<6; j++)
          for(int k=0; k<6; k++)
            E[6 * i + j] += rotator[6 * i + k] * buffer[6 * k + j];

    }
    else
    {
      printf("Error: CorotationalLinearFEM: unknown material encounterd in the mesh.\n");
      return false;
    }
  }

  return true;
}

void CorotationalLinearFEM::GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology)
{
  GetStiffnessMatrixTopology(volumetricMesh, stiffnessMatrixTopology);
}

void CorotationalLinearFEM::GetStiffnessMatrixTopology(VolumetricMesh * volumetricMesh, SparseMatrix ** stiffnessMatrixTopology)
{
  int numVertices = volumetricMesh->getNumVertices();
  SparseMatrixOutline * emptyMatrix = new SparseMatrixOutline(3 * numVertices);

  int numElements = volumetricMesh->getNumElements();
  int numElementVertices = volumetricMesh->getNumElementVertices();
  vector<int> vtxIndex(numElementVertices);
  for (int el = 0; el < numElements; el++)
  {
    for(int vtx = 0; vtx < numElementVertices; vtx++)
      vtxIndex[vtx] = volumetricMesh->getVertexIndex(el, vtx);

    for (int i = 0; i < numElementVertices; i++)
      for (int j = 0; j < numElementVertices; j++)
      {
        // add 3x3 block corresponding to pair of vertices (i,j)
        for(int k = 0; k < 3; k++)
          for(int l = 0; l < 3; l++)
            emptyMatrix->AddEntry(3 * vtxIndex[i] + k, 3 * vtxIndex[j] + l, 0.0);
      }
  }

  *stiffnessMatrixTopology = new SparseMatrix(emptyMatrix);
  delete(emptyMatrix);
}

// compute RK = R * K and RKRT = R * K * R^T (block-wise)
// input: K, R
// output: RK, RKRT
void CorotationalLinearFEM::WarpMatrix(double * K, double * R, double * RK, double * RKRT)
{
  int numElementVertices = volumetricMesh->getNumElementVertices();
  int Ksize = numElementVertices * 3;
  memset(RK, 0, sizeof(double) * Ksize * Ksize);
  memset(RKRT, 0, sizeof(double) * Ksize * Ksize);
  for(int i=0; i<numElementVertices; i++)
    for(int j=0; j<numElementVertices; j++)
    {
      // RK = R * K
      for(int k=0; k<3; k++)
         for(int l=0; l<3; l++)
           for(int m=0; m<3; m++)
             RK[Ksize * (3 * i + k) + (3 * j + l)] += R[3 * k + m] * K[Ksize * (3 * i + m) + (3 * j + l)];

      // RKRT = RK * R^T
      for(int k=0; k<3; k++)
        for(int l=0; l<3; l++)
          for(int m=0; m<3; m++)
            RKRT[Ksize * (3 * i + k) + (3 * j + l)] += RK[Ksize * (3 * i + k) + (3 * j + m)] * R[3 * l + m];
    }
}

void CorotationalLinearFEM::ComputeEnergyAndForceAndStiffnessMatrix(const double * u, double * energy, double * f, SparseMatrix * stiffnessMatrix, int warp)
{
  if (energy != NULL)
    *energy = 0.0;
  if (f != NULL) // clear f to zero
    memset(f, 0, sizeof(double) * 3 * numVertices);
  if (stiffnessMatrix != NULL) // clear stiffness matrix to zero
    stiffnessMatrix->ResetToZero();

  AddEnergyAndForceAndStiffnessMatrixOfSubmesh(u, energy, f, stiffnessMatrix, warp, 0, volumetricMesh->getNumElements());
}

void CorotationalLinearFEM::ComputeElementEnergyAndForceAndStiffnessMatrix(int el, const double * u, double * elementEnergy, 
    double * elementInternalForces, double * elementStiffnessMatrix, int warp)
{
  const int maxNumElementVertices = 8;
  const int maxNumElementDOFs = maxNumElementVertices * 3;

  const int numElementVertices = volumetricMesh->getNumElementVertices();
  const int numElementDOFs = numElementVertices * 3;
  const int elementStiffnessMatrixSpace = numElementDOFs * numElementDOFs;
  VolumetricMesh::elementType type = volumetricMesh->getElementType();

  // int vtxIndex[4];
  const int * vtxIndex = volumetricMesh->getVertexIndices(el);

  Vec3d deformedPos[maxNumElementVertices];
  for (int i = 0; i < numElementVertices; i++)
    deformedPos[i] = Vec3d(&undeformedPositions[3 * vtxIndex[i]]) + Vec3d(&u[3 * vtxIndex[i]]);

  if (warp > 0)
  {
    double P[16] ={0}; // the current world-coordinate positions (row-major)
    
    // P = [ v0   v1   v2   v3 ]
    //     [  1    1    1    1 ]
    // rows 1,2,3
    if(type == VolumetricMesh::TET)
    {
      for(int i=0; i<3; i++)
        for(int j=0; j<4; j++)
          P[4 * i + j] = deformedPos[j][i];
    }
    else if(type == VolumetricMesh::CUBIC)
    {
      for(int i=0; i<3; i++)
        for(int j=0; j<4; j++)
          P[4 * i + j] = deformedPos[cubeRotationTetIndex[j]][i];
    }
    else
      assert(0);
    // row 4
    for(int j=0; j<4; j++)
      P[12 + j] = 1.0;

    // F = P * Inverse(M)
    double F[9]; // upper-left 3x3 block
    for(int i=0; i<3; i++) 
      for(int j=0; j<3; j++) 
      {
        F[3 * i + j] = 0;
        for(int k=0; k<4; k++)
          F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
      }

    double R[9]; // rotation (row-major)
    double S[9]; // symmetric (row-major)
    double tolerance = 1E-6;
    int forceRotation = 1;
    PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

    // RK = R * K
    // KElement = R * K * R^T
    double RK[maxNumElementDOFs * maxNumElementDOFs]; // row-major
    if (elementStiffnessMatrix != NULL)
      WarpMatrix(KElementUndeformed[el], R, RK, elementStiffnessMatrix);

    double z[maxNumElementDOFs]; // z = RT x - x0
    Mat3d Rmat(R);
    Mat3d RTmat = trans(Rmat);
    for(int vtx = 0; vtx < numElementVertices; vtx++)
    {
      Vec3d x0 = Vec3d(&undeformedPositions[3 * vtxIndex[vtx]]);
      Vec3d vtxz = RTmat * deformedPos[vtx] - x0;
      vtxz.convertToArray(z + 3*vtx);
    }
    
    if (elementEnergy != NULL) // energy = 1/2 <Kz, z>, z = RT x - x0
    {
      double Kz[maxNumElementDOFs];
      memset(Kz, 0, sizeof(Kz));
      for(int i = 0; i < numElementDOFs; i++)
        for(int j = 0; j < numElementDOFs; j++)
          Kz[i] += KElementUndeformed[el][numElementDOFs*i+j] * z[j];

      double e = 0.0;
      for(int i = 0; i < numElementDOFs; i++)
        e += Kz[i] * z[i];
      *elementEnergy = 0.5 * e;
    }

    if (elementInternalForces != NULL) // f = RK (RT x - x0) = RK z
    {      
      memset(elementInternalForces, 0, sizeof(double) * numElementDOFs);
      for(int i = 0; i < numElementDOFs; i++)
        for(int j = 0; j < numElementDOFs; j++)
          elementInternalForces[i] += RK[numElementDOFs*i+j] * z[j];
    }

    // compute exact stiffness matrix
    if (warp == 2 && type == VolumetricMesh::TET && elementStiffnessMatrix != NULL)
    {
      // compute G = (tr(S) I - S) R^T
      double G[9]; 
      double tr = S[0] + S[4] + S[8];
      double temp[9];
      for(int i=0; i<9; i++)
        temp[i] = -S[i];
      temp[0] += tr;
      temp[4] += tr;
      temp[8] += tr;
      // G = temp * R^T
      MATRIX_MULTIPLY3X3ABT(temp, R, G);

      double invG[9]; // invG = G^{-1}
      inverse3x3(G, invG);

      double rhs[27]; // 3 x 9 matrix (column-major)
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
        {
          double temp[9];
          for(int k=0; k<9; k++)
            temp[k] = 0.0;
          // copy i-th row of R into column j of temp      
          for(int k=0; k<3; k++)
            temp[3 * k + j] = R[3 * i + k];
          // extract the skew-symmetric part
          SKEW_PART_NO_DIV2(temp, &rhs[3 * (3 * i + j)]);
        }

      // solve G * omega = rhs
      double omega[27]; // column-major
      for(int i=0; i<9; i++)
      {
        MATRIX_VECTOR_MULTIPLY3X3(invG, &rhs[3 * i], &omega[3 * i]);
      }

      double dRdF[81]; // each column is skew(omega) * R ; column-major
      for(int i=0; i<9; i++)
      {
        double skew[9];
        SKEW_MATRIX(&omega[3 * i], skew);
        MATRIX_MULTIPLY3X3(skew, R, &dRdF[9 * i]);
      }

      double B[3][3][9];
      // re-arrange dRdF into B, for easier dRdF * dFdx multiplication (to exploit sparsity of dFdx)
      for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
          for(int k=0; k<3; k++)
            for(int l=0; l<3; l++)
            {
              int row = 3 * i + k;
              int column = 3 * j + l;
              B[i][j][3 * k + l] = dRdF[9 * column + row];
            }

      // four pointers to a 3-vector
      double * minv[4] = { &MInverse[el][0], &MInverse[el][4], &MInverse[el][8], &MInverse[el][12] }; // the four rows of MInverse (last column ignored)

      double dRdx[108]; // derivative of the element rotation matrix with respect to the positions of the tet vertices; column-major
      for(int k=0; k<4; k++)
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
          {
            double temp[3];
            MATRIX_VECTOR_MULTIPLY3X3(B[i][j], minv[k], temp);
            int row = 3 * i;
            int column = 3 * k + j;
            VECTOR_SET3(&dRdx[9 * column + row], temp);
          }

      // add contribution of dRdx to KElement

      // term 1: \hat{dR/dxl} K (R^T x - m)

      // compute K (R^T x - m)
      double tempVec[12]; // R^T x - m
      for(int vtx=0; vtx<4; vtx++)
      {
        double pos[3];
        for(int i=0; i<3; i++)
          pos[i] = P[4 * i + vtx];
        MATRIX_VECTOR_MULTIPLY3X3T(R, pos, &tempVec[3*vtx]);
        // subtract m
        for(int i=0; i<3; i++)
          tempVec[3*vtx+i] -= undeformedPositions[3 * vtxIndex[vtx] + i];
      }
      double a[12]; // a = K * tempVec
      for (int i=0; i<12; i++)
      {
        a[i] = 0.0;
        for (int j=0; j<12; j++)
          a[i] += KElementUndeformed[el][12 * i + j] * tempVec[j];
      }

      // add [\hat{dR/dxl} K R^T x]_l, l=1 to 12
      for(int column=0; column<12; column++)
      {
        double b[12]; // b = \hat{dR/dxl} * a
        for(int j=0; j<4; j++)
        {
          MATRIX_VECTOR_MULTIPLY3X3(&dRdx[9 * column], &a[3*j], &b[3*j]);
        }
        // write b into KElement (add b to i-th column)
        for(int row=0; row<12; row++)
          elementStiffnessMatrix[12 * row + column] += b[row]; // KElement is row-major
      }

      // term 2: (R K \hat{dRdxl}^T)x

      // re-write positions into a
      for(int vtx=0; vtx<4; vtx++)
      {
        for(int i=0; i<3; i++)
          a[3 * vtx + i] = P[4 * i + vtx];
      }

      // compute [\hat{dRdxl}^T x)]_l, l=1 to 12
      for(int column=0; column<12; column++)
      {
        double b[12]; // b = \hat{dRdxl}^T * a
        for(int j=0; j<4; j++)
        {
          MATRIX_VECTOR_MULTIPLY3X3T(&dRdx[9 * column], &a[3*j], &b[3*j]);
        }

        // add RK * b to column of KElement
        int rowStart = 0;
        for (int row=0; row<12; row++)
        {
          double contrib = 0.0;
          for (int j=0; j<12; j++)
            contrib += RK[rowStart + j] * b[j];
          elementStiffnessMatrix[rowStart + column] += contrib;
          rowStart += 12;
        }
      }
    }
  }
  else // no warping, warp == 0
  {
    if (elementStiffnessMatrix != NULL)
      memcpy(elementStiffnessMatrix, KElementUndeformed[el], sizeof(double) * elementStiffnessMatrixSpace);
    // f = K u
    if (elementInternalForces != NULL || elementEnergy != NULL)
    {
      double fElementBuffer[maxNumElementDOFs];
      double * fEle = (elementInternalForces ? elementInternalForces : fElementBuffer);
      for(int i=0; i<numElementDOFs; i++)
      {
        fEle[i] = 0.0;
        for(int j=0; j<numElementVertices; j++)
        {
          fEle[i] += 
            KElementUndeformed[el][numElementDOFs * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
            KElementUndeformed[el][numElementDOFs * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
            KElementUndeformed[el][numElementDOFs * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
        }
      }
      if (elementEnergy != NULL)
      {
        double energy = 0.0;
        for(int j = 0; j < numElementVertices; j++)
          energy += fEle[3*j + 0] * u[3*vtxIndex[j] + 0] +
                    fEle[3*j + 1] * u[3*vtxIndex[j] + 1] +
                    fEle[3*j + 2] * u[3*vtxIndex[j] + 2];
        *elementEnergy = 0.5 * energy;
      }
    }
  }
}

void CorotationalLinearFEM::AddEnergyAndForceAndStiffnessMatrixOfSubmesh(const double * u, double * energy, double * f, SparseMatrix * stiffnessMatrix, int warp, int elementLo, int elementHi)
{
  const int maxNumElementVertices = 8;
  const int maxNumElementDOFs = maxNumElementVertices * 3;

  int numElementVertices = volumetricMesh->getNumElementVertices();

  int numElementDOFs = numElementVertices * 3;
  int elementStiffnessMatrixSpace = numElementDOFs * numElementDOFs;
  VolumetricMesh::elementType type = volumetricMesh->getElementType();
  for (int el=elementLo; el < elementHi; el++)
  {
    const int * vtxIndex = volumetricMesh->getVertexIndices(el);

    // element stiffness matrix, to be computed below; row-major
    double KElement[maxNumElementDOFs * maxNumElementDOFs];

    Vec3d deformedPos[maxNumElementVertices];
    for (int i = 0; i < numElementVertices; i++)
      deformedPos[i] = Vec3d(&undeformedPositions[3 * vtxIndex[i]]) + Vec3d(&u[3 * vtxIndex[i]]);

    if (warp > 0)
    {
      double P[16] ={0}; // the current world-coordinate positions (row-major)
      //
      // P = [ v0   v1   v2   v3 ]
      //     [  1    1    1    1 ]
      // rows 1,2,3
      if(type == VolumetricMesh::TET)
      {
        for(int i=0; i<3; i++)
          for(int j=0; j<4; j++)
            P[4 * i + j] = deformedPos[j][i];
      }
      else if(type == VolumetricMesh::CUBIC)
      {
        for(int i=0; i<3; i++)
          for(int j=0; j<4; j++)
            P[4 * i + j] = deformedPos[cubeRotationTetIndex[j]][i];
      }
      else
        assert(0);
      // row 4
      for(int j=0; j<4; j++)
        P[12 + j] = 1.0;

      // F = P * Inverse(M)
      double F[9]; // upper-left 3x3 block
      for(int i=0; i<3; i++) 
        for(int j=0; j<3; j++) 
        {
          F[3 * i + j] = 0;
          for(int k=0; k<4; k++)
            F[3 * i + j] += P[4 * i + k] * MInverse[el][4 * k + j];
        }

      double R[9]; // rotation (row-major)
      double S[9]; // symmetric (row-major)
      double tolerance = 1E-6;
      int forceRotation = 1;
      PolarDecomposition::Compute(F, R, S, tolerance, forceRotation);

      // RK = R * K
      // KElement = R * K * R^T
      double RK[maxNumElementDOFs * maxNumElementDOFs]; // row-major
      WarpMatrix(KElementUndeformed[el], R, RK, KElement);

      double z[maxNumElementDOFs]; // z = RT x - x0
      Mat3d Rmat(R);
      Mat3d RTmat = trans(Rmat);
      for(int vtx = 0; vtx < numElementVertices; vtx++)
      {
        Vec3d x0 = Vec3d(&undeformedPositions[3 * vtxIndex[vtx]]);
        Vec3d vtxz = RTmat * deformedPos[vtx] - x0;
        vtxz.convertToArray(z + 3*vtx);
      }
      
      if (energy != NULL) // energy = 1/2 <Kz, z>, z = RT x - x0
      {
        double Kz[maxNumElementDOFs];
        memset(Kz, 0, sizeof(Kz));
        for(int i = 0; i < numElementDOFs; i++)
          for(int j = 0; j < numElementDOFs; j++)
            Kz[i] += KElementUndeformed[el][numElementDOFs*i+j] * z[j];

        double e = 0.0;
        for(int i = 0; i < numElementDOFs; i++)
          e += Kz[i] * z[i];
        *energy += 0.5 * e;
      }

      double fElement[maxNumElementDOFs];
      if (f != NULL) // f = RK (RT x - x0) = RK z
      {      
        memset(fElement, 0, sizeof(fElement));
        for(int i = 0; i < numElementDOFs; i++)
          for(int j = 0; j < numElementDOFs; j++)
            fElement[i] += RK[numElementDOFs*i+j] * z[j];

        // write to global matrix
        for(int j=0; j<numElementVertices; j++)
          for(int l=0; l<3; l++)
            f[3 * vtxIndex[j] + l] += fElement[3 * j + l];
      }

      // compute exact stiffness matrix
      if (warp == 2 && type == VolumetricMesh::TET)
      {
        // compute G = (tr(S) I - S) R^T
        double G[9]; 
        double tr = S[0] + S[4] + S[8];
        double temp[9];
        for(int i=0; i<9; i++)
          temp[i] = -S[i];
        temp[0] += tr;
        temp[4] += tr;
        temp[8] += tr;
        // G = temp * R^T
        MATRIX_MULTIPLY3X3ABT(temp, R, G);

        double invG[9]; // invG = G^{-1}
        inverse3x3(G, invG);

        double rhs[27]; // 3 x 9 matrix (column-major)
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
          {
            double temp[9];
            for(int k=0; k<9; k++)
              temp[k] = 0.0;
            // copy i-th row of R into column j of temp      
            for(int k=0; k<3; k++)
              temp[3 * k + j] = R[3 * i + k];
            // extract the skew-symmetric part
            SKEW_PART_NO_DIV2(temp, &rhs[3 * (3 * i + j)]);
          }

        // solve G * omega = rhs
        double omega[27]; // column-major
        for(int i=0; i<9; i++)
        {
          MATRIX_VECTOR_MULTIPLY3X3(invG, &rhs[3 * i], &omega[3 * i]);
        }

        double dRdF[81]; // each column is skew(omega) * R ; column-major
        for(int i=0; i<9; i++)
        {
          double skew[9];
          SKEW_MATRIX(&omega[3 * i], skew);
          MATRIX_MULTIPLY3X3(skew, R, &dRdF[9 * i]);
        }

        double B[3][3][9];
        // re-arrange dRdF into B, for easier dRdF * dFdx multiplication (to exploit sparsity of dFdx)
        for(int i=0; i<3; i++)
          for(int j=0; j<3; j++)
            for(int k=0; k<3; k++)
              for(int l=0; l<3; l++)
              {
                int row = 3 * i + k;
                int column = 3 * j + l;
                B[i][j][3 * k + l] = dRdF[9 * column + row];
              }

        // four pointers to a 3-vector
        double * minv[4] = { &MInverse[el][0], &MInverse[el][4], &MInverse[el][8], &MInverse[el][12] }; // the four rows of MInverse (last column ignored)

        double dRdx[108]; // derivative of the element rotation matrix with respect to the positions of the tet vertices; column-major
        for(int k=0; k<4; k++)
          for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
            {
              double temp[3];
              MATRIX_VECTOR_MULTIPLY3X3(B[i][j], minv[k], temp);
              int row = 3 * i;
              int column = 3 * k + j;
              VECTOR_SET3(&dRdx[9 * column + row], temp);
            }

        // add contribution of dRdx to KElement

        // term 1: \hat{dR/dxl} K (R^T x - m)

        // compute K (R^T x - m)
        double tempVec[12]; // R^T x - m
        for(int vtx=0; vtx<4; vtx++)
        {
          double pos[3];
          for(int i=0; i<3; i++)
            pos[i] = P[4 * i + vtx];
          MATRIX_VECTOR_MULTIPLY3X3T(R, pos, &tempVec[3*vtx]);
          // subtract m
          for(int i=0; i<3; i++)
            tempVec[3*vtx+i] -= undeformedPositions[3 * vtxIndex[vtx] + i];
        }
        double a[12]; // a = K * tempVec
        for (int i=0; i<12; i++)
        {
          a[i] = 0.0;
          for (int j=0; j<12; j++)
            a[i] += KElementUndeformed[el][12 * i + j] * tempVec[j];
        }

        // add [\hat{dR/dxl} K R^T x]_l, l=1 to 12
        for(int column=0; column<12; column++)
        {
          double b[12]; // b = \hat{dR/dxl} * a
          for(int j=0; j<4; j++)
          {
            MATRIX_VECTOR_MULTIPLY3X3(&dRdx[9 * column], &a[3*j], &b[3*j]);
          }
          // write b into KElement (add b to i-th column)
          for(int row=0; row<12; row++)
            KElement[12 * row + column] += b[row]; // KElement is row-major
        }

        // term 2: (R K \hat{dRdxl}^T)x

        // re-write positions into a
        for(int vtx=0; vtx<4; vtx++)
        {
          for(int i=0; i<3; i++)
            a[3 * vtx + i] = P[4 * i + vtx];
        }

        // compute [\hat{dRdxl}^T x)]_l, l=1 to 12
        for(int column=0; column<12; column++)
        {
          double b[12]; // b = \hat{dRdxl}^T * a
          for(int j=0; j<4; j++)
          {
            MATRIX_VECTOR_MULTIPLY3X3T(&dRdx[9 * column], &a[3*j], &b[3*j]);
          }

          // add RK * b to column of KElement
          int rowStart = 0;
          for (int row=0; row<12; row++)
          {
            double contrib = 0.0;
            for (int j=0; j<12; j++)
              contrib += RK[rowStart + j] * b[j];
            KElement[rowStart + column] += contrib;
            rowStart += 12;
          }
        }
      }
    }
    else
    {
      // no warp
      memcpy(KElement, KElementUndeformed[el], sizeof(double) * elementStiffnessMatrixSpace);
      // f = K u
      double fElement[maxNumElementDOFs];
      if (f != NULL || energy != NULL)
      {
        for(int i=0; i<numElementDOFs; i++)
        {
          fElement[i] = 0.0;
          for(int j=0; j<numElementVertices; j++)
          {
            fElement[i] += 
              KElement[numElementDOFs * i + 3 * j + 0] * u[3 * vtxIndex[j] + 0] +
              KElement[numElementDOFs * i + 3 * j + 1] * u[3 * vtxIndex[j] + 1] +
              KElement[numElementDOFs * i + 3 * j + 2] * u[3 * vtxIndex[j] + 2];
          }
        }
        if (energy != NULL)
        {
          double elementEnergy = 0.0;
          for(int j = 0; j < numElementVertices; j++)
            elementEnergy += fElement[3*j + 0] * u[3*vtxIndex[j] + 0] +
                             fElement[3*j + 1] * u[3*vtxIndex[j] + 1] +
                             fElement[3*j + 2] * u[3*vtxIndex[j] + 2];
          *energy += 0.5 * elementEnergy;
        }

        // add fElement into the global f
        if (f != NULL)
          for(int j=0; j<numElementVertices; j++)
          {
            f[3 * vtxIndex[j] + 0] += fElement[3 * j + 0];
            f[3 * vtxIndex[j] + 1] += fElement[3 * j + 1];
            f[3 * vtxIndex[j] + 2] += fElement[3 * j + 2];
          }
      }
    }

    if (stiffnessMatrix != NULL)
    {
      int * rowIndex = rowIndices[el];
      int * columnIndex = columnIndices[el];

      // add KElement to the global stiffness matrix
      for (int i=0; i<numElementVertices; i++)
        for (int j=0; j<numElementVertices; j++)
          for(int k=0; k<3; k++)
            for(int l=0; l<3; l++)
              stiffnessMatrix->AddEntry(3 * rowIndex[i] + k, 3 * columnIndex[numElementVertices * i + j] + l, KElement[numElementDOFs * (3 * i + k) + 3 * j + l]);
    }
  }
}

void CorotationalLinearFEM::ClearRowColumnIndices()
{
  for (int el=0; el < volumetricMesh->getNumElements(); el++)
  {
    free(rowIndices[el]);
    free(columnIndices[el]);
  }

  free(rowIndices);
  free(columnIndices);
  rowIndices = NULL;
  columnIndices = NULL;
}

void CorotationalLinearFEM::BuildRowColumnIndices(SparseMatrix * sparseMatrix)
{
  int numElements = volumetricMesh->getNumElements();
  int numElementVertices = volumetricMesh->getNumElementVertices();
  rowIndices = (int**) malloc (sizeof(int*) * numElements);
  columnIndices = (int**) malloc (sizeof(int*) * numElements);

  for (int el=0; el < numElements; el++)
  {
    // the rows corresponding to each vertices in the element
    rowIndices[el] = (int*) malloc (sizeof(int) * numElementVertices);
    for(int i=0; i<numElementVertices; i++)
      rowIndices[el][i] = volumetricMesh->getVertexIndex(el, i);

    // the columns corresponding to all element vertices, in row of each vertex
    columnIndices[el] = (int*) malloc (sizeof(int) * numElementVertices * numElementVertices);
    // find index of vertex j in row of vertex i, and cache it
    for(int i=0; i<numElementVertices; i++)
      for(int j=0; j<numElementVertices; j++)
        columnIndices[el][numElementVertices * i + j] = sparseMatrix->GetInverseIndex(3*rowIndices[el][i], 3*rowIndices[el][j]) / 3;
  }
}

// inverse of a 3x3 matrix
// row-major format
void CorotationalLinearFEM::inverse3x3(double * A, double * AInv)
{
  // converted to C from Mathematica output   
  AInv[0] = -A[5] * A[7] + A[4] * A[8]; 
  AInv[1] = A[2] * A[7] - A[1] * A[8]; 
  AInv[2] = -A[2] * A[4] + A[1] * A[5];
  AInv[3] = A[5] * A[6] - A[3] * A[8]; 
  AInv[4] = -A[2] * A[6] + A[0] * A[8]; 
  AInv[5] = A[2] * A[3] - A[0] * A[5];
  AInv[6] = -A[4] * A[6] + A[3] * A[7]; 
  AInv[7] = A[1] * A[6] - A[0] * A[7];
  AInv[8] = -A[1] * A[3] + A[0] * A[4];

  double invDet = 1.0 / (A[0] * AInv[0] + A[1] * AInv[3] + A[2] * AInv[6]);

  for(int i=0; i<9; i++)
    AInv[i] *= invDet;
}

// inverse of a 4x4 matrix
// row-major format
void CorotationalLinearFEM::inverse4x4(double * A, double * AInv)
{
  // converted to C from Mathematica output   
  AInv[0] = -A[11] * A[14] * A[5] + A[10] * A[15] * A[5] + A[11] * A[13] * A[6] - A[10] * A[13] * A[7] - A[15] * A[6] * A[9] + A[14] * A[7] * A[9];
  AInv[1] = A[1] * A[11] * A[14] - A[1] * A[10] * A[15] - A[11] * A[13] * A[2] + A[10] * A[13] * A[3] + A[15] * A[2] * A[9] - A[14] * A[3] * A[9];
  AInv[2] = -A[15] * A[2] * A[5] + A[14] * A[3] * A[5] + A[1] * A[15] * A[6] - A[13] * A[3] * A[6] - A[1] * A[14] * A[7] + A[13] * A[2] * A[7];
  AInv[3] = A[11] * A[2] * A[5] - A[10] * A[3] * A[5] - A[1] * A[11] * A[6] + A[1] * A[10] * A[7] + A[3] * A[6] * A[9] - A[2] * A[7] * A[9];
  AInv[4] = A[11] * A[14] * A[4] - A[10] * A[15] * A[4] - A[11] * A[12] * A[6] + A[10] * A[12] * A[7] + A[15] * A[6] * A[8] - A[14] * A[7] * A[8];
  AInv[5] = -A[0] * A[11] * A[14] + A[0] * A[10] * A[15] + A[11] * A[12] * A[2] - A[10] * A[12] * A[3] - A[15] * A[2] * A[8] + A[14] * A[3] * A[8];
  AInv[6] = A[15] * A[2] * A[4] - A[14] * A[3] * A[4] - A[0] * A[15] * A[6] + A[12] * A[3] * A[6] + A[0] * A[14] * A[7] - A[12] * A[2] * A[7];
  AInv[7] = -A[11] * A[2] * A[4] + A[10] * A[3] * A[4] + A[0] * A[11] * A[6] - A[0] * A[10] * A[7] - A[3] * A[6] * A[8] + A[2] * A[7] * A[8];
  AInv[8] = -A[11] * A[13] * A[4] + A[11] * A[12] * A[5] - A[15] * A[5] * A[8] + A[13] * A[7] * A[8] + A[15] * A[4] * A[9] - A[12] * A[7] * A[9];
  AInv[9] = -A[1] * A[11] * A[12] + A[0] * A[11] * A[13] + A[1] * A[15] * A[8] - A[13] * A[3] * A[8] - A[0] * A[15] * A[9] + A[12] * A[3] * A[9];
  AInv[10] = -A[1] * A[15] * A[4] + A[13] * A[3] * A[4] + A[0] * A[15] * A[5] - A[12] * A[3] * A[5] + A[1] * A[12] * A[7] - A[0] * A[13] * A[7];
  AInv[11] = A[1] * A[11] * A[4] - A[0] * A[11] * A[5] + A[3] * A[5] * A[8] - A[1] * A[7] * A[8] - A[3] * A[4] * A[9] + A[0] * A[7] * A[9]; 
  AInv[12] = A[10] * A[13] * A[4] - A[10] * A[12] * A[5] + A[14] * A[5] * A[8] - A[13] * A[6] * A[8] - A[14] * A[4] * A[9] + A[12] * A[6] * A[9];
  AInv[13] = A[1] * A[10] * A[12] - A[0] * A[10] * A[13] - A[1] * A[14] * A[8] + A[13] * A[2] * A[8] + A[0] * A[14] * A[9] - A[12] * A[2] * A[9]; 
  AInv[14] = A[1] * A[14] * A[4] - A[13] * A[2] * A[4] - A[0] * A[14] * A[5] + A[12] * A[2] * A[5] - A[1] * A[12] * A[6] + A[0] * A[13] * A[6];
  AInv[15] = -A[1] * A[10] * A[4] + A[0] * A[10] * A[5] - A[2] * A[5] * A[8] + A[1] * A[6] * A[8] + A[2] * A[4] *A[9] - A[0] * A[6] * A[9];

  double invDet = 1.0 / (A[0] * AInv[0] + A[1] * AInv[4] + A[2] * AInv[8] + A[3] * AInv[12]);

  for(int i=0; i<16; i++)
    AInv[i] *= invDet;
}

