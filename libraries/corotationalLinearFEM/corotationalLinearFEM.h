/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "corotational linear FEM" library , Copyright (C) 2018 USC            *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Jernej Barbic, Yijing Li                                *
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

#ifndef _COROTATIONALLINEARFEM_H_
#define _COROTATIONALLINEARFEM_H_

/*
  Corotational linear FEM deformable model for tet meshes and cubic meshes.

  1. For tet meshes, this class implements the deformable model described in the following paper:

    M. Mueller, M. Gross: Interactive Virtual Materials.
    In Proc. of Graphics Interface 2004 (2004), pp. 239–246.

  In [Mueller 2004], the tangent stiffness matrix is approximate (warp=1). 
  This class can also compute the exact tangent stiffness matrix (warp=2).
  The implementation is described in:
  J. Barbic: Exact Corotational Linear FEM Stiffness Matrix, Technical Report, USC, 2012

  It is also possible to turn warping off (warp=0). This gives fast linear FEM dynamics,
  but large deformations are not well-represented.

  2. For cubic meshes, the class implements the following paper:

    Jesus Perez, Alvaro G. Perez and Miguel A. Otaduy:
    Simulation of Hyperelastic Materials Using Energy Constraints
    Proc. of Congreso Espanol de Informatica Grafica, 2013

  Note that for cubic meshes, only the warp=1 and warp=0 modes are supported.
  For cubic meshes with warp=2, the code behaves the same as for warp=1.

  3. This class supports both isotropic and orthotropic materials. By default, isotropic material is used.
  The choice of material is controlled by the material specification in the .veg file.
  For orthotropic materials, we implement the following:

     Yijing Li and Jernej Barbic: 
     Stable Orthotropic Materials, 
     Symposium on Computer Animation 2014.
     
  See also:
     http://en.wikipedia.org/wiki/Orthotropic_material
     http://www.solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
*/

#include "tetMesh.h"
#include "sparseMatrix.h"

class CorotationalLinearFEM
{
public:

  // initializes corotational linear FEM
  // input: tetMesh and cubicMesh
  CorotationalLinearFEM(VolumetricMesh * volumetricMesh);
  virtual ~CorotationalLinearFEM();

  void GetStiffnessMatrixTopology(SparseMatrix ** stiffnessMatrixTopology); // returns a zero matrix containing the locations of non-zero elements in the stiffness matrix

  // computes elastic energy, the internal forces and (warped) stiffness matrix for the entire mesh
  // vertex displacements (input) and internal forces (output) must be (pre-allocated) vectors of length 3 * numVertices
  // the internal forces are returned with the sign corresponding to f_int(x) on the left side of the equation M * x'' + f_int(x) = f_ext
  // i.e., the computed internal forces are *negatives* of the actual physical internal forces acting on the material
  // warp:
  //   0: no warping (linear FEM)
  //   1: stiffness warping (corotational linear FEM with approximate stiffness matrix) [Mueller 2004]
  //   2: corotational linear FEM with exact tangent stiffness matrix (see the technical report [Barbic 2012]), only works on tetMesh
  // if you do not want to compute the energy, internal forces or stiffness matrix (any combination), pass a NULL pointer for that argument
  virtual void ComputeEnergyAndForceAndStiffnessMatrix(const double * vertexDisplacements, double * energy, double * internalForces, SparseMatrix * stiffnessMatrix, int warp=1);

  // this routine is same as above, except that (1) it adds to the existing value, (2) only traverses elements from elementLo <= element <= elementHi - 1
  // if you do not want to compute the energy, internal forces or stiffness matrix (any combination), pass a NULL pointer for that argument
  void AddEnergyAndForceAndStiffnessMatrixOfSubmesh(const double * vertexDisplacements, double * energy, double * internalForces, SparseMatrix * stiffnessMatrix, int warp, int elementLo, int elementHi);

  // this routine is same as above, except that (1) it only computes a single element. (2) the output forces and the matrix 
  // are stored in dense format.
  void ComputeElementEnergyAndForceAndStiffnessMatrix(int elementID, const double * vertexDisplacements, double * elementEnergy, 
    double * elementInternalForces, double * elementStiffnessMatrix, int warp);

  inline VolumetricMesh * GetVolumetricMesh() { return volumetricMesh; }

  static void inverse3x3(double * A, double * AInv); // inverse of a row-major 3x3 matrix
  static void inverse4x4(double * A, double * AInv); // inverse of a row-major 4x4 matrix

protected:
  int numVertices;
  VolumetricMesh * volumetricMesh;
  double * undeformedPositions;
  double ** MInverse;
  double ** KElementUndeformed;

  void WarpMatrix(double * K, double * R, double * RK, double * RKRT);

  // acceleration indices
  int ** rowIndices;
  int ** columnIndices;
  void ClearRowColumnIndices();
  void BuildRowColumnIndices(SparseMatrix * sparseMatrix);

  bool computeElasticityStiffnessTensor(double E[36], int el);
  // build inverse of M = [ v0   v1   v2   v3 ]
  //                      [  1    1    1    1 ]
  // volumetricMesh must be TET of CUBIC. If it's CUBIC, M consists of vertices forming a tet in the center of the cube.
  static void buildMInverse(double * MInverse, VolumetricMesh * volumetricMesh);

  void clear();

  static void GetStiffnessMatrixTopology(VolumetricMesh * mesh, SparseMatrix ** stiffnessMatrixTopology);
};

#endif

