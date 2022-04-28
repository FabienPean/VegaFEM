/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "finite difference tester" utility , Copyright (C) 2018 USC           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Yijing Li                                                *
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
  This driver utility tests the correctness of Vega's material classes.
  It tests internal forces by performing finite differences on the energies,
  and it tests stiffness matrices by performing finite differences on the
  internal forces.
*/
#include "finiteDifferenceTester.h"

#include "forceModel.h"
#include "configFile.h"

#include "tetMesh.h"
#include "volumetricMesh.h"
#include "volumetricMeshLoader.h"
#include "generateMassMatrix.h"

#include "clothBWFromObjMesh.h"
#include "clothBWStencilForceModel.h"

#include "matrix.h"
#include "vectorHelper.h"

#include "neoHookeanIsotropicMaterial.h"
#include "StVKIsotropicMaterial.h"
#include "MooneyRivlinIsotropicMaterial.h"
#include "isotropicHyperelasticFEM.h"
#include "isotropicHyperelasticFEMStencilForceModel.h"

#include "corotationalLinearFEM.h"
#include "corotationalLinearFEMStencilForceModel.h"

#include "StVKFEM.h"
#include "StVKElementABCDLoader.h"
#include "StVKStencilForceModel.h"

#include "massSpringSystem.h"
#include "massSpringSystemFromObjMeshConfigFile.h"
#include "massSpringSystemFromTetMeshConfigFile.h"
#include "massSpringSystemFromCubicMeshConfigFile.h"
#include "massSpringStencilForceModel.h"

#include "stencilForceModel.h"
#include "forceModelAssembler.h"

#include "linearFEMStencilForceModel.h"

#include <iostream>
#include <cassert>
#include <numeric>

using namespace std;

#define ADD_CONFIG(v) config.addOptionOptional(#v, &v, v)

static double h = 0.0001;
int n = 0, r = 0; // n: #vtx, r = n * 3

static string volumetricMeshFilename;
static string clothMeshFilename;
static string animationFilename;

static string deformableObjectMethod;
static string invertibleMaterial;
static int numInternalForceThreads = 1;

static double surfaceDensity = 100;
static double tensileStiffness = 4e6;
static double shearStiffness = 2e6;
static double bendStiffnessU = 0.2;
static double bendStiffnessV = 0.2;
static bool computeStShForce = true;
static bool computeBendForce = true;
static bool computeStShStiffness = true;
static bool computeBendStiffness = true;
static bool useRestAngle = true;
static string customMassSpringSystem;
static string massSpringSystemObjConfigFilename;
static string massSpringSystemTetMeshConfigFilename;
static string massSpringSystemCubicMeshConfigFilename;

static int enableCompressionResistance = 0.;
static double compressionResistance = 500;
static double inversionThreshold = -1;

static int corotationalLinearFEM_warp = 1;

static ObjMesh * mesh = nullptr;
static ClothBW * cloth = nullptr;
static VolumetricMesh * volumetricMesh = nullptr;

static ForceModel * forceModel = nullptr;
static IsotropicMaterial * isotropicMaterial = nullptr;
static IsotropicHyperelasticFEM * isotropicHyperelasticFEM = nullptr;
CorotationalLinearFEM * corotationalLinearFEM = nullptr;
StVKElementABCD * stVKPrecomputedIntegrals = nullptr;
MassSpringSystem * massSpringSystem = nullptr;
StVKFEM * stvkFEM = nullptr;

static StencilForceModel *stencilForceModel = nullptr;
static ForceModelAssembler *fmAsm = nullptr;

static int numTestThreads = 1;
static bool useFivePoint = true;
static bool testStiffness = false;
static bool hApproachZero = false;

static Matrix<double> * animation = nullptr;
static int numFrames = 0;

static FiniteDifferenceTester * tester = nullptr;

static void initConfig(const char * configFilename)
{
  ConfigFile config;

  ADD_CONFIG(deformableObjectMethod);
  ADD_CONFIG(invertibleMaterial);
  ADD_CONFIG(volumetricMeshFilename);
  ADD_CONFIG(clothMeshFilename);
  ADD_CONFIG(surfaceDensity);
  ADD_CONFIG(tensileStiffness);
  ADD_CONFIG(shearStiffness);
  ADD_CONFIG(bendStiffnessU);
  ADD_CONFIG(bendStiffnessV);
  ADD_CONFIG(computeStShForce);
  ADD_CONFIG(computeBendForce);
  ADD_CONFIG(computeStShStiffness);
  ADD_CONFIG(computeBendStiffness);
  ADD_CONFIG(useRestAngle);
  ADD_CONFIG(animationFilename);
  ADD_CONFIG(h);
  ADD_CONFIG(numTestThreads);
  ADD_CONFIG(numInternalForceThreads);
  ADD_CONFIG(useFivePoint);
  ADD_CONFIG(enableCompressionResistance);
  ADD_CONFIG(compressionResistance);
  ADD_CONFIG(inversionThreshold);
  ADD_CONFIG(corotationalLinearFEM_warp);
  ADD_CONFIG(testStiffness);
  ADD_CONFIG(hApproachZero);
  ADD_CONFIG(customMassSpringSystem);
  ADD_CONFIG(massSpringSystemObjConfigFilename);
  ADD_CONFIG(massSpringSystemTetMeshConfigFilename);
  ADD_CONFIG(massSpringSystemCubicMeshConfigFilename);

  if (config.parseOptions(configFilename) != 0)
  {
    cout << "Error parsing " << configFilename << endl;
    exit(1);
  }
  config.printOptions();

}

static void initForceModel()
{

  if (deformableObjectMethod == "ClothBW")
  {
    assert(clothMeshFilename.size() > 0);
    mesh = new ObjMesh(clothMeshFilename);
    n = mesh->getNumVertices();
    r = n * 3;
    ClothBW::MaterialGroup group;
    group.tensileStiffness = tensileStiffness;
    group.shearStiffness = shearStiffness;
    group.bendStiffnessU = bendStiffnessU;
    group.bendStiffnessV = bendStiffnessV;
    cloth = ClothBWFromObjMesh::GenerateClothBW(mesh, surfaceDensity, group, 0);
    bool parameters[4];
    parameters[0] = computeStShForce;
    parameters[1] = computeBendForce;
    parameters[2] = computeStShStiffness;
    parameters[3] = computeBendStiffness;
    cloth->SetComputationMode(parameters);
    cloth->UseRestAnglesForBendingForces(useRestAngle);

    stencilForceModel = new ClothBWStencilForceModel(cloth);
  }
  else if (deformableObjectMethod == "InvertibleFEM")
  {
    TetMesh * tetMesh = new TetMesh(volumetricMeshFilename.c_str());
    n = tetMesh->getNumVertices();
    r = 3 * n;
    volumetricMesh = tetMesh;
    if (invertibleMaterial == "StVK")
      isotropicMaterial = new StVKIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
    else if (invertibleMaterial == "neoHookean")
      isotropicMaterial = new NeoHookeanIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
    else if (invertibleMaterial == "MooneyRivlin")
      isotropicMaterial = new MooneyRivlinIsotropicMaterial(tetMesh, enableCompressionResistance, compressionResistance);
    assert(isotropicMaterial);

    // create the invertible FEM deformable model
    isotropicHyperelasticFEM = new IsotropicHyperelasticFEM(tetMesh, isotropicMaterial, inversionThreshold);
    stencilForceModel = new IsotropicHyperelasticFEMStencilForceModel(isotropicHyperelasticFEM);
  }
  else if (deformableObjectMethod == "CLFEM")
  {
    TetMesh * tetMesh = new TetMesh(volumetricMeshFilename.c_str());
    n = tetMesh->getNumVertices();
    r = 3 * n;
    volumetricMesh = tetMesh;
    corotationalLinearFEM = new CorotationalLinearFEM(tetMesh);
    stencilForceModel = new CorotationalLinearFEMStencilForceModel(corotationalLinearFEM);
  }
  else if (deformableObjectMethod == "StVK" || deformableObjectMethod == "LinearFEM")
  {
    volumetricMesh = VolumetricMeshLoader::load(volumetricMeshFilename.c_str());
    assert(volumetricMesh);
    n = volumetricMesh->getNumVertices();
    r = 3 * n;

    unsigned int loadingFlag = 0; // 0 = use low-memory version, 1 = use high-memory version
    stVKPrecomputedIntegrals = StVKElementABCDLoader::load(volumetricMesh, loadingFlag);
    assert(stVKPrecomputedIntegrals);

    stvkFEM = new StVKFEM(volumetricMesh, stVKPrecomputedIntegrals);

    if (deformableObjectMethod == "StVK")
      stencilForceModel = new StVKStencilForceModel(stvkFEM);
    else if (deformableObjectMethod == "LinearFEM")
      stencilForceModel = new LinearFEMStencilForceModel(new StVKStencilForceModel(stvkFEM));
  }
  else if (deformableObjectMethod == "MassSpring")
  {
    if (massSpringSystemObjConfigFilename.size() > 0)
    {
      printf("Loading mass spring system from an obj file...\n");
      MassSpringSystemFromObjMeshConfigFile massSpringSystemFromObjMeshConfigFile;
      MassSpringSystemObjMeshConfiguration massSpringSystemObjMeshConfiguration;
      if (massSpringSystemFromObjMeshConfigFile.GenerateMassSpringSystem(massSpringSystemObjConfigFilename.c_str(), &massSpringSystem, &massSpringSystemObjMeshConfiguration) != 0)
      {
        printf("Error initializing the mass spring system.\n");
        exit(1);
      }
    }
    else if (massSpringSystemTetMeshConfigFilename.size() > 0)
    {
      printf("Loading mass spring system from a tet mesh file...\n");
      MassSpringSystemFromTetMeshConfigFile massSpringSystemFromTetMeshConfigFile;
      MassSpringSystemTetMeshConfiguration massSpringSystemTetMeshConfiguration;
      if (massSpringSystemFromTetMeshConfigFile.GenerateMassSpringSystem(massSpringSystemTetMeshConfigFilename.c_str(), &massSpringSystem, &massSpringSystemTetMeshConfiguration) != 0)
      {
        printf("Error initializing the mass spring system.\n");
        exit(1);
      }
    }
    else if (massSpringSystemCubicMeshConfigFilename.size() > 0)
    {
      printf("Loading mass spring system from a cubic mesh file...\n");
      MassSpringSystemFromCubicMeshConfigFile massSpringSystemFromCubicMeshConfigFile;
      MassSpringSystemCubicMeshConfiguration massSpringSystemCubicMeshConfiguration;
      if (massSpringSystemFromCubicMeshConfigFile.GenerateMassSpringSystem(massSpringSystemCubicMeshConfigFilename.c_str(), &massSpringSystem, &massSpringSystemCubicMeshConfiguration) != 0)
      {
        printf("Error initializing the mass spring system.\n");
        exit(1);
      }
    }
    else if (customMassSpringSystem.size() > 0)
    {
      int numParticles;
      double groupStiffness;
      sscanf(customMassSpringSystem.c_str(), "chain,%d,%lf", &numParticles, &groupStiffness);
      printf("Creating a chain mass-spring system with %d particles...\n", numParticles);

      double * masses = (double*) malloc (sizeof(double) * numParticles);
      for(int i=0; i<numParticles; i++)
        masses[i] = 1.0;

      double * restPositions = (double*) malloc (sizeof(double) * 3 * numParticles);
      for(int i=0; i<numParticles; i++)
      {
        restPositions[3*i+0] = 0;
        restPositions[3*i+1] = (numParticles == 1) ? 0.0 : 1.0 * i / (numParticles-1);
        restPositions[3*i+2] = 0;
      }
      int * edges = (int*) malloc (sizeof(int) * 2 * (numParticles - 1));
      for(int i=0; i<numParticles-1; i++)
      {
        edges[2*i+0] = i;
        edges[2*i+1] = i+1;
      }

      int * edgeGroups = (int*) malloc (sizeof(int) * (numParticles - 1));
      for(int i=0; i<numParticles-1; i++)
        edgeGroups[i] = 0;
      double groupDamping = 0;

      massSpringSystem = new MassSpringSystem(numParticles, masses, restPositions, numParticles - 1, edges, edgeGroups, 1, &groupStiffness, &groupDamping, false);

      free(edgeGroups);
      free(edges);
      free(restPositions);
      free(masses);
    }
    assert(massSpringSystem);

    n = massSpringSystem->GetNumParticles();
    r = n * 3;

    stencilForceModel = new MassSpringStencilForceModel(massSpringSystem);
  }

  assert(stencilForceModel != nullptr);
  fmAsm = new ForceModelAssembler(stencilForceModel);

  forceModel = fmAsm;
}

static void testAnimation(double & maxForceRelErr, double & maxStiffnessRelErr)
{
  maxStiffnessRelErr = 0;
//  int maxRelErrFrame = 0;

  for(int i = 0; i < numFrames; i++)
  {
    double * u = animation->GetData() + r * i;

    double relForceErr = tester->testInternalForce(u);
    if (maxForceRelErr < relForceErr)
    {
      maxForceRelErr = relForceErr;
//      maxRelErrFrame = i;
    }

    if (hApproachZero == false)
      cout << "Frame " << i << " force relerr: " << relForceErr * 100 << "%";

    if (testStiffness)
    {
      double unsymErr = 0.0;
      double relStiffErr = tester->testStiffnessMatrix(u, &unsymErr);
      if (maxStiffnessRelErr < relStiffErr)
        maxStiffnessRelErr = relStiffErr;

      if (hApproachZero == false)
        cout << ", stiffness relerr: " << relStiffErr * 100 << "% unsymErr: " << unsymErr * 100 << "%";
    } // end if testStiffness
    if (hApproachZero == false)
      cout << endl;
  }
}


int main(int argc, char ** argv)
{
  int numFixedArgs = 2;
  if (argc < numFixedArgs)
  {
    cout << "Usage: " << argv[0] << ": <config filename>" << endl;
    return 0;
  }

  initConfig(argv[1]);

  initForceModel();

  animation = new Matrix<double>(animationFilename.c_str());
  assert(animation->Getm() == r);
  numFrames = animation->Getn();
  assert(numFrames > 0);
  assert(numTestThreads >= 1);

  vector<double> approachForceRelErr, approachStiffnessRelErr;
  double maxForceRelErr = 0;
  double maxStiffnessRelErr= 0;

  FiniteDifferenceTester::Mode mode = (useFivePoint ? FiniteDifferenceTester::FIVE_POINT : FiniteDifferenceTester::TWO_POINT);
  tester = new FiniteDifferenceTester(forceModel, h, mode, numTestThreads);

  if (hApproachZero)
    for(h = 1; h >= 1e-11; h = h / 10.)
    {
      tester->setTimestep(h);
      maxForceRelErr = 0;
      maxStiffnessRelErr = 0;
      testAnimation(maxForceRelErr, maxStiffnessRelErr);

      if (hApproachZero)
      {
        cout << "h = " << h << endl;
        approachForceRelErr.push_back(maxForceRelErr);
        approachStiffnessRelErr.push_back(maxStiffnessRelErr);
      }
      cout << "max force rel error " << maxForceRelErr * 100 << "% ";
      if (testStiffness)
      {
        cout << "max stiffness rel error " << maxStiffnessRelErr * 100 << "%" << endl;
      }
    }
  else
    testAnimation(maxForceRelErr, maxStiffnessRelErr);

  if (approachForceRelErr.size() > 0)
  {
    cout << "======================================" << endl;
    for(size_t i = 0; i < approachForceRelErr.size(); i++)
      cout << approachForceRelErr[i] * 100 << "%\t";
    cout << endl;
    for(size_t i = 0; i < approachStiffnessRelErr.size(); i++)
      cout << approachStiffnessRelErr[i] * 100 << "%\t";
    cout << endl;
  }
  return 0;
}

