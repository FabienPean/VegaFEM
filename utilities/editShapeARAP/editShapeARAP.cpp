/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "Edit shape with ARAP" driver application,                            *
 *  Copyright (C) 2007 CMU, 2009 MIT, 2018 USC                           *
 *                                                                       *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code authors: Yijing Li, Hongyi Xu, Koki Nagano, Jernej Barbic        *
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
 * This utility is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This utility is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

/*  
   An OpenGL driver to our "shapeEdit" ARAP library.
   Allows the user to interactively deform a 3d model using the ARAP energy.
*/ 

#include "getopts.h"
#include "sceneObjectDeformable.h"
#include "objMesh.h"
#include "performanceCounter.h"
#include "configFile.h"
#include "volumetricMeshLoader.h"
#include <GL/glui.h>
#include "lighting.h"
#include "valueIndex.h"
#include "listIO.h"
#include "openGLHelper.h"
#include "matrix.h"
#include "arapDeformer.h"
#include "matrixMultiplyMacros.h"
#include "constrainedDOFs.h"
#include "matrixIO.h"
#ifdef USE_GLSL
  #include "glslPhong.h"
#endif
#include "saveScreenShot.h"
#include "camera.h"
#include "matrix.h"
#include "handleControl.h"
#include "tetMesh.h"
#include "generateSurfaceMesh.h"
#include "handleRender.h"
#include "sceneGroundPlane.h"
#include "averagingBuffer.h"
#include "barycentricCoordinates.h"
#include "generateTetMeshFromCubicMesh.h"
#include "simulationRecorder.h"
#include "stringHelper.h"
#include "cameraLighting.h"
#include "inputDevice.h"
#include "cameraChangeLoad.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdio>
#include <algorithm>
#include <cassert>
using namespace std;

// graphics
static char windowTitleBase[4096] = "DM";
static int windowID;
static int windowWidth = 1024;
static int windowHeight = 576;
static double zNear = 0.01, zFar = 10.0;
static double cameraRadius = 1.0;
Vec3d cameraFocusPosition(0.0);
static InputDevice id;

static int moveVerticesTogether = 0;

static int renderTextInfo = 1;

// glui
static GLUI * glui;

static double cameraLongitude, cameraLattitude;
static SphericalCamera * camera = nullptr;
string cameraPosFilename;
static int renderWireframe = 1;
static int renderEmbeddedWireframe = 0;
static int renderDeformableObject = 1;
static int renderEmbeddedObject = 1;
static int renderAxes=1;
static int renderGroundPlane = 0;
static double lineThickness = 3.0;

static ARAPDeformer * deformer = nullptr;

static int numFixedDOFs = 0;
static int * fixedDOFs = nullptr;
static Lighting * lighting = nullptr;
static CameraLighting * cameraLighting = nullptr;
double cameraLightIntensity = 1.0;
int lightOn[8] = {0, 0, 1, 1, 0, 0, 1, 1 };
float lightIntensity = 0.5;
int showSpecularLight = 0;
string lightingConfigSavingFilename = "lighting.config";

static SceneObjectDeformable * deformableMesh = nullptr, * comparedMesh = nullptr, * embeddedMesh = nullptr;
static SceneObject * extraSceneGeometry = nullptr;
static int buildExtraSceneGeometryNormals=1;
static int renderExtraSceneGeometryShadow = 1;
static int renderExtraSceneGeometry = 1;
static int enableTextures = 1;
static string groundPlaneString;
static SceneGroundPlane * groundPlane = nullptr;

#ifdef USE_GLSL
  static GLSLPhong * glslPhong = nullptr;
#endif
static int enableGLSLPhong = 1;

// config file
static string configFilename;
static string surfaceMeshFilename;
static string volumetricMeshFilename;
static string fixedVerticesFilename;
static string initialHandlesFilename;
// how many frames to deform initial handles to the input disp.
// < 0 means no pre-disp frames
static int numInitialHandlesPreDispFrames = 0;

//static char deformableObjectFilename[4096];
static char extraSceneGeometryFilename[4096];
static string groundPlaneTextureFilename;
static string lightingConfigFilename;
static string embeddedMeshFilename;
static string embeddedMeshInterpolationFilename;
static vector<double> u_embed;

static char backgroundColorString[4096] = "255 255 255";
static int numSolverThreads;
static int numOptIterations = 4;
static int saveScreenToFile = 0;
static int dumpScreenshots = 0;
static string screenshotBaseName;
static int numFramesToExit = -1; // negative means never exit automatically

// simulation
static double cpuLoad = 0;
static int graphicFrame = 0;

static PerformanceCounter titleBarCounter, explosionCounter, cpuLoadCounter;
static int numFixedVertices;
static int * fixedVertices;
static int renderFixedVertices = 1;

static PerformanceCounter worldCounter;
static int renderWorldCounter = 1;
// static int renderFPS = 1;

static int renderVertices = 0;
static int renderIKVertices = 1;

// selected manipulation vertices
typedef map<int, Vec3d> HandleMap;
typedef HandleMap::iterator HandleIter;
static HandleMap selectedIKVertices;

static int n = 0;
static AveragingBuffer fpsBuffer(5);
static int dumpDeformation = 0;

static bool interactiveDeformation = true;
static int maxIterations = 20;
static double iterationThreshold = 1e-5;

static string outputConstraintsFilename = "outputConstraints.txt";
static string initialDeformationFilename;

static vector<double> disp, comparedDisp;

static VolumetricMesh * volumetricMesh = nullptr;
static CubicMesh * cubicMesh = nullptr;
static TetMesh * tetMesh = nullptr;
static ObjMesh * volSurfaceMesh = nullptr;

static bool reverseHandle = false;
static int numHiddenHandles = 50;
static int numFrames = 1;

static int currentFrame = 0;
static bool verticalWindowSplit = true;
static bool renderPrimaryWindow = true;
static bool renderSecondaryWindow = false;
static bool compareReducedWithFull = false;
static BarycentricCoordinates * interpCoords = nullptr;
static ConfigFile configFile;

static map<int, Vec3d> initialHandleDisp;

static FrameRecorder * dispRecorder = nullptr;
static string dispRecordFilename;
static HandleControl handleControl;
static CameraChangeLoad cameraChangeLoad;
static string cameraChangeLoadFilename;
static int renderHandlesWithOffset = 0;
static int globalFrameCount = -1;
static int fullScreen = 0;
static string handleObjMeshFilename;
static ObjMesh * handleObjMesh = nullptr;
static string handleObjMeshDispRecordFilename;
static FrameRecorder * handleObjMeshRecorder = nullptr;
static bool key_h_pressed = false;

const static char * mouseControlStr[8] =
{
  "CTRL + click: and/remove a handle",
  "(on Mac, use key 'h' + click)",
  "Left click: show a nearby handle",
  "Drag handle: change handle location",
  "Middle + drag: zoom in/out",
  "Right + drag: rotate camera",
  "Keys 'e', 'w': toggle volumetric mesh",
  "Keys 'E', 'W': toggle embedded mesh"
};

#define ADD_CONFIG(v) configFile.addOptionOptional(#v, &v, v)
#define HAS_STRING(v) ((v).size() > 0)

// Dump constraints to a file. The format is the same as for reading constraints from a file.
static void saveConstraints()
{
  ofstream fout(outputConstraintsFilename.c_str(), ios::binary);
  if(!fout)
  {
    cout << "Cannot open " << outputConstraintsFilename << endl;
    return;
  }

  for(HandleIter it = selectedIKVertices.begin(); it != selectedIKVertices.end(); it++)
  {
    int vtx = it->first;
    Vec3d disp = it->second;
    fout << vtx << " " << disp[0] << " " << disp[1] << " " << disp[2] << endl;
  }
  fout.close();
}

static void ClearSelectedIKVertices()
{
  selectedIKVertices.clear();
  handleControl.clearHandleSelection();
}

static void Sync_GLUI()
{
  glui->sync_live();
}

// displayFunction is always called before idleFunction
static void idleFunction(void)
{
  cpuLoadCounter.StartCounter();
  glutSetWindow(windowID);

  globalFrameCount++;

  // Because displayFunction() is called before idleFunction(), we dump animation for the rendered shape.
  if (dumpDeformation == 1)
  {
    dumpDeformation = 0;
    static int count = 0;
    string filename = "def" + to_string(count) + ".u";
    cout << "Writing deformation to " << filename << "..." << flush;
    WriteMatrixToDisk_(filename.c_str(),3*n,1,&disp[0]);
    cout << "Done" << endl;
    count++;
  }
  else if (dumpDeformation == 2)
  {
    dumpDeformation = 0;
    static int count = 0;

    if (cubicMesh)
    {
      CubicMesh newCubicMesh(*cubicMesh);
      newCubicMesh.applyDeformation(&disp[0]);
      string filename = "cubicMesh" + to_string(count) + ".veg";
      newCubicMesh.save(filename.c_str());
      filename += ".obj";
      deformableMesh->ResetDeformationToRest();
      deformableMesh->AddVertexDeformations(&disp[0]);
      deformableMesh->GetMesh()->save(filename);
    } else if (tetMesh)
    {
      TetMesh newTetMesh(*tetMesh);
      newTetMesh.applyDeformation(&disp[0]);
      string filename = "tetMesh" + to_string(count) + ".veg";
      newTetMesh.save(filename.c_str());
      filename += ".obj";
      deformableMesh->ResetDeformationToRest();
      deformableMesh->AddVertexDeformations(&disp[0]);
      deformableMesh->GetMesh()->save(filename);
    }
    if (embeddedMesh)
    {
      string filename = "embMesh" + to_string(count) + ".obj";
      embeddedMesh->GetMesh()->save(filename);
    }
    count++;
  }
  if (dispRecorder)
  {
    dispRecorder->addFrame(&disp[0]);
  }

  if (handleObjMeshRecorder) // record the motion for the displacement on the idleFunc
  {
    vector<double> buffer;
    for(const pair<int,Vec3d> & p : selectedIKVertices)
    {
      int vtx = p.first;
      Vec3d pos = deformableMesh->GetSingleVertexRestPosition(vtx);
      pos += p.second;
      int numObjVtx = handleObjMesh->getNumVertices();
      for(int i = 0; i < numObjVtx; i++)
        for(int j = 0; j < 3; j++)
          buffer.push_back(pos[j]);
    }
    handleObjMeshRecorder->addFrame(buffer.data());
  }

  // Because displayFunction() is called before idleFunction(), we have to give a +1 below.
  if (numInitialHandlesPreDispFrames > 0 && globalFrameCount+1 <= numInitialHandlesPreDispFrames)
  {
    double fac = (globalFrameCount+1) / (double)numInitialHandlesPreDispFrames;
    for(const pair<int, Vec3d> & p : initialHandleDisp)
    {
      selectedIKVertices[p.first] = fac * p.second;
    }
  }

  if (cameraChangeLoadFilename.size() > 0)
  {
    cameraChangeLoad.controlCamera(globalFrameCount, *camera);
  }

  if (cameraLighting)
  {
    cameraLighting->SetLightIntensity(cameraLightIntensity);
  }
  else if (lighting)
  {
    for(int i = 0; i < 8; i++)
    {
      lighting->SetLightEnabled(i, lightOn[i]);
      lighting->SetLightIntensity(i, lightIntensity);
    }
    lighting->SetSpecularEnabled(showSpecularLight);
  }

  deformer->updateHandles(selectedIKVertices);

  // Take appropriate action in case the user is dragging a vertex.
  auto processDrag = [&](int vertex, Vec3d posDiff)
  {
    if (!moveVerticesTogether)
    {
      auto IKPulledVertexIter = selectedIKVertices.find(vertex);
      (IKPulledVertexIter->second) += posDiff;
      if (len2(posDiff) > 0)
        cout << IKPulledVertexIter->first << ": " << IKPulledVertexIter->second << endl;
    }
    else
    {
      for(auto iter = selectedIKVertices.begin(); iter != selectedIKVertices.end(); iter++)
        (iter->second) += posDiff;
    }
  };
  handleControl.processHandleMovement(id.getMousePosX(), id.getMousePosY(), id.shiftPressed(), processDrag);

  if(interactiveDeformation)
  {
    deformer->deformOneIter(&disp[0], &disp[0]);
  }
  else
  {
    int scale = 1;
    assert(maxIterations >= 1);
    assert(iterationThreshold >= 0);

    deformer->deform(&disp[0], &disp[0], iterationThreshold / scale, maxIterations * scale);
  }

  if (interpCoords && embeddedMesh && (renderEmbeddedObject || renderEmbeddedWireframe))
  {
    interpCoords->deform(&disp[0], &u_embed[0]);
    embeddedMesh->SetVertexDeformations(&u_embed[0]);
    embeddedMesh->BuildNormals();
  }

  // render the deformable object
  if (renderDeformableObject || renderWireframe)
  {
    deformableMesh->SetVertexDeformations(&disp[0]);
    deformableMesh->BuildNormals();
  }

  graphicFrame++;

  // update title bar at 4 Hz
  titleBarCounter.StopCounter();
  double elapsedTime = titleBarCounter.GetElapsedTime();
  if (elapsedTime >= 1.0 / 4)
  {
    titleBarCounter.StartCounter();
    double fps = graphicFrame / elapsedTime;

    fpsBuffer.addValue(fps);

    // update menu bar
    char windowTitle[4096];

    sprintf(windowTitle, "Vertices: %d | %.1f FPS | iter %d ", n, fpsBuffer.getAverage(), globalFrameCount);

    glutSetWindowTitle(windowTitle);
    graphicFrame = 0;

    //char ptext[96];
    //sprintf(ptext, "Force assembly: %G", forceAssemblyTime);
    //forceAssemblyStaticText->set_text(ptext);
    //sprintf(ptext, "System solve: %G", systemSolveTime);
    //systemSolveStaticText->set_text(ptext);

    Sync_GLUI();
  }

  cpuLoadCounter.StopCounter();
  double cpuTimePerGraphicsFrame = cpuLoadCounter.GetElapsedTime();
  cpuLoad = cpuTimePerGraphicsFrame * fpsBuffer.getAverage();

  if (saveScreenToFile || dumpScreenshots)
  {
    char defaultBaseName[] = "pic";
    const char * baseName = defaultBaseName;
    if (screenshotBaseName.size() > 0)
      baseName = screenshotBaseName.c_str();

    static int sprite = 0;
    char s[100] = "";
    sprintf(s,"%s.%04d.png",baseName,sprite);
    Screenshot::SaveScreenshot(s, ImageIO::FORMAT_PNG,windowWidth, windowHeight);
    sprite++;

    saveScreenToFile=0; // save only once
  }

  if (numFramesToExit >= 0)
  {
    if (globalFrameCount+1 >= numFramesToExit)
    {
      cout << "Exiting after " << numFramesToExit << " frames" << endl;
      exit(0);
    }
  }
  glutPostRedisplay();
}

static void initScene()
{
  // initialize the rendering mesh for the deformable object
  if(HAS_STRING(volumetricMeshFilename))
  {
    volumetricMesh = VolumetricMeshLoader::load(volumetricMeshFilename.c_str());
    assert(volumetricMesh);
    if (volumetricMesh->getElementType() == VolumetricMesh::TET)
    {
      tetMesh = dynamic_cast<TetMesh*>(volumetricMesh);
    }
    else if (volumetricMesh->getElementType() == VolumetricMesh::CUBIC)
    {
      cubicMesh = dynamic_cast<CubicMesh*>(volumetricMesh);
    }
    volSurfaceMesh = GenerateSurfaceMesh::ComputeMesh(volumetricMesh);
    volSurfaceMesh->buildFaceNormals();
    deformableMesh = new SceneObjectDeformable(volSurfaceMesh, true);
  }
  else if (HAS_STRING(surfaceMeshFilename))
  {
    deformableMesh = new SceneObjectDeformable(surfaceMeshFilename.c_str());
  }

  if (deformableMesh->HasTextures() && enableTextures)
  {
    deformableMesh->EnableTextures();
    deformableMesh->SetUpTextures();
  }
  else
    deformableMesh->DisableTextures();

  n = deformableMesh->GetNumVertices();
  //deformableMesh->EnableCustomColor();
  //deformableMesh->SetCustomColor(vector<Vec3d>(n, Vec3d(1,0,0)));

  disp.resize(3*n, 0.);
  if(HAS_STRING(initialDeformationFilename))
  {
    Matrix<double> m(initialDeformationFilename.c_str());
    assert(m.Getm() == 3*n && m.Getn() > 0);
    memcpy(disp.data(), m.GetData(), sizeof(double) * m.Getm());
  }

  // init lighting
  cameraLighting = new CameraLighting;
  deformableMesh->ResetDeformationToRest();
  deformableMesh->BuildNeighboringStructure();
  deformableMesh->BuildNormals();
  deformableMesh->SetMaterialAlpha(0.5);

  // load/create fixed vertices
  // 1-indexed notation
  numFixedVertices = 0;
  fixedVertices = nullptr;
  if(fixedVerticesFilename.size() > 0)
  {
    if (ListIO::load(fixedVerticesFilename.c_str(), &numFixedVertices, &fixedVertices, 1) != 0)
    {
      printf("Error reading fixed vertices.\n");
      exit(1);
    }

    printf("Loaded %d fixed vertices.\n", numFixedVertices);
    //ListIO::print(numFixedVertices, fixedVertices);

    // check if indices are within bounds
    for(int i=0; i<numFixedVertices; i++)
    {
      if ((fixedVertices[i] < 0) || (fixedVertices[i] >= n))
      {
        printf("Error: encountered a boundary vertex out of bounds. Vertex: %d. Valid range: [0..%d).\n", fixedVertices[i], n);
        exit(1);
      }
    }
  }

  numFixedDOFs = 3 * numFixedVertices;
  fixedDOFs = (int*) malloc (sizeof(int) * numFixedDOFs);

  for(int i=0; i<numFixedVertices; i++)
  {
    fixedDOFs[3*i+0] = 3 * (fixedVertices[i]) + 0;
    fixedDOFs[3*i+1] = 3 * (fixedVertices[i]) + 1;
    fixedDOFs[3*i+2] = 3 * (fixedVertices[i]) + 2;
  }

  cout << "Building ARAP shape deformer..." << endl;
  if (tetMesh)
    deformer = new ARAPDeformer(tetMesh, numFixedVertices, fixedVertices, numSolverThreads);
  else if (cubicMesh)
    deformer = new ARAPDeformer(cubicMesh, numFixedVertices, fixedVertices, numSolverThreads);
  else
    deformer = new ARAPDeformer(deformableMesh->GetMesh(), numFixedVertices, fixedVertices, numSolverThreads);

  if (initialHandlesFilename.size() > 0)
  {
    ifstream fin(initialHandlesFilename);
    assert(fin);
    fin >> ws;
    while(!fin.eof())
    {
      int vertex = -1;
      Vec3d d(0.);
      fin >> vertex >> d[0] >> d[1] >> d[2];
      assert(vertex >= 0 && vertex < n);
      initialHandleDisp.emplace(vertex, d);
      assert(fin.fail() == false);
      fin >> ws;
    }

    if (numInitialHandlesPreDispFrames <= 0)
    {
      selectedIKVertices.insert(initialHandleDisp.begin(), initialHandleDisp.end());
    } else
    {
      for(const pair<int, Vec3d> & p : initialHandleDisp)
        selectedIKVertices.emplace(p.first, Vec3d(0.0));
    }
  }

  // Load any external geometry file (e.g. some static scene for decoration; usually there will be none).
  printf("Loading extra scene geometry...\n");
  if (strcmp(extraSceneGeometryFilename,"__none") != 0)
  {
    extraSceneGeometry = new SceneObjectDeformable(extraSceneGeometryFilename);
    if (buildExtraSceneGeometryNormals)
    {
      extraSceneGeometry->BuildNeighboringStructure();
      extraSceneGeometry->BuildNormals();
    }
    if (enableTextures)
    {
      if (extraSceneGeometry->SetUpTextures() != 0)
        exit(1);
    }
  }
  else
    extraSceneGeometry = nullptr;

  if (groundPlaneString.size() > 0)
  {
    groundPlane = new SceneGroundPlane(groundPlaneString.c_str(), groundPlaneTextureFilename.c_str());
    renderGroundPlane = 1;
  }

  // set background color
  int colorR, colorG, colorB;
  sscanf(backgroundColorString, "%d %d %d", &colorR, &colorG, &colorB);
  glClearColor(1.0 * colorR / 255, 1.0 * colorG / 255, 1.0 * colorB / 255, 0.0);
  
  #ifdef USE_GLSL
    glslPhong = new GLSLPhong();
  #endif

  if(HAS_STRING(volumetricMeshFilename) && HAS_STRING(embeddedMeshFilename))
  {
    embeddedMesh = new SceneObjectDeformable(embeddedMeshFilename.c_str());
    if (enableTextures && embeddedMesh->HasTextures())
      embeddedMesh->SetUpTextures(SceneObject::MODULATE, SceneObject::NOMIPMAP,
          SceneObject::USEANISOTROPICFILTERING, SceneObject::USETEXTURETRANSPARENCY);
    else
      embeddedMesh->DisableTextures();
    // cout << "objRdrMesh has texture ? " << embeddedMesh->HasTextures() << endl;
    embeddedMesh->ResetDeformationToRest();
    embeddedMesh->BuildNeighboringStructure();
    embeddedMesh->BuildNormals();
    // embeddedMesh->SetMaterialAlpha(0.5);
    u_embed.assign(embeddedMesh->GetNumVertices() * 3, 0.0);

    if (HAS_STRING(embeddedMeshInterpolationFilename))
    {
      interpCoords = new BarycentricCoordinates(embeddedMeshInterpolationFilename);
      assert((int)interpCoords->getNumLocations() == embeddedMesh->GetNumVertices());
    }
    else
    {
      cout << "Error: no embeding weights provided" << endl;
      exit(1);
    }
  }

  if (dispRecordFilename.size() > 0)
  {
    dispRecorder = new FrameRecorder(n*3, dispRecordFilename);
  }

  double cameraUp[3] = {0,1,0};

  Vec3d cameraFocus;
  Vec3d bmin, bmax;
  deformableMesh->GetMesh()->computeBoundingBox();
  deformableMesh->GetMesh()->getCubicBoundingBox(1.0, &bmin, &bmax);
  Vec3d modelCenter = (bmin + bmax) / 2.0;
  double modelRadius = deformableMesh->GetMesh()->getDiameter() / 2;

  if (configFile.isOptionLoaded("cameraRadius") == false)
  {
    cameraRadius = modelRadius * 3;
  }
  if (configFile.isOptionLoaded("cameraFocusPosition"))
    cameraFocus = cameraFocusPosition;
  else
    cameraFocus = modelCenter;
  zNear = cameraRadius * 0.01;
  zFar = cameraRadius * 100.0;

  camera = new SphericalCamera(cameraRadius,
            1.0 * cameraLongitude / 360 * (2.0 * PI),
            1.0 * cameraLattitude / 360 * (2.0 * PI),
            &cameraFocus[0], cameraUp, 0.05);
  if (cameraPosFilename.size() > 0)
    camera->LoadPosition(cameraPosFilename.c_str());
  // else
  //   camera->LoadPosition("cameraPos.txt");

  // update lighting parameters
  if (cameraLighting)
  {
    cameraLighting->SetLightIntensity(cameraLightIntensity);
  }
  else if (lighting)
  {
    for(int i = 0; i < 8; i++)
    {
      lighting->SetLightEnabled(i, lightOn[i]);
      lighting->SetLightIntensity(i, lightIntensity);
    }
    lighting->SetSpecularEnabled(showSpecularLight);
  }

  if (cameraChangeLoadFilename.size() > 0)
  {
    int ret = cameraChangeLoad.load(cameraChangeLoadFilename);
    assert(ret == 0);
  }

  if (handleObjMeshFilename.size() > 0)
  {
    handleObjMesh = new ObjMesh(handleObjMeshFilename);
    cout << "Loaded handle mesh: " << handleObjMeshFilename << "." << endl;

    if (handleObjMeshDispRecordFilename.size() > 0)
    {
      int numHandles = initialHandleDisp.size();
      handleObjMeshRecorder = new FrameRecorder(numHandles * handleObjMesh->getNumVertices() * 3, handleObjMeshDispRecordFilename);
    }
  }

  Sync_GLUI();

  printf("Initialization completed.\n");
}

void initConfigurations()
{
  printf("Parsing configuration file %s...\n", configFilename.c_str());

  // Specify the entries of the config file.

  // At least one of the following three must be present:
  ADD_CONFIG(fixedVerticesFilename);
  ADD_CONFIG(groundPlaneTextureFilename);
  ADD_CONFIG(initialHandlesFilename);
  ADD_CONFIG(numInitialHandlesPreDispFrames);
  ADD_CONFIG(dumpScreenshots);
  ADD_CONFIG(screenshotBaseName);
  ADD_CONFIG(numFramesToExit);
  ADD_CONFIG(renderFixedVertices);
  ADD_CONFIG(cameraPosFilename);
  ADD_CONFIG(dispRecordFilename);
  ADD_CONFIG(cameraLightIntensity);
  ADD_CONFIG(lineThickness);
  ADD_CONFIG(renderWireframe);
  ADD_CONFIG(renderIKVertices);
  ADD_CONFIG(renderDeformableObject);
  ADD_CONFIG(renderEmbeddedWireframe);
  ADD_CONFIG(renderHandlesWithOffset);
  ADD_CONFIG(cameraChangeLoadFilename);
  ADD_CONFIG(fullScreen);
  ADD_CONFIG(handleObjMeshFilename);
  ADD_CONFIG(handleObjMeshDispRecordFilename);
  ADD_CONFIG(cameraFocusPosition);

  configFile.addOptionOptional("numSolverThreads", &numSolverThreads, 1);
  configFile.addOptionOptional("numOptIterations", &numOptIterations, 4);

  configFile.addOptionOptional("windowWidth",&windowWidth,windowWidth);
  configFile.addOptionOptional("windowHeight",&windowHeight,windowHeight);
  configFile.addOptionOptional("cameraRadius",&cameraRadius,17.5);
  configFile.addOptionOptional("cameraLongitude",&cameraLongitude,-60.0);
  configFile.addOptionOptional("cameraLattitude",&cameraLattitude,20.0);
  configFile.addOptionOptional("renderVertices",&renderVertices, renderVertices);
  configFile.addOptionOptional("renderAxes",&renderAxes,renderAxes);
  configFile.addOptionOptional("extraSceneGeometry",extraSceneGeometryFilename,"__none");
  configFile.addOptionOptional("buildExtraSceneGeometryNormals", &buildExtraSceneGeometryNormals, buildExtraSceneGeometryNormals);
  configFile.addOptionOptional("renderExtraSceneGeometryShadow", &renderExtraSceneGeometryShadow, renderExtraSceneGeometryShadow);

  configFile.addOptionOptional("enableTextures",&enableTextures,enableTextures);
  configFile.addOptionOptional("backgroundColor",backgroundColorString, backgroundColorString);

  ADD_CONFIG(groundPlaneString);
  configFile.addOptionOptional("renderWorldCounter",&renderWorldCounter, renderWorldCounter);
  configFile.addOptionOptional("enableGLSLPhong",&enableGLSLPhong, enableGLSLPhong);

  ADD_CONFIG(embeddedMeshFilename);
  ADD_CONFIG(embeddedMeshInterpolationFilename);
  ADD_CONFIG(interactiveDeformation);
  ADD_CONFIG(maxIterations);
  ADD_CONFIG(iterationThreshold);

  ADD_CONFIG(outputConstraintsFilename);
  ADD_CONFIG(initialDeformationFilename);
  ADD_CONFIG(surfaceMeshFilename);
  ADD_CONFIG(volumetricMeshFilename);
  ADD_CONFIG(lightingConfigFilename);
  ADD_CONFIG(numHiddenHandles);
  ADD_CONFIG(numFrames);
  ADD_CONFIG(compareReducedWithFull);

  // parse the configuration file
  if (configFile.parseOptions((char*)configFilename.c_str()) != 0)
  {
    printf("Error parsing options.\n");
    exit(1);
  }

  // the config variables have now been loaded with their specified values

  // informatively print the variables (with assigned values) that were just parsed
  configFile.printOptions();
}

static void renderScene(int frame, bool renderHandle)
{
  // if GPU rendering, must set lighting
  // else set lights directly via OpenGL (all done inside the following function)
  if (cameraLighting)
  {
    cameraLighting->LightScene(camera);
  }
  else if (lighting)
  {
    deformableMesh->SetLighting(lighting);
    lighting->LightScene();
  }

  glEnable(GL_LIGHTING);

  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);  //only when stencil pass and z-buffer pass, set stencil value to stencil reference
  glStencilFunc(GL_ALWAYS, 0, ~(0u)); //always pass stencil test, stencil renference value is 0

  #ifdef USE_GLSL
    if (enableGLSLPhong)
      glslPhong->Enable();
  #endif

  int allOffset = 0.0;
  if (renderHandlesWithOffset) // we offset all the faces so that the handle points appear in front of everything
    allOffset = 3;
  PolygonOffsetFillState allFaceOffset(allOffset,allOffset);

  if (renderGroundPlane && groundPlane)
  {
    glEnable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    groundPlane->render();
    glDisable(GL_TEXTURE_2D);
  }

  glStencilFunc(GL_ALWAYS, 1, ~(0u)); // always pass the stencil test, stencil renference value is 1

  if (renderEmbeddedObject && embeddedMesh)
    embeddedMesh->Render();

  if (renderDeformableObject && deformableMesh)
  {
    //  cout << "set disp at frame " << currentFrame << endl;
    bool doBlending = renderEmbeddedObject && (embeddedMesh || comparedMesh);
    if (doBlending)
    {
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_BLEND);
      // PolygonOffsetFillState offsetState(1.0, 1.0);
      glDrawBuffer(GL_NONE);
      deformableMesh->Render();
      glDrawBuffer(GL_BACK);
    }

    {
      // PolygonFillOffsetState offsetState(1.0, 1.0);
      if (volSurfaceMesh)
        deformableMesh->EnableFlatFaces();

      if (doBlending)
      {
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
      }
      // cout << "render vol mesh" << endl;
      deformableMesh->Render();
      if (doBlending)
      {
        glDisable(GL_CULL_FACE);
      }
    }

    if (renderVertices)
    {
      glDisable(GL_LIGHTING);
      glColor3f(0.5,0,0);
      glPointSize(8.0);

      deformableMesh->RenderVertices();

      glEnable(GL_LIGHTING);
    }

    if (doBlending)
      glDisable(GL_BLEND);
  }

  // render any extra scene geometry
  glStencilFunc(GL_ALWAYS, 0, ~(0u));
  if ((renderExtraSceneGeometry) && (extraSceneGeometry != nullptr))
  {
    extraSceneGeometry->Render();
  }

  glDisable(GL_TEXTURE_2D);
  #ifdef USE_GLSL
    glslPhong->Disable();
  #endif

  // render shadow
  if (renderGroundPlane)
  {
    glColor3f(0.1,0.1,0.1);
    glDisable(GL_LIGHTING);

    if (renderDeformableObject && deformableMesh && groundPlane)
      groundPlane->renderShadow(deformableMesh);

    if ((renderExtraSceneGeometry) && (extraSceneGeometry) && groundPlane)
      groundPlane->renderShadow(extraSceneGeometry);

    glEnable(GL_LIGHTING);
  }

  glDisable(GL_LIGHTING);

  glStencilFunc(GL_ALWAYS, 1, ~(0u));
  glColor3f(0,0,0);
  if (renderWireframe)
  {
    glLineWidth(lineThickness);
    deformableMesh->RenderEdges();
  }

  if (renderEmbeddedWireframe)
  {
    glColor3f(0,0,0);
    glLineWidth(lineThickness);
    if (embeddedMesh)
      embeddedMesh->RenderEdges();
  }

  glColor3f(0,0.5,0);
  glPointSize(10.0);
  if (renderVertices)
    deformableMesh->RenderVertices();

  #ifdef USE_GLSL
    glslPhong->Disable();
  #endif
  allFaceOffset.restore();

  glStencilFunc(GL_ALWAYS, 1, ~(0u));
  if(renderIKVertices)
  {
    glPointSize(17.0);

    glBegin(GL_POINTS);
    glColor3f(0,0,1);
    for(map<int,Vec3d> :: iterator iter = selectedIKVertices.begin(); iter != selectedIKVertices.end(); iter++)
    {
      double pulledVertexPos[3];
      deformableMesh->GetSingleVertexRestPosition(iter->first,
          &pulledVertexPos[0], &pulledVertexPos[1], &pulledVertexPos[2]);
      pulledVertexPos[0] += iter->second[0];
      pulledVertexPos[1] += iter->second[1];
      pulledVertexPos[2] += iter->second[2];
      glVertex3f(pulledVertexPos[0], pulledVertexPos[1], pulledVertexPos[2]);
    }
    glEnd();
  }

  // turn off stencil modifications
  glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

  glColor3f(0,0,0);

  if (renderAxes)
  {
    glLineWidth(1.0);
    RenderAxes(1.0);
  }

  double pulledVertexPos[3];

  // render the vertex currently being manipulated via IK
  if (handleControl.isHandleSelected())
  {
    deformableMesh->GetSingleVertexPositionFromBuffer(handleControl.getSelectedHandle(), &pulledVertexPos[0], &pulledVertexPos[1], &pulledVertexPos[2]);
    glColor3f(1,0,0);

    glPointSize(8.0);
    glBegin(GL_POINTS);
    glVertex3f(pulledVertexPos[0], pulledVertexPos[1], pulledVertexPos[2]);
    glEnd();

    // render the move handle
    Vec3d handlePosition;
    deformableMesh->GetSingleVertexRestPosition(handleControl.getSelectedHandle(), &handlePosition[0], &handlePosition[1], &handlePosition[2]);
    auto it = selectedIKVertices.find(handleControl.getSelectedHandle());
    assert(it != selectedIKVertices.end());
    handlePosition += it->second;

    handleControl.renderHandle(camera, handlePosition, reverseHandle);
  }

  // render fixed vertices
  if (renderFixedVertices)
  {
    glPointSize(18.0);
    glColor3f(1,0,0);
    glBegin(GL_POINTS);
    for(int i=0; i<numFixedVertices; i++)
    {
      glColor3f(1,0,0);
      double fixedVertexPos[3];
      deformableMesh->GetSingleVertexRestPosition(fixedVertices[i],
          &fixedVertexPos[0], &fixedVertexPos[1], &fixedVertexPos[2]);
      glVertex3f(fixedVertexPos[0], fixedVertexPos[1], fixedVertexPos[2]);
    }
    glEnd();
  }
}

static void displayFunction()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);

  // Setup model transformations.
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  camera->Look();

  if (renderPrimaryWindow)
  {
    if (renderSecondaryWindow)
    {
      if (verticalWindowSplit)
        glViewport(0, 0, windowWidth/2, windowHeight);
      else
        glViewport(0, windowHeight/2, windowWidth, windowHeight/2);
    }
    else
      glViewport(0,0,windowWidth,windowHeight);

    renderScene(currentFrame, true);
  }

  if (renderSecondaryWindow)
  {
    if (renderPrimaryWindow)
    {
      if (verticalWindowSplit)
        glViewport(windowWidth/2,0,windowWidth/2,windowHeight);
      else
        glViewport(0,0,windowWidth,windowHeight/2);
    }
    else
      glViewport(0,0,windowWidth,windowHeight);

    renderScene(0, false);
  }

  glutSwapBuffers();
}

void reset_buttonCallBack(int code=0)
{
  ClearSelectedIKVertices();

  Sync_GLUI();
}

void keyboardFunction(unsigned char key, int x, int y)
{
  Vec3d objectCenter;
  id.setModifiers();
  id.setMousePos(x, y);

  switch (key)
  {
  case 'd':
    dumpDeformation = 1;
    break;

  case 'D':
    dumpDeformation = 2;
    break;

  case 9:
  {
    fullScreen = !fullScreen;
    if (fullScreen)
    {
      cout << "Full screen" << endl;
      glutFullScreen();
    }
    else
    {
      cout << "Normal screen" << endl;
      glutReshapeWindow(800, 600);
      glutPositionWindow(50, 50);
    }
    break;
  }

  case 27:
    exit(0);
    break;

  case 'w':
    renderWireframe = !renderWireframe;
    break;

  case 'W':
    renderEmbeddedWireframe = !renderEmbeddedWireframe;
    break;

  case '?':
    renderTextInfo = !renderTextInfo;
    break;

  case 'i':
  {
    double cameraX,cameraY,cameraZ;
    camera->GetAbsWorldPosition(cameraX,cameraY,cameraZ);
    printf("Camera is positioned at: %G %G %G\n",cameraX,cameraY,cameraZ);
    printf("Camera radius is: %G \n",camera->GetRadius());
    printf("Camera Phi is: %G \n",180.0/M_PI*camera->GetPhi());
    printf("Camera Theta is: %G \n",180.0/M_PI*camera->GetTheta());
    double focusPos[3];
    camera->GetFocusPosition(focusPos);
    printf("Camera focus position is: %G %G %G\n", focusPos[0], focusPos[1], focusPos[2]);
  }
  break;

  case '\\':
    camera->Reset();
    break;

  case 13:
    break;

  case '0':
    reset_buttonCallBack(0); // hard reset
    break;

  case '9':
    reset_buttonCallBack(1); // soft reset
    break;

  case 'a':
    renderAxes = !renderAxes;
    Sync_GLUI();
    break;

  case '[':
    printf("Loading camera position.\n");
    camera->LoadPosition("cameraPos.txt");
    break;

  case ']':
    printf("Saving camera position.\n");
    camera->SavePosition("cameraPos.txt");
    break;

  case 'b':
    renderFixedVertices = !renderFixedVertices;
    break;

  case 'v':
    renderIKVertices = !renderIKVertices;
    break;

  case 'V':
    renderVertices = !renderVertices;
    Sync_GLUI();
    break;

  case 'e':
    renderDeformableObject = !renderDeformableObject;
    break;

  case 'E':
    renderEmbeddedObject = !renderEmbeddedObject;
    break;

  case '!':
    renderExtraSceneGeometry = !renderExtraSceneGeometry;
    break;

  case 't':
    reverseHandle = !reverseHandle;
    printf("Reverse handle is now %s .\n", reverseHandle ? "ON" : "OFF");
    break;

  case 'P':
    enableGLSLPhong = !enableGLSLPhong;
    printf("Phong shader is now %s.\n", enableGLSLPhong ? "ON" : "OFF");
    break;

  case 'g':
    renderGroundPlane = !renderGroundPlane;
    break;

  case 'x':
    saveScreenToFile = !saveScreenToFile;
    break;

  case 'y':
    moveVerticesTogether = !moveVerticesTogether;
    break;

  case 'h':
    key_h_pressed = true;
    break;

  case 'c': // dump constraints
    cout << "Dump constraints to " << outputConstraintsFilename << "." << endl;
    saveConstraints();
    break;

  case '-':
    if (currentFrame > 0)
      currentFrame--;
    cout << "current frame " << currentFrame;
    break;

  case '=':
    if (currentFrame+1 < numFrames)
      currentFrame++;
    cout << "current frame " << currentFrame;
    break;
  }
}

void specialFunction(int key, int x, int y)
{
  id.setModifiers();
  id.setMousePos(x, y);

  switch (key)
  {
  case GLUT_KEY_HOME:
    break;

  case GLUT_KEY_PAGE_UP:
    break;

  case GLUT_KEY_PAGE_DOWN:
    break;

  case GLUT_KEY_LEFT:
    camera->MoveFocusRight(+0.03 * fabs(camera->GetRadius()));
    break;

  case GLUT_KEY_RIGHT:
    camera->MoveFocusRight(-0.03 * fabs(camera->GetRadius()));
    break;

  case GLUT_KEY_DOWN:
    camera->MoveFocusUp(+0.03 * fabs(camera->GetRadius()));
    break;

  case GLUT_KEY_UP:
    camera->MoveFocusUp(-0.03 * fabs(camera->GetRadius()));
    break;

  case GLUT_KEY_END:
    break;

  case GLUT_KEY_INSERT:
    break;

  default:
    break;
  }
}

void reshape(int x, int y)
{
  windowWidth = x;
  windowHeight = y;

  glMatrixMode(GL_PROJECTION); // Select The Projection Matrix
  glLoadIdentity(); // Reset The Projection Matrix
  gluPerspective(45.0f, 1.0 * windowWidth / windowHeight, cameraRadius*0.01f, 20*cameraRadius);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void exitHandler()
{
  if (dispRecorder)
    dispRecorder->save();
  if (handleObjMeshRecorder)
    handleObjMeshRecorder->save();
}

void mouseMotionFunction(int x, int y)
{
  int mouseDeltaX = x-id.getMousePosX();
  int mouseDeltaY = y-id.getMousePosY();

  id.setMousePos(x, y);

  if (id.rightMouseButtonDown()) // right mouse button handles camera rotations
  {
    const double factor = 0.2;

    camera->MoveRight(factor * mouseDeltaX);
    camera->MoveUp(factor * mouseDeltaY);
  }

  if ((id.middleMouseButtonDown()) || (id.leftMouseButtonDown()&& id.altPressed()))// handle zoom in/out
  {
    const double factor = 0.1;
    camera->ZoomIn(cameraRadius * factor * mouseDeltaY);
  }
}
void mouseNoDrag(int x, int y)
{
  id.setMousePos(x, y);
  if (handleControl.isHandleSelected())
  {
    Vec3d worldPos(0.0);
    GLubyte stencilValue;
    float zValue;
    unprojectPointFromScreen(x,y, &worldPos[0], &stencilValue, &zValue);

    if (stencilValue == 1)
    {
      handleControl.setMousePosition(worldPos);
    }
  }
}

void mouseButtonActivityFunction(int button, int state, int x, int y)
{
  id.setButton(button, state);
  id.setMousePos(x,y);

  Vec3d worldPos(0.0);
  GLubyte stencilValue;
  float zValue;

  switch (button)
  {
  case GLUT_LEFT_BUTTON:
  {
    if (id.leftMouseButtonDown())
      glutSetCursor(GLUT_CURSOR_CROSSHAIR);
    else
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);

    unprojectPointFromScreen(x,y, &worldPos[0], &stencilValue, &zValue);

    auto getClosestHandle = [&]()
    {
      MinValueIndex vi;
      for(map<int,Vec3d> :: iterator iter = selectedIKVertices.begin(); iter != selectedIKVertices.end(); iter++)
      {
        Vec3d pos(0.0);
        deformableMesh->GetSingleVertexPositionFromBuffer(iter->first, &pos[0], &pos[1], &pos[2]);
        vi.update(len2(pos - worldPos), iter->first);
      }
      return vi.index;
    };

    auto addOrRemoveHandle = [&]()
    {
      int clickedVertex = deformableMesh->GetClosestVertex(worldPos, nullptr, nullptr);
      printf("(vertex IK handle select) Clicked on vertex: %d\n", clickedVertex);
      map<int,Vec3d> :: iterator iter = selectedIKVertices.find(clickedVertex);
      pair<int, bool> ret;
      bool createNew = false;
      if (iter == selectedIKVertices.end())
      {
        // select vertex
        createNew = true;
        double curPos[3];
        deformableMesh->GetSingleVertexPositionFromBuffer(clickedVertex, &curPos[0], &curPos[1], &curPos[2]);
        double restPos[3];
        deformableMesh->GetSingleVertexRestPosition(clickedVertex, &restPos[0], &restPos[1], &restPos[2]);

        Vec3d delta(curPos[0] - restPos[0], curPos[1] - restPos[1],curPos[2] - restPos[2]);
        selectedIKVertices.insert(make_pair(clickedVertex, delta));
        // set new constraint
        printf("Create new handle: %d\n", clickedVertex);
      }
      else
      {
        // de-select vertex
        selectedIKVertices.erase(iter);
      }
      return make_pair(clickedVertex, createNew);
    };

    handleControl.setMouseButtonActivity(id.leftMouseButtonDown(), stencilValue == 1, id.ctrlPressed() || key_h_pressed,
        worldPos, zValue, getClosestHandle, addOrRemoveHandle);

    break;
  }

  case GLUT_MIDDLE_BUTTON:

    if (id.middleMouseButtonDown())
      glutSetCursor(GLUT_CURSOR_UP_DOWN);
    else
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
    break;

  case GLUT_RIGHT_BUTTON:

    if (id.rightMouseButtonDown())
      glutSetCursor(GLUT_CURSOR_CYCLE);
    else
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
    break;
  }
}

void exit_buttonCallBack(int code)
{
  exit(0);
}

void gluiSaveLightConfigCallback(int)
{
  if (lighting)
    lighting->SaveConfig(lightingConfigSavingFilename.c_str());
}

void initGLUI()
{
  // generate the UI
  glui = GLUI_Master.create_glui("Controls", 0, windowWidth + 52, 0);

  // ******** view ********

  glui->add_checkbox("Display axes", &renderAxes);
  glui->add_checkbox("Display wireframe", &renderWireframe);
  glui->add_checkbox("Display vertices", &renderVertices);
  glui->add_separator();

  GLUI_Panel * instructions_panel =
     glui->add_panel("Mouse buttons:", GLUI_PANEL_EMBOSSED);
  instructions_panel->set_alignment(GLUI_ALIGN_LEFT);
  
  for(size_t i = 0; i < sizeof(mouseControlStr) / sizeof(const char *); i++)
    glui->add_statictext_to_panel(instructions_panel, mouseControlStr[i]);

  glui->add_separator();

  glui->add_button("Reset", 0, reset_buttonCallBack);

  glui->add_separator();
  glui->add_button("Exit program", 0, exit_buttonCallBack);

  Sync_GLUI();

  glui->set_main_gfx_window(windowID);
}

void keyboardUpFunction(unsigned char key, int x, int y)
{
  // cout << "key " << key << " up" << endl;
  switch(key)
  {
    case 'h':
      key_h_pressed = false;
      break;

    default:
      break;
  }
}

void initGLUT(int argc, char* argv[], char * windowTitle, int windowWidth, int windowHeight, int * windowID)
{
  // Initialize GLUT.
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_STENCIL | GLUT_MULTISAMPLE);
  glutInitWindowSize(windowWidth, windowHeight);
  *windowID = glutCreateWindow(windowTitle);
  glutPositionWindow(0,0);

  // Setup GLUT callbacks.
  GLUI_Master.set_glutDisplayFunc(displayFunction);
  //glutDisplayFunc(displayFunction);

  GLUI_Master.set_glutIdleFunc(idleFunction);
  GLUI_Master.set_glutKeyboardFunc(keyboardFunction);
  GLUI_Master.set_glutSpecialFunc(specialFunction);
  GLUI_Master.set_glutReshapeFunc(reshape);
  GLUI_Master.set_glutMouseFunc(mouseButtonActivityFunction);
  glutKeyboardUpFunc(keyboardUpFunction);
  glutMotionFunc(mouseMotionFunction);

  // callback for mouse movement without any buttons pressed
  glutPassiveMotionFunc(mouseNoDrag);
}

static void initGraphics(int windowWidth, int windowHeight)
{
   // clear to white
  glClearColor(256.0 / 256, 256.0 / 256, 256.0 / 256, 0.0);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_STENCIL_TEST);
  glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POLYGON_SMOOTH);
  glEnable(GL_LINE_SMOOTH);

  reshape(windowWidth,windowHeight);

  printf ("Graphics initialization complete.\n");
}

int main(int argc, char* argv[])
{
  int numFixedArgs = 2;
  if (argc < numFixedArgs)
  {
    printf("As-rigid-as-possible surface deformation editing.\n");
    printf("Usage: %s [config file]\n", argv[0]);
    return 1;
  }

  configFilename = argv[1];

  printf("Starting application.\n");

  printf("Loading scene configuration from %s.\n", configFilename.c_str());

  initConfigurations();

  atexit(exitHandler);

  printf("Opening %d x %d window.\n", windowWidth, windowHeight);
  initGLUT(argc, argv, windowTitleBase , windowWidth, windowHeight, &windowID);
  initGraphics(windowWidth, windowHeight);

  initGLUI();

  initScene();

  if (fullScreen)
    glutFullScreen();

  cout << "===================================================" << endl;
  cout << "Mouse buttons:" << endl;
  for(size_t i = 0; i < sizeof(mouseControlStr) / sizeof(const char *); i++)
    cout << "- " << mouseControlStr[i] << endl;
  
  glutMainLoop();

  return 0;
}

