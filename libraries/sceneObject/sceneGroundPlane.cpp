/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "sceneObject" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC    *
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

#include "sceneGroundPlane.h"
#include "openGL-headers.h"
#include "imageIO.h"
#include <cstring>
using namespace std;


SceneGroundPlane::SceneGroundPlane(const char * groundPlaneString, const char * groundPlaneTextureFilename)
{
  groundPlaneHeight = 0;
  lightPos = Vec3d(0,0,100);
  groundPlaneSize = 50;
  color[0] = color[1] = color[2] = 1.0f;
  ambientFactor = diffuseFactor = specularFactor = 1.0f;
  shininess = 120.0f;

  displayListGround = 0;
  displayListGroundWithTexture = 0;
  planeResolution = 100;
  groundTexName = 0;

  int ret = sscanf(groundPlaneString,"%lf,%lf,%lf,%lf,%lf,r%f,g%f,b%f,a%f,d%f,s%f,sh%f", &groundPlaneHeight,
      &lightPos[0], &lightPos[1], &lightPos[2],
      &groundPlaneSize, &color[0], &color[1], &color[2], &ambientFactor, &diffuseFactor, &specularFactor, &shininess);
  if (ret != 12)
    throw 1;

  displayListGround = glGenLists(1);
  glNewList(displayListGround, GL_COMPILE);
  buildDisplayList();
  glEndList();


  if (groundPlaneTextureFilename && strlen(groundPlaneTextureFilename) > 0)
  {
    ImageIO groundTexture;
    ImageIO::fileFormatType type;
    if (groundTexture.load(groundPlaneTextureFilename, &type) != ImageIO::OK)
    {
      cout << "Error reading image " << groundPlaneTextureFilename << "." << endl;
      throw 2;
    }

    // == create ground texture ==
    glGenTextures(1, &groundTexName);
    glBindTexture(GL_TEXTURE_2D, groundTexName);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); // repeat pattern in s
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); // repeat pattern in t
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB,
            groundTexture.getWidth(), groundTexture.getHeight(), GL_RGB,
            GL_UNSIGNED_BYTE, groundTexture.getPixels());

    //glTexImage2D(GL_TEXTURE_2D, 0 , GL_RGB,
    //groundTexture->nx, groundTexture->ny, 0,
    //GL_RGB, GL_UNSIGNED_BYTE, groundTexture->pix);

    displayListGroundWithTexture = glGenLists(1);
    glNewList(displayListGroundWithTexture, GL_COMPILE);
    buildTexturedDisplayList();
    glEndList();
  }
}

SceneGroundPlane::~SceneGroundPlane()
{
  if (displayListGround != 0)
    glDeleteLists(displayListGround, 1);
  if (displayListGroundWithTexture != 0)
  {
    glDeleteTextures(1, &groundTexName);
    glDeleteLists(displayListGroundWithTexture, 1);
  }
}

void SceneGroundPlane::buildDisplayList()
{
  glPushAttrib(GL_ENABLE_BIT|GL_LIGHTING_BIT|GL_POLYGON_BIT);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glEnable(GL_LIGHTING);
  glPolygonOffset(1.0,1.0);

  float blendFactor = 1.0f;
  float planeAmbient[4]  = {  ambientFactor * color[0],  ambientFactor * color[1],  ambientFactor * color[2], blendFactor};
  float planeDiffuse[4]  = {  diffuseFactor * color[0],  diffuseFactor * color[1],  diffuseFactor * color[2], blendFactor};
  float planeSpecular[4] = { specularFactor * color[0], specularFactor * color[1], specularFactor * color[2], blendFactor};
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, planeAmbient);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, planeDiffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, planeSpecular);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

  glNormal3f(0,1,0);

  glBegin(GL_QUADS);
  double planeIncrement = groundPlaneSize / planeResolution;
  for(int i=0; i<planeResolution; i++)
    for(int j=0; j<planeResolution; j++)
    {
      glVertex3f(-groundPlaneSize/2 + i * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + j * planeIncrement);
      glVertex3f(-groundPlaneSize/2 + i * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + (j+1) * planeIncrement);
      glVertex3f(-groundPlaneSize/2 + (i+1) * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + (j+1) * planeIncrement);
      glVertex3f(-groundPlaneSize/2 + (i+1) * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + j * planeIncrement);
    }
  glEnd();

  glPopAttrib();
}

void SceneGroundPlane::buildTexturedDisplayList()
{
  glPushAttrib(GL_ENABLE_BIT);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D, groundTexName);
  glNormal3f(0,1,0);
  const double wrapFactor = 8.0;
  double planeIncrement = groundPlaneSize / planeResolution;

  glBegin(GL_QUADS);
  for(int i=0; i<planeResolution; i++)
    for(int j=0; j<planeResolution; j++)
    {
      glTexCoord2f(wrapFactor * i / planeResolution, wrapFactor * j / planeResolution);
      glVertex3f(-groundPlaneSize/2 + i * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + j * planeIncrement);

      glTexCoord2f(wrapFactor * (i+1) / planeResolution, wrapFactor * j / planeResolution);
      glVertex3f(-groundPlaneSize/2 + (i+1) * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + j * planeIncrement);

      glTexCoord2f(wrapFactor * (i+1) / planeResolution, wrapFactor * (j+1) / planeResolution);
      glVertex3f(-groundPlaneSize/2 + (i+1) * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + (j+1) * planeIncrement);

      glTexCoord2f(wrapFactor * i / planeResolution, wrapFactor * (j+1) / planeResolution);
      glVertex3f(-groundPlaneSize/2 + i * planeIncrement, groundPlaneHeight, -groundPlaneSize/2 + (j+1) * planeIncrement);
    }
  glEnd();

  glPopAttrib();
}

void SceneGroundPlane::renderShadow(SceneObject * object)
{
  GLboolean enabled = GL_FALSE;
  glGetBooleanv(GL_TEXTURE_2D, &enabled);
  glDisable(GL_TEXTURE_2D);
//  double ground[4] = {0, 1, 0, -groundPlaneHeight - 0.01};
  double ground[4] = {0, 1, 0, -groundPlaneHeight};
  double light[4] = {lightPos[0], lightPos[1], lightPos[2], 1};
  object->RenderShadow(ground, light);
  if (enabled)
    glEnable(GL_TEXTURE_2D);
}

void SceneGroundPlane::render()
{
  if (displayListGroundWithTexture == 0)
    glCallList(displayListGround);
  else
    glCallList(displayListGroundWithTexture);
}
