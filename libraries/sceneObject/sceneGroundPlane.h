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

#ifndef SCENEGROUNDPLANE_H
#define SCENEGROUNDPLANE_H

#include <cstddef>
#include "sceneObject.h"

// render ground plane and shadow for SceneObject

class SceneGroundPlane
{
public:
  // groundPlaneString format:
  // "<ground height>,<light position X,Y,Z>,<ground size>,r<red>,g<green>,b<blue>,a<ambient factor>,d<diffuse factor>,s<specular factor>,sh<shininess>"
  // where <red>,<green>,<blue> is the ground color, ranging from 0.0 to 1.0
  //       <ambient/diffuse/specular factor> scales the ground color for ambient/diffuse/specular component, ranging from 0.0 to 1.0
  // groundPlaneTextureFilename: input texture; if it's NULL, use the color specified in groundPlaneString
  SceneGroundPlane(const char * groundPlaneString, const char * groundPlaneTextureFilename = NULL);
  virtual ~SceneGroundPlane();

  void renderShadow(SceneObject * object);

  void render();

protected:
  void buildDisplayList();
  void buildTexturedDisplayList();

  double groundPlaneHeight;
  Vec3d lightPos;
  double groundPlaneSize;
  float color[3];
  float ambientFactor, diffuseFactor, specularFactor;
  float shininess;
  GLuint displayListGround;
  GLuint displayListGroundWithTexture;
  int planeResolution;
  unsigned int groundTexName;
};



#endif
