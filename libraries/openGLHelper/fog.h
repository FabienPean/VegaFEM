/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "openGLHelper" library , Copyright (C) 2018 USC                       *
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

#ifndef FOG_H
#define FOG_H

#include "openGL-headers.h"

// a simple struct for rendering fog using openGL
struct Fog
{
  Fog() : fogStart(10.0),  fogEnd(14.0), fogDensity(0.25) {}

  // set openGL parameters for fog
  inline void setParameters() const;

  void enable() const { glEnable(GL_FOG); }

  void disable() const { glDisable(GL_FOG); }

  double fogStart, fogEnd, fogDensity;
};


inline void Fog::setParameters() const
{
  GLfloat fogColor[4] = {1.0, 1.0, 1.0, 1.0};
  glFogfv(GL_FOG_COLOR, fogColor);
  glFogf(GL_FOG_START, fogStart);
  glFogf(GL_FOG_END, fogEnd);
  glFogf(GL_FOG_DENSITY, fogDensity);
  glFogi(GL_FOG_MODE, GL_LINEAR);
}

#endif
