/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "interpolationCoordinates" library , Copyright (C) 2018 USC           *
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

#ifndef VEGAFEM_INTERPOLATIONCOORDINATES_H
#define VEGAFEM_INTERPOLATIONCOORDINATES_H

class InterpolationCoordinates
{
public:
  virtual ~InterpolationCoordinates() {}

  // deform embedded geometry
  // embedderDisp: vertex displacements of the embedder mesh
  // embeddedDisp: vertex displacements of the embedded mesh
  virtual void deform(const double * embedderDisp, double * embeddedDisp) const = 0;

protected:
  InterpolationCoordinates() {}
};

#endif

