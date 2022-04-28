/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "matrixIO" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC       *
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

/*
  Load/save multiple binary matrices from a single file.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "multiMatrixIO.h"

int MultiMatrixIO::Save(const char * filename, int numMatrices, int * m, int * n, double ** matrices)
{
  FILE * fout = fopen(filename, "wb");
  if (fout == nullptr)
  {
    printf("Error in MultiMatrixIO::Save: cannot open %s for writing.\n", filename);
    return 1;
  }

  // write the number of matrices
  fwrite(&numMatrices, sizeof(int), 1, fout);

  // write the dimensions of the matrices
  for(int i=0; i<numMatrices; i++)
  {
    fwrite(&m[i], sizeof(int), 1, fout);
    fwrite(&n[i], sizeof(int), 1, fout);
  }

  // write the matrices
  for(int i=0; i<numMatrices; i++)
  {
    fwrite(matrices[i], sizeof(double), m[i] * n[i], fout);
  }

  fclose(fout);

  return 0;
}

int MultiMatrixIO::Load(const char * filename, int * numMatrices, int ** m, int ** n, double *** matrices)
{
  *m = nullptr;
  *n = nullptr;
  *matrices = nullptr;

  FILE * fin = fopen(filename, "rb");
  if (fin == nullptr)
  {
    printf("Error in MultiMatrixIO::Load: cannot open %s for reading.\n", filename);
    return 1;
  }

  // read number of matrices
  unsigned int items = fread(numMatrices, sizeof(int), 1, fin);
  if (items != 1)
  {
    printf("Error in MultiMatrixIO::Load: mismatch in the number of matrices.\n");
    fclose(fin);
    return 1;
  }

  *m = (int*) malloc (sizeof(int) * (*numMatrices));
  *n = (int*) malloc (sizeof(int) * (*numMatrices));

  // read the matrix dimensions 
  for(int i=0; i<*numMatrices; i++)
  {
    items = fread(&((*m)[i]), sizeof(int), 1, fin);
    if (items != 1)
    {
      printf("Error in MultiMatrixIO::Load: mismatch in the number of matrices.\n");
      free(*m);
      *m = nullptr;
      free(*n);
      *n = nullptr;
      fclose(fin);
      return 1;
    }

    items = fread(&((*n)[i]), sizeof(int), 1, fin);
    if (items != 1)
    {
      printf("Error in MultiMatrixIO::Load: mismatch in the number of matrices.\n");
      free(*m);
      *m = nullptr;
      free(*n);
      *n = nullptr;
      fclose(fin);
      return 1;
    }
  }

  // allocate memory for matrices
  (*matrices) = (double **) malloc (sizeof(double*) * *numMatrices);
  for(int i=0; i<*numMatrices; i++)
  {
    int numEntries = (*m)[i] * (*n)[i];
    if (numEntries > 0)
      (*matrices)[i] = (double *) malloc(sizeof(double) * numEntries);
    else
      (*matrices)[i] = nullptr;
  }

  // read the each matrices
  for(int i=0; i<*numMatrices; i++)
  {
    unsigned int matrixSize = (*m)[i] * (*n)[i];
    items = fread((*matrices)[i], sizeof(double), matrixSize, fin);
    if (items != matrixSize)
    {
      printf("Error in MultiMatrixIO::Load: invalid number of bytes.\n");
      free(*m);
      *m = nullptr;
      free(*n);
      *n = nullptr;
      for(int i = 0; i < *numMatrices; i++) { free((*matrices)[i]); }
      free(*matrices);
      *matrices = nullptr;
      fclose(fin);
      return 1;
    }
  }
  fclose(fin);

  return 0;
}

