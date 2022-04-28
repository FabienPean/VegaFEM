/*************************************************************************
 *                                                                       *
 * Vega FEM Simulation Library Version 4.0                               *
 *                                                                       *
 * "matrix" library , Copyright (C) 2007 CMU, 2009 MIT, 2018 USC         *
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

#include <math.h>
#include "matrix.h"
#include "matrixIO.h"
#include "matrixExp.h"
#include <cassert>
#include <algorithm>
using namespace std;

template<class real>
Matrix<real>::Matrix(const char * filename)
{
  if (ReadMatrixFromDisk(filename, &m, &n, &data) != 0)
  {
    printf("Error loading matrix from %s.\n", filename);
    throw 1;
  }
  freeDataInDestructor = true;
}

template<class real>
Matrix<real>::Matrix(int m_, int n_, const real * data_, 
                     bool makeInternalDataCopy, bool freeDataInDestructor_):
  m(m_), n(n_), freeDataInDestructor(freeDataInDestructor_)
{
  if (makeInternalDataCopy)
  {
    data = (real*) malloc (sizeof(real) * m * n);
    memcpy(data, data_, sizeof(real) * m * n);
  }
  else
    data = (real*)data_;
}

template<class real>
Matrix<real>::Matrix (int m_, int n_, bool freeDataInDestructor_):
  m(m_), n(n_), freeDataInDestructor(freeDataInDestructor_)
{
  data = (real*) calloc (m * n, sizeof(real));
}

template<class real>
Matrix<real>::Matrix (int m_, int n_, real constEntry, bool freeDataInDestructor_):
  m(m_), n(n_), freeDataInDestructor(freeDataInDestructor_)
{
  data = (real*) malloc (sizeof(real) * m * n);
  for(int i=0; i<m*n; i++)
    data[i] = constEntry;
}

template<class real>
Matrix<real>::Matrix(int m_, real diagonal,
         bool freeDataInDestructor_):
  m(m_), n(m_), freeDataInDestructor(freeDataInDestructor_)
{
  data = (real*) calloc (m*m, sizeof(real));
  for(int i=0; i<m; i++)
    data[ELT(m,i,i)] = diagonal;
}

template<class real>
Matrix<real>::Matrix (int m_, const Matrix & vec, bool freeDataInDestructor_):
  m(m_), n(m_), freeDataInDestructor(freeDataInDestructor_)
{
  if (m != vec.Getm())
    throw 101;
  if (vec.Getn() != 1)
    throw 102;
  data = (real*) calloc (m*m, sizeof(real));
  for(int i=0; i<m; i++)
    data[ELT(m,i,i)] = vec(i,0);
}

template<class real>
Matrix<real>::Matrix(int m_, const real * diagonal,
         bool freeDataInDestructor_):
  m(m_), n(m_), freeDataInDestructor(freeDataInDestructor_)
{
  data = (real*) calloc (m*m, sizeof(real));
  for(int i=0; i<m; i++)
    data[ELT(m,i,i)] = diagonal[i];
}

template<class real>
Matrix<real>::Matrix(const Matrix<real> & mtx2)
{
  m = mtx2.Getm();
  n = mtx2.Getn();
  data = (real*) malloc (sizeof(real) * m * n);
  memcpy(data, mtx2.GetData(), sizeof(real) * m * n);
  freeDataInDestructor = true;
}

template<class real>
Matrix<real>::~Matrix()
{
  if (freeDataInDestructor)
    free(data);
}

template<class real>
const Matrix<real> Matrix<real>::operator+ (const Matrix<real> & mtx2) const
{
  if ((m != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator+ . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 1;
  }
  real * output = SumMatrices(m, n, data, mtx2.GetData(), (real*)NULL);
  return Matrix<real>(m,n,output,false);
}

template<class real>
const Matrix<real> Matrix<real>::operator- (const Matrix<real> & mtx2) const
{
  if ((m != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator- . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 2;
  }
  real * output = SubtractMatrices(m, n, data, mtx2.GetData(), (real*)NULL);
  return Matrix<real>(m,n,output,false);
}

template<class real>
const Matrix<real> Matrix<real>::operator* (const Matrix<real> & mtx2) const
{
  if (n != mtx2.Getm())
  {
    printf("Matrix size mismatch in Matrix::operator* . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 3;
  }
  real * output = MultiplyMatrices(m, n, mtx2.Getn(), data, mtx2.GetData(), (real*)NULL);
  return Matrix<real>(m,mtx2.Getn(),output,false);
}

// output = trans(this) * mtx2  
template<class real>
const Matrix<real> Matrix<real>::MultiplyT(const Matrix & mtx2) const
{
  if (m != mtx2.Getm())
  {
    printf("Matrix size mismatch in Matrix::MultiplyT . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 31;
  }
  real * output = MultiplyMatricesT(n, m, mtx2.Getn(), data, mtx2.GetData(), (real*)NULL);
  return Matrix<real>(n,mtx2.Getn(),output,false);
}

template<class real>
Matrix<real> & Matrix<real>::operator= (const Matrix<real> & mtx2)
{
  if ((m != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator= . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 4;
  }
 
  if (this == &mtx2) // self-assignment
    return (*this);

  memcpy(data, mtx2.GetData(), sizeof(real) * m * n);
  return (*this);
}

template<class real>
Matrix<real> & Matrix<real>::operator+= (const Matrix<real> & mtx2)
{
  if ((m != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator+= . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 5;
  }
  // works fine even if self-assignment
  SumMatrices(m, n, data, mtx2.GetData(), data);
  return (*this);
}

template<class real>
Matrix<real> & Matrix<real>::operator-= (const Matrix<real> & mtx2)
{
  if ((m != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator-= . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 6;
  }
  // works fine even if self-assignment
  SubtractMatrices(m, n, data, mtx2.GetData(), data);
  return (*this);
}

template<class real>
Matrix<real> & Matrix<real>::operator*= (const Matrix<real> & mtx2)
{
  if ((n != mtx2.Getm()) || (n != mtx2.Getn()))
  {
    printf("Matrix size mismatch in Matrix::operator*= . m = %d, n = %d, mtx2.m = %d, mtx2.n = %d\n",
      m, n, mtx2.Getm(), mtx2.Getn());
    throw 7;
  }
  // works fine even if self-assignment
  real * buffer = MultiplyMatrices(m, n, n, data, mtx2.GetData());
  memcpy(data, buffer, sizeof(real) * m * n);
  free(buffer);
  return (*this);
}

template<class real>
Matrix<real> & Matrix<real>::operator*= (const real alpha)
{
  int mn = m * n;
  for(int i=0; i<mn; i++)
    data[i] *= alpha;
  return (*this);
}

template<class real>
bool Matrix<real>::operator == (const Matrix<real> & mtx2) const
{
  if ((n != mtx2.Getm()) || (n != mtx2.Getn()))
    return false;
  int mn = m * n;
  for(int i = 0; i < mn; i++)
    if (data[i] != mtx2.data[i])
      return false;
  return true;
}

template<class real>
bool Matrix<real>::operator != (const Matrix<real> & mtx2) const
{
  return !(*this == mtx2);
}

template<class real>
real Matrix<real>::MaxAbsEntry() const
{
  real maxAbsEntry = 0.0;
  int mn = m * n;
  for(int i=0; i<mn; i++)
  {
    if (fabs(data[i]) > maxAbsEntry)
      maxAbsEntry = fabs(data[i]);
  }
  return maxAbsEntry;
}

template<class real>
real Matrix<real>::GetFrobeniusNorm() const
{
  return VectorNorm(m*n, data);
}


template<class real>
Matrix<real> & Matrix<real>::InPlaceTranspose()
{
  InPlaceTransposeMatrix(m, n, data);
  // swap m,n
  int swap = m;
  m = n;
  n = swap;
  return (*this);
}

template<class real>
void Matrix<real>::Print(int numDigits) const
{
  char formatString[6];
  if (numDigits < 0)
    strcpy(formatString,"%G ");
  else
    sprintf(formatString,"%%.%df ", numDigits);    

  for(int i=0; i<m; i++)
  {
    for(int j=0; j<n; j++)
      printf(formatString, (*this)(i,j));
    printf("\n");
  }
}

template<class real>
void Matrix<real>::SymmetricEigenDecomposition(Matrix<real> & Q, Matrix<real> & Lambda)
{    
  if (m != n)
  {
    printf("Matrix size mismatch in Matrix::SymmetricEigenDecomposition . m = %d, n = %d\n", m, n);
    throw 12;
  }    

  real * QData = (real*) malloc(sizeof(real) * m * m);
  real * LambdaData = (real*) malloc(sizeof(real) * m);
  SymmetricMatrixEigenDecomposition(m, data, QData, LambdaData);
  Q = Matrix<real>(m,m,QData,false);
  Lambda = Matrix<real>(m,1,LambdaData,false);
}

template<class real>
void Matrix<real>::SVD(Matrix<real> & U, Matrix<real> & Sigma, Matrix<real> & VT)
{    
  real * UData = (real*) malloc(sizeof(real) * m * MIN(m,n));
  real * SigmaData = (real*) malloc(sizeof(real) * MIN(m,n));
  real * VTData = (real*) malloc(sizeof(real) * MIN(m,n) * n);
  MatrixSVD(m, n, data, UData, SigmaData, VTData);
  U = Matrix<real>(m, MIN(m,n), UData, false);
  Sigma = Matrix<real>(MIN(m,n), 1, SigmaData, false);
  VT = Matrix<real>(MIN(m,n), n, VTData, false);
}

template<class real>
int Matrix<real>::EigenDecomposition(Matrix<real> & EigenVectors, Matrix<real> & LambdaRe, Matrix<real> & LambdaIm)
{    
  if (m != n)
  {
    printf("Matrix size mismatch in Matrix::EigenDecomposition . m = %d, n = %d\n", m, n);
    throw 13;
  }    

  real * EigenVectorsData = (real*) malloc(sizeof(real) * m * m);
  real * LambdaReData = (real*) malloc(sizeof(real) * m);
  real * LambdaImData = (real*) malloc(sizeof(real) * m);
  int code = MatrixEigenDecomposition(m, data, 
    EigenVectorsData, LambdaReData, LambdaImData);

  EigenVectors = Matrix<real>(m,m,EigenVectorsData,false);
  LambdaRe = Matrix<real>(m,1,LambdaReData,false);
  LambdaIm = Matrix<real>(m,1,LambdaImData,false);

  return code;
}

template<class real>
int Matrix<real>::LUSolve(Matrix<real> & x, const Matrix<real> & rhs)
{
  if ((Getm() != Getn()) || (Getm() != x.Getm()) || (x.Getm() != rhs.Getm()) || (x.Getn() != rhs.Getn()) )
  {    
    printf("Matrix size mismatch in Matrix::LUSolve .\n");
    throw 12;
  }
  
  int exitCode = 0;
  try
  {
    MatrixLUSolve(Getn(), rhs.Getn(), GetData(), x.GetData(), rhs.GetData());
  }
  catch(int exceptionCode)
  {
    exitCode = exceptionCode;
  }

  return exitCode;
}

#ifdef USE_EXPOKIT

template<class real>
void Matrix<real>::MExpv(real t, const Matrix<real> & v, Matrix<real> & w)
{
  if ((Getm() != Getn()) || (Getm() != v.Getm()) || (Getm() != w.Getm()) || (v.Getn() != 1) || (w.Getn() != 1))
  {
    printf("Matrix size mismatch in Matrix::MExpv .\n");
    throw 13;
  }

  MatrixExpv(Getm(), GetData(), t, v.GetData(), w.GetData());
}

#else

template<class real>
void Matrix<real>::MExpv(real t, const Matrix<real> & v, Matrix<real> & w)
{
  printf("Error: MExpv is not enabled.\n");
}

#endif



template<class real>
int Matrix<real>::Save(const char * filename) const
{
  if (WriteMatrixToDisk(filename, m, n, data) != 0)
  {
    printf("Error writing matrix to %s.\n", filename);
    return 1;
  }
  return 0;
}

template<class real>
int Matrix<real>::Load(const char * filename)
{
  if (freeDataInDestructor)
    free(data);
  data = NULL;
  m = n = 0;
  freeDataInDestructor = true;
  if (ReadMatrixFromDisk(filename, &m, &n, &data) != 0)
  {
    printf("Error loading matrix from %s.\n", filename);
    m = n = 0;
    data = NULL;
    return 1;
  }
  return 0;
}


template<class real>
void Matrix<real>::SetSubmatrix(int I, int J, const Matrix<real> & submatrix)
{
  int subm = submatrix.Getm();
  int subn = submatrix.Getn();
  const real * subdata = submatrix.GetData();

  if ((I < 0) || (J < 0) || (I + subm > m) || (J + subn > n))
  {
    printf("Error: matrix index out of bounds.\n");
    throw 21;
  }

  for(int j=0; j<subn; j++)
    for(int i=0; i<subm; i++)
      data[ELT(m,I+i,J+j)] = subdata[ELT(subm,i,j)];
}

template<class real>
void Matrix<real>::GetSubmatrix(int I, int J, Matrix<real> & submatrix) const
{
  int subm = submatrix.Getm();
  int subn = submatrix.Getn();
  real * subdata = submatrix.GetData();

  if ((I < 0) || (J < 0) || (I + subm > m) || (J + subn > n))
  {
    printf("Error: matrix index out of bounds.\n");
    throw 21;
  }

  for(int j=0; j<subn; j++)
    for(int i=0; i<subm; i++)
      subdata[ELT(subm,i,j)] = data[ELT(m,I+i,J+j)];
}

template<class real>
void Matrix<real>::RemoveColumns(int columnStart, int columnEnd)
{
  // write everything to the right of columnEnd to the left
  int stride = columnEnd - columnStart;
  for(int column=columnEnd; column<n; column++)
  {
    // write column to column-stride
    memcpy(&data[ELT(m, 0, column-stride)], &data[ELT(m, 0, column)], sizeof(real) * m);
  }

  // free the space
  data = (real*) realloc (data, sizeof(real) * m * (n - stride));
 
  n = n - stride;
}

template<class real>
void Matrix<real>::RemoveRows(int rowStart, int rowEnd)
{
  assert(0 <= rowStart && rowStart <= rowEnd && rowEnd <= m);
  int mNew = m - (rowEnd - rowStart);
  // for i-th column, move consecutive data block from entry (rowEnd, i) to (rowStart, i+1) to the 
  // target location in the resultant matrix
  // to avoid potential memory block overlap, we use memmove instead of memcpy; the latter is not safe
  for(int i = 0, num = n-1; i < num; i++)
    memmove(data + i*mNew + rowStart, data + i*m + rowEnd, sizeof(real) * mNew);

  // the last column is a special case because the consecutive data block is from (rowEnd, n-1) to (m-1, n-1)
  if (n >= 1)
    memmove(data + (n-1)*mNew + rowStart, data + (n-1)*m + rowEnd, sizeof(real) * (m - rowEnd));
  data = (real*) realloc (data, sizeof(real) * mNew * n);
  m = mNew;
}

template<class real>
void Matrix<real>::RemoveRowsColumns(int start, int end)
{
  assert(0 <= start && start <= end);
  assert(end <= m && end <= n);

  int stride = end - start;
  int mNew = m - stride;
  int nNew = n - stride;

  // similar method as in RemoveRows()
  // for i-th column, move consecutive data block from entry (rowEnd, i) to (rowStart, i+1) to the 
  // target location in the resultant matrix, until the (start-1) column
  for(int i = 0, num = start-1; i < num; i++)
    memmove(data + i*mNew + start, data + i*m + end, sizeof(real) * mNew);

  // the (start-1) column is a special case because the consecutive data block is from (rowEnd, start-1)
  // to (m-1, start-1)
  if (start >= 1)
    memmove(data + (start-1)*mNew + start, data + (start-1)*m + end, sizeof(real) * (m - end));

  // we skip columns from start to end-1
  // if we have columns after end
  if (end < n)
  {
    // move data from (end, 0) to (end, start-1) to the tagert location
    memmove(data + start*mNew, data + end*m, sizeof(real) * start);

    // fot the i-th column after the end column, move the consecutive data block as well
    for(int i = start, num = nNew-1; i < num; i++)
      memmove(data + i*mNew + start, data + (i+stride)*m + end, sizeof(real) * mNew);
  
    // process the special case at the last column
    if (nNew >= 1 && n >= 1)
      memmove(data + (nNew-1)*mNew+start, data + (n-1)*m + end, sizeof(real) * (m - end));
  }

  data = (real*) realloc (data, sizeof(real) * mNew * nNew);
  m = mNew;
  n = nNew;

//  RemoveColumns(start, end);
//  RemoveRows(start, end);
}

template<class real>
void Matrix<real>::AppendColumns(const Matrix<real> & columns)
{
  if (columns.Getm() != m)
  {
    printf("Error: mismatch in number of rows in AppendColumns.\n");
    throw 41;
  }

  data = (real*) realloc (data, sizeof(real) * (m * (n + columns.Getn())));
  memcpy(&data[ELT(m, 0, n)], columns.GetData(), sizeof(real) * m * columns.Getn());
  n += columns.Getn();
}

template<class real>
void Matrix<real>::AppendRows(const Matrix<real> & rows)
{
  if (rows.n != n)
  {
    printf("Error: mismatch in number of columns in AppendRows.\n");
    throw 42;
  }

  int oldm = m;
  Resize(0,0, m + rows.m, n);
  for(int i = 0; i < n; i++)
    memcpy(data + i*m + oldm, rows.data + i*rows.m, sizeof(real) * rows.m);
}

template<class real>
void Matrix<real>::AppendRowsColumns(const Matrix<real> & bottomLeftBlock, const Matrix<real> & topRightBlock, const Matrix<real> & bottomRightBlock)
{
  if (topRightBlock.Getm() != m)
  {
    printf("Error: mismatch in AppendRowsColumns.\n");
    throw 43;
  }

  if (bottomLeftBlock.Getn() != n)
  {
    printf("Error: mismatch in AppendRowsColumns.\n");
    throw 44;
  }

  if (bottomRightBlock.Getm() != bottomLeftBlock.Getm())
  {
    printf("Error: mismatch in AppendRowsColumns.\n");
    throw 45;
  }

  if (bottomRightBlock.Getn() != topRightBlock.Getn())
  {
    printf("Error: mismatch in AppendRowsColumns.\n");
    throw 46;
  }

  AppendColumns(topRightBlock);
 
  Matrix<real> blockMatrix(bottomLeftBlock.Getm(), n);
  real * blockMatrixData = blockMatrix.GetData();
  memcpy(blockMatrixData, bottomLeftBlock.GetData(), sizeof(real) * (bottomLeftBlock.Getm() * bottomLeftBlock.Getn()));
  memcpy(&blockMatrixData[ELT(bottomLeftBlock.Getm(), 0, bottomLeftBlock.Getn())], bottomRightBlock.GetData(), sizeof(real) * (bottomRightBlock.Getm() * bottomRightBlock.Getn()));
  AppendRows(blockMatrix);
}

template<class real>
void Matrix<real>::AppendRowsColumns(int numAddedRows, int numAddedCols)
{
  assert(numAddedRows >= 0 && numAddedCols >= 0);
  int newRows = m + numAddedRows;
  int newCols = n + numAddedCols;

  real * newdata = (real*) calloc(newRows * newCols, sizeof(real));

  for(int col = 0; col < n; col++)
    memcpy(&newdata[newRows*col], &data[m*col], sizeof(real) * m);

  free(data);
  data = newdata;
  m = newRows;
  n = newCols;
}

template<class real>
void Matrix<real>::Resize(int newRows, int newCols)
{
  assert(newRows >= 0 && newCols >= 0);

  real * newdata = (real*) calloc(newRows * newCols, sizeof(real));
  free(data);
  data = newdata;
  m = newRows;
  n = newCols;
}

template<class real>
void Matrix<real>::Resize(int i, int j, int newRows, int newCols)
{
  assert(newRows >= 0 && newCols >= 0);

  real * newdata = (real*) calloc(newRows * newCols, sizeof(real));
  int colStartToCopy = max(j, 0);
  int numColsToCopy = min(j+newCols, n) - colStartToCopy;
  int rowStartToCopy = max(i, 0);
  int numRowsToCopy = min(i+newRows, m) - rowStartToCopy;

  for(int col = 0; col < numColsToCopy; col++)
    memcpy(&newdata[newRows*col], &data[m*col + colStartToCopy + rowStartToCopy], sizeof(real) * numRowsToCopy);

  free(data);
  data = newdata;
  m = newRows;
  n = newCols;
}

template class Matrix<float>;
template class Matrix<double>;
