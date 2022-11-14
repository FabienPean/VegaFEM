/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

#ifndef VEGAFEM_EIG3_H
#define VEGAFEM_EIG3_H
namespace vegafem
{
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */
void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
}//namespace vegafem
#endif
