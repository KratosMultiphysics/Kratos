#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "globheads.h"
#include "protos.h"

#define TOL 1.e-17
#define max(a,b) (((a)>(b))?(a):(b))

int invGauss(int nn, double *A) {
  /* *-------------------- inversion by svd
     This calls lapack routines for inverting a dense matrix.
     dgetrf and dgetri

     ON ENTRY
     ** A = square matrix of size n x n -- dimensioned with
     ** leading dimension =  n

     ON RETURN A contains the inverse  of the input matrix.
  */
  int lWk, info;

  double *Wk;
  int *ipiv;

  lWk = 10*nn;

  /*-------------------- trivial case nn = 1                     */
  if (nn == 1) {
    if (A[0] == 0.0)
      return 1;
    else {
      A[0] = 1.0 / A[0];
      return 0;
    }
  }
  /*-------------------- general case                              */

  Wk  = (double *) malloc(lWk*sizeof(double));
  ipiv = (int *) malloc(nn*sizeof(int));
  if (Wk == NULL || ipiv == NULL)
    return -1;
  /*-------------------- get LU factorization with pivoting         */
  DGETRF (nn, nn, A, nn, ipiv, &info);
  if (info !=0 ) return info;
  /*-------------------- compute inverse                            */
  DGETRI (nn, A, nn, ipiv, Wk, lWk, &info);

  free(Wk);
  free(ipiv);
  return info;
}

int invSVD(int nn, double *A) {
  /* *-------------------- inversion by svd
     This calls lapack routine dgesvd --
     ON ENTRY
     ** A = square matrix of size n x n -- dimensioned with
     ** leading dimension =  n
     ON RETURN A contains the truncated SVD inverse of input matrix.
     ** tolerance set for truncation is TOL and can be changed in
     ** above define statement
     **--------------------
     */
  int lWk, i, info;

  double *U, *VT, *S, *Wk;
  double tmp, nrm, one=1.0, zero=0.0;

  double tol=TOL;

  lWk = 5*nn;

  U  = (double *) malloc(nn*nn*sizeof(double));
  VT = (double *) malloc(nn*nn*sizeof(double));
  S  = (double *) malloc(nn*sizeof(double));
  Wk  = (double *) malloc(lWk*sizeof(double));

  if (U == NULL || VT == NULL || S == NULL || Wk == NULL)
    return -1;
  /*-------------------- trivial case nn = 1                     */
  if (nn == 1) {
    if (A[0] == 0.0)
      return 1;
    else {
      A[0] = one / A[0];
      return 0;
    }
  }
  /*-------------------- general case                              */
  dgesvd ("A","A", &nn, &nn, A, &nn, S, U, &nn, VT, &nn, Wk, &lWk,
	  &info) ;
  if (S[0] == 0.0)
    return 1;
  nrm = S[0]*tol;
  /*-------------------- compute S\inv * VT                        */
  for (i=0; i<nn; i++) {
    tmp = one / max(S[i],nrm) ;
    DSCAL(nn, tmp, &VT[i], nn);
  }
  /*-------------------- do [V^T S\inv ] * U^T                     */
  dgemm("t","t",&nn,&nn,&nn, &one, VT, &nn, U, &nn, &zero, A, &nn);
  /*-------------------- Done -------------------------------------*/
  free(U);
  free(VT);
  free(S);
  free(Wk);
  return 0;
}

