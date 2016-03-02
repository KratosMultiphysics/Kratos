// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2009-2010 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "common.h"

int EIGEN_BLAS_FUNC(gemv)(char *opa, int *m, int *n, RealScalar *palpha, RealScalar *pa, int *lda, RealScalar *pb, int *incb, RealScalar *pbeta, RealScalar *pc, int *incc)
{
  typedef void (*functype)(int, int, const Scalar *, int, const Scalar *, int , Scalar *, int, Scalar);
  static functype func[4];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<4; ++k)
      func[k] = 0;

    func[NOTR] = (internal::general_matrix_vector_product<int,Scalar,ColMajor,false,Scalar,false>::run);
    func[TR  ] = (internal::general_matrix_vector_product<int,Scalar,RowMajor,false,Scalar,false>::run);
    func[ADJ ] = (internal::general_matrix_vector_product<int,Scalar,RowMajor,Conj, Scalar,false>::run);

    init = true;
  }

  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* b = reinterpret_cast<Scalar*>(pb);
  Scalar* c = reinterpret_cast<Scalar*>(pc);
  Scalar alpha  = *reinterpret_cast<Scalar*>(palpha);
  Scalar beta   = *reinterpret_cast<Scalar*>(pbeta);

  // check arguments
  int info = 0;
  if(OP(*opa)==INVALID)           info = 1;
  else if(*m<0)                   info = 2;
  else if(*n<0)                   info = 3;
  else if(*lda<std::max(1,*m))    info = 6;
  else if(*incb==0)               info = 8;
  else if(*incc==0)               info = 11;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"GEMV ",&info,6);

  if(*m==0 || *n==0 || (alpha==Scalar(0) && beta==Scalar(1)))
    return 0;

  int actual_m = *m;
  int actual_n = *n;
  int code = OP(*opa);
  if(code!=NOTR)
    std::swap(actual_m,actual_n);

  Scalar* actual_b = get_compact_vector(b,actual_n,*incb);
  Scalar* actual_c = get_compact_vector(c,actual_m,*incc);

  if(beta!=Scalar(1))
  {
    if(beta==Scalar(0)) vector(actual_c, actual_m).setZero();
    else                vector(actual_c, actual_m) *= beta;
  }

  if(code>=4 || func[code]==0)
    return 0;

  func[code](actual_m, actual_n, a, *lda, actual_b, 1, actual_c, 1, alpha);

  if(actual_b!=b) delete[] actual_b;
  if(actual_c!=c) delete[] copy_back(actual_c,c,actual_m,*incc);

  return 1;
}

int EIGEN_BLAS_FUNC(trsv)(char *uplo, char *opa, char *diag, int *n, RealScalar *pa, int *lda, RealScalar *pb, int *incb)
{
  typedef void (*functype)(int, const Scalar *, int, Scalar *);
  static functype func[16];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<16; ++k)
      func[k] = 0;

    func[NOTR  | (UP << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       false,ColMajor>::run);
    func[TR    | (UP << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       Conj, RowMajor>::run);

    func[NOTR  | (LO << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       false,ColMajor>::run);
    func[TR    | (LO << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (NUNIT << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       Conj, RowMajor>::run);

    func[NOTR  | (UP << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,false,ColMajor>::run);
    func[TR    | (UP << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,Conj, RowMajor>::run);

    func[NOTR  | (LO << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,false,ColMajor>::run);
    func[TR    | (LO << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (UNIT  << 3)] = (internal::triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,Conj, RowMajor>::run);

    init = true;
  }

  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* b = reinterpret_cast<Scalar*>(pb);

  int info = 0;
  if(UPLO(*uplo)==INVALID)                                            info = 1;
  else if(OP(*opa)==INVALID)                                          info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*lda<std::max(1,*n))                                        info = 6;
  else if(*incb==0)                                                   info = 8;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TRSV ",&info,6);

  Scalar* actual_b = get_compact_vector(b,*n,*incb);

  int code = OP(*opa) | (UPLO(*uplo) << 2) | (DIAG(*diag) << 3);
  func[code](*n, a, *lda, actual_b);

  if(actual_b!=b) delete[] copy_back(actual_b,b,*n,*incb);

  return 0;
}



int EIGEN_BLAS_FUNC(trmv)(char *uplo, char *opa, char *diag, int *n, RealScalar *pa, int *lda, RealScalar *pb, int *incb)
{
  typedef void (*functype)(int, int, const Scalar *, int, const Scalar *, int, Scalar *, int, const Scalar&);
  static functype func[16];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<16; ++k)
      func[k] = 0;

    func[NOTR  | (UP << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Upper|0,       Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (UP << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Lower|0,       Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Lower|0,       Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (LO << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Lower|0,       Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (LO << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Upper|0,       Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (NUNIT << 3)] = (internal::triangular_matrix_vector_product<int,Upper|0,       Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (UP << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (UP << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (LO << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (LO << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (UNIT  << 3)] = (internal::triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,Conj, Scalar,false,RowMajor>::run);

    init = true;
  }

  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* b = reinterpret_cast<Scalar*>(pb);

  int info = 0;
  if(UPLO(*uplo)==INVALID)                                            info = 1;
  else if(OP(*opa)==INVALID)                                          info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*lda<std::max(1,*n))                                        info = 6;
  else if(*incb==0)                                                   info = 8;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TRMV ",&info,6);

  if(*n==0)
    return 1;

  Scalar* actual_b = get_compact_vector(b,*n,*incb);
  Matrix<Scalar,Dynamic,1> res(*n);
  res.setZero();

  int code = OP(*opa) | (UPLO(*uplo) << 2) | (DIAG(*diag) << 3);
  if(code>=16 || func[code]==0)
    return 0;

  func[code](*n, *n, a, *lda, actual_b, 1, res.data(), 1, Scalar(1));

  copy_back(res.data(),b,*n,*incb);
  if(actual_b!=b) delete[] actual_b;

  return 1;
}

/**  GBMV  performs one of the matrix-vector operations
  *
  *     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
  *
  *  where alpha and beta are scalars, x and y are vectors and A is an
  *  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
  */
int EIGEN_BLAS_FUNC(gbmv)(char *trans, int *m, int *n, int *kl, int *ku, RealScalar *palpha, RealScalar *pa, int *lda,
                          RealScalar *px, int *incx, RealScalar *pbeta, RealScalar *py, int *incy)
{
  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* x = reinterpret_cast<Scalar*>(px);
  Scalar* y = reinterpret_cast<Scalar*>(py);
  Scalar alpha = *reinterpret_cast<Scalar*>(palpha);
  Scalar beta = *reinterpret_cast<Scalar*>(pbeta);
  int coeff_rows = *kl+*ku+1;
  
  int info = 0;
       if(OP(*trans)==INVALID)                                        info = 1;
  else if(*m<0)                                                       info = 2;
  else if(*n<0)                                                       info = 3;
  else if(*kl<0)                                                      info = 4;
  else if(*ku<0)                                                      info = 5;
  else if(*lda<coeff_rows)                                            info = 8;
  else if(*incx==0)                                                   info = 10;
  else if(*incy==0)                                                   info = 13;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"GBMV ",&info,6);
  
  if(*m==0 || *n==0 || (alpha==Scalar(0) && beta==Scalar(1)))
    return 0;
  
  int actual_m = *m;
  int actual_n = *n;
  if(OP(*trans)!=NOTR)
    std::swap(actual_m,actual_n);
  
  Scalar* actual_x = get_compact_vector(x,actual_n,*incx);
  Scalar* actual_y = get_compact_vector(y,actual_m,*incy);
  
  if(beta!=Scalar(1))
  {
    if(beta==Scalar(0)) vector(actual_y, actual_m).setZero();
    else                vector(actual_y, actual_m) *= beta;
  }
  
  MatrixType mat_coeffs(a,coeff_rows,*n,*lda);
  
  int nb = std::min(*n,(*m)+(*ku));
  for(int j=0; j<nb; ++j)
  {
    int start = std::max(0,j - *ku);
    int end = std::min((*m)-1,j + *kl);
    int len = end - start + 1;
    int offset = (*ku) - j + start;
    if(OP(*trans)==NOTR)
      vector(actual_y+start,len) += (alpha*actual_x[j]) * mat_coeffs.col(j).segment(offset,len);
    else if(OP(*trans)==TR)
      actual_y[j] += alpha * ( mat_coeffs.col(j).segment(offset,len).transpose() * vector(actual_x+start,len) ).value();
    else
      actual_y[j] += alpha * ( mat_coeffs.col(j).segment(offset,len).adjoint()   * vector(actual_x+start,len) ).value();
  }    
  
  if(actual_x!=x) delete[] actual_x;
  if(actual_y!=y) delete[] copy_back(actual_y,y,actual_m,*incy);
  
  return 0;
}

#if 0
/**  TBMV  performs one of the matrix-vector operations
  *
  *     x := A*x,   or   x := A'*x,
  *
  *  where x is an n element vector and  A is an n by n unit, or non-unit,
  *  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
  */
int EIGEN_BLAS_FUNC(tbmv)(char *uplo, char *opa, char *diag, int *n, int *k, RealScalar *pa, int *lda, RealScalar *px, int *incx)
{
  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* x = reinterpret_cast<Scalar*>(px);
  int coeff_rows = *k + 1;
  
  int info = 0;
       if(UPLO(*uplo)==INVALID)                                       info = 1;
  else if(OP(*opa)==INVALID)                                          info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*k<0)                                                       info = 5;
  else if(*lda<coeff_rows)                                            info = 7;
  else if(*incx==0)                                                   info = 9;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TBMV ",&info,6);
  
  if(*n==0)
    return 0;
  
  int actual_n = *n;
  
  Scalar* actual_x = get_compact_vector(x,actual_n,*incx);
  
  MatrixType mat_coeffs(a,coeff_rows,*n,*lda);
  
  int ku = UPLO(*uplo)==UPPER ? *k : 0;
  int kl = UPLO(*uplo)==LOWER ? *k : 0;
  
  for(int j=0; j<*n; ++j)
  {
    int start = std::max(0,j - ku);
    int end = std::min((*m)-1,j + kl);
    int len = end - start + 1;
    int offset = (ku) - j + start;
    
    if(OP(*trans)==NOTR)
      vector(actual_y+start,len) += (alpha*actual_x[j]) * mat_coeffs.col(j).segment(offset,len);
    else if(OP(*trans)==TR)
      actual_y[j] += alpha * ( mat_coeffs.col(j).segment(offset,len).transpose() * vector(actual_x+start,len) ).value();
    else
      actual_y[j] += alpha * ( mat_coeffs.col(j).segment(offset,len).adjoint()   * vector(actual_x+start,len) ).value();
  }    
  
  if(actual_x!=x) delete[] actual_x;
  if(actual_y!=y) delete[] copy_back(actual_y,y,actual_m,*incy);
  
  return 0;
}
#endif

/**  DTBSV  solves one of the systems of equations
  *
  *     A*x = b,   or   A'*x = b,
  *
  *  where b and x are n element vectors and A is an n by n unit, or
  *  non-unit, upper or lower triangular band matrix, with ( k + 1 )
  *  diagonals.
  *
  *  No test for singularity or near-singularity is included in this
  *  routine. Such tests must be performed before calling this routine.
  */
int EIGEN_BLAS_FUNC(tbsv)(char *uplo, char *op, char *diag, int *n, int *k, RealScalar *pa, int *lda, RealScalar *px, int *incx)
{
  typedef void (*functype)(int, int, const Scalar *, int, Scalar *);
  static functype func[16];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<16; ++k)
      func[k] = 0;

    func[NOTR  | (UP << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Upper|0,       Scalar,false,Scalar,ColMajor>::run);
    func[TR    | (UP << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Lower|0,       Scalar,false,Scalar,RowMajor>::run);
    func[ADJ   | (UP << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Lower|0,       Scalar,Conj, Scalar,RowMajor>::run);

    func[NOTR  | (LO << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Lower|0,       Scalar,false,Scalar,ColMajor>::run);
    func[TR    | (LO << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Upper|0,       Scalar,false,Scalar,RowMajor>::run);
    func[ADJ   | (LO << 2) | (NUNIT << 3)] = (internal::band_solve_triangular_selector<int,Upper|0,       Scalar,Conj, Scalar,RowMajor>::run);

    func[NOTR  | (UP << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Upper|UnitDiag,Scalar,false,Scalar,ColMajor>::run);
    func[TR    | (UP << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Lower|UnitDiag,Scalar,false,Scalar,RowMajor>::run);
    func[ADJ   | (UP << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Lower|UnitDiag,Scalar,Conj, Scalar,RowMajor>::run);

    func[NOTR  | (LO << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Lower|UnitDiag,Scalar,false,Scalar,ColMajor>::run);
    func[TR    | (LO << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Upper|UnitDiag,Scalar,false,Scalar,RowMajor>::run);
    func[ADJ   | (LO << 2) | (UNIT  << 3)] = (internal::band_solve_triangular_selector<int,Upper|UnitDiag,Scalar,Conj, Scalar,RowMajor>::run);

    init = true;
  }

  Scalar* a = reinterpret_cast<Scalar*>(pa);
  Scalar* x = reinterpret_cast<Scalar*>(px);
  int coeff_rows = *k+1;
  
  int info = 0;
       if(UPLO(*uplo)==INVALID)                                       info = 1;
  else if(OP(*op)==INVALID)                                           info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*k<0)                                                       info = 5;
  else if(*lda<coeff_rows)                                            info = 7;
  else if(*incx==0)                                                   info = 9;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TBSV ",&info,6);
  
  if(*n==0 || (*k==0 && DIAG(*diag)==UNIT))
    return 0;
  
  int actual_n = *n;
 
  Scalar* actual_x = get_compact_vector(x,actual_n,*incx);
  
  int code = OP(*op) | (UPLO(*uplo) << 2) | (DIAG(*diag) << 3);
  if(code>=16 || func[code]==0)
    return 0;

  func[code](*n, *k, a, *lda, actual_x);
  
  if(actual_x!=x) delete[] copy_back(actual_x,x,actual_n,*incx);
  
  return 0;
}

/**  DTPMV  performs one of the matrix-vector operations
  *
  *     x := A*x,   or   x := A'*x,
  *
  *  where x is an n element vector and  A is an n by n unit, or non-unit,
  *  upper or lower triangular matrix, supplied in packed form.
  */
int EIGEN_BLAS_FUNC(tpmv)(char *uplo, char *opa, char *diag, int *n, RealScalar *pap, RealScalar *px, int *incx)
{
  typedef void (*functype)(int, const Scalar*, const Scalar*, Scalar*, Scalar);
  static functype func[16];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<16; ++k)
      func[k] = 0;

    func[NOTR  | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|0,       Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|0,       Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|0,       Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|0,       Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|0,       Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|0,       Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,Conj, Scalar,false,RowMajor>::run);

    func[NOTR  | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Lower|UnitDiag,Scalar,false,Scalar,false,ColMajor>::run);
    func[TR    | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,false,Scalar,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_matrix_vector_product<int,Upper|UnitDiag,Scalar,Conj, Scalar,false,RowMajor>::run);

    init = true;
  }

  Scalar* ap = reinterpret_cast<Scalar*>(pap);
  Scalar* x = reinterpret_cast<Scalar*>(px);

  int info = 0;
  if(UPLO(*uplo)==INVALID)                                            info = 1;
  else if(OP(*opa)==INVALID)                                          info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*incx==0)                                                   info = 7;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TPMV ",&info,6);

  if(*n==0)
    return 1;

  Scalar* actual_x = get_compact_vector(x,*n,*incx);
  Matrix<Scalar,Dynamic,1> res(*n);
  res.setZero();

  int code = OP(*opa) | (UPLO(*uplo) << 2) | (DIAG(*diag) << 3);
  if(code>=16 || func[code]==0)
    return 0;

  func[code](*n, ap, actual_x, res.data(), Scalar(1));

  copy_back(res.data(),x,*n,*incx);
  if(actual_x!=x) delete[] actual_x;

  return 1;
}

/**  DTPSV  solves one of the systems of equations
  *
  *     A*x = b,   or   A'*x = b,
  *
  *  where b and x are n element vectors and A is an n by n unit, or
  *  non-unit, upper or lower triangular matrix, supplied in packed form.
  *
  *  No test for singularity or near-singularity is included in this
  *  routine. Such tests must be performed before calling this routine.
  */
int EIGEN_BLAS_FUNC(tpsv)(char *uplo, char *opa, char *diag, int *n, RealScalar *pap, RealScalar *px, int *incx)
{
  typedef void (*functype)(int, const Scalar*, Scalar*);
  static functype func[16];

  static bool init = false;
  if(!init)
  {
    for(int k=0; k<16; ++k)
      func[k] = 0;

    func[NOTR  | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       false,ColMajor>::run);
    func[TR    | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       Conj, RowMajor>::run);

    func[NOTR  | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|0,       false,ColMajor>::run);
    func[TR    | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (NUNIT << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|0,       Conj, RowMajor>::run);

    func[NOTR  | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,false,ColMajor>::run);
    func[TR    | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,false,RowMajor>::run);
    func[ADJ   | (UP << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,Conj, RowMajor>::run);

    func[NOTR  | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Lower|UnitDiag,false,ColMajor>::run);
    func[TR    | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,false,RowMajor>::run);
    func[ADJ   | (LO << 2) | (UNIT  << 3)] = (internal::packed_triangular_solve_vector<Scalar,Scalar,int,OnTheLeft, Upper|UnitDiag,Conj, RowMajor>::run);

    init = true;
  }

  Scalar* ap = reinterpret_cast<Scalar*>(pap);
  Scalar* x = reinterpret_cast<Scalar*>(px);

  int info = 0;
  if(UPLO(*uplo)==INVALID)                                            info = 1;
  else if(OP(*opa)==INVALID)                                          info = 2;
  else if(DIAG(*diag)==INVALID)                                       info = 3;
  else if(*n<0)                                                       info = 4;
  else if(*incx==0)                                                   info = 7;
  if(info)
    return xerbla_(SCALAR_SUFFIX_UP"TPSV ",&info,6);

  Scalar* actual_x = get_compact_vector(x,*n,*incx);

  int code = OP(*opa) | (UPLO(*uplo) << 2) | (DIAG(*diag) << 3);
  func[code](*n, ap, actual_x);

  if(actual_x!=x) delete[] copy_back(actual_x,x,*n,*incx);

  return 1;
}

