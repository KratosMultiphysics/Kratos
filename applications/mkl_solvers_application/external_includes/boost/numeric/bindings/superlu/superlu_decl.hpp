// (-*- C++ -*- header)

/*
 *
 * Copyright (c) Kresimir Fresl 2003
 * Copyright (c) Georg Baum 2004
 *
 * Permission to copy, modify, use and distribute this software
 * for any non-commercial or commercial purpose is granted provided
 * that this license appear on all copies of the software source code.
 *
 * Author assumes no responsibility whatsoever for its use and makes
 * no guarantees about its quality, correctness or reliability.
 *
 * Author acknowledges the support of the Faculty of Civil Engineering,
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_SUPERLU_DECL_H
#define BOOST_NUMERIC_BINDINGS_SUPERLU_DECL_H

#include <boost/numeric/bindings/superlu/superlu_inc.hpp>

extern "C" {

  // float versions (from superlu/ssp_defs.h)

  void sgssv (superlu_options_t *, SuperMatrix* A, int* perm_c, int* perm_r,
              SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t *, int* info);

  void sCreate_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               float* nzval, int* rowind, int* colptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void sCreate_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               float* nzval, int* colind, int* rowptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void sCreate_Dense_Matrix (SuperMatrix* A, int m, int n,
                             float* vals, int ld_vals,
                             Stype_t stype, Dtype_t dtype, Mtype_t mtype);

  void sgstrf (superlu_options_t*, SuperMatrix*, float, int, int, int*, void*, int, int*, int*,
               SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int*);
  void sgstrs (trans_t, SuperMatrix*, SuperMatrix*, int*, int*,
               SuperMatrix*, SuperLUStat_t*, int*);

  void sPrint_CompCol_Matrix (char* what, SuperMatrix* A);
  void sPrint_CompRow_Matrix (char* what, SuperMatrix* A);
  void sPrint_SuperNode_Matrix (char* what, SuperMatrix* A);
  void sPrint_Dense_Matrix (char* what, SuperMatrix* A);


  // double versions (from superlu/dsp_defs.h)

  void dgssv (superlu_options_t *, SuperMatrix* A, int* perm_c, int* perm_r,
              SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t *, int* info);

  void dCreate_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               double* nzval, int* rowind, int* colptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void dCreate_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               double* nzval, int* colind, int* rowptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void dCreate_Dense_Matrix (SuperMatrix* A, int m, int n,
                             double* vals, int ld_vals,
                             Stype_t stype, Dtype_t dtype, Mtype_t mtype);

  void dgstrf (superlu_options_t*, SuperMatrix*, double, int, int, int*, void*, int, int*, int*,
               SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int*);
  void dgstrs (trans_t, SuperMatrix*, SuperMatrix*, int*, int*,
               SuperMatrix*, SuperLUStat_t*, int*);

  void dPrint_CompCol_Matrix (char* what, SuperMatrix* A);
  void dPrint_CompRow_Matrix (char* what, SuperMatrix* A);
  void dPrint_SuperNode_Matrix (char* what, SuperMatrix* A);
  void dPrint_Dense_Matrix (char* what, SuperMatrix* A);


  // complex float versions (from superlu/csp_defs.h)

  typedef struct { float r, i; } flt_complex;

  void cgssv (superlu_options_t *, SuperMatrix* A, int* perm_c, int* perm_r,
              SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t *, int* info);

  void cCreate_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               flt_complex* nzval, int* rowind, int* colptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void cCreate_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               flt_complex* nzval, int* colind, int* rowptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void cCreate_Dense_Matrix (SuperMatrix* A, int m, int n,
                             flt_complex* vals, int ld_vals,
                             Stype_t stype, Dtype_t dtype, Mtype_t mtype);

  void cgstrf (superlu_options_t*, SuperMatrix*, float, int, int, int*, void*, int, int*, int*,
               SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int*);
  void cgstrs (trans_t, SuperMatrix*, SuperMatrix*, int*, int*,
               SuperMatrix*, SuperLUStat_t*, int *);

  void cPrint_CompCol_Matrix (char* what, SuperMatrix* A);
  void cPrint_CompRow_Matrix (char* what, SuperMatrix* A);
  void cPrint_SuperNode_Matrix (char* what, SuperMatrix* A);
  void cPrint_Dense_Matrix (char* what, SuperMatrix* A);


  // complex float versions (from superlu/zsp_defs.h)

  typedef struct { double r, i; } dbl_complex;

  void zgssv (superlu_options_t *, SuperMatrix* A, int* perm_c, int* perm_r,
              SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t *, int* info);

  void zCreate_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               dbl_complex* nzval, int* rowind, int* colptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void zCreate_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                               dbl_complex* nzval, int* colind, int* rowptr,
                               Stype_t stype, Dtype_t dtype, Mtype_t mtype);
  void zCreate_Dense_Matrix (SuperMatrix* A, int m, int n,
                             dbl_complex* vals, int ld_vals,
                             Stype_t stype, Dtype_t dtype, Mtype_t mtype);

  void zgstrf (superlu_options_t*, SuperMatrix*, double, int, int, int*, void*, int, int*, int*,
               SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int*);
  void zgstrs (trans_t, SuperMatrix*, SuperMatrix*, int*, int*,
               SuperMatrix*, SuperLUStat_t*, int*);

  void zPrint_CompCol_Matrix (char* what, SuperMatrix* A);
  void zPrint_CompRow_Matrix (char* what, SuperMatrix* A);
  void zPrint_SuperNode_Matrix (char* what, SuperMatrix* A);
  void zPrint_Dense_Matrix (char* what, SuperMatrix* A);


  // all types

  void get_perm_c (int permc_spec, SuperMatrix* A, int* perm_c);

  void Destroy_SuperMatrix_Store (SuperMatrix* A);
  void Destroy_CompCol_Matrix (SuperMatrix* A);
  void Destroy_CompRow_Matrix (SuperMatrix* A);
  void Destroy_SuperNode_Matrix (SuperMatrix* A);
  void Destroy_CompCol_Permuted (SuperMatrix* A);
  void Destroy_Dense_Matrix (SuperMatrix* A);

}

#endif
