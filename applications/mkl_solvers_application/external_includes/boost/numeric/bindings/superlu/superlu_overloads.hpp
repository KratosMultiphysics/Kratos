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

#ifndef BOOST_NUMERIC_BINDINGS_SUPERLU_OVERLOADS_H
#define BOOST_NUMERIC_BINDINGS_SUPERLU_OVERLOADS_H

#include <cassert>
#include <complex>
#include <boost/numeric/bindings/superlu/superlu_decl.hpp>

namespace boost { namespace numeric { namespace bindings {

  namespace superlu { namespace detail {

    void gssv (float, superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r,
               SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      sgssv (options, A, perm_c, perm_r, L, U, B, stat, info);
    }

    void gssv (double, superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r,
               SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      dgssv (options, A, perm_c, perm_r, L, U, B, stat, info);
    }

    void gssv (std::complex<float>, superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r,
               SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      cgssv (options, A, perm_c, perm_r, L, U, B, stat, info);
    }

    void gssv (std::complex<double>, superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r,
               SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      zgssv (options, A, perm_c, perm_r, L, U, B, stat, info);
    }


    ////////////////////////////////

    inline
    void Create_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                float* nzval, int* rowind, int* colptr,
                                Stype_t stype = SLU_NC,
                                Dtype_t dtype = SLU_S,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NC);
      assert (dtype == SLU_S);
      assert (mtype == SLU_GE);
      sCreate_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr,
                              SLU_NC, SLU_S, SLU_GE);
    }

    inline
    void Create_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                float* nzval, int* colind, int* rowptr,
                                Stype_t stype = SLU_NR,
                                Dtype_t dtype = SLU_S,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NR);
      assert (dtype == SLU_S);
      assert (mtype == SLU_GE);
      sCreate_CompRow_Matrix (A, m, n, nnz, nzval, colind, rowptr,
                              SLU_NR, SLU_S, SLU_GE);
    }

    inline
    void Create_Dense_Matrix (SuperMatrix* A, int m, int n,
                              float* vals, int ld_vals,
                              Stype_t stype = SLU_DN,
                              Dtype_t dtype = SLU_S,
                              Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_DN);
      assert (dtype == SLU_S);
      assert (mtype == SLU_GE);
      sCreate_Dense_Matrix (A, m, n, vals, ld_vals, SLU_DN, SLU_S, SLU_GE);
    }


    inline
    void Create_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                double* nzval, int* rowind, int* colptr,
                                Stype_t stype = SLU_NC,
                                Dtype_t dtype = SLU_D,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NC);
      assert (dtype == SLU_D);
      assert (mtype == SLU_GE);
      dCreate_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr,
                              SLU_NC, SLU_D, SLU_GE);
    }

    inline
    void Create_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                double* nzval, int* colind, int* rowptr,
                                Stype_t stype = SLU_NR,
                                Dtype_t dtype = SLU_D,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NR);
      assert (dtype == SLU_D);
      assert (mtype == SLU_GE);
      dCreate_CompRow_Matrix (A, m, n, nnz, nzval, colind, rowptr,
                              SLU_NR, SLU_D, SLU_GE);
    }

    inline
    void Create_Dense_Matrix (SuperMatrix* A, int m, int n,
                              double* vals, int ld_vals,
                              Stype_t stype = SLU_DN,
                              Dtype_t dtype = SLU_D,
                              Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_DN);
      assert (dtype == SLU_D);
      assert (mtype == SLU_GE);
      dCreate_Dense_Matrix (A, m, n, vals, ld_vals, SLU_DN, SLU_D, SLU_GE);
    }


    inline
    void Create_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                std::complex<float>* nzval,
                                int* rowind, int* colptr,
                                Stype_t stype = SLU_NC,
                                Dtype_t dtype = SLU_C,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NC);
      assert (dtype == SLU_C);
      assert (mtype == SLU_GE);
      cCreate_CompCol_Matrix (A, m, n, nnz,
                              reinterpret_cast<flt_complex*> (nzval),
                              rowind, colptr,
                              SLU_NC, SLU_C, SLU_GE);
    }

    inline
    void Create_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                std::complex<float>* nzval,
                                int* colind, int* rowptr,
                                Stype_t stype = SLU_NR,
                                Dtype_t dtype = SLU_C,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NR);
      assert (dtype == SLU_C);
      assert (mtype == SLU_GE);
      cCreate_CompRow_Matrix (A, m, n, nnz,
                              reinterpret_cast<flt_complex*> (nzval),
                              colind, rowptr,
                              SLU_NR, SLU_C, SLU_GE);
    }

    inline
    void Create_Dense_Matrix (SuperMatrix* A, int m, int n,
                              std::complex<float>* vals, int ld_vals,
                              Stype_t stype = SLU_DN,
                              Dtype_t dtype = SLU_C,
                              Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_DN);
      assert (dtype == SLU_C);
      assert (mtype == SLU_GE);
      cCreate_Dense_Matrix (A, m, n,
                            reinterpret_cast<flt_complex*> (vals),
                            ld_vals, SLU_DN, SLU_C, SLU_GE);
    }


    inline
    void Create_CompCol_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                std::complex<double>* nzval,
                                int* rowind, int* colptr,
                                Stype_t stype = SLU_NC,
                                Dtype_t dtype = SLU_Z,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NC);
      assert (dtype == SLU_Z);
      assert (mtype == SLU_GE);
      zCreate_CompCol_Matrix (A, m, n, nnz,
                              reinterpret_cast<dbl_complex*> (nzval),
                              rowind, colptr,
                              SLU_NC, SLU_Z, SLU_GE);
    }

    inline
    void Create_CompRow_Matrix (SuperMatrix* A, int m, int n, int nnz,
                                std::complex<double>* nzval,
                                int* colind, int* rowptr,
                                Stype_t stype = SLU_NR,
                                Dtype_t dtype = SLU_Z,
                                Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_NR);
      assert (dtype == SLU_Z);
      assert (mtype == SLU_GE);
      zCreate_CompRow_Matrix (A, m, n, nnz,
                              reinterpret_cast<dbl_complex*> (nzval),
                              colind, rowptr,
                              SLU_NR, SLU_Z, SLU_GE);
    }

    inline
    void Create_Dense_Matrix (SuperMatrix* A, int m, int n,
                              std::complex<double>* vals, int ld_vals,
                              Stype_t stype = SLU_DN,
                              Dtype_t dtype = SLU_Z,
                              Mtype_t mtype = SLU_GE)
    {
      assert (stype == SLU_DN);
      assert (dtype == SLU_Z);
      assert (mtype == SLU_GE);
      zCreate_Dense_Matrix (A, m, n,
                            reinterpret_cast<dbl_complex*> (vals),
                            ld_vals, SLU_DN, SLU_Z, SLU_GE);
    }

    template<typename T> inline
    void Create_Comp_Matrix(traits::column_major_t, SuperMatrix* A, int m, int n,
                            int nnz, T* nzval, int* colptr, int* rowind)
    {
      Create_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr);
    }

    template<typename T> inline
    void Create_Comp_Matrix(traits::row_major_t, SuperMatrix* A, int m, int n,
                            int nnz, T* nzval, int* rowptr, int* colind)
    {
      Create_CompRow_Matrix (A, m, n, nnz, nzval, colind, rowptr);
    }

    inline
    void gstrf (float, superlu_options_t* options, SuperMatrix* A, float drop_tol,
                int relax, int panel_size, int* etree, void* work, int lwork,
                int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U,
                SuperLUStat_t *stat, int *info)
    {
      sgstrf(options, A, drop_tol, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    inline
    void gstrf (double, superlu_options_t* options, SuperMatrix* A, double drop_tol,
                int relax, int panel_size, int* etree, void* work, int lwork,
                int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U,
                SuperLUStat_t *stat, int *info)
    {
      dgstrf(options, A, drop_tol, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    inline
    void gstrf (std::complex<float>, superlu_options_t* options, SuperMatrix* A, float drop_tol,
                int relax, int panel_size, int* etree, void* work, int lwork,
                int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U,
                SuperLUStat_t *stat, int *info)
    {
      cgstrf(options, A, drop_tol, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, stat, info);
    }

    inline
    void gstrf (std::complex<double>, superlu_options_t* options, SuperMatrix* A, double drop_tol,
                int relax, int panel_size, int* etree, void* work, int lwork,
                int* perm_c, int* perm_r, SuperMatrix* L, SuperMatrix* U,
                SuperLUStat_t *stat, int *info)
    {
      zgstrf(options, A, drop_tol, relax, panel_size, etree, work, lwork, perm_c, perm_r, L, U, stat, info);
    }


    inline
    void gstrs (float, trans_t trans, SuperMatrix* L, SuperMatrix* U, int* perm_c, int* perm_r,
                SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      sgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
    }

    inline
    void gstrs (double, trans_t trans, SuperMatrix* L, SuperMatrix* U, int* perm_c, int* perm_r,
                SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      dgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
    }

    inline
    void gstrs (std::complex<float>, trans_t trans, SuperMatrix* L, SuperMatrix* U, int* perm_c, int* perm_r,
                SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      cgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
    }

    inline
    void gstrs (std::complex<double>, trans_t trans, SuperMatrix* L, SuperMatrix* U, int* perm_c, int* perm_r,
                SuperMatrix* B, SuperLUStat_t* stat, int* info)
    {
      zgstrs(trans, L, U, perm_c, perm_r, B, stat, info);
    }

    ////////////////////////////////

    void Print_CompCol_Matrix (float, char* what, SuperMatrix* A) {
      sPrint_CompCol_Matrix (what, A);
    }
    void Print_Dense_Matrix (float, char* what, SuperMatrix* A) {
      sPrint_Dense_Matrix (what, A);
    }
    void Print_SuperNode_Matrix (float, char* what, SuperMatrix* A) {
      sPrint_SuperNode_Matrix (what, A);
    }

    void Print_CompCol_Matrix (double, char* what, SuperMatrix* A) {
      dPrint_CompCol_Matrix (what, A);
    }
    void Print_Dense_Matrix (double, char* what, SuperMatrix* A) {
      dPrint_Dense_Matrix (what, A);
    }
    void Print_SuperNode_Matrix (double, char* what, SuperMatrix* A) {
      dPrint_SuperNode_Matrix (what, A);
    }

    void Print_CompCol_Matrix (std::complex<float>,
                               char* what, SuperMatrix* A)
    {
      cPrint_CompCol_Matrix (what, A);
    }
    void Print_Dense_Matrix (std::complex<float>,
                             char* what, SuperMatrix* A)
    {
      cPrint_Dense_Matrix (what, A);
    }
    void Print_SuperNode_Matrix (std::complex<float>,
                                 char* what, SuperMatrix* A)
    {
      cPrint_SuperNode_Matrix (what, A);
    }

    void Print_CompCol_Matrix (std::complex<double>,
                               char* what, SuperMatrix* A)
    {
      zPrint_CompCol_Matrix (what, A);
    }
    void Print_Dense_Matrix (std::complex<double>,
                             char* what, SuperMatrix* A)
    {
      zPrint_Dense_Matrix (what, A);
    }
    void Print_SuperNode_Matrix (std::complex<double>,
                                 char* what, SuperMatrix* A)
    {
      zPrint_SuperNode_Matrix (what, A);
    }
    template<typename T> inline
    void Print_Comp_Matrix (traits::column_major_t, T t,
                            char* what, SuperMatrix* A)
    {
      Print_CompCol_Matrix (t, what, A);
    }
    template<typename T>
    void Print_CompRow_Matrix(T, char *what, SuperMatrix *A)
    {
      NRformat     *Astore;
      register int i,n;
      T            *dp;

      printf("\nCompRow matrix %s:\n", what);
      printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
      n = A->nrow;
      Astore = (NRformat *) A->Store;
      dp = (T *) Astore->nzval;
      printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
      printf("nzval: ");
      for (i = 0; i < Astore->rowptr[n]; ++i) printf("%f  ", dp[i]);
      printf("\ncolind: ");
      for (i = 0; i < Astore->rowptr[n]; ++i) printf("%d  ", Astore->colind[i]);
      printf("\nrowptr: ");
      for (i = 0; i <= n; ++i) printf("%d  ", Astore->rowptr[i]);
      printf("\n");
      fflush(stdout);
    }
    template<typename T>
    void Print_CompRow_Matrix (std::complex<T>, char* what, SuperMatrix* A)
    {
      NRformat     *Astore;
      register int i,n;
      T            *dp;

      printf("\nCompRow matrix %s:\n", what);
      printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
      n = A->nrow;
      Astore = (NRformat *) A->Store;
      dp = (T *) Astore->nzval;
      printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
      printf("nzval: ");
      for (i = 0; i < 2*Astore->rowptr[n]; ++i) printf("%f  ", dp[i]);
      printf("\ncolind: ");
      for (i = 0; i < Astore->rowptr[n]; ++i) printf("%d  ", Astore->colind[i]);
      printf("\nrowptr: ");
      for (i = 0; i <= n; ++i) printf("%d  ", Astore->rowptr[i]);
      printf("\n");
      fflush(stdout);
    }
    template<typename T> inline
    void Print_Comp_Matrix (traits::row_major_t, T t,
                            char* what, SuperMatrix* A)
    {
      Print_CompRow_Matrix (t, what, A);
    }
  }}

}}}


#endif
