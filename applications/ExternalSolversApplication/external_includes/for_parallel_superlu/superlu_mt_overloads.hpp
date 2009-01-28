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

#ifndef BOOST_NUMERIC_BINDINGS_SUPERLU_MT_OVERLOADS_H
#define BOOST_NUMERIC_BINDINGS_SUPERLU_MT_OVERLOADS_H

#include <cassert>
#include "external_includes/superlu/pdsp_defs.h"

namespace boost { namespace numeric { namespace bindings {

  namespace superlu { namespace detail {
      
      void pgssv( int nproc, SuperMatrix* A, int* perm_c, int* perm_r,
                  SuperMatrix* L, SuperMatrix* U, SuperMatrix* B, int* info )
      {
            std::cout << "calling pdgssv..." << std::endl;
          pdgssv(nproc, A, perm_c, perm_r, L, U, B, info);
          std::cout << "done." << std::endl;
      }

    ////////////////////////////////

    
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


    template<typename T> inline
    void Create_Comp_Matrix(traits::column_major_t, SuperMatrix* A, int m, int n,
                            int nnz, T* nzval, int* colptr, int* rowind)
    {
      Create_CompCol_Matrix (A, m, n, nnz, nzval, rowind, colptr);
    }
    
    template<typename T> inline
    void Destroy_Comp_Matrix(traits::column_major_t, SuperMatrix* A)
    {
        Destroy_CompCol_Matrix(A);
    }

  }}

}}}


#endif
