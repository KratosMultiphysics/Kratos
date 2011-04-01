/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory
======================================================================= */

#ifndef _VIENNACL_TOOLS_MATRIX_PROD_KERNEL_CLASS_DEDUCER_HPP_
#define _VIENNACL_TOOLS_MATRIX_PROD_KERNEL_CLASS_DEDUCER_HPP_

/** @file matrix_prod_kernel_class_deducer.hpp
    @brief Implementation of a helper meta class for deducing the correct kernels for matrix-matrix products
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/linalg/kernels/matrix_prod_col_col_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_col_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_row_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_row_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_col_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_col_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_row_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_row_row_kernels.h"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {
    /** @brief deduces kernel type for C=A*B, where A, B, C are MatrixType1, MatrixType2 and MatrixType3 respectively */
    template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER
    {};
    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_row_row_row<SCALARTYPE, ALIGNMENT>     ResultType;
    };

    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_row_row_col<SCALARTYPE, ALIGNMENT>     ResultType;
    };
    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_row_col_row<SCALARTYPE, ALIGNMENT>     ResultType;
    };

    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_row_col_col<SCALARTYPE, ALIGNMENT>     ResultType;
    };

    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_col_row_row<SCALARTYPE, ALIGNMENT>     ResultType;
    };

    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_col_row_col<SCALARTYPE, ALIGNMENT>     ResultType;
    };
    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_col_col_row<SCALARTYPE, ALIGNMENT>     ResultType;
    };

    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT>,
                                             viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_prod_col_col_col<SCALARTYPE, ALIGNMENT>     ResultType;
    };
    
  }

}

#endif
