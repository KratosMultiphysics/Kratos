#ifndef VIENNACL_TOOLS_MATRIX_PROD_KERNEL_CLASS_DEDUCER_HPP_
#define VIENNACL_TOOLS_MATRIX_PROD_KERNEL_CLASS_DEDUCER_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

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
    namespace detail
    {
      template <typename MatrixType>
      struct extract_matrix
      {
        typedef typename MatrixType::ERROR_UNKNOWN_MATRIX_TYPE_PROVIDED   error_type;
      };
      
      template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
      struct extract_matrix < viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
      {
        typedef viennacl::matrix<SCALARTYPE, F, ALIGNMENT>   type;
      };

      template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
      struct extract_matrix < const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
      {
        typedef viennacl::matrix<SCALARTYPE, F, ALIGNMENT>   type;
      };

      
      template <typename MatrixType>
      struct extract_matrix < viennacl::matrix_range<MatrixType> >
      {
        typedef typename extract_matrix<MatrixType>::type   type;
      };

      template <typename MatrixType>
      struct extract_matrix < const viennacl::matrix_range<MatrixType> >
      {
        typedef typename extract_matrix<MatrixType>::type   type;
      };
      
      
    }
    
    
    
    /** @brief deduces kernel type for C=A*B, where A, B, C are MatrixType1, MatrixType2 and MatrixType3 respectively */
    template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
    struct MATRIX_PROD_KERNEL_CLASS_DEDUCER
    {
      typedef typename MATRIX_PROD_KERNEL_CLASS_DEDUCER< typename detail::extract_matrix<MatrixType1>::type,
                                                         typename detail::extract_matrix<MatrixType2>::type,
                                                         typename detail::extract_matrix<MatrixType3>::type>::ResultType   ResultType;
    };
    
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
