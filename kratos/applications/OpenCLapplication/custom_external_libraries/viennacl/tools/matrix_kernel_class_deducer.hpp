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

#ifndef _VIENNACL_TOOLS_MATRIX_KERNEL_CLASS_DEDUCER_HPP_
#define _VIENNACL_TOOLS_MATRIX_KERNEL_CLASS_DEDUCER_HPP_

/** @file matrix_kernel_class_deducer.hpp
    @brief Implementation of a helper meta class for deducing the correct kernels for the supplied matrix
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/linalg/kernels/matrix_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_row_kernels.h"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {
    /**     @brief Implementation of a helper meta class for deducing the correct kernels for the supplied matrix */
    template <typename MatrixType1>
    struct MATRIX_KERNEL_CLASS_DEDUCER
    {};
    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::row_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_row<SCALARTYPE, ALIGNMENT>     ResultType;
    };
    
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    struct MATRIX_KERNEL_CLASS_DEDUCER< viennacl::matrix<SCALARTYPE, viennacl::column_major, ALIGNMENT> >
    {
      typedef viennacl::linalg::kernels::matrix_col<SCALARTYPE, ALIGNMENT>     ResultType;
    };

  }

}

#endif
