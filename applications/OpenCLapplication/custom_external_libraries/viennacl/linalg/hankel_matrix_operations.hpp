#ifndef VIENNACL_LINALG_HANKEL_MATRIX_OPERATIONS_HPP_
#define VIENNACL_LINALG_HANKEL_MATRIX_OPERATIONS_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/hankel_matrix_operations.hpp
    @brief Implementations of operations using hankel_matrix. Experimental.
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/fft.hpp"
#include "viennacl/linalg/toeplitz_matrix_operations.hpp"

namespace viennacl
{
  namespace linalg
  {

    // A * x

    /** @brief Carries out matrix-vector multiplication with a hankel_matrix
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void prod_impl(const viennacl::hankel_matrix<SCALARTYPE, ALIGNMENT> & mat,
                   const viennacl::vector_base<SCALARTYPE> & vec,
                         viennacl::vector_base<SCALARTYPE> & result)
    {
      assert(mat.size1() == result.size());
      assert(mat.size2() == vec.size());

      prod_impl(mat.elements(), vec, result);
      viennacl::detail::fft::reverse(result);
    }

  } //namespace linalg


} //namespace viennacl


#endif
