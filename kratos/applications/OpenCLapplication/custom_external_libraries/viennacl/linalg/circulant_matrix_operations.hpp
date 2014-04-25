#ifndef VIENNACL_LINALG_CIRCULANT_MATRIX_OPERATIONS_HPP_
#define VIENNACL_LINALG_CIRCULANT_MATRIX_OPERATIONS_HPP_

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

/** @file viennacl/linalg/circulant_matrix_operations.hpp
    @brief Implementations of operations using circulant_matrix. Experimental.
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/fft.hpp"
//#include "viennacl/linalg/kernels/coordinate_matrix_kernels.h"

namespace viennacl
{
  namespace linalg
  {

    // A * x

    /** @brief Carries out matrix-vector multiplication with a circulant_matrix
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    */
      template<class SCALARTYPE, unsigned int ALIGNMENT>
      void prod_impl(const viennacl::circulant_matrix<SCALARTYPE, ALIGNMENT> & mat,
                     const viennacl::vector_base<SCALARTYPE> & vec,
                           viennacl::vector_base<SCALARTYPE> & result)
      {
        assert(mat.size1() == result.size());
        assert(mat.size2() == vec.size());
        //result.clear();

        //std::cout << "prod(circulant_matrix" << ALIGNMENT << ", vector) called with internal_nnz=" << mat.internal_nnz() << std::endl;

        viennacl::vector<SCALARTYPE> circ(mat.elements().size() * 2);
        viennacl::detail::fft::real_to_complex(mat.elements(), circ, mat.elements().size());

        viennacl::vector<SCALARTYPE> tmp(vec.size() * 2);
        viennacl::vector<SCALARTYPE> tmp2(vec.size() * 2);

        viennacl::detail::fft::real_to_complex(vec, tmp, vec.size());
        viennacl::linalg::convolve(circ, tmp, tmp2);
        viennacl::detail::fft::complex_to_real(tmp2, result, vec.size());

      }

  } //namespace linalg


} //namespace viennacl


#endif
