#ifndef VIENNACL_LINALG_OPENCL_VANDERMONDE_MATRIX_OPERATIONS_HPP_
#define VIENNACL_LINALG_OPENCL_VANDERMONDE_MATRIX_OPERATIONS_HPP_

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

/** @file viennacl/linalg/opencl/vandermonde_matrix_operations.hpp
    @brief Implementations of operations using vandermonde_matrix
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
    namespace opencl
    {

      /** @brief Carries out matrix-vector multiplication with a vandermonde_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
        template<class SCALARTYPE, unsigned int ALIGNMENT>
        void prod_impl(const viennacl::vandermonde_matrix<SCALARTYPE, ALIGNMENT> & mat,
                      const viennacl::vector_base<SCALARTYPE> & vec,
                            viennacl::vector_base<SCALARTYPE> & result)
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(mat).context());
          viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::init(ctx);

          viennacl::ocl::kernel & kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::fft<SCALARTYPE>::program_name(), "vandermonde_prod");
          viennacl::ocl::enqueue(kernel(viennacl::traits::opencl_handle(mat),
                                        viennacl::traits::opencl_handle(vec),
                                        viennacl::traits::opencl_handle(result),
                                        static_cast<cl_uint>(mat.size1())));
        }

    } //namespace opencl
  } //namespace linalg
} //namespace viennacl


#endif
