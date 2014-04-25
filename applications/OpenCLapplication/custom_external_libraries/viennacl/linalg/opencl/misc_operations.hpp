#ifndef VIENNACL_LINALG_OPENCL_MISC_OPERATIONS_HPP_
#define VIENNACL_LINALG_OPENCL_MISC_OPERATIONS_HPP_

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

/** @file viennacl/linalg/opencl/misc_operations.hpp
    @brief Implementations of operations using compressed_matrix and OpenCL
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/opencl/kernels/ilu.hpp"


namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {

      namespace detail
      {

        template <typename ScalarType>
        void level_scheduling_substitute(vector<ScalarType> & vec,
                                     viennacl::backend::mem_handle const & row_index_array,
                                     viennacl::backend::mem_handle const & row_buffer,
                                     viennacl::backend::mem_handle const & col_buffer,
                                     viennacl::backend::mem_handle const & element_buffer,
                                     vcl_size_t num_rows
                                    )
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(vec).context());

          viennacl::linalg::opencl::kernels::ilu<ScalarType>::init(ctx);
          viennacl::ocl::kernel & k = ctx.get_kernel(viennacl::linalg::opencl::kernels::ilu<ScalarType>::program_name(), "level_scheduling_substitute");

          viennacl::ocl::enqueue(k(row_index_array.opencl_handle(), row_buffer.opencl_handle(), col_buffer.opencl_handle(), element_buffer.opencl_handle(),
                                   vec,
                                   static_cast<cl_uint>(num_rows)));
        }

      } //namespace detail


    } // namespace opencl
  } //namespace linalg
} //namespace viennacl


#endif
