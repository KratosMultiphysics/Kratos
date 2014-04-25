#ifndef VIENNACL_LINALG_OPENCL_SCALAR_OPERATIONS_HPP_
#define VIENNACL_LINALG_OPENCL_SCALAR_OPERATIONS_HPP_

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

/** @file viennacl/linalg/opencl/scalar_operations.hpp
    @brief Implementations of scalar operations using OpenCL
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/opencl/kernels/scalar.hpp"
#include "viennacl/linalg/opencl/common.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/result_of.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/handle.hpp"
#include "viennacl/traits/stride.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      template <typename S1,
                typename S2, typename ScalarType1>
      typename viennacl::enable_if< viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                    && viennacl::is_any_scalar<ScalarType1>::value
                                  >::type
      as(S1 & s1,
         S2 const & s2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha)
      {
        assert( &viennacl::traits::opencl_handle(s1).context() == &viennacl::traits::opencl_handle(s2).context() && bool("Operands not in the same OpenCL context!"));

        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(s1).context());
        viennacl::linalg::opencl::kernels::scalar<value_type>::init(ctx);

        cl_uint options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);

        viennacl::ocl::kernel & k = ctx.get_kernel(viennacl::linalg::opencl::kernels::scalar<value_type>::program_name(),
                                                   (viennacl::is_cpu_scalar<ScalarType1>::value ? "as_cpu" : "as_gpu"));
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(s1),
                                 viennacl::traits::opencl_handle(viennacl::tools::promote_if_host_scalar<value_type>(alpha)),
                                 options_alpha,
                                 viennacl::traits::opencl_handle(s2) )
                              );
      }


      template <typename S1,
                typename S2, typename ScalarType1,
                typename S3, typename ScalarType2>
      typename viennacl::enable_if< viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                    && viennacl::is_scalar<S3>::value
                                    && viennacl::is_any_scalar<ScalarType1>::value
                                    && viennacl::is_any_scalar<ScalarType2>::value
                                  >::type
      asbs(S1 & s1,
           S2 const & s2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha,
           S3 const & s3, ScalarType2 const & beta,  vcl_size_t len_beta,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        assert( &viennacl::traits::opencl_handle(s1).context() == &viennacl::traits::opencl_handle(s2).context() && bool("Operands not in the same OpenCL context!"));
        assert( &viennacl::traits::opencl_handle(s2).context() == &viennacl::traits::opencl_handle(s3).context() && bool("Operands not in the same OpenCL context!"));

        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(s1).context());
        viennacl::linalg::opencl::kernels::scalar<value_type>::init(ctx);

        std::string kernel_name;
        if (viennacl::is_cpu_scalar<ScalarType1>::value && viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_cpu_cpu";
        else if (viennacl::is_cpu_scalar<ScalarType1>::value && !viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_cpu_gpu";
        else if (!viennacl::is_cpu_scalar<ScalarType1>::value && viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_gpu_cpu";
        else
          kernel_name = "asbs_gpu_gpu";

        cl_uint options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);
        cl_uint options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        viennacl::ocl::kernel & k = ctx.get_kernel(viennacl::linalg::opencl::kernels::scalar<value_type>::program_name(), kernel_name);
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(s1),
                                 viennacl::traits::opencl_handle(viennacl::tools::promote_if_host_scalar<value_type>(alpha)),
                                 options_alpha,
                                 viennacl::traits::opencl_handle(s2),
                                 viennacl::traits::opencl_handle(viennacl::tools::promote_if_host_scalar<value_type>(beta)),
                                 options_beta,
                                 viennacl::traits::opencl_handle(s3) )
                              );
      }


      template <typename S1,
                typename S2, typename ScalarType1,
                typename S3, typename ScalarType2>
      typename viennacl::enable_if< viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                    && viennacl::is_scalar<S3>::value
                                    && viennacl::is_any_scalar<ScalarType1>::value
                                    && viennacl::is_any_scalar<ScalarType2>::value
                                  >::type
      asbs_s(S1 & s1,
             S2 const & s2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha,
             S3 const & s3, ScalarType2 const & beta,  vcl_size_t len_beta,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        assert( &viennacl::traits::opencl_handle(s1).context() == &viennacl::traits::opencl_handle(s2).context() && bool("Operands not in the same OpenCL context!"));
        assert( &viennacl::traits::opencl_handle(s2).context() == &viennacl::traits::opencl_handle(s3).context() && bool("Operands not in the same OpenCL context!"));

        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(s1).context());
        viennacl::linalg::opencl::kernels::scalar<value_type>::init(ctx);

        std::string kernel_name;
        if (viennacl::is_cpu_scalar<ScalarType1>::value && viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_s_cpu_cpu";
        else if (viennacl::is_cpu_scalar<ScalarType1>::value && !viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_s_cpu_gpu";
        else if (!viennacl::is_cpu_scalar<ScalarType1>::value && viennacl::is_cpu_scalar<ScalarType2>::value)
          kernel_name = "asbs_s_gpu_cpu";
        else
          kernel_name = "asbs_s_gpu_gpu";

        cl_uint options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);
        cl_uint options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        viennacl::ocl::kernel & k = ctx.get_kernel(viennacl::linalg::opencl::kernels::scalar<value_type>::program_name(), kernel_name);
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(s1),
                                 viennacl::traits::opencl_handle(viennacl::tools::promote_if_host_scalar<value_type>(alpha)),
                                 options_alpha,
                                 viennacl::traits::opencl_handle(s2),
                                 viennacl::traits::opencl_handle(viennacl::tools::promote_if_host_scalar<value_type>(beta)),
                                 options_beta,
                                 viennacl::traits::opencl_handle(s3) )
                              );
      }


      /** @brief Swaps the contents of two scalars, data is copied
      *
      * @param s1   The first scalar
      * @param s2   The second scalar
      */
      template <typename S1, typename S2>
      typename viennacl::enable_if<    viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                  >::type
      swap(S1 & s1, S2 & s2)
      {
        assert( &viennacl::traits::opencl_handle(s1).context() == &viennacl::traits::opencl_handle(s2).context() && bool("Operands not in the same OpenCL context!"));

        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(s1).context());
        viennacl::linalg::opencl::kernels::scalar<value_type>::init(ctx);

        viennacl::ocl::kernel & k = ctx.get_kernel(viennacl::linalg::opencl::kernels::scalar<value_type>::program_name(), "swap");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(s1),
                                 viennacl::traits::opencl_handle(s2))
                              );
      }



    } //namespace opencl
  } //namespace linalg
} //namespace viennacl


#endif
