#ifndef VIENNACL_LINALG_CUDA_SCALAR_OPERATIONS_HPP_
#define VIENNACL_LINALG_CUDA_SCALAR_OPERATIONS_HPP_

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

/** @file viennacl/linalg/cuda/scalar_operations.hpp
    @brief Implementations of scalar operations using CUDA
*/

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/stride.hpp"
#include "viennacl/linalg/cuda/common.hpp"

// includes CUDA
#include <cuda_runtime.h>


namespace viennacl
{
  namespace linalg
  {
    namespace cuda
    {

      namespace detail
      {

      }

      /////////////////// as /////////////////////////////

      template <typename T>
      __global__ void as_kernel(T * s1, const T * fac2, unsigned int options2, const T * s2)
      {
          T alpha = *fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          *s1 = *s2 * alpha;
      }

      template <typename T>
      __global__ void as_kernel(T * s1, T fac2, unsigned int options2, const T * s2)
      {
          T alpha = fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          *s1 = *s2 * alpha;
      }

      template <typename S1,
                typename S2, typename ScalarType1>
      typename viennacl::enable_if< viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                    && viennacl::is_any_scalar<ScalarType1>::value
                                  >::type
      as(S1 & s1,
         S2 const & s2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha)
      {
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        as_kernel<<<1, 1>>>(detail::cuda_arg<value_type>(s1),
                            detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                            options_alpha,
                            detail::cuda_arg<value_type>(s2));
        VIENNACL_CUDA_LAST_ERROR_CHECK("as_kernel");
      }

      //////////////////// asbs ////////////////////////////

      // alpha and beta on GPU
      template <typename T>
      __global__ void asbs_kernel(T * s1,
                                  const T * fac2, unsigned int options2, const T * s2,
                                  const T * fac3, unsigned int options3, const T * s3)
      {
          T alpha = *fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = *fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 = *s2 * alpha + *s3 * beta;
      }

      // alpha on CPU, beta on GPU
      template <typename T>
      __global__ void asbs_kernel(T * s1,
                                  T fac2, unsigned int options2, const T * s2,
                                  const T * fac3, unsigned int options3, const T * s3)
      {
          T alpha = fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = *fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 = *s2 * alpha + *s3 * beta;
      }

      // alpha on GPU, beta on CPU
      template <typename T>
      __global__ void asbs_kernel(T * s1,
                                  const T * fac2, unsigned int options2, const T * s2,
                                  T fac3, unsigned int options3, const T * s3)
      {
          T alpha = *fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 = *s2 * alpha + *s3 * beta;
      }

      // alpha and beta on CPU
      template <typename T>
      __global__ void asbs_kernel(T * s1,
                                  T fac2, unsigned int options2, const T * s2,
                                  T fac3, unsigned int options3, const T * s3)
      {
          T alpha = fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 = *s2 * alpha + *s3 * beta;
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
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);
        unsigned int options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        value_type temporary_beta = 0;
        if (viennacl::is_cpu_scalar<ScalarType2>::value)
          temporary_beta = beta;

        asbs_kernel<<<1, 1>>>(detail::cuda_arg<value_type>(s1),
                              detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                              options_alpha,
                              detail::cuda_arg<value_type>(s2),
                              detail::cuda_arg<value_type>(detail::arg_reference(beta, temporary_beta)),
                              options_beta,
                              detail::cuda_arg<value_type>(s3) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("asbs_kernel");
      }

      //////////////////// asbs_s ////////////////////

      // alpha and beta on GPU
      template <typename T>
      __global__ void asbs_s_kernel(T * s1,
                                    const T * fac2, unsigned int options2, const T * s2,
                                    const T * fac3, unsigned int options3, const T * s3)
      {
          T alpha = *fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = *fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 += *s2 * alpha + *s3 * beta;
      }

      // alpha on CPU, beta on GPU
      template <typename T>
      __global__ void asbs_s_kernel(T * s1,
                                    T fac2, unsigned int options2, const T * s2,
                                    const T * fac3, unsigned int options3, const T * s3)
      {
          T alpha = fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = *fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 += *s2 * alpha + *s3 * beta;
      }

      // alpha on GPU, beta on CPU
      template <typename T>
      __global__ void asbs_s_kernel(T * s1,
                                    const T * fac2, unsigned int options2, const T * s2,
                                    T fac3, unsigned int options3, const T * s3)
      {
          T alpha = *fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 += *s2 * alpha + *s3 * beta;
      }

      // alpha and beta on CPU
      template <typename T>
      __global__ void asbs_s_kernel(T * s1,
                                    T fac2, unsigned int options2, const T * s2,
                                    T fac3, unsigned int options3, const T * s3)
      {
          T alpha = fac2;
          if (options2 & (1 << 0))
            alpha = -alpha;
          if (options2 & (1 << 1))
            alpha = ((T)(1)) / alpha;

          T beta = fac3;
          if (options3 & (1 << 0))
            beta = -beta;
          if (options3 & (1 << 1))
            beta = ((T)(1)) / beta;

          *s1 += *s2 * alpha + *s3 * beta;
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
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);
        unsigned int options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        value_type temporary_beta = 0;
        if (viennacl::is_cpu_scalar<ScalarType2>::value)
          temporary_beta = beta;

        std::cout << "Launching asbs_s_kernel..." << std::endl;
        asbs_s_kernel<<<1, 1>>>(detail::cuda_arg<value_type>(s1),
                                detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                                options_alpha,
                                detail::cuda_arg<value_type>(s2),
                                detail::cuda_arg<value_type>(detail::arg_reference(beta, temporary_beta)),
                                options_beta,
                                detail::cuda_arg<value_type>(s3) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("asbs_s_kernel");
      }

      ///////////////// swap //////////////////

      template <typename T>
      __global__ void scalar_swap_kernel(T * s1, T * s2)
      {
        T tmp = *s2;
        *s2 = *s1;
        *s1 = tmp;
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
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        scalar_swap_kernel<<<1, 1>>>(detail::cuda_arg<value_type>(s1),detail::cuda_arg<value_type>(s2));
      }



    } //namespace single_threaded
  } //namespace linalg
} //namespace viennacl


#endif
