#ifndef VIENNACL_LINALG_CUDA_VECTOR_OPERATIONS_HPP_
#define VIENNACL_LINALG_CUDA_VECTOR_OPERATIONS_HPP_

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

/** @file viennacl/linalg/cuda/vector_operations.hpp
    @brief Implementations of vector operations using a plain single-threaded execution on CPU
*/

#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/scalar.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/stride.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace cuda
    {

      //
      // Introductory note: By convention, all dimensions are already checked in the dispatcher frontend. No need to double-check again in here!
      //


      //////////////////////// av /////////////////////////////

      // gpu scalar
      template <typename T>
      __global__ void av_kernel(T * vec1,
                                unsigned int start1,
                                unsigned int inc1,
                                unsigned int size1,

                                const T * fac2,
                                unsigned int options2,
                                const T * vec2,
                                unsigned int start2,
                                unsigned int inc2)
      {
        T alpha = *fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        if (options2 & (1 << 1))
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
            vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha;
        }
        else
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
            vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha;
        }
      }

      // cpu scalar
      template <typename T>
      __global__ void av_kernel(T * vec1,
                                unsigned int start1,
                                unsigned int inc1,
                                unsigned int size1,

                                T fac2,
                                unsigned int options2,
                                const T * vec2,
                                unsigned int start2,
                                unsigned int inc2)
      {
        T alpha = fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        if (options2 & (1 << 1))
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
            vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha;
        }
        else
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
            vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha;
        }
      }



      template <typename T, typename ScalarType1>
      void av(vector_base<T> & vec1,
              vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha)
      {
        typedef T        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        av_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                                options_alpha,
                                detail::cuda_arg<value_type>(vec2),
                                static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                static_cast<unsigned int>(viennacl::traits::stride(vec2)) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("av_kernel");
      }


      ///////////////////// avbv //////////////////////////////////

      // alpha and beta on GPU
      template <typename T>
      __global__ void avbv_kernel(T * vec1,
                                  unsigned int start1,
                                  unsigned int inc1,
                                  unsigned int size1,

                                  const T * fac2,
                                  unsigned int options2,
                                  const T * vec2,
                                  unsigned int start2,
                                  unsigned int inc2,

                                  const T * fac3,
                                  unsigned int options3,
                                  const T * vec3,
                                  unsigned int start3,
                                  unsigned int inc3)
      {
        T alpha = *fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = *fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha on CPU, beta on GPU
      template <typename T>
      __global__ void avbv_kernel(T * vec1,
                                  unsigned int start1,
                                  unsigned int inc1,
                                  unsigned int size1,

                                  T fac2,
                                  unsigned int options2,
                                  const T * vec2,
                                  unsigned int start2,
                                  unsigned int inc2,

                                  const T * fac3,
                                  unsigned int options3,
                                  const T * vec3,
                                  unsigned int start3,
                                  unsigned int inc3)
      {
        T alpha = fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = *fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha on GPU, beta on CPU
      template <typename T>
      __global__ void avbv_kernel(T * vec1,
                                  unsigned int start1,
                                  unsigned int inc1,
                                  unsigned int size1,

                                  const T * fac2,
                                  unsigned int options2,
                                  const T * vec2,
                                  unsigned int start2,
                                  unsigned int inc2,

                                  T fac3,
                                  unsigned int options3,
                                  const T * vec3,
                                  unsigned int start3,
                                  unsigned int inc3)
      {
        T alpha = *fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha and beta on CPU
      template <typename T>
      __global__ void avbv_kernel(T * vec1,
                                  unsigned int start1,
                                  unsigned int inc1,
                                  unsigned int size1,

                                  T fac2,
                                  unsigned int options2,
                                  const T * vec2,
                                  unsigned int start2,
                                  unsigned int inc2,

                                  T fac3,
                                  unsigned int options3,
                                  const T * vec3,
                                  unsigned int start3,
                                  unsigned int inc3)
      {
        T alpha = fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] = vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }




      template <typename T, typename ScalarType1, typename ScalarType2>
      void avbv(vector_base<T> & vec1,
                vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha,
                vector_base<T> const & vec3, ScalarType2 const & beta,  vcl_size_t len_beta,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef T        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        unsigned int options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        value_type temporary_beta = 0;
        if (viennacl::is_cpu_scalar<ScalarType2>::value)
          temporary_beta = beta;


        avbv_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                  static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                  static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                  static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                  detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                                  options_alpha,
                                  detail::cuda_arg<value_type>(vec2),
                                  static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                  static_cast<unsigned int>(viennacl::traits::stride(vec2)),

                                  detail::cuda_arg<value_type>(detail::arg_reference(beta, temporary_beta)),
                                  options_beta,
                                  detail::cuda_arg<value_type>(vec3),
                                  static_cast<unsigned int>(viennacl::traits::start(vec3)),
                                  static_cast<unsigned int>(viennacl::traits::stride(vec3)) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("avbv_kernel");
      }


      ////////////////////////// avbv_v //////////////////////////////////////


      // alpha and beta on GPU
      template <typename T>
      __global__ void avbv_v_kernel(T * vec1,
                                    unsigned int start1,
                                    unsigned int inc1,
                                    unsigned int size1,

                                    const T * fac2,
                                    unsigned int options2,
                                    const T * vec2,
                                    unsigned int start2,
                                    unsigned int inc2,

                                    const T * fac3,
                                    unsigned int options3,
                                    const T * vec3,
                                    unsigned int start3,
                                    unsigned int inc3)
      {
        T alpha = *fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = *fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha on CPU, beta on GPU
      template <typename T>
      __global__ void avbv_v_kernel(T * vec1,
                                    unsigned int start1,
                                    unsigned int inc1,
                                    unsigned int size1,

                                    T fac2,
                                    unsigned int options2,
                                    const T * vec2,
                                    unsigned int start2,
                                    unsigned int inc2,

                                    const T * fac3,
                                    unsigned int options3,
                                    const T * vec3,
                                    unsigned int start3,
                                    unsigned int inc3)
      {
        T alpha = fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = *fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha on GPU, beta on CPU
      template <typename T>
      __global__ void avbv_v_kernel(T * vec1,
                                    unsigned int start1,
                                    unsigned int inc1,
                                    unsigned int size1,

                                    const T * fac2,
                                    unsigned int options2,
                                    const T * vec2,
                                    unsigned int start2,
                                    unsigned int inc2,

                                    T fac3,
                                    unsigned int options3,
                                    const T * vec3,
                                    unsigned int start3,
                                    unsigned int inc3)
      {
        T alpha = *fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }

      // alpha and beta on CPU
      template <typename T>
      __global__ void avbv_v_kernel(T * vec1,
                                    unsigned int start1,
                                    unsigned int inc1,
                                    unsigned int size1,

                                    T fac2,
                                    unsigned int options2,
                                    const T * vec2,
                                    unsigned int start2,
                                    unsigned int inc2,

                                    T fac3,
                                    unsigned int options3,
                                    const T * vec3,
                                    unsigned int start3,
                                    unsigned int inc3)
      {
        T alpha = fac2;
        if (options2 & (1 << 0))
          alpha = -alpha;

        T beta = fac3;
        if (options3 & (1 << 0))
          beta = -beta;

        if (options2 & (1 << 1))
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] / alpha + vec3[i*inc3+start3] * beta;
          }
        }
        else
        {
          if (options3 & (1 << 1))
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] / beta;
          }
          else
          {
            for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                              i < size1;
                              i += gridDim.x * blockDim.x)
              vec1[i*inc1+start1] += vec2[i*inc2+start2] * alpha + vec3[i*inc3+start3] * beta;
          }
        }
      }


      template <typename T, typename ScalarType1, typename ScalarType2>
      void avbv_v(vector_base<T> & vec1,
                  vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t len_alpha, bool reciprocal_alpha, bool flip_sign_alpha,
                  vector_base<T> const & vec3, ScalarType2 const & beta,  vcl_size_t len_beta,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef T        value_type;

        unsigned int options_alpha = detail::make_options(len_alpha, reciprocal_alpha, flip_sign_alpha);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<ScalarType1>::value)
          temporary_alpha = alpha;

        unsigned int options_beta  = detail::make_options(len_beta,  reciprocal_beta,  flip_sign_beta);

        value_type temporary_beta = 0;
        if (viennacl::is_cpu_scalar<ScalarType2>::value)
          temporary_beta = beta;


        avbv_v_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                    static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                    static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                    static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                    detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                                    options_alpha,
                                    detail::cuda_arg<value_type>(vec2),
                                    static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                    static_cast<unsigned int>(viennacl::traits::stride(vec2)),

                                    detail::cuda_arg<value_type>(detail::arg_reference(beta, temporary_beta)),
                                    options_beta,
                                    detail::cuda_arg<value_type>(vec3),
                                    static_cast<unsigned int>(viennacl::traits::start(vec3)),
                                    static_cast<unsigned int>(viennacl::traits::stride(vec3)) );
      }


      //////////////////////////

      template <typename T>
      __global__ void vector_assign_kernel(T * vec1,
                                           unsigned int start1,
                                           unsigned int inc1,
                                           unsigned int size1,
                                           unsigned int internal_size1,

                                           T alpha)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                          i < size1;
                          i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] =  (i < size1) ? alpha : 0;
      }

      /** @brief Assign a constant value to a vector (-range/-slice)
      *
      * @param vec1   The vector to which the value should be assigned
      * @param alpha  The value to be assigned
      * @param up_to_internal_size  Specifies whether alpha should also be written to padded memory (mostly used for clearing the whole buffer).
      */
      template <typename T, typename S1>
      void vector_assign(vector_base<T> & vec1, const S1 & alpha, bool up_to_internal_size = false)
      {
        typedef T        value_type;

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<S1>::value)
          temporary_alpha = alpha;

        unsigned int size = up_to_internal_size ? static_cast<unsigned int>(vec1.internal_size()) : static_cast<unsigned int>(viennacl::traits::size(vec1));

        vector_assign_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                           static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                           static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                           size,
                                           static_cast<unsigned int>(vec1.internal_size()),  //Note: Do NOT use traits::internal_size() here, because vector proxies don't require padding.

                                           detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("avbv_v_kernel");
      }

      //////////////////////////

      template <typename T>
      __global__ void vector_swap_kernel(T * vec1,
                                         unsigned int start1,
                                         unsigned int inc1,
                                         unsigned int size1,

                                         T * vec2,
                                         unsigned int start2,
                                         unsigned int inc2)
      {
        T tmp;
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                          i < size1;
                          i += gridDim.x * blockDim.x)
        {
          tmp = vec2[i*inc2+start2];
          vec2[i*inc2+start2] = vec1[i*inc1+start1];
          vec1[i*inc1+start1] = tmp;
        }
      }


      /** @brief Swaps the contents of two vectors, data is copied
      *
      * @param vec1   The first vector (or -range, or -slice)
      * @param vec2   The second vector (or -range, or -slice)
      */
      template <typename T>
      void vector_swap(vector_base<T> & vec1, vector_base<T> & vec2)
      {
        typedef T      value_type;

        vector_swap_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                         static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                         static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                         static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                         detail::cuda_arg<value_type>(vec2),
                                         static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                         static_cast<unsigned int>(viennacl::traits::stride(vec2)) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vector_swap_kernel");
      }

      ///////////////////////// Binary Elementwise operations /////////////

      template <typename T>
      __global__ void element_op_kernel(T * vec1,
                                         unsigned int start1,
                                         unsigned int inc1,
                                         unsigned int size1,

                                         T const * vec2,
                                         unsigned int start2,
                                         unsigned int inc2,

                                         T const * vec3,
                                         unsigned int start3,
                                         unsigned int inc3,

                                         unsigned int op_type
                                       )
      {
        if (op_type == 2)
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
          {
            vec1[i*inc1+start1] = pow(vec2[i*inc2+start2], vec3[i*inc3+start3]);
          }
        }
        else if (op_type == 1)
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
          {
            vec1[i*inc1+start1] = vec2[i*inc2+start2] / vec3[i*inc3+start3];
          }
        }
        else if (op_type == 0)
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
          {
            vec1[i*inc1+start1] = vec2[i*inc2+start2] * vec3[i*inc3+start3];
          }
        }
      }

      template <typename T>
      __global__ void element_op_int_kernel(T * vec1,
                                         unsigned int start1,
                                         unsigned int inc1,
                                         unsigned int size1,

                                         T const * vec2,
                                         unsigned int start2,
                                         unsigned int inc2,

                                         T const * vec3,
                                         unsigned int start3,
                                         unsigned int inc3,

                                         unsigned int op_type
                                       )
      {
        if (op_type == 1)
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
          {
            vec1[i*inc1+start1] = vec2[i*inc2+start2] / vec3[i*inc3+start3];
          }
        }
        else if (op_type == 0)
        {
          for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x;
                            i < size1;
                            i += gridDim.x * blockDim.x)
          {
            vec1[i*inc1+start1] = vec2[i*inc2+start2] * vec3[i*inc3+start3];
          }
        }
      }

      /** @brief Implementation of the element-wise operation v1 = v2 .* v3 and v1 = v2 ./ v3    (using MATLAB syntax)
      *
      * @param vec1   The result vector (or -range, or -slice)
      * @param proxy  The proxy object holding v2, v3 and the operation
      */
      template <typename T, typename OP>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> > const & proxy)
      {
        typedef T        value_type;

        unsigned int op_type = 2; //0: product, 1: division, 2: power
        if (viennacl::is_division<OP>::value)
          op_type = 1;
        else if (viennacl::is_product<OP>::value)
          op_type = 0;

        element_op_int_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                        static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                        detail::cuda_arg<value_type>(proxy.lhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs())),

                                        detail::cuda_arg<value_type>(proxy.rhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.rhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.rhs())),

                                        op_type
                                       );
        VIENNACL_CUDA_LAST_ERROR_CHECK("element_op_kernel");
      }

      template <typename OP>
      void element_op(vector_base<float> & vec1,
                      vector_expression<const vector_base<float>, const vector_base<float>, op_element_binary<OP> > const & proxy)
      {
        typedef float        value_type;

        unsigned int op_type = 2; //0: product, 1: division, 2: power
        if (viennacl::is_division<OP>::value)
          op_type = 1;
        else if (viennacl::is_product<OP>::value)
          op_type = 0;

        element_op_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                        static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                        detail::cuda_arg<value_type>(proxy.lhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs())),

                                        detail::cuda_arg<value_type>(proxy.rhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.rhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.rhs())),

                                        op_type
                                       );
        VIENNACL_CUDA_LAST_ERROR_CHECK("element_op_kernel");
      }

      template <typename OP>
      void element_op(vector_base<double> & vec1,
                      vector_expression<const vector_base<double>, const vector_base<double>, op_element_binary<OP> > const & proxy)
      {
        typedef double        value_type;

        unsigned int op_type = 2; //0: product, 1: division, 2: power
        if (viennacl::is_division<OP>::value)
          op_type = 1;
        else if (viennacl::is_product<OP>::value)
          op_type = 0;

        element_op_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                        static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec1)),

                                        detail::cuda_arg<value_type>(proxy.lhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs())),

                                        detail::cuda_arg<value_type>(proxy.rhs()),
                                        static_cast<unsigned int>(viennacl::traits::start(proxy.rhs())),
                                        static_cast<unsigned int>(viennacl::traits::stride(proxy.rhs())),

                                        op_type
                                       );
        VIENNACL_CUDA_LAST_ERROR_CHECK("element_op_kernel");
      }

      ///////////////////////// Unary Elementwise operations /////////////

// Note: Trying to automate things with macros or template metaprogramming failed (preprocessor with nvcc did not work as expected), so this is terribly hand-rolled code
// Question (Karl Rupp): Why is CUDA code always such a hassle when trying to use it in a library context?

      // acos
      template <typename T> __global__ void vec_element_acos_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = acos(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_acos> > const & proxy)
      {
        typedef T        value_type;

        vec_element_acos_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_acos_kernel");
      }

      // asin
      template <typename T> __global__ void vec_element_asin_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = asin(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_asin> > const & proxy)
      {
        typedef T        value_type;

        vec_element_asin_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_asin_kernel");
      }


      // atan
      template <typename T> __global__ void vec_element_atan_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = atan(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_atan> > const & proxy)
      {
        typedef T        value_type;

        vec_element_atan_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_atan_kernel");
      }


      // ceil
      template <typename T> __global__ void vec_element_ceil_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = ceil(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_ceil> > const & proxy)
      {
        typedef T        value_type;

        vec_element_ceil_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_ceil_kernel");
      }


      // cos
      template <typename T> __global__ void vec_element_cos_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = cos(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_cos> > const & proxy)
      {
        typedef T        value_type;

        vec_element_cos_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_cos_kernel");
      }


      // cosh
      template <typename T> __global__ void vec_element_cosh_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = cosh(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_cosh> > const & proxy)
      {
        typedef T        value_type;

        vec_element_cosh_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_cosh_kernel");
      }


      // exp
      template <typename T> __global__ void vec_element_exp_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = exp(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_exp> > const & proxy)
      {
        typedef T        value_type;

        vec_element_exp_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_exp_kernel");
      }


      // fabs
      template <typename T> __global__ void vec_element_fabs_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = fabs(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_fabs> > const & proxy)
      {
        typedef T        value_type;

        vec_element_fabs_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_fabs_kernel");
      }

      // abs
      template <typename T> __global__ void vec_element_abs_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = abs(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_abs> > const & proxy)
      {
        typedef T        value_type;

        vec_element_abs_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                             static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                             static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                             static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                             detail::cuda_arg<value_type>(proxy.lhs()),
                                             static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                             static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                            );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_abs_kernel");
      }



      // floor
      template <typename T> __global__ void vec_element_floor_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = floor(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_floor> > const & proxy)
      {
        typedef T        value_type;

        vec_element_floor_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_floor_kernel");
      }


      // log
      template <typename T> __global__ void vec_element_log_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = log(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_log> > const & proxy)
      {
        typedef T        value_type;

        vec_element_log_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_log_kernel");
      }


      // log10
      template <typename T> __global__ void vec_element_log10_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = log10(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_log10> > const & proxy)
      {
        typedef T        value_type;

        vec_element_log10_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_log10_kernel");
      }


      // sin
      template <typename T> __global__ void vec_element_sin_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = sin(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_sin> > const & proxy)
      {
        typedef T        value_type;

        vec_element_sin_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_sin_kernel");
      }


      // sinh
      template <typename T> __global__ void vec_element_sinh_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = sinh(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_sinh> > const & proxy)
      {
        typedef T        value_type;

        vec_element_sinh_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_sinh_kernel");
      }


      // sqrt
      template <typename T> __global__ void vec_element_sqrt_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = sqrt(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_sqrt> > const & proxy)
      {
        typedef T        value_type;

        vec_element_sqrt_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_sqrt_kernel");
      }


      // tan
      template <typename T> __global__ void vec_element_tan_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = tan(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_tan> > const & proxy)
      {
        typedef T        value_type;

        vec_element_tan_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_tan_kernel");
      }


      // tanh
      template <typename T> __global__ void vec_element_tanh_kernel(
          T       * vec1, unsigned int start1, unsigned int inc1, unsigned int size1,
          T const * vec2, unsigned int start2, unsigned int inc2)
      {
        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += gridDim.x * blockDim.x)
          vec1[i*inc1+start1] = tanh(vec2[i*inc2+start2]);
      }

      template <typename T>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<op_tanh> > const & proxy)
      {
        typedef T        value_type;

        vec_element_tanh_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                              static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                              static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                              detail::cuda_arg<value_type>(proxy.lhs()),
                                              static_cast<unsigned int>(viennacl::traits::start(proxy.lhs())),
                                              static_cast<unsigned int>(viennacl::traits::stride(proxy.lhs()))
                                             );
        VIENNACL_CUDA_LAST_ERROR_CHECK("vec_element_tanh_kernel");
      }



      ///////////////////////// Norms and inner product ///////////////////


      template <typename T>
      __global__ void inner_prod_kernel(const T * vec1,
                                        unsigned int start1,
                                        unsigned int inc1,
                                        unsigned int size1,
                                        const T * vec2,
                                        unsigned int start2,
                                        unsigned int inc2,
                                        unsigned int size2,
                                        T * group_buffer)
      {
        __shared__ T tmp_buffer[128];
        unsigned int group_start1 = (blockIdx.x * size1) / (gridDim.x) * inc1 + start1;
        unsigned int group_start2 = (blockIdx.x * size2) / (gridDim.x) * inc2 + start2;

        unsigned int group_size1 = ((blockIdx.x + 1) * size1) / (gridDim.x)
                                     - (  blockIdx.x * size1) / (gridDim.x);


        T tmp = 0;
        for (unsigned int i = threadIdx.x; i < group_size1; i += blockDim.x)
          tmp += vec1[i*inc1+group_start1] * vec2[i*inc2+group_start2];
        tmp_buffer[threadIdx.x] = tmp;

        // parallel reduction
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
            tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x+stride];
        }

        if (threadIdx.x == 0)
          group_buffer[blockIdx.x] = tmp_buffer[0];

      }



      // sums the array 'vec1' and writes to result. Makes use of a single work-group only.
      template <typename T>
      __global__ void vector_sum_kernel_floats(
                const T * vec1,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int option, //0: use fmax, 1: just sum, 2: sum and return sqrt of sum
                T * result)
      {
        __shared__ T tmp_buffer[128];
        T thread_sum = 0;
        for (unsigned int i = threadIdx.x; i<size1; i += blockDim.x)
        {
          if (option > 0)
            thread_sum += vec1[i*inc1+start1];
          else
            thread_sum = fmax(thread_sum, fabs(vec1[i*inc1+start1]));
        }

        tmp_buffer[threadIdx.x] = thread_sum;

        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
          {
            if (option > 0)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x + stride];
            else
              tmp_buffer[threadIdx.x] = fmax(tmp_buffer[threadIdx.x], tmp_buffer[threadIdx.x + stride]);
          }
        }

        if (threadIdx.x == 0)
        {
          if (option == 2)
            *result = sqrt(tmp_buffer[0]);
          else
            *result = tmp_buffer[0];
        }
      }

      template <typename T>
      __global__ void vector_sum_kernel_integers(
                const T * vec1,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int option, //0: use max, 1: just sum
                T * result)
      {
        __shared__ T tmp_buffer[128];
        T thread_sum = 0;
        for (unsigned int i = threadIdx.x; i<size1; i += blockDim.x)
        {
          if (option > 0)
            thread_sum += vec1[i*inc1+start1];
          else
            thread_sum = thread_sum > abs(vec1[i*inc1+start1]) ? thread_sum : abs(vec1[i*inc1+start1]);
        }

        tmp_buffer[threadIdx.x] = thread_sum;

        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
          {
            if (option > 0)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x + stride];
            else
              tmp_buffer[threadIdx.x] = tmp_buffer[threadIdx.x] > tmp_buffer[threadIdx.x + stride] ? tmp_buffer[threadIdx.x] : tmp_buffer[threadIdx.x + stride];
          }
        }

        if (threadIdx.x == 0)
          *result = tmp_buffer[0];
      }

      template <typename T>
      __global__ void vector_sum_kernel_unsigned_integers(
                const T * vec1,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int option, //0: use max, 1: just sum
                T * result)
      {
        __shared__ T tmp_buffer[128];
        T thread_sum = 0;
        for (unsigned int i = threadIdx.x; i<size1; i += blockDim.x)
        {
          if (option > 0)
            thread_sum += vec1[i*inc1+start1];
          else
            thread_sum = (thread_sum > vec1[i*inc1+start1]) ? thread_sum : vec1[i*inc1+start1];
        }

        tmp_buffer[threadIdx.x] = thread_sum;

        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
          {
            if (option > 0)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x + stride];
            else
              tmp_buffer[threadIdx.x] = tmp_buffer[threadIdx.x] > tmp_buffer[threadIdx.x + stride] ? tmp_buffer[threadIdx.x] : tmp_buffer[threadIdx.x + stride];
          }
        }

        if (threadIdx.x == 0)
          *result = tmp_buffer[0];
      }

      namespace detail
      {
        /** \cond */
        struct vector_sum_kernel_launcher_integers
        {
          template <typename T, typename S3>
          static void apply(vector_base<T> const & temp,
                            unsigned int option,
                            S3 & result)
          {
            typedef T        value_type;
            vector_sum_kernel_integers<<<1, 128>>>(detail::cuda_arg<value_type>(temp),
                                                  static_cast<unsigned int>(viennacl::traits::start(temp)),
                                                  static_cast<unsigned int>(viennacl::traits::stride(temp)),
                                                  static_cast<unsigned int>(viennacl::traits::size(temp)),
                                                  static_cast<unsigned int>(option),
                                                  detail::cuda_arg<value_type>(result) );
            VIENNACL_CUDA_LAST_ERROR_CHECK("vector_sum_kernel");
          }
        };

        struct vector_sum_kernel_launcher_unsigned_integers
        {
          template <typename T, typename S3>
          static void apply(vector_base<T> const & temp,
                            unsigned int option,
                            S3 & result)
          {
            typedef T        value_type;
            vector_sum_kernel_unsigned_integers<<<1, 128>>>(detail::cuda_arg<value_type>(temp),
                                                            static_cast<unsigned int>(viennacl::traits::start(temp)),
                                                            static_cast<unsigned int>(viennacl::traits::stride(temp)),
                                                            static_cast<unsigned int>(viennacl::traits::size(temp)),
                                                            static_cast<unsigned int>(option),
                                                            detail::cuda_arg<value_type>(result) );
            VIENNACL_CUDA_LAST_ERROR_CHECK("vector_sum_kernel");
          }
        };

        struct vector_sum_kernel_launcher_floats
        {
          template <typename T, typename S3>
          static void apply(vector_base<T> const & temp,
                            unsigned int option,
                            S3 & result)
          {
            typedef T        value_type;
            vector_sum_kernel_floats<<<1, 128>>>(detail::cuda_arg<value_type>(temp),
                                                  static_cast<unsigned int>(viennacl::traits::start(temp)),
                                                  static_cast<unsigned int>(viennacl::traits::stride(temp)),
                                                  static_cast<unsigned int>(viennacl::traits::size(temp)),
                                                  static_cast<unsigned int>(option),
                                                  detail::cuda_arg<value_type>(result) );
            VIENNACL_CUDA_LAST_ERROR_CHECK("vector_sum_kernel");
          }
        };

        template <typename T>
        struct vector_sum_kernel_launcher : public vector_sum_kernel_launcher_integers {};

        template <>
        struct vector_sum_kernel_launcher<unsigned char>  : public vector_sum_kernel_launcher_unsigned_integers {};

        template <>
        struct vector_sum_kernel_launcher<unsigned short>  : public vector_sum_kernel_launcher_unsigned_integers {};

        template <>
        struct vector_sum_kernel_launcher<unsigned int>  : public vector_sum_kernel_launcher_unsigned_integers {};

        template <>
        struct vector_sum_kernel_launcher<unsigned long>  : public vector_sum_kernel_launcher_unsigned_integers {};

        template <>
        struct vector_sum_kernel_launcher<float>  : public vector_sum_kernel_launcher_floats {};

        template <>
        struct vector_sum_kernel_launcher<double> : public vector_sum_kernel_launcher_floats {};

        /** \endcond */
      }


      //implementation of inner product:
      //namespace {
      /** @brief Computes the inner product of two vectors - implementation. Library users should call inner_prod(vec1, vec2).
      *
      * @param vec1 The first vector
      * @param vec2 The second vector
      * @param result The result scalar (on the gpu)
      */
      template <typename T, typename S3>
      void inner_prod_impl(vector_base<T> const & vec1,
                           vector_base<T> const & vec2,
                           S3 & result)
      {
        typedef T        value_type;

        static const unsigned int work_groups = 128;
        static viennacl::vector<value_type> temp(work_groups);

        inner_prod_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                        static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                        detail::cuda_arg<value_type>(vec2),
                                        static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec2)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec2)),
                                        detail::cuda_arg<value_type>(temp)
                                       );
        VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_kernel");

        detail::vector_sum_kernel_launcher<T>::apply(temp, 1, result);
      }


      /** @brief Computes the inner product of two vectors - implementation. Library users should call inner_prod(vec1, vec2).
      *
      * @param vec1 The first vector
      * @param vec2 The second vector
      * @param result The result scalar (on the host)
      */
      template <typename T>
      void inner_prod_cpu(vector_base<T> const & vec1,
                          vector_base<T> const & vec2,
                          T & result)
      {
        typedef T        value_type;

        const unsigned int work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        inner_prod_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                        static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                        detail::cuda_arg<value_type>(vec2),
                                        static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                        static_cast<unsigned int>(viennacl::traits::stride(vec2)),
                                        static_cast<unsigned int>(viennacl::traits::size(vec2)),
                                        detail::cuda_arg<value_type>(temp)
                                       );
        VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_kernel");

        // Now copy partial results from GPU back to CPU and run reduction there:
        std::vector<value_type> temp_cpu(work_groups);
        viennacl::fast_copy(temp.begin(), temp.end(), temp_cpu.begin());

        result = 0;
        for (typename std::vector<value_type>::const_iterator it = temp_cpu.begin(); it != temp_cpu.end(); ++it)
          result += *it;
      }

      ///////////////////////////////////

#define VIENNACL_MDOT_WORKGROUP_SIZE  128
#define VIENNACL_MDOT_WORKGROUP_NUM   128
      // M = 2:
      template <typename NumericT>
      __global__ void inner_prod_2_kernel(const NumericT *x,  unsigned int startx, unsigned int stridex, unsigned int sizex,
                                          const NumericT *y0, unsigned int start0, unsigned int stride0,
                                          const NumericT *y1, unsigned int start1, unsigned int stride1,
                                          NumericT *group_results)
      {
        __shared__ NumericT tmp_buffer[2*VIENNACL_MDOT_WORKGROUP_SIZE];
        unsigned int entries_per_thread = (sizex - 1) / (blockDim.x * gridDim.x) + 1;
        unsigned int vec_start_index = blockIdx.x * blockDim.x * entries_per_thread;
        unsigned int vec_stop_index  = min((blockIdx.x + 1) * blockDim.x * entries_per_thread, sizex); // don't go beyond size of x

        NumericT entry_x    = 0;
        NumericT group_sum0 = 0;
        NumericT group_sum1 = 0;
        for (unsigned int i = vec_start_index + threadIdx.x; i < vec_stop_index; i += blockDim.x) {
          entry_x     = x[i * stridex + startx];   // load only once from global memory!
          group_sum0 += entry_x * y0[i * stride0 + start0];
          group_sum1 += entry_x * y1[i * stride1 + start1];
        }
        tmp_buffer[threadIdx.x]              = group_sum0;
        tmp_buffer[threadIdx.x + blockDim.x] = group_sum1;

        // parallel reduction
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2) {
          __syncthreads();
          if (threadIdx.x < stride) {
            tmp_buffer[threadIdx.x             ] += tmp_buffer[threadIdx.x+stride             ];
            tmp_buffer[threadIdx.x + blockDim.x] += tmp_buffer[threadIdx.x+stride + blockDim.x];
          }
        }

        // write result of group to group_results
        if (threadIdx.x == 0) {
          group_results[blockIdx.x]             = tmp_buffer[0];
          group_results[blockIdx.x + gridDim.x] = tmp_buffer[blockDim.x];
        }
      }

      // M = 3:
      template <typename NumericT>
      __global__ void inner_prod_3_kernel(const NumericT *x,  unsigned int startx, unsigned int stridex, unsigned int sizex,
                                          const NumericT *y0, unsigned int start0, unsigned int stride0,
                                          const NumericT *y1, unsigned int start1, unsigned int stride1,
                                          const NumericT *y2, unsigned int start2, unsigned int stride2,
                                          NumericT *group_results)
      {
        __shared__ NumericT tmp_buffer[3*VIENNACL_MDOT_WORKGROUP_SIZE];
        unsigned int entries_per_thread = (sizex - 1) / (blockDim.x * gridDim.x) + 1;
        unsigned int vec_start_index = blockIdx.x * blockDim.x * entries_per_thread;
        unsigned int vec_stop_index  = min((blockIdx.x + 1) * blockDim.x * entries_per_thread, sizex); // don't go beyond vec size

        NumericT entry_x    = 0;
        NumericT group_sum0 = 0;
        NumericT group_sum1 = 0;
        NumericT group_sum2 = 0;
        for (unsigned int i = vec_start_index + threadIdx.x; i < vec_stop_index; i += blockDim.x) {
          entry_x     = x[i * stridex + startx];   // load only once from global memory!
          group_sum0 += entry_x * y0[i * stride0 + start0];
          group_sum1 += entry_x * y1[i * stride1 + start1];
          group_sum2 += entry_x * y2[i * stride2 + start2];
        }
        tmp_buffer[threadIdx.x]                  = group_sum0;
        tmp_buffer[threadIdx.x +     blockDim.x] = group_sum1;
        tmp_buffer[threadIdx.x + 2 * blockDim.x] = group_sum2;

        // parallel reduction
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2) {
          __syncthreads();
          if (threadIdx.x < stride) {
            tmp_buffer[threadIdx.x                 ] += tmp_buffer[threadIdx.x+stride                 ];
            tmp_buffer[threadIdx.x +     blockDim.x] += tmp_buffer[threadIdx.x+stride +     blockDim.x];
            tmp_buffer[threadIdx.x + 2 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 2 * blockDim.x];
          }
        }

        // write result of group to group_results
        if (threadIdx.x == 0) {
          group_results[blockIdx.x                ] = tmp_buffer[0];
          group_results[blockIdx.x +     gridDim.x] = tmp_buffer[    blockDim.x];
          group_results[blockIdx.x + 2 * gridDim.x] = tmp_buffer[2 * blockDim.x];
        }
      }

      // M = 4:
      template <typename NumericT>
      __global__ void inner_prod_4_kernel(const NumericT *x,  unsigned int startx, unsigned int stridex, unsigned int sizex,
                                          const NumericT *y0, unsigned int start0, unsigned int stride0,
                                          const NumericT *y1, unsigned int start1, unsigned int stride1,
                                          const NumericT *y2, unsigned int start2, unsigned int stride2,
                                          const NumericT *y3, unsigned int start3, unsigned int stride3,
                                          NumericT *group_results)
      {
        __shared__ NumericT tmp_buffer[4*VIENNACL_MDOT_WORKGROUP_SIZE];
        unsigned int entries_per_thread = (sizex - 1) / (blockDim.x * gridDim.x) + 1;
        unsigned int vec_start_index = blockIdx.x * blockDim.x * entries_per_thread;
        unsigned int vec_stop_index  = min((blockIdx.x + 1) * blockDim.x * entries_per_thread, sizex); // don't go beyond vec size

        NumericT entry_x    = 0;
        NumericT group_sum0 = 0;
        NumericT group_sum1 = 0;
        NumericT group_sum2 = 0;
        NumericT group_sum3 = 0;
        for (unsigned int i = vec_start_index + threadIdx.x; i < vec_stop_index; i += blockDim.x) {
          entry_x     = x[i * stridex + startx];   // load only once from global memory!
          group_sum0 += entry_x * y0[i * stride0 + start0];
          group_sum1 += entry_x * y1[i * stride1 + start1];
          group_sum2 += entry_x * y2[i * stride2 + start2];
          group_sum3 += entry_x * y3[i * stride3 + start3];
        }
        tmp_buffer[threadIdx.x]                  = group_sum0;
        tmp_buffer[threadIdx.x +     blockDim.x] = group_sum1;
        tmp_buffer[threadIdx.x + 2 * blockDim.x] = group_sum2;
        tmp_buffer[threadIdx.x + 3 * blockDim.x] = group_sum3;

        // parallel reduction
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2) {
          __syncthreads();
          if (threadIdx.x < stride) {
            tmp_buffer[threadIdx.x                 ] += tmp_buffer[threadIdx.x+stride                 ];
            tmp_buffer[threadIdx.x +     blockDim.x] += tmp_buffer[threadIdx.x+stride +     blockDim.x];
            tmp_buffer[threadIdx.x + 2 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 2 * blockDim.x];
            tmp_buffer[threadIdx.x + 3 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 3 * blockDim.x];
          }
        }

        // write result of group to group_results
        if (threadIdx.x == 0) {
          group_results[blockIdx.x                ] = tmp_buffer[0];
          group_results[blockIdx.x +     gridDim.x] = tmp_buffer[    blockDim.x];
          group_results[blockIdx.x + 2 * gridDim.x] = tmp_buffer[2 * blockDim.x];
          group_results[blockIdx.x + 3 * gridDim.x] = tmp_buffer[3 * blockDim.x];
        }
      }

      // M = 8:
      template <typename NumericT>
      __global__ void inner_prod_8_kernel(const NumericT *x,  unsigned int startx, unsigned int stridex, unsigned int sizex,
                                          const NumericT *y0, unsigned int start0, unsigned int stride0,
                                          const NumericT *y1, unsigned int start1, unsigned int stride1,
                                          const NumericT *y2, unsigned int start2, unsigned int stride2,
                                          const NumericT *y3, unsigned int start3, unsigned int stride3,
                                          const NumericT *y4, unsigned int start4, unsigned int stride4,
                                          const NumericT *y5, unsigned int start5, unsigned int stride5,
                                          const NumericT *y6, unsigned int start6, unsigned int stride6,
                                          const NumericT *y7, unsigned int start7, unsigned int stride7,
                                          NumericT *group_results)
      {
        __shared__ NumericT tmp_buffer[8*VIENNACL_MDOT_WORKGROUP_SIZE];
        unsigned int entries_per_thread = (sizex - 1) / (blockDim.x * gridDim.x) + 1;
        unsigned int vec_start_index = blockIdx.x * blockDim.x * entries_per_thread;
        unsigned int vec_stop_index  = min((blockIdx.x + 1) * blockDim.x * entries_per_thread, sizex); // don't go beyond vec size

        NumericT entry_x    = 0;
        NumericT group_sum0 = 0;
        NumericT group_sum1 = 0;
        NumericT group_sum2 = 0;
        NumericT group_sum3 = 0;
        NumericT group_sum4 = 0;
        NumericT group_sum5 = 0;
        NumericT group_sum6 = 0;
        NumericT group_sum7 = 0;
        for (unsigned int i = vec_start_index + threadIdx.x; i < vec_stop_index; i += blockDim.x) {
          entry_x     = x[i * stridex + startx];   // load only once from global memory!
          group_sum0 += entry_x * y0[i * stride0 + start0];
          group_sum1 += entry_x * y1[i * stride1 + start1];
          group_sum2 += entry_x * y2[i * stride2 + start2];
          group_sum3 += entry_x * y3[i * stride3 + start3];
          group_sum4 += entry_x * y4[i * stride4 + start4];
          group_sum5 += entry_x * y5[i * stride5 + start5];
          group_sum6 += entry_x * y6[i * stride6 + start6];
          group_sum7 += entry_x * y7[i * stride7 + start7];
        }
        tmp_buffer[threadIdx.x]                  = group_sum0;
        tmp_buffer[threadIdx.x +     blockDim.x] = group_sum1;
        tmp_buffer[threadIdx.x + 2 * blockDim.x] = group_sum2;
        tmp_buffer[threadIdx.x + 3 * blockDim.x] = group_sum3;
        tmp_buffer[threadIdx.x + 4 * blockDim.x] = group_sum4;
        tmp_buffer[threadIdx.x + 5 * blockDim.x] = group_sum5;
        tmp_buffer[threadIdx.x + 6 * blockDim.x] = group_sum6;
        tmp_buffer[threadIdx.x + 7 * blockDim.x] = group_sum7;

        // parallel reduction
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2) {
          __syncthreads();
          if (threadIdx.x < stride) {
            tmp_buffer[threadIdx.x                 ] += tmp_buffer[threadIdx.x+stride                 ];
            tmp_buffer[threadIdx.x +     blockDim.x] += tmp_buffer[threadIdx.x+stride +     blockDim.x];
            tmp_buffer[threadIdx.x + 2 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 2 * blockDim.x];
            tmp_buffer[threadIdx.x + 3 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 3 * blockDim.x];
            tmp_buffer[threadIdx.x + 4 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 4 * blockDim.x];
            tmp_buffer[threadIdx.x + 5 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 5 * blockDim.x];
            tmp_buffer[threadIdx.x + 6 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 6 * blockDim.x];
            tmp_buffer[threadIdx.x + 7 * blockDim.x] += tmp_buffer[threadIdx.x+stride + 7 * blockDim.x];
          }
        }

        // write result of group to group_results
        if (threadIdx.x == 0) {
          group_results[blockIdx.x                ] = tmp_buffer[0];
          group_results[blockIdx.x +     gridDim.x] = tmp_buffer[    blockDim.x];
          group_results[blockIdx.x + 2 * gridDim.x] = tmp_buffer[2 * blockDim.x];
          group_results[blockIdx.x + 3 * gridDim.x] = tmp_buffer[3 * blockDim.x];
          group_results[blockIdx.x + 4 * gridDim.x] = tmp_buffer[4 * blockDim.x];
          group_results[blockIdx.x + 5 * gridDim.x] = tmp_buffer[5 * blockDim.x];
          group_results[blockIdx.x + 6 * gridDim.x] = tmp_buffer[6 * blockDim.x];
          group_results[blockIdx.x + 7 * gridDim.x] = tmp_buffer[7 * blockDim.x];
        }
      }

      // sums the array 'vec1' and writes to result. Makes use of a single work-group only.
      template <typename T>
      __global__ void vector_multi_sum_kernel(
                T const * vec1,
                T * result,
                unsigned int start_result,
                unsigned int inc_result)
      {
        __shared__ T tmp_buffer[VIENNACL_MDOT_WORKGROUP_SIZE];

        tmp_buffer[threadIdx.x] = vec1[threadIdx.x + blockIdx.x * VIENNACL_MDOT_WORKGROUP_SIZE];

        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
            tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x + stride];
        }

        if (threadIdx.x == 0)
          result[start_result + inc_result * blockIdx.x] = tmp_buffer[0];
      }

      template <typename T>
      void inner_prod_impl(vector_base<T> const & x,
                           vector_tuple<T> const & vec_tuple,
                           vector_base<T> & result)
      {
        typedef T        value_type;

        static viennacl::vector<value_type> temp(8 * VIENNACL_MDOT_WORKGROUP_NUM);

        vcl_size_t current_index = 0;
        while (vec_tuple.const_size() > current_index)
        {
          switch (vec_tuple.const_size() - current_index)
          {
            case 7:
            case 6:
            case 5:
            case 4:
            {
              vector_base<T> const & y0 = vec_tuple.const_at(current_index);
              vector_base<T> const & y1 = vec_tuple.const_at(current_index + 1);
              vector_base<T> const & y2 = vec_tuple.const_at(current_index + 2);
              vector_base<T> const & y3 = vec_tuple.const_at(current_index + 3);

              inner_prod_4_kernel<<<VIENNACL_MDOT_WORKGROUP_NUM,
                                    VIENNACL_MDOT_WORKGROUP_SIZE>>>( detail::cuda_arg<value_type>(x),
                                                                     static_cast<unsigned int>(viennacl::traits::start(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::size(x)),
                                                                     detail::cuda_arg<value_type>(y0),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y0)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y0)),
                                                                     detail::cuda_arg<value_type>(y1),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y1)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y1)),
                                                                     detail::cuda_arg<value_type>(y2),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y2)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y2)),
                                                                     detail::cuda_arg<value_type>(y3),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y3)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y3)),
                                                                     detail::cuda_arg<value_type>(temp)
                                                                    );
              VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_4_kernel");
              vector_multi_sum_kernel<<<4, VIENNACL_MDOT_WORKGROUP_NUM>>>(detail::cuda_arg<value_type>(temp),
                                                                          detail::cuda_arg<value_type>(result),
                                                                          static_cast<unsigned int>(viennacl::traits::start(result) + viennacl::traits::stride(result) * current_index),
                                                                          static_cast<unsigned int>(viennacl::traits::stride(result))
                                                                         );
              VIENNACL_CUDA_LAST_ERROR_CHECK("vector_multi_sum_kernel");
            }
              current_index += 4;
              break;
            case 3:
            {
              vector_base<T> const & y0 = vec_tuple.const_at(current_index);
              vector_base<T> const & y1 = vec_tuple.const_at(current_index + 1);
              vector_base<T> const & y2 = vec_tuple.const_at(current_index + 2);

              inner_prod_3_kernel<<<VIENNACL_MDOT_WORKGROUP_NUM,
                                    VIENNACL_MDOT_WORKGROUP_SIZE>>>( detail::cuda_arg<value_type>(x),
                                                                     static_cast<unsigned int>(viennacl::traits::start(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::size(x)),
                                                                     detail::cuda_arg<value_type>(y0),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y0)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y0)),
                                                                     detail::cuda_arg<value_type>(y1),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y1)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y1)),
                                                                     detail::cuda_arg<value_type>(y2),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y2)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y2)),
                                                                     detail::cuda_arg<value_type>(temp)
                                                                    );
              VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_3_kernel");
              vector_multi_sum_kernel<<<3, VIENNACL_MDOT_WORKGROUP_NUM>>>(detail::cuda_arg<value_type>(temp),
                                                                          detail::cuda_arg<value_type>(result),
                                                                          static_cast<unsigned int>(viennacl::traits::start(result) + viennacl::traits::stride(result) * current_index),
                                                                          static_cast<unsigned int>(viennacl::traits::stride(result))
                                                                         );
              VIENNACL_CUDA_LAST_ERROR_CHECK("vector_multi_sum_kernel");
            }
              current_index += 3;
              break;
            case 2:
            {
              vector_base<T> const & y0 = vec_tuple.const_at(current_index);
              vector_base<T> const & y1 = vec_tuple.const_at(current_index + 1);

              inner_prod_2_kernel<<<VIENNACL_MDOT_WORKGROUP_NUM,
                                    VIENNACL_MDOT_WORKGROUP_SIZE>>>( detail::cuda_arg<value_type>(x),
                                                                     static_cast<unsigned int>(viennacl::traits::start(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::size(x)),
                                                                     detail::cuda_arg<value_type>(y0),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y0)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y0)),
                                                                     detail::cuda_arg<value_type>(y1),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y1)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y1)),
                                                                     detail::cuda_arg<value_type>(temp)
                                                                    );
              VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_2_kernel");
              vector_multi_sum_kernel<<<2, VIENNACL_MDOT_WORKGROUP_NUM>>>(detail::cuda_arg<value_type>(temp),
                                                                          detail::cuda_arg<value_type>(result),
                                                                          static_cast<unsigned int>(viennacl::traits::start(result) + viennacl::traits::stride(result) * current_index),
                                                                          static_cast<unsigned int>(viennacl::traits::stride(result))
                                                                         );
              VIENNACL_CUDA_LAST_ERROR_CHECK("vector_multi_sum_kernel");
            }
              current_index += 2;
              break;
            case 1:
            {
              vector_base<T> const & y0 = vec_tuple.const_at(current_index);
              inner_prod_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(x),
                                              static_cast<unsigned int>(viennacl::traits::start(x)),
                                              static_cast<unsigned int>(viennacl::traits::stride(x)),
                                              static_cast<unsigned int>(viennacl::traits::size(x)),
                                              detail::cuda_arg<value_type>(y0),
                                              static_cast<unsigned int>(viennacl::traits::start(y0)),
                                              static_cast<unsigned int>(viennacl::traits::stride(y0)),
                                              static_cast<unsigned int>(viennacl::traits::size(y0)),
                                              detail::cuda_arg<value_type>(temp)
                                             );
              VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_kernel");

              vector_multi_sum_kernel<<<1, 128>>>(detail::cuda_arg<value_type>(temp),
                                                  detail::cuda_arg<value_type>(result),
                                                  static_cast<unsigned int>(viennacl::traits::start(result) + viennacl::traits::stride(result) * current_index),
                                                  static_cast<unsigned int>(viennacl::traits::stride(result))
                                                 );
              VIENNACL_CUDA_LAST_ERROR_CHECK("vector_multi_sum_kernel");
            }
              current_index += 1;
              break;

            default:
            {
              vector_base<T> const & y0 = vec_tuple.const_at(current_index);
              vector_base<T> const & y1 = vec_tuple.const_at(current_index + 1);
              vector_base<T> const & y2 = vec_tuple.const_at(current_index + 2);
              vector_base<T> const & y3 = vec_tuple.const_at(current_index + 3);
              vector_base<T> const & y4 = vec_tuple.const_at(current_index + 4);
              vector_base<T> const & y5 = vec_tuple.const_at(current_index + 5);
              vector_base<T> const & y6 = vec_tuple.const_at(current_index + 6);
              vector_base<T> const & y7 = vec_tuple.const_at(current_index + 7);

              inner_prod_8_kernel<<<VIENNACL_MDOT_WORKGROUP_NUM,
                                    VIENNACL_MDOT_WORKGROUP_SIZE>>>( detail::cuda_arg<value_type>(x),
                                                                     static_cast<unsigned int>(viennacl::traits::start(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(x)),
                                                                     static_cast<unsigned int>(viennacl::traits::size(x)),
                                                                     detail::cuda_arg<value_type>(y0),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y0)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y0)),
                                                                     detail::cuda_arg<value_type>(y1),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y1)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y1)),
                                                                     detail::cuda_arg<value_type>(y2),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y2)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y2)),
                                                                     detail::cuda_arg<value_type>(y3),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y3)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y3)),
                                                                     detail::cuda_arg<value_type>(y4),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y4)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y4)),
                                                                     detail::cuda_arg<value_type>(y5),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y5)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y5)),
                                                                     detail::cuda_arg<value_type>(y6),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y6)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y6)),
                                                                     detail::cuda_arg<value_type>(y7),
                                                                     static_cast<unsigned int>(viennacl::traits::start(y7)),
                                                                     static_cast<unsigned int>(viennacl::traits::stride(y7)),
                                                                     detail::cuda_arg<value_type>(temp)
                                                                    );
              VIENNACL_CUDA_LAST_ERROR_CHECK("inner_prod_8_kernel");
              vector_multi_sum_kernel<<<8, VIENNACL_MDOT_WORKGROUP_NUM>>>(detail::cuda_arg<value_type>(temp),
                                                                          detail::cuda_arg<value_type>(result),
                                                                          static_cast<unsigned int>(viennacl::traits::start(result) + viennacl::traits::stride(result) * current_index),
                                                                          static_cast<unsigned int>(viennacl::traits::stride(result))
                                                                         );
              VIENNACL_CUDA_LAST_ERROR_CHECK("vector_multi_sum_kernel");
            }
              current_index += 8;
              break;
          }
        }
      }

#undef VIENNACL_MDOT_WORKGROUP_NUM
#undef VIENNACL_MDOT_WORKGROUP_SIZE

      ///////////////////////////////////

      template <typename T>
      __global__ void norm_kernel_floats(
                 const T * vec,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int norm_selector,
                T * group_buffer)
      {
        __shared__ T tmp_buffer[128];

        T tmp = 0;
        unsigned int work_per_thread = (size1 - 1) / (gridDim.x * blockDim.x) + 1;
        unsigned int group_start = blockIdx.x * work_per_thread * blockDim.x;
        unsigned int group_stop  = (blockIdx.x + 1) * work_per_thread * blockDim.x;
        group_stop = (group_stop > size1) ? size1 : group_stop;

        if (norm_selector == 1) //norm_1
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp += fabs(vec[i*inc1 + start1]);
        }
        else if (norm_selector == 2) //norm_2
        {
          T vec_entry = 0;
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
          {
            vec_entry = vec[i*inc1 + start1];
            tmp += vec_entry * vec_entry;
          }
        }
        else if (norm_selector == 0) //norm_inf
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp = fmax(fabs(vec[i*inc1 + start1]), tmp);
        }

        tmp_buffer[threadIdx.x] = tmp;

        if (norm_selector > 0) //parallel reduction for norm_1 or norm_2:
        {
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x+stride];
          }
        }
        else
        {
          //norm_inf:
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] = fmax(tmp_buffer[threadIdx.x], tmp_buffer[threadIdx.x+stride]);
          }
        }

        if (threadIdx.x == 0)
          group_buffer[blockIdx.x] = tmp_buffer[0];
      }

      template <typename T>
      __global__ void norm_kernel_integers(
                 const T * vec,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int norm_selector,
                T * group_buffer)
      {
        __shared__ T tmp_buffer[128];

        T tmp = 0;
        unsigned int work_per_thread = (size1 - 1) / (gridDim.x * blockDim.x) + 1;
        unsigned int group_start = blockIdx.x * work_per_thread * blockDim.x;
        unsigned int group_stop  = (blockIdx.x + 1) * work_per_thread * blockDim.x;
        group_stop = (group_stop > size1) ? size1 : group_stop;

        if (norm_selector == 1) //norm_1
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp += abs(vec[i*inc1 + start1]);
        }
        else if (norm_selector == 0) //norm_inf
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp = (tmp > abs(vec[i*inc1 + start1])) ? tmp : abs(vec[i*inc1 + start1]);
        }

        tmp_buffer[threadIdx.x] = tmp;

        if (norm_selector > 0) //parallel reduction for norm_1 or norm_2:
        {
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x+stride];
          }
        }
        else
        {
          //norm_inf:
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] = (tmp_buffer[threadIdx.x] > tmp_buffer[threadIdx.x+stride]) ? tmp_buffer[threadIdx.x] : tmp_buffer[threadIdx.x+stride];
          }
        }

        if (threadIdx.x == 0)
          group_buffer[blockIdx.x] = tmp_buffer[0];
      }

      template <typename T>
      __global__ void norm_kernel_unsigned_integers(
                 const T * vec,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                unsigned int norm_selector,
                T * group_buffer)
      {
        __shared__ T tmp_buffer[128];

        T tmp = 0;
        unsigned int work_per_thread = (size1 - 1) / (gridDim.x * blockDim.x) + 1;
        unsigned int group_start = blockIdx.x * work_per_thread * blockDim.x;
        unsigned int group_stop  = (blockIdx.x + 1) * work_per_thread * blockDim.x;
        group_stop = (group_stop > size1) ? size1 : group_stop;

        if (norm_selector == 1) //norm_1
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp += vec[i*inc1 + start1];
        }
        else if (norm_selector == 0) //norm_inf
        {
          for (unsigned int i = group_start + threadIdx.x; i < group_stop; i += blockDim.x)
            tmp = (tmp > vec[i*inc1 + start1]) ? tmp : vec[i*inc1 + start1];
        }

        tmp_buffer[threadIdx.x] = tmp;

        if (norm_selector > 0) //parallel reduction for norm_1 or norm_2:
        {
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] += tmp_buffer[threadIdx.x+stride];
          }
        }
        else
        {
          //norm_inf:
          for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
          {
            __syncthreads();
            if (threadIdx.x < stride)
              tmp_buffer[threadIdx.x] = (tmp_buffer[threadIdx.x] > tmp_buffer[threadIdx.x+stride]) ? tmp_buffer[threadIdx.x] : tmp_buffer[threadIdx.x+stride];
          }
        }

        if (threadIdx.x == 0)
          group_buffer[blockIdx.x] = tmp_buffer[0];
      }

      /** \cond */
      namespace detail
      {
        struct norm_kernel_launcher_integers
        {
          template <typename T>
          static void apply(vector_base<T> const & vec1,
                            vector_base<T> & temp,
                            unsigned int option)
          {
            typedef T        value_type;
            norm_kernel_integers<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                               static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                               static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                               static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                               static_cast<unsigned int>(option),
                                               detail::cuda_arg<value_type>(temp)
                                              );
            VIENNACL_CUDA_LAST_ERROR_CHECK("norm_kernel");
          }
        };

        struct norm_kernel_launcher_unsigned_integers
        {
          template <typename T>
          static void apply(vector_base<T> const & vec1,
                            vector_base<T> & temp,
                            unsigned int option)
          {
            typedef T        value_type;
            norm_kernel_unsigned_integers<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                                       static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                                       static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                                       static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                                       static_cast<unsigned int>(option),
                                                       detail::cuda_arg<value_type>(temp)
                                                      );
            VIENNACL_CUDA_LAST_ERROR_CHECK("norm_kernel");
          }
        };


        struct norm_kernel_launcher_floats
        {
          template <typename T>
          static void apply(vector_base<T> const & vec1,
                            vector_base<T> & temp,
                            unsigned int option)
          {
            typedef T        value_type;
            norm_kernel_floats<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                             static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                             static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                             static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                             static_cast<unsigned int>(option),
                                             detail::cuda_arg<value_type>(temp)
                                            );
            VIENNACL_CUDA_LAST_ERROR_CHECK("norm_kernel");
          }
        };

        template <typename T>
        struct norm_kernel_launcher : public norm_kernel_launcher_integers {};

        template <>
        struct norm_kernel_launcher<unsigned char>  : public norm_kernel_launcher_unsigned_integers {};

        template <>
        struct norm_kernel_launcher<unsigned short>  : public norm_kernel_launcher_unsigned_integers {};

        template <>
        struct norm_kernel_launcher<unsigned int>  : public norm_kernel_launcher_unsigned_integers {};

        template <>
        struct norm_kernel_launcher<unsigned long>  : public norm_kernel_launcher_unsigned_integers {};

        template <>
        struct norm_kernel_launcher<float>  : public norm_kernel_launcher_floats {};

        template <>
        struct norm_kernel_launcher<double> : public norm_kernel_launcher_floats {};

      }
      /** \endcond */


      /** @brief Computes the l^1-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_1_impl(vector_base<T> const & vec1,
                       scalar<T> & result)
      {
        typedef T        value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 1);
        detail::vector_sum_kernel_launcher<T>::apply(temp, 1, result);
      }

      /** @brief Computes the l^1-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_1_cpu(vector_base<T> const & vec1,
                      T & result)
      {
        typedef T        value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 1);

        // Now copy partial results from GPU back to CPU and run reduction there:
        std::vector<value_type> temp_cpu(work_groups);
        viennacl::fast_copy(temp.begin(), temp.end(), temp_cpu.begin());

        result = 0;
        for (typename std::vector<value_type>::const_iterator it = temp_cpu.begin(); it != temp_cpu.end(); ++it)
          result += *it;
      }

      ///// norm_2

      /** @brief Computes the l^2-norm of a vector - implementation
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_2_impl(vector_base<T> const & vec1,
                       scalar<T> & result)
      {
        typedef T       value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 2);

        detail::vector_sum_kernel_launcher<T>::apply(temp, 2, result);
      }

      /** @brief Computes the l^2-norm of a vector - implementation
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_2_cpu(vector_base<T> const & vec1,
                      T & result)
      {
        typedef T        value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 2);

        std::vector<value_type> temp_cpu(work_groups);
        viennacl::fast_copy(temp.begin(), temp.end(), temp_cpu.begin());

        result = 0;
        for (typename std::vector<value_type>::const_iterator it = temp_cpu.begin(); it != temp_cpu.end(); ++it)
          result += *it;
        result = std::sqrt(result);
      }


      ////// norm_inf

      /** @brief Computes the supremum-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_inf_impl(vector_base<T> const & vec1,
                         scalar<T> & result)
      {
        typedef T      value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 0);
        detail::vector_sum_kernel_launcher<T>::apply(temp, 0, result);
      }



      /** @brief Computes the supremum-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T>
      void norm_inf_cpu(vector_base<T> const & vec1,
                        T & result)
      {
        typedef T        value_type;

        vcl_size_t work_groups = 128;
        viennacl::vector<value_type> temp(work_groups);

        detail::norm_kernel_launcher<T>::apply(vec1, temp, 0);

        std::vector<value_type> temp_cpu(work_groups);
        viennacl::fast_copy(temp.begin(), temp.end(), temp_cpu.begin());

        result = 0;
        for (typename std::vector<value_type>::const_iterator it = temp_cpu.begin(); it != temp_cpu.end(); ++it)
          result = std::max(result, *it);
      }


      //////////////////////////////////////



      //index_norm_inf:

      // fixes the problem of not having (f)abs available in a consistent manner
      template <typename T>
      __device__ T              cuda_abs(T val) { return (val < 0) ? -val : val; }
      __device__ inline unsigned long  cuda_abs(unsigned long  val) { return val; }
      __device__ inline unsigned int   cuda_abs(unsigned int   val) { return val; }
      __device__ inline unsigned short cuda_abs(unsigned short val) { return val; }
      __device__ inline unsigned char  cuda_abs(unsigned char  val) { return val; }

      template <typename T>
      __global__ void index_norm_inf_kernel(const T * vec,
                                            unsigned int start1,
                                            unsigned int inc1,
                                            unsigned int size1,
                                            unsigned int * result)
      {
        __shared__ T float_buffer[128];
        __shared__ unsigned int index_buffer[128];

        float_buffer[threadIdx.x] = 0;
        index_buffer[threadIdx.x] = 0;

        //step 1: fill buffer:
        T cur_max = (T)0;
        T tmp;
        for (unsigned int i = threadIdx.x; i < size1; i += blockDim.x)
        {
          tmp = vec[i*inc1+start1];
          tmp = cuda_abs(tmp);
          if (cur_max < tmp)
          {
            float_buffer[threadIdx.x] = tmp;
            index_buffer[threadIdx.x] = i;
            cur_max = tmp;
          }
        }

        //step 2: parallel reduction:
        for (unsigned int stride = blockDim.x/2; stride > 0; stride /= 2)
        {
          __syncthreads();
          if (threadIdx.x < stride)
          {
            //find the first occurring index
            if (float_buffer[threadIdx.x] < float_buffer[threadIdx.x+stride])
            {
              index_buffer[threadIdx.x] = index_buffer[threadIdx.x+stride];
              float_buffer[threadIdx.x] = float_buffer[threadIdx.x+stride];
            }
          }
        }

        if (threadIdx.x == 0)
          *result = index_buffer[0];
      }

      //This function should return a CPU scalar, otherwise statements like
      // vcl_rhs[index_norm_inf(vcl_rhs)]
      // are ambiguous
      /** @brief Computes the index of the first entry that is equal to the supremum-norm in modulus.
      *
      * @param vec1 The vector
      * @return The result. Note that the result must be a CPU scalar (unsigned int), since gpu scalars are floating point types.
      */
      template <typename T>
      vcl_size_t index_norm_inf(vector_base<T> const & vec1)
      {
        typedef T       value_type;

        viennacl::backend::mem_handle h;
        viennacl::backend::memory_create(h, sizeof(unsigned int), viennacl::traits::context(vec1));

        index_norm_inf_kernel<<<1, 128>>>(detail::cuda_arg<value_type>(vec1),
                                          static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                          static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                          static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                          //detail::cuda_arg<unsigned int>(h.cuda_handle())
                                          reinterpret_cast<unsigned int *>(h.cuda_handle().get())
                                        );
        VIENNACL_CUDA_LAST_ERROR_CHECK("index_norm_inf_kernel");

        unsigned int ret = 0;
        viennacl::backend::memory_read(h, 0, sizeof(unsigned int), &ret);
        return static_cast<vcl_size_t>(ret);
      }

      ///////////////////////////////////////////

      template <typename T>
      __global__ void plane_rotation_kernel(
                T * vec1,
                unsigned int start1,
                unsigned int inc1,
                unsigned int size1,
                T * vec2,
                unsigned int start2,
                unsigned int inc2,
                unsigned int size2,
                T alpha,
                T beta)
      {
        T tmp1 = 0;
        T tmp2 = 0;

        for (unsigned int i = blockDim.x * blockIdx.x + threadIdx.x; i < size1; i += blockDim.x * gridDim.x)
        {
          tmp1 = vec1[i*inc1+start1];
          tmp2 = vec2[i*inc2+start2];

          vec1[i*inc1+start1] = alpha * tmp1 + beta * tmp2;
          vec2[i*inc2+start2] = alpha * tmp2 - beta * tmp1;
        }

      }

      /** @brief Computes a plane rotation of two vectors.
      *
      * Computes (x,y) <- (alpha * x + beta * y, -beta * x + alpha * y)
      *
      * @param vec1   The first vector
      * @param vec2   The second vector
      * @param alpha  The first transformation coefficient
      * @param beta   The second transformation coefficient
      */
      template <typename T>
      void plane_rotation(vector_base<T> & vec1,
                          vector_base<T> & vec2,
                          T alpha, T beta)
      {
        typedef T     value_type;

        value_type temporary_alpha = 0;
        if (viennacl::is_cpu_scalar<value_type>::value)
          temporary_alpha = alpha;

        value_type temporary_beta = 0;
        if (viennacl::is_cpu_scalar<value_type>::value)
          temporary_beta = beta;

        plane_rotation_kernel<<<128, 128>>>(detail::cuda_arg<value_type>(vec1),
                                            static_cast<unsigned int>(viennacl::traits::start(vec1)),
                                            static_cast<unsigned int>(viennacl::traits::stride(vec1)),
                                            static_cast<unsigned int>(viennacl::traits::size(vec1)),
                                            detail::cuda_arg<value_type>(vec2),
                                            static_cast<unsigned int>(viennacl::traits::start(vec2)),
                                            static_cast<unsigned int>(viennacl::traits::stride(vec2)),
                                            static_cast<unsigned int>(viennacl::traits::size(vec2)),
                                            detail::cuda_arg<value_type>(detail::arg_reference(alpha, temporary_alpha)),
                                            detail::cuda_arg<value_type>(detail::arg_reference(beta, temporary_beta)) );
        VIENNACL_CUDA_LAST_ERROR_CHECK("plane_rotation_kernel");
      }

    } //namespace opencl
  } //namespace linalg
} //namespace viennacl


#endif
