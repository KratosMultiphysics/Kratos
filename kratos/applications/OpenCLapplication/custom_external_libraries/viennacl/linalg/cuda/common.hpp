#ifndef VIENNACL_LINALG_CUDA_COMMON_HPP_
#define VIENNACL_LINALG_CUDA_COMMON_HPP_

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

/** @file viennacl/linalg/cuda/common.hpp
    @brief Common routines for CUDA execution
*/

#include "viennacl/traits/handle.hpp"

#define VIENNACL_CUDA_LAST_ERROR_CHECK(message)  detail::cuda_last_error_check (message, __FILE__, __LINE__)

namespace viennacl
{
  namespace linalg
  {
    namespace cuda
    {
      namespace detail
      {
        inline unsigned int make_options(vcl_size_t length, bool reciprocal, bool flip_sign)
        {
          return static_cast<unsigned int>( ((length > 1) ? (static_cast<unsigned int>(length) << 2) : 0) + (reciprocal ? 2 : 0) + (flip_sign ? 1 : 0) );
        }

        inline void cuda_last_error_check(const char * message, const char * file, const int line )
        {
          cudaError_t error_code = cudaGetLastError();

          if(cudaSuccess != error_code)
          {
            std::cerr << file << "(" << line << "): " << ": getLastCudaError() CUDA error " << error_code << ": " << cudaGetErrorString( error_code ) << " @ " << message << std::endl;
            throw "CUDA error";
          }
        }

        template <typename T, typename U>
        T * cuda_arg(vector_base<U> & obj)
        {
          return reinterpret_cast<T *>(viennacl::traits::handle(obj).cuda_handle().get());
        }

        template <typename T, typename U>
        const T * cuda_arg(vector_base<U> const & obj)
        {
          return reinterpret_cast<const T *>(viennacl::traits::handle(obj).cuda_handle().get());
        }

        template <typename NumericT, typename F>
        NumericT * cuda_arg(matrix_base<NumericT, F> & obj)
        {
          return reinterpret_cast<NumericT *>(viennacl::traits::handle(obj).cuda_handle().get());
        }

        template <typename NumericT, typename F>
        const NumericT * cuda_arg(matrix_base<NumericT, F> const & obj)
        {
          return reinterpret_cast<const NumericT *>(viennacl::traits::handle(obj).cuda_handle().get());
        }


        template <typename ScalarType, typename T>
        typename viennacl::enable_if< viennacl::is_scalar<T>::value,
                                      ScalarType *>::type
        cuda_arg(T & obj)
        {
          return reinterpret_cast<ScalarType *>(viennacl::traits::handle(obj).cuda_handle().get());
        }

        template <typename ScalarType, typename T>
        typename viennacl::enable_if< viennacl::is_scalar<T>::value,
                                      const ScalarType *>::type
        cuda_arg(T const & obj)
        {
          return reinterpret_cast<const ScalarType *>(viennacl::traits::handle(obj).cuda_handle().get());
        }

        template <typename ScalarType>
        ScalarType *  cuda_arg(viennacl::backend::mem_handle::cuda_handle_type & h)
        {
          return reinterpret_cast<ScalarType *>(h.get());
        }

        template <typename ScalarType>
        ScalarType const *  cuda_arg(viennacl::backend::mem_handle::cuda_handle_type const & h)
        {
          return reinterpret_cast<const ScalarType *>(h.get());
        }

        //template <typename ScalarType>
        //ScalarType cuda_arg(ScalarType const & val)  { return val; }

        inline unsigned int cuda_arg(unsigned int val)  { return val; }

        template <typename T> char           cuda_arg(char val)           { return val; }
        template <typename T> unsigned char  cuda_arg(unsigned char val)  { return val; }

        template <typename T> short          cuda_arg(short val)          { return val; }
        template <typename T> unsigned short cuda_arg(unsigned short val) { return val; }

        template <typename T> int            cuda_arg(int val)            { return val; }
        template <typename T> unsigned int   cuda_arg(unsigned int val)   { return val; }

        template <typename T> long           cuda_arg(long val)           { return val; }
        template <typename T> unsigned long  cuda_arg(unsigned long val)  { return val; }

        template <typename T> float          cuda_arg(float val)          { return val; }
        template <typename T> double         cuda_arg(double val)         { return val; }

        template <typename T, typename U>
        typename viennacl::backend::mem_handle::cuda_handle_type & arg_reference(viennacl::scalar<T> & s, U) { return s.handle().cuda_handle(); }

        template <typename T, typename U>
        typename viennacl::backend::mem_handle::cuda_handle_type const & arg_reference(viennacl::scalar<T> const & s, U) { return s.handle().cuda_handle(); }

        // all other cases where T is not a ViennaCL scalar
        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      char const &>::type
        arg_reference(T, char const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      unsigned char const &>::type
        arg_reference(T, unsigned char const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      short const &>::type
        arg_reference(T, short const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      unsigned short const &>::type
        arg_reference(T, unsigned short const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      int const &>::type
        arg_reference(T, int const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      unsigned int const &>::type
        arg_reference(T, unsigned int const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      long const &>::type
        arg_reference(T, long const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      unsigned long const &>::type
        arg_reference(T, unsigned long const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      float const &>::type
        arg_reference(T, float const & val)  { return val; }

        template <typename T>
        typename viennacl::enable_if< viennacl::is_cpu_scalar<T>::value,
                                      double const &>::type
        arg_reference(T, double const & val)  { return val; }
      } //namespace detail

    } //namespace cuda
  } //namespace linalg
} //namespace viennacl


#endif
