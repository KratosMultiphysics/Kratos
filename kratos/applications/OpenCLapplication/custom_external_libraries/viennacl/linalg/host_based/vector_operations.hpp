#ifndef VIENNACL_LINALG_HOST_BASED_VECTOR_OPERATIONS_HPP_
#define VIENNACL_LINALG_HOST_BASED_VECTOR_OPERATIONS_HPP_

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

/** @file viennacl/linalg/host_based/vector_operations.hpp
    @brief Implementations of vector operations using a plain single-threaded or OpenMP-enabled execution on CPU
*/

#include <cmath>
#include <algorithm>  //for std::max and std::min

#include "viennacl/forwards.h"
#include "viennacl/scalar.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/linalg/host_based/common.hpp"
#include "viennacl/linalg/detail/op_applier.hpp"
#include "viennacl/traits/stride.hpp"


// Minimum vector size for using OpenMP on vector operations:
#ifndef VIENNACL_OPENMP_VECTOR_MIN_SIZE
  #define VIENNACL_OPENMP_VECTOR_MIN_SIZE  5000
#endif

namespace viennacl
{
  namespace linalg
  {
    namespace host_based
    {
      namespace detail
      {
        template <typename NumericT>
        NumericT flip_sign(NumericT val) { return -val; }
        inline unsigned long  flip_sign(unsigned long  val) { return val; }
        inline unsigned int   flip_sign(unsigned int   val) { return val; }
        inline unsigned short flip_sign(unsigned short val) { return val; }
        inline unsigned char  flip_sign(unsigned char  val) { return val; }
      }

      //
      // Introductory note: By convention, all dimensions are already checked in the dispatcher frontend. No need to double-check again in here!
      //

      template <typename T, typename ScalarType1>
      void av(vector_base<T> & vec1,
              vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t /*len_alpha*/, bool reciprocal_alpha, bool flip_sign_alpha)
      {
        typedef T        value_type;

        value_type       * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = detail::flip_sign(data_alpha);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

        if (reciprocal_alpha)
        {
#ifdef VIENNACL_WITH_OPENMP
          #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
          for (long i = 0; i < static_cast<long>(size1); ++i)
            data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] / data_alpha;
        }
        else
        {
#ifdef VIENNACL_WITH_OPENMP
          #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
          for (long i = 0; i < static_cast<long>(size1); ++i)
            data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] * data_alpha;
        }
      }


      template <typename T, typename ScalarType1, typename ScalarType2>
      void avbv(vector_base<T> & vec1,
                vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t /* len_alpha */, bool reciprocal_alpha, bool flip_sign_alpha,
                vector_base<T> const & vec3, ScalarType2 const & beta,  vcl_size_t /* len_beta */,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef T        value_type;

        value_type       * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);
        value_type const * data_vec3 = detail::extract_raw_pointer<value_type>(vec3);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = detail::flip_sign(data_alpha);

        value_type data_beta = beta;
        if (flip_sign_beta)
          data_beta = detail::flip_sign(data_beta);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

        vcl_size_t start3 = viennacl::traits::start(vec3);
        vcl_size_t inc3   = viennacl::traits::stride(vec3);

        if (reciprocal_alpha)
        {
          if (reciprocal_beta)
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] / data_alpha + data_vec3[i*inc3+start3] / data_beta;
          }
          else
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] / data_alpha + data_vec3[i*inc3+start3] * data_beta;
          }
        }
        else
        {
          if (reciprocal_beta)
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] * data_alpha + data_vec3[i*inc3+start3] / data_beta;
          }
          else
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] = data_vec2[i*inc2+start2] * data_alpha + data_vec3[i*inc3+start3] * data_beta;
          }
        }
      }


      template <typename T, typename ScalarType1, typename ScalarType2>
      void avbv_v(vector_base<T> & vec1,
                  vector_base<T> const & vec2, ScalarType1 const & alpha, vcl_size_t /*len_alpha*/, bool reciprocal_alpha, bool flip_sign_alpha,
                  vector_base<T> const & vec3, ScalarType2 const & beta,  vcl_size_t /*len_beta*/,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef T        value_type;

        value_type       * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);
        value_type const * data_vec3 = detail::extract_raw_pointer<value_type>(vec3);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = detail::flip_sign(data_alpha);

        value_type data_beta = beta;
        if (flip_sign_beta)
          data_beta = detail::flip_sign(data_beta);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

        vcl_size_t start3 = viennacl::traits::start(vec3);
        vcl_size_t inc3   = viennacl::traits::stride(vec3);

        if (reciprocal_alpha)
        {
          if (reciprocal_beta)
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] += data_vec2[i*inc2+start2] / data_alpha + data_vec3[i*inc3+start3] / data_beta;
          }
          else
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] += data_vec2[i*inc2+start2] / data_alpha + data_vec3[i*inc3+start3] * data_beta;
          }
        }
        else
        {
          if (reciprocal_beta)
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] += data_vec2[i*inc2+start2] * data_alpha + data_vec3[i*inc3+start3] / data_beta;
          }
          else
          {
#ifdef VIENNACL_WITH_OPENMP
            #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
            for (long i = 0; i < static_cast<long>(size1); ++i)
              data_vec1[i*inc1+start1] += data_vec2[i*inc2+start2] * data_alpha + data_vec3[i*inc3+start3] * data_beta;
          }
        }
      }




      /** @brief Assign a constant value to a vector (-range/-slice)
      *
      * @param vec1   The vector to which the value should be assigned
      * @param alpha  The value to be assigned
      * @param up_to_internal_size  Specifies whether alpha should also be written to padded memory (mostly used for clearing the whole buffer).
      */
      template <typename T>
      void vector_assign(vector_base<T> & vec1, const T & alpha, bool up_to_internal_size = false)
      {
        typedef T        value_type;

        value_type * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);
        vcl_size_t loop_bound  = up_to_internal_size ? vec1.internal_size() : size1;  //Note: Do NOT use traits::internal_size() here, because vector proxies don't require padding.

        value_type data_alpha = static_cast<value_type>(alpha);

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for if (loop_bound > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(loop_bound); ++i)
          data_vec1[i*inc1+start1] = data_alpha;
      }


      /** @brief Swaps the contents of two vectors, data is copied
      *
      * @param vec1   The first vector (or -range, or -slice)
      * @param vec2   The second vector (or -range, or -slice)
      */
      template <typename T>
      void vector_swap(vector_base<T> & vec1, vector_base<T> & vec2)
      {
        typedef T        value_type;

        value_type * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
        {
          value_type temp = data_vec2[i*inc2+start2];
          data_vec2[i*inc2+start2] = data_vec1[i*inc1+start1];
          data_vec1[i*inc1+start1] = temp;
        }
      }


      ///////////////////////// Elementwise operations /////////////

      /** @brief Implementation of the element-wise operation v1 = v2 .* v3 and v1 = v2 ./ v3    (using MATLAB syntax)
      *
      * @param vec1   The result vector (or -range, or -slice)
      * @param proxy  The proxy object holding v2, v3 and the operation
      */
      template <typename T, typename OP>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> > const & proxy)
      {
        typedef T                                              value_type;
        typedef viennacl::linalg::detail::op_applier<op_element_binary<OP> >    OpFunctor;

        value_type       * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(proxy.lhs());
        value_type const * data_vec3 = detail::extract_raw_pointer<value_type>(proxy.rhs());

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(proxy.lhs());
        vcl_size_t inc2   = viennacl::traits::stride(proxy.lhs());

        vcl_size_t start3 = viennacl::traits::start(proxy.rhs());
        vcl_size_t inc3   = viennacl::traits::stride(proxy.rhs());

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
          OpFunctor::apply(data_vec1[i*inc1+start1], data_vec2[i*inc2+start2], data_vec3[i*inc3+start3]);
      }

      /** @brief Implementation of the element-wise operation v1 = v2 .* v3 and v1 = v2 ./ v3    (using MATLAB syntax)
      *
      * @param vec1   The result vector (or -range, or -slice)
      * @param proxy  The proxy object holding v2, v3 and the operation
      */
      template <typename T, typename OP>
      void element_op(vector_base<T> & vec1,
                      vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> > const & proxy)
      {
        typedef T        value_type;
        typedef viennacl::linalg::detail::op_applier<op_element_unary<OP> >    OpFunctor;

        value_type       * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(proxy.lhs());

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(proxy.lhs());
        vcl_size_t inc2   = viennacl::traits::stride(proxy.lhs());

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
          OpFunctor::apply(data_vec1[i*inc1+start1], data_vec2[i*inc2+start2]);
      }


      ///////////////////////// Norms and inner product ///////////////////


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

        value_type const * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type const * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

        value_type temp = 0;

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for reduction(+: temp) if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
          temp += data_vec1[i*inc1+start1] * data_vec2[i*inc2+start2];

        result = temp;  //Note: Assignment to result might be expensive, thus 'temp' is used for accumulation
      }

      template <typename T>
      void inner_prod_impl(vector_base<T> const & x,
                           vector_tuple<T> const & vec_tuple,
                           vector_base<T> & result)
      {
        typedef T        value_type;

        value_type const * data_x = detail::extract_raw_pointer<value_type>(x);

        vcl_size_t start_x = viennacl::traits::start(x);
        vcl_size_t inc_x   = viennacl::traits::stride(x);
        vcl_size_t size_x  = viennacl::traits::size(x);

        std::vector<value_type> temp(vec_tuple.const_size());
        std::vector<value_type const *> data_y(vec_tuple.const_size());
        std::vector<vcl_size_t> start_y(vec_tuple.const_size());
        std::vector<vcl_size_t> stride_y(vec_tuple.const_size());

        for (vcl_size_t j=0; j<vec_tuple.const_size(); ++j)
        {
          data_y[j] = detail::extract_raw_pointer<value_type>(vec_tuple.const_at(j));
          start_y[j] = viennacl::traits::start(vec_tuple.const_at(j));
          stride_y[j] = viennacl::traits::stride(vec_tuple.const_at(j));
        }

        // Note: No OpenMP here because it cannot perform a reduction on temp-array. Savings in memory bandwidth are expected to still justify this approach...
        for (vcl_size_t i = 0; i < size_x; ++i)
        {
          value_type entry_x = data_x[i*inc_x+start_x];
          for (vcl_size_t j=0; j < vec_tuple.const_size(); ++j)
            temp[j] += entry_x * data_y[j][i*stride_y[j]+start_y[j]];
        }

        for (vcl_size_t j=0; j < vec_tuple.const_size(); ++j)
          result[j] = temp[j];  //Note: Assignment to result might be expensive, thus 'temp' is used for accumulation
      }


      /** @brief Computes the l^1-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T, typename S2>
      void norm_1_impl(vector_base<T> const & vec1,
                       S2 & result)
      {
        typedef T        value_type;

        value_type const * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        value_type temp = 0;

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for reduction(+: temp) if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
          temp += static_cast<value_type>(std::fabs(static_cast<double>(data_vec1[i*inc1+start1])));  //casting to double in order to avoid problems if T is an integer type

        result = temp;  //Note: Assignment to result might be expensive, thus 'temp' is used for accumulation
      }

      /** @brief Computes the l^2-norm of a vector - implementation
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T, typename S2>
      void norm_2_impl(vector_base<T> const & vec1,
                       S2 & result)
      {
        typedef T        value_type;

        value_type const * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        value_type temp = 0;
        value_type data = 0;

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for reduction(+: temp) private(data) if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
        {
          data = data_vec1[i*inc1+start1];
          temp += data * data;
        }

        result = std::sqrt(temp);  //Note: Assignment to result might be expensive, thus 'temp' is used for accumulation
      }

      /** @brief Computes the supremum-norm of a vector
      *
      * @param vec1 The vector
      * @param result The result scalar
      */
      template <typename T, typename S2>
      void norm_inf_impl(vector_base<T> const & vec1,
                         S2 & result)
      {
        typedef T        value_type;

        value_type const * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        value_type temp = 0;

        // Note: No max() reduction in OpenMP yet
        for (vcl_size_t i = 0; i < size1; ++i)
          temp = std::max<value_type>(temp, static_cast<value_type>(std::fabs(static_cast<double>(data_vec1[i*inc1+start1]))));  //casting to double in order to avoid problems if T is an integer type

        result = temp;  //Note: Assignment to result might be expensive, thus 'temp' is used for accumulation
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
        typedef T        value_type;

        value_type const * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        value_type temp = 0;
        value_type data;
        vcl_size_t index = start1;

        // Note: No suitable reduction in OpenMP yet
        for (vcl_size_t i = 0; i < size1; ++i)
        {
          data = static_cast<value_type>(std::fabs(static_cast<double>(data_vec1[i*inc1+start1])));  //casting to double in order to avoid problems if T is an integer type
          if (data > temp)
          {
            index = i;
            temp = data;
          }
        }

        return index;
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
        typedef T   value_type;

        value_type * data_vec1 = detail::extract_raw_pointer<value_type>(vec1);
        value_type * data_vec2 = detail::extract_raw_pointer<value_type>(vec2);

        vcl_size_t start1 = viennacl::traits::start(vec1);
        vcl_size_t inc1   = viennacl::traits::stride(vec1);
        vcl_size_t size1  = viennacl::traits::size(vec1);

        vcl_size_t start2 = viennacl::traits::start(vec2);
        vcl_size_t inc2   = viennacl::traits::stride(vec2);

        value_type temp1 = 0;
        value_type temp2 = 0;
        value_type data_alpha = alpha;
        value_type data_beta  = beta;

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for private(temp1, temp2) if (size1 > VIENNACL_OPENMP_VECTOR_MIN_SIZE)
#endif
        for (long i = 0; i < static_cast<long>(size1); ++i)
        {
          temp1 = data_vec1[i*inc1+start1];
          temp2 = data_vec2[i*inc2+start2];

          data_vec1[i*inc1+start1] = data_alpha * temp1 + data_beta * temp2;
          data_vec2[i*inc2+start2] = data_alpha * temp2 - data_beta * temp1;
        }
      }

    } //namespace host_based
  } //namespace linalg
} //namespace viennacl


#endif
