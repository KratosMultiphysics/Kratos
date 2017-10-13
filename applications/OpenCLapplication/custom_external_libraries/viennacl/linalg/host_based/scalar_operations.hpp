#ifndef VIENNACL_LINALG_HOST_BASED_SCALAR_OPERATIONS_HPP_
#define VIENNACL_LINALG_HOST_BASED_SCALAR_OPERATIONS_HPP_

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

/** @file viennacl/linalg/host_based/scalar_operations.hpp
    @brief Implementations of scalar operations using a plain single-threaded or OpenMP-enabled execution on CPU
*/

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/stride.hpp"
#include "viennacl/linalg/host_based/common.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace host_based
    {
      template <typename S1,
                typename S2, typename ScalarType1>
      typename viennacl::enable_if< viennacl::is_scalar<S1>::value
                                    && viennacl::is_scalar<S2>::value
                                    && viennacl::is_any_scalar<ScalarType1>::value
                                  >::type
      as(S1 & s1,
         S2 const & s2, ScalarType1 const & alpha, vcl_size_t /*len_alpha*/, bool reciprocal_alpha, bool flip_sign_alpha)
      {
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        value_type       * data_s1 = detail::extract_raw_pointer<value_type>(s1);
        value_type const * data_s2 = detail::extract_raw_pointer<value_type>(s2);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        *data_s1 = *data_s2 * data_alpha;
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
           S2 const & s2, ScalarType1 const & alpha, vcl_size_t /*len_alpha*/, bool reciprocal_alpha, bool flip_sign_alpha,
           S3 const & s3, ScalarType2 const & beta,  vcl_size_t /*len_beta*/,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        value_type       * data_s1 = detail::extract_raw_pointer<value_type>(s1);
        value_type const * data_s2 = detail::extract_raw_pointer<value_type>(s2);
        value_type const * data_s3 = detail::extract_raw_pointer<value_type>(s3);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        value_type data_beta = beta;
        if (flip_sign_beta)
          data_beta = -data_beta;
        if (reciprocal_beta)
          data_beta = static_cast<value_type>(1) / data_beta;

        *data_s1 = *data_s2 * data_alpha + *data_s3 * data_beta;
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
             S2 const & s2, ScalarType1 const & alpha, vcl_size_t /*len_alpha*/, bool reciprocal_alpha, bool flip_sign_alpha,
             S3 const & s3, ScalarType2 const & beta,  vcl_size_t /*len_beta*/,  bool reciprocal_beta,  bool flip_sign_beta)
      {
        typedef typename viennacl::result_of::cpu_value_type<S1>::type        value_type;

        value_type       * data_s1 = detail::extract_raw_pointer<value_type>(s1);
        value_type const * data_s2 = detail::extract_raw_pointer<value_type>(s2);
        value_type const * data_s3 = detail::extract_raw_pointer<value_type>(s3);

        value_type data_alpha = alpha;
        if (flip_sign_alpha)
          data_alpha = -data_alpha;
        if (reciprocal_alpha)
          data_alpha = static_cast<value_type>(1) / data_alpha;

        value_type data_beta = beta;
        if (flip_sign_beta)
          data_beta = -data_beta;
        if (reciprocal_beta)
          data_beta = static_cast<value_type>(1) / data_beta;

        *data_s1 += *data_s2 * data_alpha + *data_s3 * data_beta;
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

        value_type * data_s1 = detail::extract_raw_pointer<value_type>(s1);
        value_type * data_s2 = detail::extract_raw_pointer<value_type>(s2);

        value_type temp = *data_s2;
        *data_s2 = *data_s1;
        *data_s1 = temp;
      }



    } //namespace host_based
  } //namespace linalg
} //namespace viennacl


#endif
