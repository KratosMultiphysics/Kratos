#ifndef VIENNACL_LINALG_BISECT_HPP_
#define VIENNACL_LINALG_BISECT_HPP_

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

/** @file viennacl/linalg/bisect.hpp
*   @brief Implementation of the algorithm for finding eigenvalues of a tridiagonal matrix.
*
*   Contributed by Guenther Mader and Astrid Rupp.
*/

#include <vector>
#include <cmath>
#include <limits>
#include <cstddef>
#include "viennacl/meta/result_of.hpp"

namespace viennacl
{
  namespace linalg
  {

    namespace detail
    {
      /**
      *    @brief overloaded function for copying vectors
      */
      template <typename T, typename OtherVectorType>
      void copy_vec_to_vec(viennacl::vector<T> const & src, OtherVectorType & dest)
      {
        viennacl::copy(src, dest);
      }

      template <typename OtherVectorType, typename T>
      void copy_vec_to_vec(OtherVectorType const & src, viennacl::vector<T> & dest)
      {
        viennacl::copy(src, dest);
      }

      template <typename VectorType1, typename VectorType2>
      void copy_vec_to_vec(VectorType1 const & src, VectorType2 & dest)
      {
        for (vcl_size_t i=0; i<src.size(); ++i)
          dest[i] = src[i];
      }
    }

    /**
    *   @brief Implementation of the bisect-algorithm for the calculation of the eigenvalues of a tridiagonal matrix. Experimental - interface might change.
    *
    *   @param alphas       Elements of the main diagonal
    *   @param betas        Elements of the secondary diagonal
    *   @return             Returns the eigenvalues of the tridiagonal matrix defined by alpha and beta
    */
    template< typename VectorT >
    std::vector<
            typename viennacl::result_of::cpu_value_type<typename VectorT::value_type>::type
            >
    bisect(VectorT const & alphas, VectorT const & betas)
    {
      typedef typename viennacl::result_of::value_type<VectorT>::type           ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

      vcl_size_t size = betas.size();
      std::vector<CPU_ScalarType>  x_temp(size);


      std::vector<CPU_ScalarType> beta_bisect;
      std::vector<CPU_ScalarType> wu;

      double rel_error = std::numeric_limits<CPU_ScalarType>::epsilon();
      beta_bisect.push_back(0);

      for(vcl_size_t i = 1; i < size; i++){
              beta_bisect.push_back(betas[i] * betas[i]);
      }

      double xmin = alphas[size - 1] - std::fabs(betas[size - 1]);
      double xmax = alphas[size - 1] + std::fabs(betas[size - 1]);

      for(vcl_size_t i = 0; i < size - 1; i++)
      {
        double h = std::fabs(betas[i]) + std::fabs(betas[i + 1]);
        if (alphas[i] + h > xmax)
          xmax = alphas[i] + h;
        if (alphas[i] - h < xmin)
          xmin = alphas[i] - h;
      }


      double eps1 = 1e-6;
      /*double eps2 = (xmin + xmax > 0) ? (rel_error * xmax) : (-rel_error * xmin);
      if(eps1 <= 0)
        eps1 = eps2;
      else
        eps2 = 0.5 * eps1 + 7.0 * eps2; */

      double x0 = xmax;

      for(vcl_size_t i = 0; i < size; i++)
      {
        x_temp[i] = xmax;
        wu.push_back(xmin);
      }

      for(long k = static_cast<long>(size) - 1; k >= 0; --k)
      {
        double xu = xmin;
        for(long i = k; i >= 0; --i)
        {
          if(xu < wu[k-i])
          {
            xu = wu[i];
            break;
          }
        }

        if(x0 > x_temp[k])
          x0 = x_temp[k];

        double x1 = (xu + x0) / 2.0;
        while (x0 - xu > 2.0 * rel_error * (std::fabs(xu) + std::fabs(x0)) + eps1)
        {
          vcl_size_t a = 0;
          double q = 1;
          for(vcl_size_t i = 0; i < size; i++)
          {
            if(q != 0)
              q = alphas[i] - x1 - beta_bisect[i] / q;
            else
              q = alphas[i] - x1 - std::fabs(betas[i] / rel_error);

            if(q < 0)
              a++;
          }

          if (a <= static_cast<vcl_size_t>(k))
          {
            xu = x1;
            if(a < 1)
              wu[0] = x1;
            else
            {
              wu[a] = x1;
              if(x_temp[a - 1] > x1)
                  x_temp[a - 1] = x1;
            }
          }
          else
            x0 = x1;

          x1 = (xu + x0) / 2.0;
        }
        x_temp[k] = x1;
      }
      return x_temp;
    }

  } // end namespace linalg
} // end namespace viennacl
#endif
