#ifndef VIENNACL_LINALG_POWER_ITER_HPP_
#define VIENNACL_LINALG_POWER_ITER_HPP_

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

/** @file viennacl/linalg/power_iter.hpp
    @brief Defines a tag for the configuration of the power iteration method.

    Contributed by Astrid Rupp.
*/

#include <cmath>
#include <vector>
#include "viennacl/linalg/bisect.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"

namespace viennacl
{
  namespace linalg
  {
    /** @brief A tag for the power iteration algorithm. */
    class power_iter_tag
    {
      public:

        /** @brief The constructor
        *
        * @param tfac      If the eigenvalue does not change more than this termination factor, the algorithm stops
        * @param max_iters Maximum number of iterations for the power iteration
        */
        power_iter_tag(double tfac = 1e-8, vcl_size_t max_iters = 50000) : termination_factor_(tfac), max_iterations_(max_iters) {}

        /** @brief Sets the factor for termination */
        void factor(double fct){ termination_factor_ = fct; }

          /** @brief Returns the factor for termination */
        double factor() const { return termination_factor_; }

        vcl_size_t max_iterations() const { return max_iterations_; }
        void max_iterations(vcl_size_t new_max) { max_iterations_ = new_max; }

      private:
        double termination_factor_;
        vcl_size_t max_iterations_;

    };

   /**
    *   @brief Implementation of the calculation of eigenvalues using poweriteration
    *
    *   @param matrix        The system matrix
    *   @param tag           Tag with termination factor
    *   @return              Returns the largest eigenvalue computed by the power iteration method
    */
    template< typename MatrixT >
    typename viennacl::result_of::cpu_value_type<typename MatrixT::value_type>::type
    eig(MatrixT const& matrix, power_iter_tag const & tag)
    {

      typedef typename viennacl::result_of::value_type<MatrixT>::type           ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      typedef typename viennacl::result_of::vector_for_matrix<MatrixT>::type    VectorT;

      CPU_ScalarType eigenvalue;
      vcl_size_t matrix_size = matrix.size1();
      VectorT r(matrix_size);
      VectorT r2(matrix_size);
      std::vector<CPU_ScalarType> s(matrix_size);

      for(vcl_size_t i=0; i<s.size(); ++i)
        s[i] = (i % 3) * CPU_ScalarType(0.1234) - CPU_ScalarType(0.5);   //'random' starting vector

      detail::copy_vec_to_vec(s,r);

      //std::cout << s << std::endl;

      double epsilon = tag.factor();
      CPU_ScalarType norm = norm_2(r);
      CPU_ScalarType norm_prev = 0;
      long numiter = 0;

      for (vcl_size_t i=0; i<tag.max_iterations(); ++i)
      {
        if (std::fabs(norm - norm_prev) / std::fabs(norm) < epsilon)
          break;

        r /= norm;
        r2 = viennacl::linalg::prod(matrix, r);  //using helper vector r2 for the computation of r <- A * r in order to avoid the repeated creation of temporaries
        r = r2;
        norm_prev = norm;
        norm = norm_2(r);
        numiter++;
      }

      eigenvalue = norm;
      return eigenvalue;
    }


  } // end namespace linalg
} // end namespace viennacl
#endif
