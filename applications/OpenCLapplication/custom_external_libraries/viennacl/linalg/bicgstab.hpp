#ifndef VIENNACL_LINALG_BICGSTAB_HPP_
#define VIENNACL_LINALG_BICGSTAB_HPP_

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

/** @file bicgstab.hpp
    @brief The stabilized bi-conjugate gradient method is implemented here
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/traits/clear.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/meta/result_of.hpp"

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for the stabilized Bi-conjugate gradient solver. Used for supplying solver parameters and for dispatching the solve() function
    */
    class bicgstab_tag
    {
      public:
        /** @brief The constructor
        *
        * @param tol              Relative tolerance for the residual (solver quits if ||r|| < tol * ||r_initial||)
        * @param max_iters        The maximum number of iterations
        * @param max_iters_before_restart   The maximum number of iterations before BiCGStab is reinitialized (to avoid accumulation of round-off errors)
        */
        bicgstab_tag(double tol = 1e-8, vcl_size_t max_iters = 400, vcl_size_t max_iters_before_restart = 200)
          : tol_(tol), iterations_(max_iters), iterations_before_restart_(max_iters_before_restart) {}

        /** @brief Returns the relative tolerance */
        double tolerance() const { return tol_; }
        /** @brief Returns the maximum number of iterations */
        vcl_size_t max_iterations() const { return iterations_; }
        /** @brief Returns the maximum number of iterations before a restart*/
        vcl_size_t max_iterations_before_restart() const { return iterations_before_restart_; }

        /** @brief Return the number of solver iterations: */
        vcl_size_t iters() const { return iters_taken_; }
        void iters(vcl_size_t i) const { iters_taken_ = i; }

        /** @brief Returns the estimated relative error at the end of the solver run */
        double error() const { return last_error_; }
        /** @brief Sets the estimated relative error at the end of the solver run */
        void error(double e) const { last_error_ = e; }

      private:
        double tol_;
        vcl_size_t iterations_;
        vcl_size_t iterations_before_restart_;

        //return values from solver
        mutable vcl_size_t iters_taken_;
        mutable double last_error_;
    };


    /** @brief Implementation of the stabilized Bi-conjugate gradient solver
    *
    * Following the description in "Iterative Methods for Sparse Linear Systems" by Y. Saad
    *
    * @param matrix     The system matrix
    * @param rhs        The load vector
    * @param tag        Solver configuration tag
    * @return The result vector
    */
    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, bicgstab_tag const & tag)
    {
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      VectorType result = rhs;
      viennacl::traits::clear(result);

      VectorType residual = rhs;
      VectorType p = rhs;
      VectorType r0star = rhs;
      VectorType tmp0 = rhs;
      VectorType tmp1 = rhs;
      VectorType s = rhs;

      CPU_ScalarType norm_rhs_host = viennacl::linalg::norm_2(residual);
      CPU_ScalarType ip_rr0star = norm_rhs_host * norm_rhs_host;
      CPU_ScalarType beta;
      CPU_ScalarType alpha;
      CPU_ScalarType omega;
      //ScalarType inner_prod_temp; //temporary variable for inner product computation
      CPU_ScalarType new_ip_rr0star = 0;
      CPU_ScalarType residual_norm = norm_rhs_host;

      if (norm_rhs_host == 0) //solution is zero if RHS norm is zero
        return result;

      bool restart_flag = true;
      vcl_size_t last_restart = 0;
      for (vcl_size_t i = 0; i < tag.max_iterations(); ++i)
      {
        if (restart_flag)
        {
          residual = rhs;
          residual -= viennacl::linalg::prod(matrix, result);
          p = residual;
          r0star = residual;
          ip_rr0star = viennacl::linalg::norm_2(residual);
          ip_rr0star *= ip_rr0star;
          restart_flag = false;
          last_restart = i;
        }

        tag.iters(i+1);
        tmp0 = viennacl::linalg::prod(matrix, p);
        alpha = ip_rr0star / viennacl::linalg::inner_prod(tmp0, r0star);

        s = residual - alpha*tmp0;

        tmp1 = viennacl::linalg::prod(matrix, s);
        CPU_ScalarType norm_tmp1 = viennacl::linalg::norm_2(tmp1);
        omega = viennacl::linalg::inner_prod(tmp1, s) / (norm_tmp1 * norm_tmp1);

        result += alpha * p + omega * s;
        residual = s - omega * tmp1;

        new_ip_rr0star = viennacl::linalg::inner_prod(residual, r0star);
        residual_norm = viennacl::linalg::norm_2(residual);
        if (std::fabs(residual_norm / norm_rhs_host) < tag.tolerance())
          break;

        beta = new_ip_rr0star / ip_rr0star * alpha/omega;
        ip_rr0star = new_ip_rr0star;

        if (ip_rr0star == 0 || omega == 0 || i - last_restart > tag.max_iterations_before_restart()) //search direction degenerate. A restart might help
          restart_flag = true;

        // Execution of
        //  p = residual + beta * (p - omega*tmp0);
        // without introducing temporary vectors:
        p -= omega * tmp0;
        p = residual + beta * p;
      }

      //store last error estimate:
      tag.error(residual_norm / norm_rhs_host);

      return result;
    }

    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, bicgstab_tag const & tag, viennacl::linalg::no_precond)
    {
      return solve(matrix, rhs, tag);
    }

    /** @brief Implementation of the preconditioned stabilized Bi-conjugate gradient solver
    *
    * Following the description of the unpreconditioned case in "Iterative Methods for Sparse Linear Systems" by Y. Saad
    *
    * @param matrix     The system matrix
    * @param rhs        The load vector
    * @param tag        Solver configuration tag
    * @param precond    A preconditioner. Precondition operation is done via member function apply()
    * @return The result vector
    */
    template <typename MatrixType, typename VectorType, typename PreconditionerType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, bicgstab_tag const & tag, PreconditionerType const & precond)
    {
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      VectorType result = rhs;
      viennacl::traits::clear(result);

      VectorType residual = rhs;
      VectorType r0star = residual;  //can be chosen arbitrarily in fact
      VectorType tmp0 = rhs;
      VectorType tmp1 = rhs;
      VectorType s = rhs;

      VectorType p = residual;

      CPU_ScalarType ip_rr0star = viennacl::linalg::norm_2(residual);
      CPU_ScalarType norm_rhs_host = viennacl::linalg::norm_2(residual);
      CPU_ScalarType beta;
      CPU_ScalarType alpha;
      CPU_ScalarType omega;
      CPU_ScalarType new_ip_rr0star = 0;
      CPU_ScalarType residual_norm = norm_rhs_host;

      if (norm_rhs_host == 0) //solution is zero if RHS norm is zero
        return result;

      bool restart_flag = true;
      vcl_size_t last_restart = 0;
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        if (restart_flag)
        {
          residual = rhs;
          residual -= viennacl::linalg::prod(matrix, result);
          precond.apply(residual);
          p = residual;
          r0star = residual;
          ip_rr0star = viennacl::linalg::norm_2(residual);
          ip_rr0star *= ip_rr0star;
          restart_flag = false;
          last_restart = i;
        }

        tag.iters(i+1);
        tmp0 = viennacl::linalg::prod(matrix, p);
        precond.apply(tmp0);
        alpha = ip_rr0star / viennacl::linalg::inner_prod(tmp0, r0star);

        s = residual - alpha*tmp0;

        tmp1 = viennacl::linalg::prod(matrix, s);
        precond.apply(tmp1);
        CPU_ScalarType norm_tmp1 = viennacl::linalg::norm_2(tmp1);
        omega = viennacl::linalg::inner_prod(tmp1, s) / (norm_tmp1 * norm_tmp1);

        result += alpha * p + omega * s;
        residual = s - omega * tmp1;

        residual_norm = viennacl::linalg::norm_2(residual);
        if (residual_norm / norm_rhs_host < tag.tolerance())
          break;

        new_ip_rr0star = viennacl::linalg::inner_prod(residual, r0star);

        beta = new_ip_rr0star / ip_rr0star * alpha/omega;
        ip_rr0star = new_ip_rr0star;

        if (ip_rr0star == 0 || omega == 0 || i - last_restart > tag.max_iterations_before_restart()) //search direction degenerate. A restart might help
          restart_flag = true;

        // Execution of
        //  p = residual + beta * (p - omega*tmp0);
        // without introducing temporary vectors:
        p -= omega * tmp0;
        p = residual + beta * p;

        //std::cout << "Rel. Residual in current step: " << std::sqrt(std::fabs(viennacl::linalg::inner_prod(residual, residual) / norm_rhs_host)) << std::endl;
      }

      //store last error estimate:
      tag.error(residual_norm / norm_rhs_host);

      return result;
    }

  }
}

#endif
