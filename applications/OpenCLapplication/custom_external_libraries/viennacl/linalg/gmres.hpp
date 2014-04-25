#ifndef VIENNACL_GMRES_HPP_
#define VIENNACL_GMRES_HPP_

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

/** @file gmres.hpp
    @brief Implementations of the generalized minimum residual method are in this file.
*/

#include <vector>
#include <cmath>
#include <limits>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/traits/clear.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/meta/result_of.hpp"

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for the solver GMRES. Used for supplying solver parameters and for dispatching the solve() function
    */
    class gmres_tag       //generalized minimum residual
    {
      public:
        /** @brief The constructor
        *
        * @param tol            Relative tolerance for the residual (solver quits if ||r|| < tol * ||r_initial||)
        * @param max_iterations The maximum number of iterations (including restarts
        * @param krylov_dim     The maximum dimension of the Krylov space before restart (number of restarts is found by max_iterations / krylov_dim)
        */
        gmres_tag(double tol = 1e-10, unsigned int max_iterations = 300, unsigned int krylov_dim = 20)
         : tol_(tol), iterations_(max_iterations), krylov_dim_(krylov_dim), iters_taken_(0) {}

        /** @brief Returns the relative tolerance */
        double tolerance() const { return tol_; }
        /** @brief Returns the maximum number of iterations */
        unsigned int max_iterations() const { return iterations_; }
        /** @brief Returns the maximum dimension of the Krylov space before restart */
        unsigned int krylov_dim() const { return krylov_dim_; }
        /** @brief Returns the maximum number of GMRES restarts */
        unsigned int max_restarts() const
        {
          unsigned int ret = iterations_ / krylov_dim_;
          if (ret > 0 && (ret * krylov_dim_ == iterations_) )
            return ret - 1;
          return ret;
        }

        /** @brief Return the number of solver iterations: */
        unsigned int iters() const { return iters_taken_; }
        /** @brief Set the number of solver iterations (should only be modified by the solver) */
        void iters(unsigned int i) const { iters_taken_ = i; }

        /** @brief Returns the estimated relative error at the end of the solver run */
        double error() const { return last_error_; }
        /** @brief Sets the estimated relative error at the end of the solver run */
        void error(double e) const { last_error_ = e; }

      private:
        double tol_;
        unsigned int iterations_;
        unsigned int krylov_dim_;

        //return values from solver
        mutable unsigned int iters_taken_;
        mutable double last_error_;
    };

    namespace detail
    {

      template <typename SRC_VECTOR, typename DEST_VECTOR>
      void gmres_copy_helper(SRC_VECTOR const & src, DEST_VECTOR & dest, vcl_size_t len, vcl_size_t start = 0)
      {
        for (vcl_size_t i=0; i<len; ++i)
          dest[start+i] = src[start+i];
      }

      template <typename ScalarType, typename DEST_VECTOR>
      void gmres_copy_helper(viennacl::vector<ScalarType> const & src, DEST_VECTOR & dest, vcl_size_t len, vcl_size_t start = 0)
      {
        typedef typename viennacl::vector<ScalarType>::difference_type   difference_type;
        viennacl::copy( src.begin() + static_cast<difference_type>(start),
                        src.begin() + static_cast<difference_type>(start + len),
                       dest.begin() + static_cast<difference_type>(start));
      }

      /** @brief Computes the householder vector 'hh_vec' which rotates 'input_vec' such that all entries below the j-th entry of 'v' become zero.
        *
        * @param input_vec       The input vector
        * @param hh_vec          The householder vector defining the relection (I - beta * hh_vec * hh_vec^T)
        * @param beta            The coefficient beta in (I - beta  * hh_vec * hh_vec^T)
        * @param mu              The norm of the input vector part relevant for the reflection: norm_2(input_vec[j:size])
        * @param j               Index of the last nonzero index in 'input_vec' after applying the reflection
      */
      template <typename VectorType, typename ScalarType>
      void gmres_setup_householder_vector(VectorType const & input_vec, VectorType & hh_vec, ScalarType & beta, ScalarType & mu, vcl_size_t j)
      {
        ScalarType input_j = input_vec(j);

        // copy entries from input vector to householder vector:
        detail::gmres_copy_helper(input_vec, hh_vec, viennacl::traits::size(hh_vec) - (j+1), j+1);

        ScalarType sigma = viennacl::linalg::norm_2(hh_vec);
        sigma *= sigma;

        if (sigma == 0)
        {
          beta = 0;
          mu = input_j;
        }
        else
        {
          mu = std::sqrt(sigma + input_j*input_j);

          ScalarType hh_vec_0 = (input_j <= 0) ? (input_j - mu) : (-sigma / (input_j + mu));

          beta = ScalarType(2) * hh_vec_0 * hh_vec_0 / (sigma + hh_vec_0 * hh_vec_0);

          //divide hh_vec by its diagonal element hh_vec_0
          hh_vec /= hh_vec_0;
          hh_vec[j] = ScalarType(1);
        }
      }

      // Apply (I - beta h h^T) to x (Householder reflection with Householder vector h)
      template <typename VectorType, typename ScalarType>
      void gmres_householder_reflect(VectorType & x, VectorType const & h, ScalarType beta)
      {
        ScalarType hT_in_x = viennacl::linalg::inner_prod(h, x);
        x -= (beta * hT_in_x) * h;
      }

    }

    /** @brief Implementation of the GMRES solver.
    *
    * Following the algorithm proposed by Walker in "A Simpler GMRES"
    *
    * @param matrix     The system matrix
    * @param rhs        The load vector
    * @param tag        Solver configuration tag
    * @param precond    A preconditioner. Precondition operation is done via member function apply()
    * @return The result vector
    */
    template <typename MatrixType, typename VectorType, typename PreconditionerType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, gmres_tag const & tag, PreconditionerType const & precond)
    {
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      unsigned int problem_size = static_cast<unsigned int>(viennacl::traits::size(rhs));
      VectorType result = rhs;
      viennacl::traits::clear(result);

      unsigned int krylov_dim = tag.krylov_dim();
      if (problem_size < tag.krylov_dim())
        krylov_dim = problem_size; //A Krylov space larger than the matrix would lead to seg-faults (mathematically, error is certain to be zero already)

      VectorType res = rhs;
      VectorType v_k_tilde = rhs;
      VectorType v_k_tilde_temp = rhs;

      std::vector< std::vector<CPU_ScalarType> > R(krylov_dim, std::vector<CPU_ScalarType>(tag.krylov_dim()));
      std::vector<CPU_ScalarType> projection_rhs(krylov_dim);

      std::vector<VectorType>      householder_reflectors(krylov_dim, rhs);
      std::vector<CPU_ScalarType>  betas(krylov_dim);

      CPU_ScalarType norm_rhs = viennacl::linalg::norm_2(rhs);

      if (norm_rhs == 0) //solution is zero if RHS norm is zero
        return result;

      tag.iters(0);

      for (unsigned int it = 0; it <= tag.max_restarts(); ++it)
      {
        //
        // (Re-)Initialize residual: r = b - A*x (without temporary for the result of A*x)
        //
        res = rhs;
        res -= viennacl::linalg::prod(matrix, result);  //initial guess zero
        precond.apply(res);

        CPU_ScalarType rho_0 = viennacl::linalg::norm_2(res);

        //
        // Check for premature convergence
        //
        if (rho_0 / norm_rhs < tag.tolerance() ) // norm_rhs is known to be nonzero here
        {
          tag.error(rho_0 / norm_rhs);
          return result;
        }

        //
        // Normalize residual and set 'rho' to 1 as requested in 'A Simpler GMRES' by Walker and Zhou.
        //
        res /= rho_0;
        CPU_ScalarType rho = static_cast<CPU_ScalarType>(1.0);


        //
        // Iterate up until maximal Krylove space dimension is reached:
        //
        unsigned int k = 0;
        for (k = 0; k < krylov_dim; ++k)
        {
          tag.iters( tag.iters() + 1 ); //increase iteration counter

          // prepare storage:
          viennacl::traits::clear(R[k]);
          viennacl::traits::clear(householder_reflectors[k]);

          //compute v_k = A * v_{k-1} via Householder matrices
          if (k == 0)
          {
            v_k_tilde = viennacl::linalg::prod(matrix, res);
            precond.apply(v_k_tilde);
          }
          else
          {
            viennacl::traits::clear(v_k_tilde);
            v_k_tilde[k-1] = CPU_ScalarType(1);

            //Householder rotations, part 1: Compute P_1 * P_2 * ... * P_{k-1} * e_{k-1}
            for (int i = k-1; i > -1; --i)
              detail::gmres_householder_reflect(v_k_tilde, householder_reflectors[i], betas[i]);

            v_k_tilde_temp = viennacl::linalg::prod(matrix, v_k_tilde);
            precond.apply(v_k_tilde_temp);
            v_k_tilde = v_k_tilde_temp;

            //Householder rotations, part 2: Compute P_{k-1} * ... * P_{1} * v_k_tilde
            for (unsigned int i = 0; i < k; ++i)
              detail::gmres_householder_reflect(v_k_tilde, householder_reflectors[i], betas[i]);
          }

          //
          // Compute Householder reflection for v_k_tilde such that all entries below k-th entry are zero:
          //
          CPU_ScalarType rho_k_k = 0;
          detail::gmres_setup_householder_vector(v_k_tilde, householder_reflectors[k], betas[k], rho_k_k, k);

          //
          // copy first k entries from v_k_tilde to R[k] in order to fill k-th column with result of
          // P_k * v_k_tilde = (v[0], ... , v[k-1], norm(v), 0, 0, ...) =: (rho_{1,k}, rho_{2,k}, ..., rho_{k,k}, 0, ..., 0);
          //
          detail::gmres_copy_helper(v_k_tilde, R[k], k);
          R[k][k] = rho_k_k;

          //
          // Update residual: r = P_k r
          // Set zeta_k = r[k] including machine precision considerations: mathematically we have |r[k]| <= rho
          // Set rho *= sin(acos(r[k] / rho))
          //
          detail::gmres_householder_reflect(res, householder_reflectors[k], betas[k]);

          if (res[k] > rho) //machine precision reached
            res[k] = rho;
          if (res[k] < -rho) //machine precision reached
            res[k] = -rho;
          projection_rhs[k] = res[k];

          rho *= std::sin( std::acos(projection_rhs[k] / rho) );

          if (std::fabs(rho * rho_0 / norm_rhs) < tag.tolerance())  // Residual is sufficiently reduced, stop here
          {
            tag.error( std::fabs(rho*rho_0 / norm_rhs) );
            ++k;
            break;
          }
        } // for k

        //
        // Triangular solver stage:
        //

        for (int i=k-1; i>-1; --i)
        {
          for (unsigned int j=i+1; j<k; ++j)
            projection_rhs[i] -= R[j][i] * projection_rhs[j];     //R is transposed

          projection_rhs[i] /= R[i][i];
        }

        //
        // Note: 'projection_rhs' now holds the solution (eta_1, ..., eta_k)
        //

        res *= projection_rhs[0];

        if (k > 0)
        {
          for (unsigned int i = 0; i < k-1; ++i)
            res[i] += projection_rhs[i+1];
        }

        //
        // Form z inplace in 'res' by applying P_1 * ... * P_{k}
        //
        for (int i=k-1; i>=0; --i)
          detail::gmres_householder_reflect(res, householder_reflectors[i], betas[i]);

        res *= rho_0;
        result += res;  // x += rho_0 * z    in the paper

        //
        // Check for convergence:
        //
        tag.error(std::fabs(rho*rho_0 / norm_rhs));
        if ( tag.error() < tag.tolerance() )
          return result;
      }

      return result;
    }

    /** @brief Convenience overload of the solve() function using GMRES. Per default, no preconditioner is used
    */
    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, gmres_tag const & tag)
    {
      return solve(matrix, rhs, tag, no_precond());
    }


  }
}

#endif
