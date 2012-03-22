#ifndef VIENNACL_CG_HPP_
#define VIENNACL_CG_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file cg.hpp
    @brief The conjugate gradient method is implemented here
*/

#include <vector>
#include <map>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/traits/clear.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/meta/result_of.hpp"

namespace viennacl
{
  namespace linalg
  {
    
    /** @brief A tag for the conjugate gradient Used for supplying solver parameters and for dispatching the solve() function
    */
    class cg_tag
    {
      public:
        /** @brief The constructor
        *
        * @param tol              Relative tolerance for the residual (solver quits if ||r|| < tol * ||r_initial||)
        * @param max_iterations   The maximum number of iterations
        */
        cg_tag(double tol = 1e-8, unsigned int max_iterations = 300) : _tol(tol), _iterations(max_iterations) {};
      
        /** @brief Returns the relative tolerance */
        double tolerance() const { return _tol; }
        /** @brief Returns the maximum number of iterations */
        unsigned int max_iterations() const { return _iterations; }
        
        /** @brief Return the number of solver iterations: */
        unsigned int iters() const { return iters_taken_; }
        void iters(unsigned int i) const { iters_taken_ = i; }
        
        /** @brief Returns the estimated relative error at the end of the solver run */
        double error() const { return last_error_; }
        /** @brief Sets the estimated relative error at the end of the solver run */
        void error(double e) const { last_error_ = e; }
        
        
      private:
        double _tol;
        unsigned int _iterations;
        
        //return values from solver
        mutable unsigned int iters_taken_;
        mutable double last_error_;
    };
    

    /** @brief Implementation of the conjugate gradient solver without preconditioner
    *
    * Following the algorithm in the book by Y. Saad "Iterative Methods for sparse linear systems"
    *
    * @param matrix     The system matrix
    * @param rhs        The load vector
    * @param tag        Solver configuration tag
    * @return The result vector
    */
    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, cg_tag const & tag)
    {
      //typedef typename VectorType::value_type      ScalarType;
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      //std::cout << "Starting CG" << std::endl;
      std::size_t problem_size = viennacl::traits::size(rhs);
      VectorType result(problem_size);
      viennacl::traits::clear(result);

      VectorType residual = rhs;
      VectorType p = rhs;
      VectorType tmp(problem_size);

      ScalarType tmp_in_p;
      ScalarType residual_norm_squared;
      CPU_ScalarType ip_rr = viennacl::linalg::inner_prod(rhs,rhs);
      CPU_ScalarType alpha;
      CPU_ScalarType new_ip_rr = 0;
      CPU_ScalarType beta;
      CPU_ScalarType norm_rhs_squared = ip_rr;
      
      //std::cout << "Starting CG solver iterations... " << std::endl;
      
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);
        tmp = viennacl::linalg::prod(matrix, p);

        //alpha = ip_rr / viennacl::linalg::inner_prod(tmp, p);
        tmp_in_p = viennacl::linalg::inner_prod(tmp, p);
        alpha = ip_rr / static_cast<CPU_ScalarType>(tmp_in_p);
        
        result += alpha * p;
        residual -= alpha * tmp;
        
        residual_norm_squared = viennacl::linalg::inner_prod(residual,residual);
        new_ip_rr = static_cast<CPU_ScalarType>(residual_norm_squared);
        if (new_ip_rr / norm_rhs_squared < tag.tolerance() *  tag.tolerance())//squared norms involved here
          break;
        
        beta = new_ip_rr / ip_rr;
        ip_rr = new_ip_rr;

        //p = residual + beta*p;
        p *= beta;
        p += residual;
      } 
      
      //store last error estimate:
      tag.error(sqrt(new_ip_rr / norm_rhs_squared));
      
      return result;
    }

    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, cg_tag const & tag, viennacl::linalg::no_precond)
    {
      return solve(matrix, rhs, tag);
    }

    /** @brief Implementation of the preconditioned conjugate gradient solver
    *
    * Following Algorithm 9.1 in "Iterative Methods for Sparse Linear Systems" by Y. Saad
    *
    * @param matrix     The system matrix
    * @param rhs        The load vector
    * @param tag        Solver configuration tag
    * @param precond    A preconditioner. Precondition operation is done via member function apply()
    * @return The result vector
    */
    template <typename MatrixType, typename VectorType, typename PreconditionerType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, cg_tag const & tag, PreconditionerType const & precond)
    {
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      unsigned int problem_size = viennacl::traits::size(rhs);
      
      VectorType result(problem_size);
      result.clear();
      
      VectorType residual = rhs;
      VectorType tmp(problem_size);
      VectorType z = rhs;

      precond.apply(z);
      VectorType p = z;

      ScalarType tmp_in_p;
      ScalarType residual_in_z;
      CPU_ScalarType ip_rr = viennacl::linalg::inner_prod(residual, z);
      CPU_ScalarType alpha;
      CPU_ScalarType new_ip_rr = 0;
      CPU_ScalarType beta;
      CPU_ScalarType norm_rhs_squared = ip_rr;
      CPU_ScalarType new_ipp_rr_over_norm_rhs;
      
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);
        tmp = viennacl::linalg::prod(matrix, p);
        
        tmp_in_p = viennacl::linalg::inner_prod(tmp, p);
        alpha = ip_rr / static_cast<CPU_ScalarType>(tmp_in_p);
        
        result += alpha * p;
        residual -= alpha * tmp;
        z = residual;
        precond.apply(z);
        
        residual_in_z = viennacl::linalg::inner_prod(residual, z);
        new_ip_rr = static_cast<CPU_ScalarType>(residual_in_z);
        new_ipp_rr_over_norm_rhs = new_ip_rr / norm_rhs_squared;
        if (std::fabs(new_ipp_rr_over_norm_rhs) < tag.tolerance() *  tag.tolerance())    //squared norms involved here
          break;
        
        beta = new_ip_rr / ip_rr;
        ip_rr = new_ip_rr;
        
        //p = z + beta*p;
        p *= beta;
        p += z;
      } 

      //store last error estimate:
      tag.error(sqrt(std::fabs(new_ip_rr / norm_rhs_squared)));

      return result;
    }

  }
}

#endif
