/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory
======================================================================= */

#ifndef _VIENNACL_CG_HPP_
#define _VIENNACL_CG_HPP_

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
      typedef typename viennacl::tools::result_of::value_type<VectorType>::type    ScalarType;
      unsigned int problem_size = viennacl::tools::traits::size(rhs);
      VectorType result(problem_size);
      viennacl::tools::traits::clear(result);

      VectorType residual = rhs;
      VectorType p = rhs;
      VectorType tmp(problem_size);

      ScalarType ip_rr = viennacl::linalg::inner_prod(rhs,rhs);
      ScalarType alpha;
      ScalarType new_ip_rr = 0;
      ScalarType beta;
      ScalarType norm_rhs = ip_rr;
      
      //std::cout << "Starting CG solver... " << std::endl;
      
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);
        tmp = viennacl::linalg::prod(matrix, p);

        alpha = ip_rr / viennacl::linalg::inner_prod(tmp, p);
        result += alpha * p;
        residual -= alpha * tmp;
        
        new_ip_rr = viennacl::linalg::inner_prod(residual,residual);
        if (new_ip_rr / norm_rhs < tag.tolerance() *  tag.tolerance())//squared norms involved here
          break;
        
        beta = new_ip_rr / ip_rr;
        ip_rr = new_ip_rr;

        p = residual + beta*p;
      } 
      
      //store last error estimate:
      tag.error(sqrt(new_ip_rr / norm_rhs));
      
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
      typedef typename VectorType::value_type      ScalarType;
      typedef typename viennacl::tools::CPU_SCALAR_TYPE_DEDUCER<ScalarType>::ResultType    CPU_ScalarType;
      unsigned int problem_size = viennacl::tools::traits::size(rhs);
      
      VectorType result(problem_size);
      result.clear();
      
      VectorType residual = rhs;
      VectorType tmp(problem_size);
      VectorType z = rhs;

      precond.apply(z);
      VectorType p = z;

      ScalarType ip_rr = viennacl::linalg::inner_prod(residual, z);
      ScalarType alpha;
      ScalarType new_ip_rr = 0;
      ScalarType beta;
      ScalarType norm_rhs = ip_rr;
      ScalarType new_ipp_rr_over_norm_rhs;
      
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);
        tmp = viennacl::linalg::prod(matrix, p);
        alpha = ip_rr / viennacl::linalg::inner_prod(tmp, p);
        result += alpha * p;
        residual -= alpha * tmp;
        z = residual;
        precond.apply(z);
        
        new_ip_rr = viennacl::linalg::inner_prod(residual, z);
        new_ipp_rr_over_norm_rhs = new_ip_rr / norm_rhs;
        if (std::fabs(new_ipp_rr_over_norm_rhs) < tag.tolerance() *  tag.tolerance())    //squared norms involved here
          break;
        
        beta = new_ip_rr / ip_rr;
        ip_rr = new_ip_rr;
        p = z + beta*p;
      } 

      //store last error estimate:
      tag.error(sqrt(std::fabs(new_ip_rr / norm_rhs)));

      return result;
    }

  }
}

#endif
