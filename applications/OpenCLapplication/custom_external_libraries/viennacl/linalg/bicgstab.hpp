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

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_BICGSTAB_HPP_
#define _VIENNACL_BICGSTAB_HPP_

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"

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
        * @param max_iterations   The maximum number of iterations
        */
        bicgstab_tag(double tol = 1e-8, unsigned int max_iterations = 300) : _tol(tol), _iterations(max_iterations) {};
      
        double tolerance() const { return _tol; }
        unsigned int iterations() const { return _iterations; }
        
        //number of solver iterations:
        unsigned int iters() const { return iters_taken_; }
        void iters(unsigned int i) const { iters_taken_ = i; }
        
        //esimated error :
        double error() const { return last_error_; }
        void error(double e) const { last_error_ = e; }
        
      private:
        double _tol;
        unsigned int _iterations;

        //return values from solver
        mutable unsigned int iters_taken_;
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
      typedef typename VectorType::value_type      ScalarType;
      VectorType result(rhs.size());
      result.clear();

      VectorType residual = rhs;
      VectorType p = rhs;
      VectorType r0star = rhs;
      VectorType tmp0(rhs.size());
      VectorType tmp1(rhs.size());
      VectorType s(rhs.size());

      ScalarType ip_rr0star = viennacl::linalg::inner_prod(rhs,r0star);
      ScalarType alpha;
      ScalarType omega;
      ScalarType new_ip_rr0star = 0;
      ScalarType beta;
      ScalarType norm_rhs = ip_rr0star;
      
      for (unsigned int i = 0; i < tag.iterations(); ++i)
      {
        tag.iters(i+1);
        tmp0 = viennacl::linalg::prod(matrix, p);
        alpha = ip_rr0star / viennacl::linalg::inner_prod(tmp0, r0star);

        s = residual - alpha*tmp0;
        
        tmp1 = viennacl::linalg::prod(matrix, s);
        omega = viennacl::linalg::inner_prod(tmp1, s) / viennacl::linalg::inner_prod(tmp1, tmp1);
        
        result += alpha * p + omega * s;
        residual = s - omega*tmp1;
        
        new_ip_rr0star = viennacl::linalg::inner_prod(residual,r0star);
        if (std::fabs(new_ip_rr0star / norm_rhs) < tag.tolerance() * tag.tolerance())
          break;
        
        beta = new_ip_rr0star / ip_rr0star * alpha/omega;
        ip_rr0star = new_ip_rr0star;

        p = residual + beta * (p - omega*tmp0);
      }
      
      //store last error estimate:
      tag.error(std::sqrt(std::fabs(new_ip_rr0star / norm_rhs)));
      
      return result;
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
      typedef typename VectorType::value_type      ScalarType;
      VectorType result(rhs.size());
      result.clear();

      VectorType residual = rhs;
      precond.apply(residual);
      VectorType r0star = residual;  //can be chosen arbitrarily in fact
      VectorType tmp0(rhs.size());
      VectorType tmp1(rhs.size());
      VectorType s(rhs.size());
      
      VectorType p = residual;

      ScalarType ip_rr0star = viennacl::linalg::inner_prod(residual,r0star);
      ScalarType alpha;
      ScalarType omega;
      ScalarType new_ip_rr0star = 0;
      ScalarType beta;
      ScalarType norm_rhs = ip_rr0star;
      
      for (unsigned int i = 0; i < tag.iterations(); ++i)
      {
        tag.iters(i+1);
        tmp0 = viennacl::linalg::prod(matrix, p);
        precond.apply(tmp0);
        alpha = ip_rr0star / viennacl::linalg::inner_prod(tmp0, r0star);

        s = residual - alpha*tmp0;

        tmp1 = viennacl::linalg::prod(matrix, s);
        precond.apply(tmp1);
        omega = viennacl::linalg::inner_prod(tmp1, s) / viennacl::linalg::inner_prod(tmp1, tmp1);
        
        result += alpha * p + omega * s;
        residual = s - omega*tmp1;
        
        new_ip_rr0star = viennacl::linalg::inner_prod(residual,r0star);
        if (std::fabs(new_ip_rr0star / norm_rhs) < tag.tolerance() * tag.tolerance() )
          break;
        
        beta = new_ip_rr0star / ip_rr0star * alpha/omega;
        ip_rr0star = new_ip_rr0star;

        p = residual + beta * (p - omega*tmp0);
      }
      
      //store last error estimate:
      tag.error(std::sqrt(std::fabs(new_ip_rr0star / norm_rhs)));
      
      return result;
    }

  }
}

#endif
