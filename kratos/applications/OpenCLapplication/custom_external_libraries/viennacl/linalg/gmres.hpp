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

#ifndef _VIENNACL_GMRES_HPP_
#define _VIENNACL_GMRES_HPP_

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"

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
        * @param krylov_dim     The maximum dimension of the Krylov space before restart
        * @param max_restarts   The maximum number of restarts 
        */
        gmres_tag(double tol = 1e-8, unsigned int krylov_dim = 20, unsigned int max_restarts = 5) 
         : _tol(tol), _max_restarts(max_restarts), _krylov_dim(krylov_dim), iters_taken_(0) {};
        
        double tolerance() const { return _tol; }
        unsigned int max_restarts() const { return _max_restarts; }
        unsigned int krylov_dim() const { return _krylov_dim; }
        
        //number of solver iterations:
        unsigned int iters() const { return iters_taken_; }
        void iters(unsigned int i) const { iters_taken_ = i; }
        
        //esimated error :
        double error() const { return last_error_; }
        void error(double e) const { last_error_ = e; }
        
      private:
        double _tol;
        unsigned int _max_restarts;
        unsigned int _krylov_dim;
        
        //return values from solver
        mutable unsigned int iters_taken_;
        mutable double last_error_;
    };
    
    namespace
    {
      
      template <typename SRC_VECTOR, typename DEST_VECTOR>
      void gmres_copy_helper(SRC_VECTOR const & src, DEST_VECTOR & dest, unsigned int len)
      {
        for (unsigned int i=0; i<len; ++i)
          dest[i] = src[i];
      }

      template <typename ScalarType, typename DEST_VECTOR>
      void gmres_copy_helper(viennacl::vector<ScalarType> const & src, DEST_VECTOR & dest, unsigned int len)
      {
        viennacl::copy(src.begin(), src.begin() + len, dest.begin());
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
      typedef typename MatrixType::value_type                             ScalarType;
      typedef typename viennacl::tools::CPU_SCALAR_TYPE_DEDUCER<ScalarType>::ResultType    CPU_ScalarType;
      VectorType result(rhs.size());
      result.clear();
      
      VectorType res(rhs.size());
      VectorType v_k_tilde(rhs.size());
      VectorType v_k_tilde_temp(rhs.size());
      
      std::vector< std::vector<CPU_ScalarType> > R(tag.krylov_dim());
      std::vector<CPU_ScalarType> projection_rhs(tag.krylov_dim());
      std::vector<VectorType> U(tag.krylov_dim());
      std::vector<CPU_ScalarType> cpu_temp(rhs.size());

      const CPU_ScalarType gpu_scalar_minus_1 = static_cast<CPU_ScalarType>(-1);    //representing the scalar '-1' on the GPU. Prevents blocking write operations
      const CPU_ScalarType gpu_scalar_1 = static_cast<CPU_ScalarType>(1);    //representing the scalar '1' on the GPU. Prevents blocking write operations
      const CPU_ScalarType gpu_scalar_2 = static_cast<CPU_ScalarType>(2);    //representing the scalar '2' on the GPU. Prevents blocking write operations
      
      CPU_ScalarType norm_rhs = viennacl::linalg::norm_2(rhs);
      
      unsigned int k;
      for (k = 0; k < tag.krylov_dim(); ++k)
      {
        R[k].resize(tag.krylov_dim()); 
        U[k].resize(rhs.size());
      }

      //std::cout << "Starting GMRES..." << std::endl;
      tag.iters(0);

      for (unsigned int it = 0; it < tag.max_restarts(); ++it)
      {
        res = rhs;
        res -= viennacl::linalg::prod(matrix, result);  //initial guess zero
        precond.apply(res);
        
        CPU_ScalarType rho_0 = viennacl::linalg::norm_2(res); 
        CPU_ScalarType rho = 1.0;

        if (rho_0 / norm_rhs < tag.tolerance())
        {
          //std::cout << "Allowed Error reached at begin of loop" << std::endl;
          tag.error(rho_0 / norm_rhs);
          return result;
        }

        res /= rho_0;
        for (k=0; k<tag.krylov_dim(); ++k)
        {
          R[k].clear();
          U[k].clear();
          R[k].resize(tag.krylov_dim()); 
          U[k].resize(rhs.size());
        }

        for (k = 0; k < tag.krylov_dim(); ++k)
        {
          tag.iters( tag.iters() + 1 ); //increase iteration counter

          //compute v_k = A * v_{k-1} via Householder matrices
          if (k == 0)
          {
            v_k_tilde = viennacl::linalg::prod(matrix, res);
            precond.apply(v_k_tilde);
          }
          else
          {
            v_k_tilde.clear();
            v_k_tilde[k-1] = gpu_scalar_1;
            //Householder rotations part 1
            for (int i = k-1; i > -1; --i)
              v_k_tilde -= U[i] * (viennacl::linalg::inner_prod(U[i], v_k_tilde) * gpu_scalar_2);

            v_k_tilde_temp = viennacl::linalg::prod(matrix, v_k_tilde);
            precond.apply(v_k_tilde_temp);
            v_k_tilde = v_k_tilde_temp;

            //Householder rotations part 2
            for (unsigned int i = 0; i < k; ++i)
              v_k_tilde -= U[i] * (viennacl::linalg::inner_prod(U[i], v_k_tilde) * gpu_scalar_2);
          }

          U[k].clear();
          U[k].resize(rhs.size());
          for (unsigned int i = 0; i < k; ++i)
            U[k][i] = v_k_tilde[i];

          U[k][k] = std::sqrt( viennacl::linalg::inner_prod(v_k_tilde, v_k_tilde) - viennacl::linalg::inner_prod(U[k], U[k]) );

          //generic copy from U to R (use of plain copy() leads to problems under Visual Studio 2008)
          gmres_copy_helper(U[k], R[k], k+1);
          
          U[k] -= v_k_tilde;
          U[k] *= gpu_scalar_minus_1 / viennacl::linalg::norm_2( U[k] );

          res -= U[k] * (viennacl::linalg::inner_prod( U[k], res ) * gpu_scalar_2);

          projection_rhs[k] = res[k];
          rho *= std::sin( std::acos(projection_rhs[k] / rho) );

          if (std::fabs(rho * rho_0 / norm_rhs) < tag.tolerance())
          {
            //std::cout << "Krylov space big enough" << endl;
            tag.error( std::fabs(rho*rho_0 / norm_rhs) );
            break;
          }
          
        }
        
        //inplace solution of the upper triangular matrix:
        for (int i=k-1; i>-1; --i)
        {
          for (unsigned int j=i+1; j<k; ++j)
            //temp_rhs[i] -= R[i][j] * temp_rhs[j];   //if R is not transposed
            projection_rhs[i] -= R[j][i] * projection_rhs[j];     //R is transposed
            
          projection_rhs[i] /= R[i][i];
        }
        
        res *= projection_rhs[0];
        
        for (unsigned int i = 0; i < k-1; ++i)
          res[i] += projection_rhs[i+1];

        for (int i = k-1; i > -1; --i)
          res -= U[i] * (viennacl::linalg::inner_prod(U[i], res) * gpu_scalar_2);

        res *= rho_0;
        result += res;

        if ( std::fabs(rho*rho_0 / norm_rhs) < tag.tolerance() )
        {
          //std::cout << "Allowed Error reached at end of loop" << std::endl;
          tag.error(std::fabs(rho*rho_0 / norm_rhs));
          return result;
        }

        res = rhs;
        res -= viennacl::linalg::prod(matrix, result);
        //std::cout << "norm_2(r)=" << norm_2(r) << std::endl;
        //std::cout << "std::abs(rho*rho_0)=" << std::abs(rho*rho_0) << std::endl;
        //std::cout << r << std::endl; 

        tag.error(std::fabs(rho*rho_0));
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
