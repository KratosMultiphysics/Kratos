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

#ifndef VIENNACL_BICGSTAB_TUNED_HPP
#define VIENNACL_BICGSTAB_TUNED_HPP

/** @file bicgstab.hpp
    @brief The stabilized bi-conjugate gradient method is implemented here
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"

namespace viennacl
{
  namespace linalg
  {
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
    VectorType solve_tuned(const MatrixType & matrix, VectorType const & rhs, bicgstab_tag const & tag)
    {
      typedef typename viennacl::tools::result_of::value_type<VectorType>::type            ScalarType;
      typedef typename viennacl::tools::CPU_SCALAR_TYPE_DEDUCER<ScalarType>::ResultType    CPU_ScalarType;
      unsigned int problem_size = viennacl::tools::traits::size(rhs);
      VectorType result(problem_size);
      viennacl::tools::traits::clear(result);

      VectorType residual = rhs;
      VectorType p = rhs;
      VectorType r0star = rhs;
      VectorType tmp0(problem_size);
      VectorType tmp1(problem_size);
      VectorType s(problem_size);
      VectorType s2(problem_size);

      ScalarType ip_rr0star = viennacl::linalg::inner_prod(rhs,r0star);
      ScalarType error_estimate = 0;
      CPU_ScalarType norm_rhs_host = ip_rr0star;
      //CPU_ScalarType beta;
      ScalarType alpha;
      //CPU_ScalarType omega;
      ScalarType inner_prod_temp; //temporary variable for inner product computation
      ScalarType new_ip_rr0star = 0;

      
      size_t num_threads = 256; //also give 512 a try
      viennacl::ocl::kernel & k1 = viennacl::ocl::get_kernel(viennacl::linalg::kernels::compressed_matrix<CPU_ScalarType, 1>::program_name(),
                                                             "bicgstab_kernel1");
      k1.local_work_size(0, num_threads);
      k1.global_work_size(0, k1.local_work_size(0));

      
      viennacl::ocl::kernel & k2 = viennacl::ocl::get_kernel(viennacl::linalg::kernels::compressed_matrix<CPU_ScalarType, 1>::program_name(),
                                                             "bicgstab_kernel2");
      k2.local_work_size(0, num_threads);
      k2.global_work_size(0, k2.local_work_size(0));
      
      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);
        tmp0 = viennacl::linalg::prod(matrix, p);
        
        //// kernel 1 start ////
        
        //
        // what it does:
        //
        
        //alpha = ip_rr0star / viennacl::linalg::inner_prod(tmp0, r0star);
        //s = residual - alpha*tmp0;
        
        viennacl::ocl::enqueue(k1(tmp0, r0star, residual, s, alpha, ip_rr0star, 
                                  viennacl::ocl::local_mem(sizeof(CPU_ScalarType) * k1.local_work_size(0)),
                                  static_cast<cl_uint>(tmp0.size()))
                              );        
        
        //// kernel 1 end ////
        
        tmp1 = viennacl::linalg::prod(matrix, s);
        
        //// kernel 2 start ////
        
        //
        // what it does:
        //
        
        //omega = viennacl::linalg::inner_prod(tmp1, s) / viennacl::linalg::inner_prod(tmp1, tmp1);
        //result += alpha * p + omega * s;
        //residual = s - omega * tmp1;
        
        //new_ip_rr0star = viennacl::linalg::inner_prod(residual,r0star);
        //beta = new_ip_rr0star / ip_rr0star * alpha/omega;
        //p = residual + beta * (p - omega*tmp0);
        //ip_rr0star = new_ip_rr0star;

        viennacl::ocl::enqueue(k2(tmp0, tmp1, r0star, s, p, result, residual,
                                  alpha, ip_rr0star, error_estimate,
                                  viennacl::ocl::local_mem(sizeof(CPU_ScalarType) * k2.local_work_size(0)),
                                  static_cast<cl_uint>(tmp1.size()))
                              );        
        
        //// kernel 2 end ////
        
        //termination criterion:
        if (std::fabs( static_cast<CPU_ScalarType>(error_estimate) / norm_rhs_host) < tag.tolerance() * tag.tolerance())
          break;
      }
      
      //store last error estimate:
      tag.error(std::sqrt(std::fabs(viennacl::linalg::inner_prod(residual, residual) / norm_rhs_host)));
      
      return result;
    }

    template <typename MatrixType, typename VectorType>
    VectorType solve_tuned(const MatrixType & matrix, VectorType const & rhs, bicgstab_tag const & tag, viennacl::linalg::no_precond)
    {
      return solve_tuned(matrix, rhs, tag);
    }


  }
}

#endif
