#ifndef VIENNACL_LINALG_NMF_HPP
#define VIENNACL_LINALG_NMF_HPP

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

/** @file viennacl/linalg/nmf.hpp
    @brief Provides a nonnegative matrix factorization implementation.  Experimental in 1.3.x.
    
    Contributed by Volodymyr Kysenko.
*/


#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/linalg/kernels/nmf_kernels.h"

namespace viennacl
{
  namespace linalg
  {
    
    class nmf_config
    {
      public:
        nmf_config(double val_epsilon = 1e-4,
                   double val_epsilon_stagnation = 1e-5,
                   std::size_t num_max_iters = 10000,
                   std::size_t num_check_iters = 100)
         : eps_(val_epsilon), stagnation_eps_(val_epsilon_stagnation),
           max_iters_(num_max_iters),
           check_after_steps_( (num_check_iters > 0) ? num_check_iters : 1),
           iters_(0) {}
         
        /** @brief Returns the relative tolerance for convergence */
        double tolerance() const { return eps_; }
        
        /** @brief Sets the relative tolerance for convergence, i.e. norm(V - W * H) / norm(V - W_init * H_init) */
        void tolerance(double e) { eps_ = e; }
        
        /** @brief Relative tolerance for the stagnation check */
        double stagnation_tolerance() const { return stagnation_eps_; }
        
        /** @brief Sets the tolerance for the stagnation check (i.e. the minimum required relative change of the residual between two iterations) */
        void stagnation_tolerance(double e) { stagnation_eps_ = e; }
        
        /** @brief Returns the maximum number of iterations for the NMF algorithm */
        std::size_t max_iterations() const { return max_iters_; }
        /** @brief Sets the maximum number of iterations for the NMF algorithm */
        void max_iterations(std::size_t m) { max_iters_ = m; }
        
        
        /** @brief Returns the number of iterations of the last NMF run using this configuration object */
        std::size_t iters() const { return iters_; }
        
        
        /** @brief Number of steps after which the convergence of NMF should be checked (again) */
        std::size_t check_after_steps() const { return check_after_steps_; }
        
        /** @brief Set the number of steps after which the convergence of NMF should be checked (again) */
        void check_after_steps(std::size_t c) { if (c > 0) check_after_steps_ = c; }
        
        
        template <typename ScalarType>
        friend void nmf(viennacl::matrix<ScalarType> const & V,
                        viennacl::matrix<ScalarType> & W,
                        viennacl::matrix<ScalarType> & H,
                        nmf_config const & conf);
        
      private:
        double eps_;
        double stagnation_eps_;
        std::size_t max_iters_;
        std::size_t check_after_steps_;
        mutable std::size_t iters_;
    };

    
    /** @brief The nonnegative matrix factorization (approximation) algorithm as suggested by Lee and Seung. Factorizes a matrix V with nonnegative entries into matrices W and H such that ||V - W*H|| is minimized.
     * 
     * @param V     Input matrix 
     * @param W     First factor
     * @param H     Second factor
     * @param conf  A configuration object holding tolerances and the like 
     */
    template <typename ScalarType>
    void nmf(viennacl::matrix<ScalarType> const & V,
             viennacl::matrix<ScalarType> & W,
             viennacl::matrix<ScalarType> & H,
             nmf_config const & conf)
    {
      const std::string NMF_MUL_DIV_KERNEL = "el_wise_mul_div";
      const std::string NMF_SUB_KERNEL = "sub_wise";
      
      viennacl::linalg::kernels::nmf<ScalarType, 1>::init();
      
      assert(V.size1() == W.size1() && V.size2() == H.size2() && "Dimensions of W and H don't allow for V = W * H");
      assert(W.size2() == H.size1() && "Dimensions of W and H don't match, prod(W, H) impossible");

      std::size_t k = W.size2();
      conf.iters_ = 0;
      
      viennacl::matrix<ScalarType> wn(V.size1(), k);
      viennacl::matrix<ScalarType> wd(V.size1(), k);
      viennacl::matrix<ScalarType> wtmp(V.size1(), V.size2());

      viennacl::matrix<ScalarType> hn(k, V.size2());
      viennacl::matrix<ScalarType> hd(k, V.size2());
      viennacl::matrix<ScalarType> htmp(k, k);

      viennacl::matrix<ScalarType> appr(V.size1(), V.size2());
      viennacl::vector<ScalarType> diff(V.size1() * V.size2());

      ScalarType last_diff = 0;
      ScalarType diff_init = 0;
      bool stagnation_flag = false;

      
      for (std::size_t i = 0; i < conf.max_iterations(); i++)
      {
        conf.iters_ = i + 1;
        {
          hn   = viennacl::linalg::prod(trans(W), V);
          htmp = viennacl::linalg::prod(trans(W), W);
          hd   = viennacl::linalg::prod(htmp, H);

          viennacl::ocl::kernel & mul_div_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::nmf<ScalarType, 1>::program_name(), 
                                                                             NMF_MUL_DIV_KERNEL);
          viennacl::ocl::enqueue(mul_div_kernel(H, hn, hd, cl_uint(H.internal_size1() * H.internal_size2())));
        }
        {
          wn   = viennacl::linalg::prod(V, trans(H));
          wtmp = viennacl::linalg::prod(W, H);
          wd   = viennacl::linalg::prod(wtmp, trans(H));

          viennacl::ocl::kernel & mul_div_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::nmf<ScalarType, 1>::program_name(), 
                                                                             NMF_MUL_DIV_KERNEL);
          
          viennacl::ocl::enqueue(mul_div_kernel(W, wn, wd, cl_uint(W.internal_size1() * W.internal_size2())));
        }

        if (i % conf.check_after_steps() == 0)  //check for convergence
        {
          appr = viennacl::linalg::prod(W, H);

          viennacl::ocl::kernel & sub_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::nmf<ScalarType, 1>::program_name(), 
                                                                         NMF_SUB_KERNEL);
          //this is a cheat: save difference of two matrix into vector to get norm_2
          viennacl::ocl::enqueue(sub_kernel(appr, V, diff, cl_uint(V.size1() * V.size2())));
          ScalarType diff_val = viennacl::linalg::norm_2(diff);
          
          if (i == 0)
            diff_init = diff_val;
          
          std::cout << diff_val / diff_init << std::endl;

          // Approximation check
          if (diff_val / diff_init < conf.tolerance())
            break;

          // Stagnation check
          if (fabs(diff_val - last_diff) / (diff_val * conf.check_after_steps()) < conf.stagnation_tolerance()) //avoid situations where convergence stagnates
          {
            if (stagnation_flag)       // iteration stagnates (two iterates with no notable progress)
              break;
            else                       // record stagnation in this iteration
              stagnation_flag = true;
          }
          else                         // good progress in this iteration, so unset stagnation flag
            stagnation_flag = false;

          // prepare for next iterate:
          last_diff = diff_val;
        }
      }
      
      
    }
  }
}

#endif