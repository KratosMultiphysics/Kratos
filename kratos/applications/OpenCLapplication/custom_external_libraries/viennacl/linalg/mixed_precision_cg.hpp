#ifndef VIENNACL_LINALG_MIXED_PRECISION_CG_HPP_
#define VIENNACL_LINALG_MIXED_PRECISION_CG_HPP_

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

/** @file viennacl/linalg/mixed_precision_cg.hpp
    @brief The conjugate gradient method using mixed precision is implemented here. Experimental.
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
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/backend/memory.hpp"

#include "viennacl/vector_proxy.hpp"

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for the conjugate gradient Used for supplying solver parameters and for dispatching the solve() function
    */
    class mixed_precision_cg_tag
    {
      public:
        /** @brief The constructor
        *
        * @param tol              Relative tolerance for the residual (solver quits if ||r|| < tol * ||r_initial||)
        * @param max_iterations   The maximum number of iterations
        * @param inner_tol        Inner tolerance for the low-precision iterations
        */
        mixed_precision_cg_tag(double tol = 1e-8, unsigned int max_iterations = 300, float inner_tol = 1e-2f) : tol_(tol), iterations_(max_iterations), inner_tol_(inner_tol) {}

        /** @brief Returns the relative tolerance */
        double tolerance() const { return tol_; }
        /** @brief Returns the relative tolerance */
        float inner_tolerance() const { return inner_tol_; }
        /** @brief Returns the maximum number of iterations */
        unsigned int max_iterations() const { return iterations_; }

        /** @brief Return the number of solver iterations: */
        unsigned int iters() const { return iters_taken_; }
        void iters(unsigned int i) const { iters_taken_ = i; }

        /** @brief Returns the estimated relative error at the end of the solver run */
        double error() const { return last_error_; }
        /** @brief Sets the estimated relative error at the end of the solver run */
        void error(double e) const { last_error_ = e; }


      private:
        double tol_;
        unsigned int iterations_;
        float inner_tol_;

        //return values from solver
        mutable unsigned int iters_taken_;
        mutable double last_error_;
    };


    const char * double_float_conversion_program =
    "#if defined(cl_khr_fp64)\n"
    "#  pragma OPENCL EXTENSION cl_khr_fp64: enable\n"
    "#elif defined(cl_amd_fp64)\n"
    "#  pragma OPENCL EXTENSION cl_amd_fp64: enable\n"
    "#endif\n"
    "__kernel void assign_double_to_float(\n"
    "          __global float * vec1,\n"
    "          __global const double * vec2, \n"
    "          unsigned int size) \n"
    "{ \n"
    "  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
    "    vec1[i] = (float)(vec2[i]);\n"
    "};\n\n"
    "__kernel void inplace_add_float_to_double(\n"
    "          __global double * vec1,\n"
    "          __global const float * vec2, \n"
    "          unsigned int size) \n"
    "{ \n"
    "  for (unsigned int i = get_global_id(0); i < size; i += get_global_size(0))\n"
    "    vec1[i] += (double)(vec2[i]);\n"
    "};\n";


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
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, mixed_precision_cg_tag const & tag)
    {
      //typedef typename VectorType::value_type      ScalarType;
      typedef typename viennacl::result_of::value_type<VectorType>::type        ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

      //TODO: Assert CPU_ScalarType == double

      //std::cout << "Starting CG" << std::endl;
      vcl_size_t problem_size = viennacl::traits::size(rhs);
      VectorType result(rhs);
      viennacl::traits::clear(result);

      VectorType residual = rhs;

      CPU_ScalarType ip_rr = viennacl::linalg::inner_prod(rhs, rhs);
      CPU_ScalarType new_ip_rr = 0;
      CPU_ScalarType norm_rhs_squared = ip_rr;

      if (norm_rhs_squared == 0) //solution is zero if RHS norm is zero
        return result;

      static bool first = true;

      viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(matrix).context());
      if (first)
      {
        ctx.add_program(double_float_conversion_program, "double_float_conversion_program");
      }

      viennacl::vector<float> residual_low_precision(problem_size, viennacl::traits::context(rhs));
      viennacl::vector<float> result_low_precision(problem_size, viennacl::traits::context(rhs));
      viennacl::vector<float> p_low_precision(problem_size, viennacl::traits::context(rhs));
      viennacl::vector<float> tmp_low_precision(problem_size, viennacl::traits::context(rhs));
      float inner_ip_rr = static_cast<float>(ip_rr);
      float new_inner_ip_rr = 0;
      float initial_inner_rhs_norm_squared = static_cast<float>(ip_rr);
      float alpha;
      float beta;

      viennacl::ocl::kernel & assign_double_to_float      = ctx.get_kernel("double_float_conversion_program", "assign_double_to_float");
      viennacl::ocl::kernel & inplace_add_float_to_double = ctx.get_kernel("double_float_conversion_program", "inplace_add_float_to_double");

      // transfer rhs to single precision:
      viennacl::ocl::enqueue( assign_double_to_float(p_low_precision.handle().opencl_handle(),
                                                     rhs.handle().opencl_handle(),
                                                     cl_uint(rhs.size())
                                                    ) );
      //std::cout << "copying p_low_precision..." << std::endl;
      //assign_double_to_float(p_low_precision.handle(), residual.handle(), residual.size());
      residual_low_precision = p_low_precision;

      // transfer matrix to single precision:
      viennacl::compressed_matrix<float> matrix_low_precision(matrix.size1(), matrix.size2(), matrix.nnz(), viennacl::traits::context(rhs));
      viennacl::backend::memory_copy(matrix.handle1(), const_cast<viennacl::backend::mem_handle &>(matrix_low_precision.handle1()), 0, 0, sizeof(cl_uint) * (matrix.size1() + 1) );
      viennacl::backend::memory_copy(matrix.handle2(), const_cast<viennacl::backend::mem_handle &>(matrix_low_precision.handle2()), 0, 0, sizeof(cl_uint) * (matrix.nnz()) );

      viennacl::ocl::enqueue( assign_double_to_float(matrix_low_precision.handle().opencl_handle(),
                                                     matrix.handle().opencl_handle(),
                                                     cl_uint(matrix.nnz())
                                                    ) );
      //std::cout << "copying matrix_low_precision..." << std::endl;
      //assign_double_to_float(const_cast<viennacl::backend::mem_handle &>(matrix_low_precision.handle()), matrix.handle(), matrix.nnz());

      //std::cout << "Starting CG solver iterations... " << std::endl;


      for (unsigned int i = 0; i < tag.max_iterations(); ++i)
      {
        tag.iters(i+1);

        // lower precision 'inner iteration'
        tmp_low_precision = viennacl::linalg::prod(matrix_low_precision, p_low_precision);

        alpha = inner_ip_rr / viennacl::linalg::inner_prod(tmp_low_precision, p_low_precision);
        result_low_precision += alpha * p_low_precision;
        residual_low_precision -= alpha * tmp_low_precision;

        new_inner_ip_rr = viennacl::linalg::inner_prod(residual_low_precision, residual_low_precision);

        beta = new_inner_ip_rr / inner_ip_rr;
        inner_ip_rr = new_inner_ip_rr;

        p_low_precision = residual_low_precision + beta * p_low_precision;



        if (new_inner_ip_rr < tag.inner_tolerance() * initial_inner_rhs_norm_squared || i == tag.max_iterations()-1)
        {
          //std::cout << "outer correction at i=" << i << std::endl;
          //result += result_low_precision;
          viennacl::ocl::enqueue( inplace_add_float_to_double(result.handle().opencl_handle(),
                                                              result_low_precision.handle().opencl_handle(),
                                                              cl_uint(result.size())
                                                             ) );

          // residual = b - Ax  (without introducing a temporary)
          residual = viennacl::linalg::prod(matrix, result);
          residual = rhs - residual;

          new_ip_rr = viennacl::linalg::inner_prod(residual, residual);
          if (new_ip_rr / norm_rhs_squared < tag.tolerance() *  tag.tolerance())//squared norms involved here
            break;

          // p_low_precision = residual;
          viennacl::ocl::enqueue( assign_double_to_float(p_low_precision.handle().opencl_handle(),
                                                         residual.handle().opencl_handle(),
                                                         cl_uint(residual.size())
                                                        ) );
          result_low_precision.clear();
          residual_low_precision = p_low_precision;
          initial_inner_rhs_norm_squared = static_cast<float>(new_ip_rr);
          inner_ip_rr = static_cast<float>(new_ip_rr);
        }
      }

      //store last error estimate:
      tag.error(std::sqrt(new_ip_rr / norm_rhs_squared));

      return result;
    }

    template <typename MatrixType, typename VectorType>
    VectorType solve(const MatrixType & matrix, VectorType const & rhs, mixed_precision_cg_tag const & tag, viennacl::linalg::no_precond)
    {
      return solve(matrix, rhs, tag);
    }


  }
}

#endif
