#ifndef VIENNACL_LINALG_OPENCL_DIRECT_SOLVE_HPP
#define VIENNACL_LINALG_OPENCL_DIRECT_SOLVE_HPP

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

/** @file viennacl/linalg/opencl/direct_solve.hpp
    @brief Implementations of dense direct solvers are found here.
*/

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/linalg/opencl/kernels/matrix_solve.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace opencl
    {
      namespace detail
      {
        inline cl_uint get_option_for_solver_tag(viennacl::linalg::upper_tag)      { return 0; }
        inline cl_uint get_option_for_solver_tag(viennacl::linalg::unit_upper_tag) { return (1 << 0); }
        inline cl_uint get_option_for_solver_tag(viennacl::linalg::lower_tag)      { return (1 << 2); }
        inline cl_uint get_option_for_solver_tag(viennacl::linalg::unit_lower_tag) { return (1 << 2) | (1 << 0); }

        template <typename M1, typename M2, typename KernelType>
        void inplace_solve_impl(M1 const & A, M2 & B, KernelType & k)
        {
          viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(A),
                                   cl_uint(viennacl::traits::start1(A)),         cl_uint(viennacl::traits::start2(A)),
                                   cl_uint(viennacl::traits::stride1(A)),        cl_uint(viennacl::traits::stride2(A)),
                                   cl_uint(viennacl::traits::size1(A)),          cl_uint(viennacl::traits::size2(A)),
                                   cl_uint(viennacl::traits::internal_size1(A)), cl_uint(viennacl::traits::internal_size2(A)),
                                   viennacl::traits::opencl_handle(B),
                                   cl_uint(viennacl::traits::start1(B)),         cl_uint(viennacl::traits::start2(B)),
                                   cl_uint(viennacl::traits::stride1(B)),        cl_uint(viennacl::traits::stride2(B)),
                                   cl_uint(viennacl::traits::size1(B)),          cl_uint(viennacl::traits::size2(B)),
                                   cl_uint(viennacl::traits::internal_size1(B)), cl_uint(viennacl::traits::internal_size2(B))
                                  )
                                );
        }
      }


      //
      // Note: By convention, all size checks are performed in the calling frontend. No need to double-check here.
      //

      ////////////////// upper triangular solver (upper_tag) //////////////////////////////////////
      /** @brief Direct inplace solver for dense triangular systems. Matlab notation: A \ B
      *
      * @param A    The system matrix
      * @param B    The matrix of row vectors, where the solution is directly written to
      */
      template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
      void inplace_solve(const matrix_base<NumericT, F1> & A, matrix_base<NumericT, F2> & B, SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

        typedef viennacl::linalg::opencl::kernels::matrix_solve<NumericT, F1, F2>    KernelClass;
        KernelClass::init(ctx);

        std::stringstream ss;
        ss << SOLVERTAG::name() << "_solve";
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), ss.str());

        k.global_work_size(0, B.size2() * k.local_work_size());
        detail::inplace_solve_impl(A, B, k);
      }

      /** @brief Direct inplace solver for dense triangular systems with transposed right hand side
      *
      * @param A       The system matrix
      * @param proxy_B The transposed matrix of row vectors, where the solution is directly written to
      */
      template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
      void inplace_solve(const matrix_base<NumericT, F1> & A,
                         matrix_expression< const matrix_base<NumericT, F2>, const matrix_base<NumericT, F2>, op_trans> proxy_B,
                         SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

        typedef viennacl::linalg::opencl::kernels::matrix_solve<NumericT, F1, F2>    KernelClass;
        KernelClass::init(ctx);

        std::stringstream ss;
        ss << SOLVERTAG::name() << "_trans_solve";
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), ss.str());

        k.global_work_size(0, proxy_B.lhs().size1() * k.local_work_size());
        detail::inplace_solve_impl(A, proxy_B.lhs(), k);
      }

      //upper triangular solver for transposed lower triangular matrices
      /** @brief Direct inplace solver for dense triangular systems that stem from transposed triangular systems
      *
      * @param proxy_A  The system matrix proxy
      * @param B        The matrix holding the load vectors, where the solution is directly written to
      */
      template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
      void inplace_solve(const matrix_expression< const matrix_base<NumericT, F1>, const matrix_base<NumericT, F1>, op_trans> & proxy_A,
                         matrix_base<NumericT, F2> & B,
                         SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(B).context());

        typedef viennacl::linalg::opencl::kernels::matrix_solve<NumericT, F1, F2>    KernelClass;
        KernelClass::init(ctx);

        std::stringstream ss;
        ss << "trans_" << SOLVERTAG::name() << "_solve";
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), ss.str());

        k.global_work_size(0, B.size2() * k.local_work_size());
        detail::inplace_solve_impl(proxy_A.lhs(), B, k);
      }

      /** @brief Direct inplace solver for dense transposed triangular systems with transposed right hand side. Matlab notation: A' \ B'
      *
      * @param proxy_A  The system matrix proxy
      * @param proxy_B  The matrix holding the load vectors, where the solution is directly written to
      */
      template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
      void inplace_solve(const matrix_expression< const matrix_base<NumericT, F1>, const matrix_base<NumericT, F1>, op_trans> & proxy_A,
                               matrix_expression< const matrix_base<NumericT, F2>, const matrix_base<NumericT, F2>, op_trans>   proxy_B,
                         SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(proxy_A.lhs()).context());

        typedef viennacl::linalg::opencl::kernels::matrix_solve<NumericT, F1, F2>    KernelClass;
        KernelClass::init(ctx);

        std::stringstream ss;
        ss << "trans_" << SOLVERTAG::name() << "_trans_solve";
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), ss.str());

        k.global_work_size(0, proxy_B.lhs().size1() * k.local_work_size());
        detail::inplace_solve_impl(proxy_A.lhs(), proxy_B.lhs(), k);
      }



      //
      //  Solve on vector
      //

      template <typename NumericT, typename F, typename SOLVERTAG>
      void inplace_solve(const matrix_base<NumericT, F> & mat,
                               vector_base<NumericT> & vec,
                         SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(mat).context());

        typedef viennacl::linalg::opencl::kernels::matrix<NumericT, F>  KernelClass;
        KernelClass::init(ctx);

        cl_uint options = detail::get_option_for_solver_tag(SOLVERTAG());
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), "triangular_substitute_inplace");

        k.global_work_size(0, k.local_work_size());
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(mat),
                                 cl_uint(viennacl::traits::start1(mat)),         cl_uint(viennacl::traits::start2(mat)),
                                 cl_uint(viennacl::traits::stride1(mat)),        cl_uint(viennacl::traits::stride2(mat)),
                                 cl_uint(viennacl::traits::size1(mat)),          cl_uint(viennacl::traits::size2(mat)),
                                 cl_uint(viennacl::traits::internal_size1(mat)), cl_uint(viennacl::traits::internal_size2(mat)),
                                 viennacl::traits::opencl_handle(vec),
                                 cl_uint(viennacl::traits::start(vec)),
                                 cl_uint(viennacl::traits::stride(vec)),
                                 cl_uint(viennacl::traits::size(vec)),
                                 options
                                )
                              );
      }

      /** @brief Direct inplace solver for dense upper triangular systems that stem from transposed lower triangular systems
      *
      * @param proxy    The system matrix proxy
      * @param vec    The load vector, where the solution is directly written to
      */
      template <typename NumericT, typename F, typename SOLVERTAG>
      void inplace_solve(const matrix_expression< const matrix_base<NumericT, F>, const matrix_base<NumericT, F>, op_trans> & proxy,
                         vector_base<NumericT> & vec,
                         SOLVERTAG)
      {
        viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(vec).context());

        typedef viennacl::linalg::opencl::kernels::matrix<NumericT, F>  KernelClass;
        KernelClass::init(ctx);

        cl_uint options = detail::get_option_for_solver_tag(SOLVERTAG()) | 0x02;  //add transpose-flag
        viennacl::ocl::kernel & k = ctx.get_kernel(KernelClass::program_name(), "triangular_substitute_inplace");

        k.global_work_size(0, k.local_work_size());
        viennacl::ocl::enqueue(k(viennacl::traits::opencl_handle(proxy.lhs()),
                                 cl_uint(viennacl::traits::start1(proxy.lhs())),         cl_uint(viennacl::traits::start2(proxy.lhs())),
                                 cl_uint(viennacl::traits::stride1(proxy.lhs())),        cl_uint(viennacl::traits::stride2(proxy.lhs())),
                                 cl_uint(viennacl::traits::size1(proxy.lhs())),          cl_uint(viennacl::traits::size2(proxy.lhs())),
                                 cl_uint(viennacl::traits::internal_size1(proxy.lhs())), cl_uint(viennacl::traits::internal_size2(proxy.lhs())),
                                 viennacl::traits::opencl_handle(vec),
                                 cl_uint(viennacl::traits::start(vec)),
                                 cl_uint(viennacl::traits::stride(vec)),
                                 cl_uint(viennacl::traits::size(vec)),
                                 options
                                )
                              );
      }


    }
  }
}

#endif
