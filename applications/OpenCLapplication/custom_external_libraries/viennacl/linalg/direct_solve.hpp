#ifndef VIENNACL_LINALG_DIRECT_SOLVE_HPP_
#define VIENNACL_LINALG_DIRECT_SOLVE_HPP_

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

/** @file viennacl/linalg/direct_solve.hpp
    @brief Implementations of dense direct solvers are found here.
*/

#include "viennacl/forwards.h"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/host_based/direct_solve.hpp"

#ifdef VIENNACL_WITH_OPENCL
  #include "viennacl/linalg/opencl/direct_solve.hpp"
#endif

#ifdef VIENNACL_WITH_CUDA
  #include "viennacl/linalg/cuda/direct_solve.hpp"
#endif

namespace viennacl
{
  namespace linalg
  {

    //
    // A \ B:
    //

    /** @brief Direct inplace solver for dense triangular systems. Matlab notation: A \ B
    *
    * @param A    The system matrix
    * @param B    The matrix of row vectors, where the solution is directly written to
    */
    template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
    void inplace_solve(const matrix_base<NumericT, F1> & A, matrix_base<NumericT, F2> & B, SOLVERTAG)
    {
      assert( (viennacl::traits::size1(A) == viennacl::traits::size2(A)) && bool("Size check failed in inplace_solve(): size1(A) != size2(A)"));
      assert( (viennacl::traits::size1(A) == viennacl::traits::size1(B)) && bool("Size check failed in inplace_solve(): size1(A) != size1(B)"));

      switch (viennacl::traits::handle(A).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(A, B, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(A, B, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(A, B, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
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
      assert( (viennacl::traits::size1(A) == viennacl::traits::size2(A))       && bool("Size check failed in inplace_solve(): size1(A) != size2(A)"));
      assert( (viennacl::traits::size1(A) == viennacl::traits::size1(proxy_B)) && bool("Size check failed in inplace_solve(): size1(A) != size1(B^T)"));

      switch (viennacl::traits::handle(A).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(A, proxy_B, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(A, proxy_B, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(A, proxy_B, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
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
      assert( (viennacl::traits::size1(proxy_A) == viennacl::traits::size2(proxy_A)) && bool("Size check failed in inplace_solve(): size1(A) != size2(A)"));
      assert( (viennacl::traits::size1(proxy_A) == viennacl::traits::size1(B))       && bool("Size check failed in inplace_solve(): size1(A^T) != size1(B)"));

      switch (viennacl::traits::handle(proxy_A.lhs()).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(proxy_A, B, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(proxy_A, B, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(proxy_A, B, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
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
      assert( (viennacl::traits::size1(proxy_A) == viennacl::traits::size2(proxy_A)) && bool("Size check failed in inplace_solve(): size1(A) != size2(A)"));
      assert( (viennacl::traits::size1(proxy_A) == viennacl::traits::size1(proxy_B)) && bool("Size check failed in inplace_solve(): size1(A^T) != size1(B^T)"));

      switch (viennacl::traits::handle(proxy_A.lhs()).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(proxy_A, proxy_B, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(proxy_A, proxy_B, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(proxy_A, proxy_B, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
    }

    //
    // A \ b
    //

    template <typename NumericT, typename F, typename SOLVERTAG>
    void inplace_solve(const matrix_base<NumericT, F> & mat,
                             vector_base<NumericT> & vec,
                       SOLVERTAG)
    {
      assert( (mat.size1() == vec.size()) && bool("Size check failed in inplace_solve(): size1(A) != size(b)"));
      assert( (mat.size2() == vec.size()) && bool("Size check failed in inplace_solve(): size2(A) != size(b)"));

      switch (viennacl::traits::handle(mat).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(mat, vec, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(mat, vec, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(mat, vec, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
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
      assert( (proxy.lhs().size1() == vec.size()) && bool("Size check failed in inplace_solve(): size1(A) != size(b)"));
      assert( (proxy.lhs().size2() == vec.size()) && bool("Size check failed in inplace_solve(): size2(A) != size(b)"));

      switch (viennacl::traits::handle(proxy.lhs()).get_active_handle_id())
      {
        case viennacl::MAIN_MEMORY:
          viennacl::linalg::host_based::inplace_solve(proxy, vec, SOLVERTAG());
          break;
#ifdef VIENNACL_WITH_OPENCL
        case viennacl::OPENCL_MEMORY:
          viennacl::linalg::opencl::inplace_solve(proxy, vec, SOLVERTAG());
          break;
#endif
#ifdef VIENNACL_WITH_CUDA
        case viennacl::CUDA_MEMORY:
          viennacl::linalg::cuda::inplace_solve(proxy, vec, SOLVERTAG());
          break;
#endif
        case viennacl::MEMORY_NOT_INITIALIZED:
          throw memory_exception("not initialised!");
        default:
          throw memory_exception("not implemented");
      }
    }

    /////////////////// general wrappers for non-inplace solution //////////////////////


    /** @brief Convenience functions for C = solve(A, B, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param A    The system matrix
    * @param B    The matrix of load vectors
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
    matrix<NumericT, F2> solve(const matrix_base<NumericT, F1> & A,
                               const matrix_base<NumericT, F2> & B,
                               SOLVERTAG tag)
    {
      // do an inplace solve on the result vector:
      matrix<NumericT, F2> result(B);

      inplace_solve(A, result, tag);

      return result;
    }


    //////////

    /** @brief Convenience functions for C = solve(A, B^T, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param A    The system matrix
    * @param proxy  The transposed load vector
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
    matrix<NumericT, F2> solve(const matrix_base<NumericT, F1> & A,
                               const matrix_expression< const matrix_base<NumericT, F2>, const matrix_base<NumericT, F2>, op_trans> & proxy,
                               SOLVERTAG tag)
    {
      // do an inplace solve on the result vector:
      matrix<NumericT, F2> result(proxy);

      inplace_solve(A, result, tag);

      return result;
    }

    /** @brief Convenience functions for result = solve(mat, vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param mat    The system matrix
    * @param vec    The load vector
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename SOLVERTAG>
    vector<NumericT> solve(const matrix_base<NumericT, F1> & mat,
                           const vector_base<NumericT> & vec,
                           SOLVERTAG const & tag)
    {
      // do an inplace solve on the result vector:
      vector<NumericT> result(vec);

      inplace_solve(mat, result, tag);

      return result;
    }


    ///////////// transposed system matrix:
    /** @brief Convenience functions for result = solve(trans(mat), B, some_tag()); Creates a temporary result matrix and forwards the request to inplace_solve()
    *
    * @param proxy  The transposed system matrix proxy
    * @param B      The matrix of load vectors
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
    matrix<NumericT, F2> solve(const matrix_expression< const matrix_base<NumericT, F1>, const matrix_base<NumericT, F1>, op_trans> & proxy,
                               const matrix_base<NumericT, F2> & B,
                               SOLVERTAG tag)
    {
      // do an inplace solve on the result vector:
      matrix<NumericT, F2> result(B);

      inplace_solve(proxy, result, tag);

      return result;
    }


    /** @brief Convenience functions for result = solve(trans(mat), vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param proxy_A  The transposed system matrix proxy
    * @param proxy_B  The transposed matrix of load vectors, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename F2, typename SOLVERTAG>
    matrix<NumericT, F2> solve(const matrix_expression< const matrix_base<NumericT, F1>, const matrix_base<NumericT, F1>, op_trans> & proxy_A,
                               const matrix_expression< const matrix_base<NumericT, F2>, const matrix_base<NumericT, F2>, op_trans> & proxy_B,
                               SOLVERTAG tag)
    {
      // do an inplace solve on the result vector:
      matrix<NumericT, F2> result(proxy_B);

      inplace_solve(proxy_A, result, tag);

      return result;
    }

    /** @brief Convenience functions for result = solve(trans(mat), vec, some_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param proxy  The transposed system matrix proxy
    * @param vec    The load vector, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template <typename NumericT, typename F1, typename SOLVERTAG>
    vector<NumericT> solve(const matrix_expression< const matrix_base<NumericT, F1>, const matrix_base<NumericT, F1>, op_trans> & proxy,
                           const vector_base<NumericT> & vec,
                           SOLVERTAG const & tag)
    {
      // do an inplace solve on the result vector:
      vector<NumericT> result(vec);

      inplace_solve(proxy, result, tag);

      return result;
    }


  }
}

#endif
