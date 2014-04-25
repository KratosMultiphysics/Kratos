#ifndef VIENNACL_LINALG_LU_HPP
#define VIENNACL_LINALG_LU_HPP

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

/** @file viennacl/linalg/lu.hpp
    @brief Implementations of LU factorization for row-major and column-major dense matrices.
*/

#include <algorithm>    //for std::min

#include "viennacl/matrix.hpp"
#include "viennacl/matrix_proxy.hpp"

#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/direct_solve.hpp"

namespace viennacl
{
  namespace linalg
  {
    /** @brief LU factorization of a row-major dense matrix.
    *
    * @param A    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    */
    template<typename SCALARTYPE>
    void lu_factorize(matrix<SCALARTYPE, viennacl::row_major> & A)
    {
      typedef matrix<SCALARTYPE, viennacl::row_major>  MatrixType;

      vcl_size_t max_block_size = 32;
      vcl_size_t num_blocks = (A.size2() - 1) / max_block_size + 1;
      std::vector<SCALARTYPE> temp_buffer(A.internal_size2() * max_block_size);

      // Iterate over panels
      for (vcl_size_t panel_id = 0; panel_id < num_blocks; ++panel_id)
      {
        vcl_size_t row_start = panel_id * max_block_size;
        vcl_size_t current_block_size = std::min<vcl_size_t>(A.size1() - row_start, max_block_size);

        viennacl::range     block_range(row_start, row_start + current_block_size);
        viennacl::range remainder_range(row_start + current_block_size, A.size1());

        //
        // Perform LU factorization on panel:
        //


        // Read from matrix to buffer:
        viennacl::backend::memory_read(A.handle(),
                                       sizeof(SCALARTYPE) * row_start          * A.internal_size2(),
                                       sizeof(SCALARTYPE) * current_block_size * A.internal_size2(),
                                       &(temp_buffer[0]));

        // Factorize (kij-version):
        for (vcl_size_t k=0; k < current_block_size - 1; ++k)
        {
          for (vcl_size_t i=k+1; i < current_block_size; ++i)
          {
            temp_buffer[row_start + i * A.internal_size2() + k] /= temp_buffer[row_start + k * A.internal_size2() + k];  // write l_ik

            SCALARTYPE l_ik = temp_buffer[row_start + i * A.internal_size2() + k];

            for (vcl_size_t j = row_start + k + 1; j < A.size1(); ++j)
              temp_buffer[i * A.internal_size2() + j] -= l_ik * temp_buffer[k * A.internal_size2() + j];  // l_ik * a_kj
          }
        }

        // Write back:
        viennacl::backend::memory_write(A.handle(),
                                        sizeof(SCALARTYPE) * row_start          * A.internal_size2(),
                                        sizeof(SCALARTYPE) * current_block_size * A.internal_size2(),
                                        &(temp_buffer[0]));

        if (remainder_range.size() > 0)
        {
          //
          // Compute L_12 = [ (U_11)^{T}^{-1} A_{21}^T ]^T
          //
          viennacl::matrix_range<MatrixType> U_11(A, block_range,     block_range);
          viennacl::matrix_range<MatrixType> A_21(A, remainder_range, block_range);
          viennacl::linalg::inplace_solve(trans(U_11), trans(A_21), viennacl::linalg::lower_tag());

          //
          // Update remainder of A
          //
          viennacl::matrix_range<MatrixType> L_21(A, remainder_range, block_range);
          viennacl::matrix_range<MatrixType> U_12(A, block_range,     remainder_range);
          viennacl::matrix_range<MatrixType> A_22(A, remainder_range, remainder_range);

          A_22 -= viennacl::linalg::prod(L_21, U_12);
        }
      }

    }


    /** @brief LU factorization of a column-major dense matrix.
    *
    * @param A    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    */
    template<typename SCALARTYPE>
    void lu_factorize(matrix<SCALARTYPE, viennacl::column_major> & A)
    {
      typedef matrix<SCALARTYPE, viennacl::column_major>  MatrixType;

      vcl_size_t max_block_size = 32;
      vcl_size_t num_blocks = (A.size1() - 1) / max_block_size + 1;
      std::vector<SCALARTYPE> temp_buffer(A.internal_size1() * max_block_size);

      // Iterate over panels
      for (vcl_size_t panel_id = 0; panel_id < num_blocks; ++panel_id)
      {
        vcl_size_t col_start = panel_id * max_block_size;
        vcl_size_t current_block_size = std::min<vcl_size_t>(A.size1() - col_start, max_block_size);

        viennacl::range     block_range(col_start, col_start + current_block_size);
        viennacl::range remainder_range(col_start + current_block_size, A.size1());

        //
        // Perform LU factorization on panel:
        //


        // Read from matrix to buffer:
        viennacl::backend::memory_read(A.handle(),
                                       sizeof(SCALARTYPE) * col_start          * A.internal_size1(),
                                       sizeof(SCALARTYPE) * current_block_size * A.internal_size1(),
                                       &(temp_buffer[0]));

        // Factorize (kji-version):
        for (vcl_size_t k=0; k < current_block_size; ++k)
        {
          SCALARTYPE a_kk = temp_buffer[col_start + k + k * A.internal_size1()];
          for (vcl_size_t i=col_start+k+1; i < A.size1(); ++i)
            temp_buffer[i + k * A.internal_size1()] /= a_kk;  // write l_ik

          for (vcl_size_t j=k+1; j < current_block_size; ++j)
          {
            SCALARTYPE a_kj = temp_buffer[col_start + k + j * A.internal_size1()];
            for (vcl_size_t i=col_start+k+1; i < A.size1(); ++i)
              temp_buffer[i + j * A.internal_size1()] -= temp_buffer[i + k * A.internal_size1()] * a_kj;  // l_ik * a_kj
          }
        }

        // Write back:
        viennacl::backend::memory_write(A.handle(),
                                        sizeof(SCALARTYPE) * col_start          * A.internal_size1(),
                                        sizeof(SCALARTYPE) * current_block_size * A.internal_size1(),
                                        &(temp_buffer[0]));

        if (remainder_range.size() > 0)
        {
          //
          // Compute U_12:
          //
          viennacl::matrix_range<MatrixType> L_11(A, block_range,     block_range);
          viennacl::matrix_range<MatrixType> A_12(A, block_range, remainder_range);
          viennacl::linalg::inplace_solve(L_11, A_12, viennacl::linalg::unit_lower_tag());

          //
          // Update remainder of A
          //
          viennacl::matrix_range<MatrixType> L_21(A, remainder_range, block_range);
          viennacl::matrix_range<MatrixType> U_12(A, block_range,     remainder_range);
          viennacl::matrix_range<MatrixType> A_22(A, remainder_range, remainder_range);

          A_22 -= viennacl::linalg::prod(L_21, U_12);
        }

      }

    }


    //
    // Convenience layer:
    //

    /** @brief LU substitution for the system LU = rhs.
    *
    * @param A    The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    * @param B    The matrix of load vectors, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F1, typename F2, unsigned int ALIGNMENT_A, unsigned int ALIGNMENT_B>
    void lu_substitute(matrix<SCALARTYPE, F1, ALIGNMENT_A> const & A,
                       matrix<SCALARTYPE, F2, ALIGNMENT_B> & B)
    {
      assert(A.size1() == A.size2() && bool("Matrix must be square"));
      assert(A.size1() == B.size1() && bool("Matrix must be square"));
      inplace_solve(A, B, unit_lower_tag());
      inplace_solve(A, B, upper_tag());
    }

    /** @brief LU substitution for the system LU = rhs.
    *
    * @param A      The system matrix, where the LU matrices are directly written to. The implicit unit diagonal of L is not written.
    * @param vec    The load vector, where the solution is directly written to
    */
    template<typename SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void lu_substitute(matrix<SCALARTYPE, F, ALIGNMENT> const & A,
                       vector<SCALARTYPE, VEC_ALIGNMENT> & vec)
    {
      assert(A.size1() == A.size2() && bool("Matrix must be square"));
      inplace_solve(A, vec, unit_lower_tag());
      inplace_solve(A, vec, upper_tag());
    }

  }
}

#endif
