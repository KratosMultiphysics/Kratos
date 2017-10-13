#ifndef VIENNACL_LINALG_HOST_BASED_SPARSE_MATRIX_OPERATIONS_HPP_
#define VIENNACL_LINALG_HOST_BASED_SPARSE_MATRIX_OPERATIONS_HPP_

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

/** @file viennacl/linalg/host_based/sparse_matrix_operations.hpp
    @brief Implementations of operations using sparse matrices on the CPU using a single thread or OpenMP.
*/

#include <list>

#include "viennacl/forwards.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/host_based/common.hpp"
#include "viennacl/linalg/host_based/vector_operations.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace host_based
    {
      //
      // Compressed matrix
      //

      namespace detail
      {
        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void row_info(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & mat,
                      vector_base<ScalarType> & vec,
                      viennacl::linalg::detail::row_info_types info_selector)
        {
          ScalarType         * result_buf = detail::extract_raw_pointer<ScalarType>(vec.handle());
          ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(mat.handle());
          unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle1());
          unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle2());

          for (vcl_size_t row = 0; row < mat.size1(); ++row)
          {
            ScalarType value = 0;
            unsigned int row_end = row_buffer[row+1];

            switch (info_selector)
            {
              case viennacl::linalg::detail::SPARSE_ROW_NORM_INF: //inf-norm
                for (unsigned int i = row_buffer[row]; i < row_end; ++i)
                  value = std::max<ScalarType>(value, std::fabs(elements[i]));
                break;

              case viennacl::linalg::detail::SPARSE_ROW_NORM_1: //1-norm
                for (unsigned int i = row_buffer[row]; i < row_end; ++i)
                  value += std::fabs(elements[i]);
                break;

              case viennacl::linalg::detail::SPARSE_ROW_NORM_2: //2-norm
                for (unsigned int i = row_buffer[row]; i < row_end; ++i)
                  value += elements[i] * elements[i];
                value = std::sqrt(value);
                break;

              case viennacl::linalg::detail::SPARSE_ROW_DIAGONAL: //diagonal entry
                for (unsigned int i = row_buffer[row]; i < row_end; ++i)
                {
                  if (col_buffer[i] == row)
                  {
                    value = elements[i];
                    break;
                  }
                }
                break;

              default:
                break;
            }
            result_buf[row] = value;
          }
        }
      }


      /** @brief Carries out matrix-vector multiplication with a compressed_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
      template<class ScalarType, unsigned int ALIGNMENT>
      void prod_impl(const viennacl::compressed_matrix<ScalarType, ALIGNMENT> & mat,
                     const viennacl::vector_base<ScalarType> & vec,
                           viennacl::vector_base<ScalarType> & result)
      {
        ScalarType         * result_buf = detail::extract_raw_pointer<ScalarType>(result.handle());
        ScalarType   const * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(mat.handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle2());

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
        for (long row = 0; row < static_cast<long>(mat.size1()); ++row)
        {
          ScalarType dot_prod = 0;
          vcl_size_t row_end = row_buffer[row+1];
          for (vcl_size_t i = row_buffer[row]; i < row_end; ++i)
            dot_prod += elements[i] * vec_buf[col_buffer[i] * vec.stride() + vec.start()];
          result_buf[row * result.stride() + result.start()] = dot_prod;
        }

      }

      /** @brief Carries out sparse_matrix-matrix multiplication first matrix being compressed
      *
      * Implementation of the convenience expression result = prod(sp_mat, d_mat);
      *
      * @param sp_mat     The sparse matrix
      * @param d_mat      The dense matrix
      * @param result     The result matrix
      */
      template< class ScalarType, typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::compressed_matrix<ScalarType, ALIGNMENT> & sp_mat,
                     const viennacl::matrix_base<NumericT, F1> & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result) {

        ScalarType   const * sp_mat_elements   = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_row_buffer = detail::extract_raw_pointer<unsigned int>(sp_mat.handle1());
        unsigned int const * sp_mat_col_buffer = detail::extract_raw_pointer<unsigned int>(sp_mat.handle2());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat);
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat);
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat);
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat);
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat);
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat);
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat);

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row) {
            vcl_size_t row_start = sp_mat_row_buffer[row];
            vcl_size_t row_end = sp_mat_row_buffer[row+1];
            for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {
              NumericT temp = 0;
              for (vcl_size_t k = row_start; k < row_end; ++k) {
                temp += sp_mat_elements[k] * d_mat_wrapper(sp_mat_col_buffer[k], col);
              }
              result_wrapper(row, col) = temp;
            }
          }
        }
        else {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col) {
            for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row) {
              vcl_size_t row_start = sp_mat_row_buffer[row];
              vcl_size_t row_end = sp_mat_row_buffer[row+1];
              NumericT temp = 0;
              for (vcl_size_t k = row_start; k < row_end; ++k) {
                temp += sp_mat_elements[k] * d_mat_wrapper(sp_mat_col_buffer[k], col);
              }
              result_wrapper(row, col) = temp;
            }
          }
        }

      }

      /** @brief Carries out matrix-trans(matrix) multiplication first matrix being compressed
      *          and the second transposed
      *
      * Implementation of the convenience expression result = prod(sp_mat, trans(d_mat));
      *
      * @param sp_mat             The sparse matrix
      * @param d_mat              The transposed dense matrix
      * @param result             The result matrix
      */
      template< class ScalarType, typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::compressed_matrix<ScalarType, ALIGNMENT> & sp_mat,
                const viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                                   const viennacl::matrix_base<NumericT, F1>,
                                                   viennacl::op_trans > & d_mat,
                      viennacl::matrix_base<NumericT, F2> & result) {

        ScalarType   const * sp_mat_elements   = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_row_buffer = detail::extract_raw_pointer<unsigned int>(sp_mat.handle1());
        unsigned int const * sp_mat_col_buffer = detail::extract_raw_pointer<unsigned int>(sp_mat.handle2());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat.lhs());
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat.lhs());
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat.lhs());
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat.lhs());
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat.lhs());
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat.lhs());
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat.lhs());

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row) {
            vcl_size_t row_start = sp_mat_row_buffer[row];
            vcl_size_t row_end = sp_mat_row_buffer[row+1];
            for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {
              NumericT temp = 0;
              for (vcl_size_t k = row_start; k < row_end; ++k) {
                temp += sp_mat_elements[k] * d_mat_wrapper(col, sp_mat_col_buffer[k]);
              }
              result_wrapper(row, col) = temp;
            }
          }
        }
        else {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col) {
            for (vcl_size_t row = 0; row < sp_mat.size1(); ++row) {
              vcl_size_t row_start = sp_mat_row_buffer[row];
              vcl_size_t row_end = sp_mat_row_buffer[row+1];
              NumericT temp = 0;
              for (vcl_size_t k = row_start; k < row_end; ++k) {
                temp += sp_mat_elements[k] * d_mat_wrapper(col, sp_mat_col_buffer[k]);
              }
              result_wrapper(row, col) = temp;
            }
          }
        }

      }


      //
      // Triangular solve for compressed_matrix, A \ b
      //
      namespace detail
      {
        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_inplace_solve(SizeTypeArray const & row_buffer,
                               SizeTypeArray const & col_buffer,
                               ConstScalarTypeArray const & element_buffer,
                               ScalarTypeArray & vec_buffer,
                               vcl_size_t num_cols,
                               viennacl::linalg::unit_lower_tag)
        {
          vcl_size_t row_begin = row_buffer[1];
          for (vcl_size_t row = 1; row < num_cols; ++row)
          {
            NumericT vec_entry = vec_buffer[row];
            vcl_size_t row_end = row_buffer[row+1];
            for (vcl_size_t i = row_begin; i < row_end; ++i)
            {
              vcl_size_t col_index = col_buffer[i];
              if (col_index < row)
                vec_entry -= vec_buffer[col_index] * element_buffer[i];
            }
            vec_buffer[row] = vec_entry;
            row_begin = row_end;
          }
        }

        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_inplace_solve(SizeTypeArray const & row_buffer,
                               SizeTypeArray const & col_buffer,
                               ConstScalarTypeArray const & element_buffer,
                               ScalarTypeArray & vec_buffer,
                               vcl_size_t num_cols,
                               viennacl::linalg::lower_tag)
        {
          vcl_size_t row_begin = row_buffer[0];
          for (vcl_size_t row = 0; row < num_cols; ++row)
          {
            NumericT vec_entry = vec_buffer[row];

            // substitute and remember diagonal entry
            vcl_size_t row_end = row_buffer[row+1];
            NumericT diagonal_entry = 0;
            for (vcl_size_t i = row_begin; i < row_end; ++i)
            {
              vcl_size_t col_index = col_buffer[i];
              if (col_index < row)
                vec_entry -= vec_buffer[col_index] * element_buffer[i];
              else if (col_index == row)
                diagonal_entry = element_buffer[i];
            }

            vec_buffer[row] = vec_entry / diagonal_entry;
            row_begin = row_end;
          }
        }


        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_inplace_solve(SizeTypeArray const & row_buffer,
                               SizeTypeArray const & col_buffer,
                               ConstScalarTypeArray const & element_buffer,
                               ScalarTypeArray & vec_buffer,
                               vcl_size_t num_cols,
                               viennacl::linalg::unit_upper_tag)
        {
          for (vcl_size_t row2 = 1; row2 < num_cols; ++row2)
          {
            vcl_size_t row = (num_cols - row2) - 1;
            NumericT vec_entry = vec_buffer[row];
            vcl_size_t row_begin = row_buffer[row];
            vcl_size_t row_end   = row_buffer[row+1];
            for (vcl_size_t i = row_begin; i < row_end; ++i)
            {
              vcl_size_t col_index = col_buffer[i];
              if (col_index > row)
                vec_entry -= vec_buffer[col_index] * element_buffer[i];
            }
            vec_buffer[row] = vec_entry;
          }
        }

        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_inplace_solve(SizeTypeArray const & row_buffer,
                               SizeTypeArray const & col_buffer,
                               ConstScalarTypeArray const & element_buffer,
                               ScalarTypeArray & vec_buffer,
                               vcl_size_t num_cols,
                               viennacl::linalg::upper_tag)
        {
          for (vcl_size_t row2 = 0; row2 < num_cols; ++row2)
          {
            vcl_size_t row = (num_cols - row2) - 1;
            NumericT vec_entry = vec_buffer[row];

            // substitute and remember diagonal entry
            vcl_size_t row_begin = row_buffer[row];
            vcl_size_t row_end   = row_buffer[row+1];
            NumericT diagonal_entry = 0;
            for (vcl_size_t i = row_begin; i < row_end; ++i)
            {
              vcl_size_t col_index = col_buffer[i];
              if (col_index > row)
                vec_entry -= vec_buffer[col_index] * element_buffer[i];
              else if (col_index == row)
                diagonal_entry = element_buffer[i];
            }

            vec_buffer[row] = vec_entry / diagonal_entry;
          }
        }

      } //namespace detail



      /** @brief Inplace solution of a lower triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
      *
      * @param L    The matrix
      * @param vec  The vector holding the right hand side. Is overwritten by the solution.
      * @param tag  The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & L,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::unit_lower_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(L.handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(L.handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(L.handle2());

        detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, L.size2(), tag);
      }

      /** @brief Inplace solution of a lower triangular compressed_matrix. Typically used for LU substitutions
      *
      * @param L    The matrix
      * @param vec  The vector holding the right hand side. Is overwritten by the solution.
      * @param tag  The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & L,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::lower_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(L.handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(L.handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(L.handle2());

        detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, L.size2(), tag);
      }


      /** @brief Inplace solution of a upper triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
      *
      * @param U    The matrix
      * @param vec  The vector holding the right hand side. Is overwritten by the solution.
      * @param tag  The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & U,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::unit_upper_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(U.handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(U.handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(U.handle2());

        detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, U.size2(), tag);
      }

      /** @brief Inplace solution of a upper triangular compressed_matrix. Typically used for LU substitutions
      *
      * @param U    The matrix
      * @param vec  The vector holding the right hand side. Is overwritten by the solution.
      * @param tag  The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & U,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::upper_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(U.handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(U.handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(U.handle2());

        detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, U.size2(), tag);
      }







      //
      // Triangular solve for compressed_matrix, A^T \ b
      //

      namespace detail
      {
        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_trans_inplace_solve(SizeTypeArray const & row_buffer,
                                     SizeTypeArray const & col_buffer,
                                     ConstScalarTypeArray const & element_buffer,
                                     ScalarTypeArray & vec_buffer,
                                     vcl_size_t num_cols,
                                     viennacl::linalg::unit_lower_tag)
        {
          vcl_size_t col_begin = row_buffer[0];
          for (vcl_size_t col = 0; col < num_cols; ++col)
          {
            NumericT vec_entry = vec_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              unsigned int row_index = col_buffer[i];
              if (row_index > col)
                vec_buffer[row_index] -= vec_entry * element_buffer[i];
            }
            col_begin = col_end;
          }
        }

        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_trans_inplace_solve(SizeTypeArray const & row_buffer,
                                     SizeTypeArray const & col_buffer,
                                     ConstScalarTypeArray const & element_buffer,
                                     ScalarTypeArray & vec_buffer,
                                     vcl_size_t num_cols,
                                     viennacl::linalg::lower_tag)
        {
          vcl_size_t col_begin = row_buffer[0];
          for (vcl_size_t col = 0; col < num_cols; ++col)
          {
            vcl_size_t col_end = row_buffer[col+1];

            // Stage 1: Find diagonal entry:
            NumericT diagonal_entry = 0;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index == col)
              {
                diagonal_entry = element_buffer[i];
                break;
              }
            }

            // Stage 2: Substitute
            NumericT vec_entry = vec_buffer[col] / diagonal_entry;
            vec_buffer[col] = vec_entry;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index > col)
                vec_buffer[row_index] -= vec_entry * element_buffer[i];
            }
            col_begin = col_end;
          }
        }

        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_trans_inplace_solve(SizeTypeArray const & row_buffer,
                                     SizeTypeArray const & col_buffer,
                                     ConstScalarTypeArray const & element_buffer,
                                     ScalarTypeArray & vec_buffer,
                                     vcl_size_t num_cols,
                                     viennacl::linalg::unit_upper_tag)
        {
          for (vcl_size_t col2 = 0; col2 < num_cols; ++col2)
          {
            vcl_size_t col = (num_cols - col2) - 1;

            NumericT vec_entry = vec_buffer[col];
            vcl_size_t col_begin = row_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index < col)
                vec_buffer[row_index] -= vec_entry * element_buffer[i];
            }

          }
        }

        template <typename NumericT, typename ConstScalarTypeArray, typename ScalarTypeArray, typename SizeTypeArray>
        void csr_trans_inplace_solve(SizeTypeArray const & row_buffer,
                                     SizeTypeArray const & col_buffer,
                                     ConstScalarTypeArray const & element_buffer,
                                     ScalarTypeArray & vec_buffer,
                                     vcl_size_t num_cols,
                                     viennacl::linalg::upper_tag)
        {
          for (vcl_size_t col2 = 0; col2 < num_cols; ++col2)
          {
            vcl_size_t col = (num_cols - col2) - 1;
            vcl_size_t col_begin = row_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];

            // Stage 1: Find diagonal entry:
            NumericT diagonal_entry = 0;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index == col)
              {
                diagonal_entry = element_buffer[i];
                break;
              }
            }

            // Stage 2: Substitute
            NumericT vec_entry = vec_buffer[col] / diagonal_entry;
            vec_buffer[col] = vec_entry;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index < col)
                vec_buffer[row_index] -= vec_entry * element_buffer[i];
            }
          }
        }


        //
        // block solves
        //
        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void block_inplace_solve(const matrix_expression<const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         op_trans> & L,
                                 viennacl::backend::mem_handle const & /* block_indices */, vcl_size_t /* num_blocks */,
                                 vector_base<ScalarType> const & /* L_diagonal */,  //ignored
                                 vector_base<ScalarType> & vec,
                                 viennacl::linalg::unit_lower_tag)
        {
          // Note: The following could be implemented more efficiently using the block structure and possibly OpenMP.

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(L.lhs().handle1());
          unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(L.lhs().handle2());
          ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(L.lhs().handle());
          ScalarType         * vec_buffer = detail::extract_raw_pointer<ScalarType>(vec.handle());

          vcl_size_t col_begin = row_buffer[0];
          for (vcl_size_t col = 0; col < L.lhs().size1(); ++col)
          {
            ScalarType vec_entry = vec_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              unsigned int row_index = col_buffer[i];
              if (row_index > col)
                vec_buffer[row_index] -= vec_entry * elements[i];
            }
            col_begin = col_end;
          }
        }

        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void block_inplace_solve(const matrix_expression<const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         op_trans> & L,
                                 viennacl::backend::mem_handle const & /*block_indices*/, vcl_size_t /* num_blocks */,
                                 vector_base<ScalarType> const & L_diagonal,
                                 vector_base<ScalarType> & vec,
                                 viennacl::linalg::lower_tag)
        {
          // Note: The following could be implemented more efficiently using the block structure and possibly OpenMP.

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(L.lhs().handle1());
          unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(L.lhs().handle2());
          ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(L.lhs().handle());
          ScalarType   const * diagonal_buffer = detail::extract_raw_pointer<ScalarType>(L_diagonal.handle());
          ScalarType         * vec_buffer = detail::extract_raw_pointer<ScalarType>(vec.handle());

          vcl_size_t col_begin = row_buffer[0];
          for (vcl_size_t col = 0; col < L.lhs().size1(); ++col)
          {
            vcl_size_t col_end = row_buffer[col+1];

            ScalarType vec_entry = vec_buffer[col] / diagonal_buffer[col];
            vec_buffer[col] = vec_entry;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index > col)
                vec_buffer[row_index] -= vec_entry * elements[i];
            }
            col_begin = col_end;
          }
        }



        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void block_inplace_solve(const matrix_expression<const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         op_trans> & U,
                                 viennacl::backend::mem_handle const & /*block_indices*/, vcl_size_t /* num_blocks */,
                                 vector_base<ScalarType> const & /* U_diagonal */, //ignored
                                 vector_base<ScalarType> & vec,
                                 viennacl::linalg::unit_upper_tag)
        {
          // Note: The following could be implemented more efficiently using the block structure and possibly OpenMP.

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(U.lhs().handle1());
          unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(U.lhs().handle2());
          ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(U.lhs().handle());
          ScalarType         * vec_buffer = detail::extract_raw_pointer<ScalarType>(vec.handle());

          for (vcl_size_t col2 = 0; col2 < U.lhs().size1(); ++col2)
          {
            vcl_size_t col = (U.lhs().size1() - col2) - 1;

            ScalarType vec_entry = vec_buffer[col];
            vcl_size_t col_begin = row_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index < col)
                vec_buffer[row_index] -= vec_entry * elements[i];
            }

          }
        }

        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void block_inplace_solve(const matrix_expression<const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                                         op_trans> & U,
                                 viennacl::backend::mem_handle const & /* block_indices */, vcl_size_t /* num_blocks */,
                                 vector_base<ScalarType> const & U_diagonal,
                                 vector_base<ScalarType> & vec,
                                 viennacl::linalg::upper_tag)
        {
          // Note: The following could be implemented more efficiently using the block structure and possibly OpenMP.

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(U.lhs().handle1());
          unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(U.lhs().handle2());
          ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(U.lhs().handle());
          ScalarType   const * diagonal_buffer = detail::extract_raw_pointer<ScalarType>(U_diagonal.handle());
          ScalarType         * vec_buffer = detail::extract_raw_pointer<ScalarType>(vec.handle());

          for (vcl_size_t col2 = 0; col2 < U.lhs().size1(); ++col2)
          {
            vcl_size_t col = (U.lhs().size1() - col2) - 1;
            vcl_size_t col_begin = row_buffer[col];
            vcl_size_t col_end = row_buffer[col+1];

            // Stage 2: Substitute
            ScalarType vec_entry = vec_buffer[col] / diagonal_buffer[col];
            vec_buffer[col] = vec_entry;
            for (vcl_size_t i = col_begin; i < col_end; ++i)
            {
              vcl_size_t row_index = col_buffer[i];
              if (row_index < col)
                vec_buffer[row_index] -= vec_entry * elements[i];
            }
          }
        }


      } //namespace detail

      /** @brief Inplace solution of a lower triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
      *
      * @param proxy  Proxy object for a transposed CSR-matrix
      * @param vec    The right hand side vector
      * @param tag    The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(matrix_expression< const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            op_trans> const & proxy,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::unit_lower_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(proxy.lhs().handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle2());

        detail::csr_trans_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, proxy.lhs().size1(), tag);
      }

      /** @brief Inplace solution of a lower triangular compressed_matrix. Typically used for LU substitutions
      *
      * @param proxy  Proxy object for a transposed CSR-matrix
      * @param vec    The right hand side vector
      * @param tag    The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(matrix_expression< const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            op_trans> const & proxy,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::lower_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(proxy.lhs().handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle2());

        detail::csr_trans_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, proxy.lhs().size1(), tag);
      }


      /** @brief Inplace solution of a upper triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
      *
      * @param proxy  Proxy object for a transposed CSR-matrix
      * @param vec    The right hand side vector
      * @param tag    The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(matrix_expression< const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            op_trans> const & proxy,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::unit_upper_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(proxy.lhs().handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle2());

        detail::csr_trans_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, proxy.lhs().size1(), tag);
      }


      /** @brief Inplace solution of a upper triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
      *
      * @param proxy  Proxy object for a transposed CSR-matrix
      * @param vec    The right hand side vector
      * @param tag    The solver tag identifying the respective triangular solver
      */
      template<typename ScalarType, unsigned int MAT_ALIGNMENT>
      void inplace_solve(matrix_expression< const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            const compressed_matrix<ScalarType, MAT_ALIGNMENT>,
                                            op_trans> const & proxy,
                         vector_base<ScalarType> & vec,
                         viennacl::linalg::upper_tag tag)
      {
        ScalarType         * vec_buf    = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements   = detail::extract_raw_pointer<ScalarType>(proxy.lhs().handle());
        unsigned int const * row_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle1());
        unsigned int const * col_buffer = detail::extract_raw_pointer<unsigned int>(proxy.lhs().handle2());

        detail::csr_trans_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_buf, proxy.lhs().size1(), tag);
      }



      //
      // Compressed Compressed Matrix
      //

      /** @brief Carries out matrix-vector multiplication with a compressed_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
      template<class ScalarType>
      void prod_impl(const viennacl::compressed_compressed_matrix<ScalarType> & mat,
                     const viennacl::vector_base<ScalarType> & vec,
                           viennacl::vector_base<ScalarType> & result)
      {
        ScalarType         * result_buf  = detail::extract_raw_pointer<ScalarType>(result.handle());
        ScalarType   const * vec_buf     = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements    = detail::extract_raw_pointer<ScalarType>(mat.handle());
        unsigned int const * row_buffer  = detail::extract_raw_pointer<unsigned int>(mat.handle1());
        unsigned int const * row_indices = detail::extract_raw_pointer<unsigned int>(mat.handle3());
        unsigned int const * col_buffer  = detail::extract_raw_pointer<unsigned int>(mat.handle2());

        vector_assign(result, ScalarType(0));

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
        for (long i = 0; i < static_cast<long>(mat.nnz1()); ++i)
        {
          ScalarType dot_prod = 0;
          vcl_size_t row_end = row_buffer[i+1];
          for (vcl_size_t j = row_buffer[i]; j < row_end; ++j)
            dot_prod += elements[j] * vec_buf[col_buffer[j] * vec.stride() + vec.start()];
          result_buf[row_indices[i] * result.stride() + result.start()] = dot_prod;
        }

      }



      //
      // Coordinate Matrix
      //

      namespace detail
      {
        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void row_info(coordinate_matrix<ScalarType, MAT_ALIGNMENT> const & mat,
                      vector_base<ScalarType> & vec,
                      viennacl::linalg::detail::row_info_types info_selector)
        {
          ScalarType         * result_buf   = detail::extract_raw_pointer<ScalarType>(vec.handle());
          ScalarType   const * elements     = detail::extract_raw_pointer<ScalarType>(mat.handle());
          unsigned int const * coord_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle12());

          ScalarType value = 0;
          unsigned int last_row = 0;

          for (vcl_size_t i = 0; i < mat.nnz(); ++i)
          {
            unsigned int current_row = coord_buffer[2*i];

            if (current_row != last_row)
            {
              if (info_selector == viennacl::linalg::detail::SPARSE_ROW_NORM_2)
                value = std::sqrt(value);

              result_buf[last_row] = value;
              value = 0;
              last_row = current_row;
            }

            switch (info_selector)
            {
              case viennacl::linalg::detail::SPARSE_ROW_NORM_INF: //inf-norm
                value = std::max<ScalarType>(value, std::fabs(elements[i]));
                break;

              case viennacl::linalg::detail::SPARSE_ROW_NORM_1: //1-norm
                value += std::fabs(elements[i]);
                break;

              case viennacl::linalg::detail::SPARSE_ROW_NORM_2: //2-norm
                value += elements[i] * elements[i];
                break;

              case viennacl::linalg::detail::SPARSE_ROW_DIAGONAL: //diagonal entry
                if (coord_buffer[2*i+1] == current_row)
                  value = elements[i];
                break;

              default:
                break;
            }
          }

          if (info_selector == viennacl::linalg::detail::SPARSE_ROW_NORM_2)
            value = std::sqrt(value);

          result_buf[last_row] = value;
        }
      }

      /** @brief Carries out matrix-vector multiplication with a coordinate_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
      template<class ScalarType, unsigned int ALIGNMENT>
      void prod_impl(const viennacl::coordinate_matrix<ScalarType, ALIGNMENT> & mat,
                     const viennacl::vector_base<ScalarType> & vec,
                           viennacl::vector_base<ScalarType> & result)
      {
        ScalarType         * result_buf   = detail::extract_raw_pointer<ScalarType>(result.handle());
        ScalarType   const * vec_buf      = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements     = detail::extract_raw_pointer<ScalarType>(mat.handle());
        unsigned int const * coord_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle12());

        for (vcl_size_t i = 0; i< result.size(); ++i)
          result_buf[i * result.stride() + result.start()] = 0;

        for (vcl_size_t i = 0; i < mat.nnz(); ++i)
          result_buf[coord_buffer[2*i] * result.stride() + result.start()]
            += elements[i] * vec_buf[coord_buffer[2*i+1] * vec.stride() + vec.start()];
      }

      /** @brief Carries out Compressed Matrix(COO)-Dense Matrix multiplication
      *
      * Implementation of the convenience expression result = prod(sp_mat, d_mat);
      *
      * @param sp_mat     The Sparse Matrix (Coordinate format)
      * @param d_mat      The Dense Matrix
      * @param result     The Result Matrix
      */
      template<class ScalarType, unsigned int ALIGNMENT, class NumericT, typename F1, typename F2>
      void prod_impl(const viennacl::coordinate_matrix<ScalarType, ALIGNMENT> & sp_mat,
                     const viennacl::matrix_base<NumericT, F1> & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result) {

        ScalarType   const * sp_mat_elements     = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_coords       = detail::extract_raw_pointer<unsigned int>(sp_mat.handle12());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat);
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat);
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat);
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat);
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat);
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat);
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat);

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
                for (vcl_size_t col = 0; col < d_mat.size2(); ++col)
                  result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long i = 0; i < static_cast<long>(sp_mat.nnz()); ++i) {
            NumericT x = static_cast<NumericT>(sp_mat_elements[i]);
            unsigned int r = sp_mat_coords[2*i];
            unsigned int c = sp_mat_coords[2*i+1];
            for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {
              NumericT y = d_mat_wrapper( c, col);
              result_wrapper(r, col) += x * y;
            }
          }
        }

        else {

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col)
            for (vcl_size_t row = 0; row < sp_mat.size1(); ++row)
                result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col) {

            for (vcl_size_t i = 0; i < sp_mat.nnz(); ++i) {

              NumericT x = static_cast<NumericT>(sp_mat_elements[i]);
              unsigned int r = sp_mat_coords[2*i];
              unsigned int c = sp_mat_coords[2*i+1];
              NumericT y = d_mat_wrapper( c, col);

              result_wrapper( r, col) += x*y;
            }

          }
        }

      }


      /** @brief Carries out Compressed Matrix(COO)-Dense Transposed Matrix multiplication
      *
      * Implementation of the convenience expression result = prod(sp_mat, trans(d_mat));
      *
      * @param sp_mat     The Sparse Matrix (Coordinate format)
      * @param d_mat      The Dense Transposed Matrix
      * @param result     The Result Matrix
      */
      template<class ScalarType, unsigned int ALIGNMENT, class NumericT, typename F1, typename F2>
      void prod_impl(const viennacl::coordinate_matrix<ScalarType, ALIGNMENT> & sp_mat,
                     const viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                                        const viennacl::matrix_base<NumericT, F1>,
                                                        viennacl::op_trans > & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result) {

        ScalarType   const * sp_mat_elements     = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_coords       = detail::extract_raw_pointer<unsigned int>(sp_mat.handle12());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat.lhs());
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat.lhs());
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat.lhs());
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat.lhs());
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat.lhs());
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat.lhs());
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat.lhs());

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
            for (vcl_size_t col = 0; col < d_mat.size2(); ++col)
              result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */
        }
        else {

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col)
            for (vcl_size_t row = 0; row < sp_mat.size1(); ++row)
              result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */
        }

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
        for (long i = 0; i < static_cast<long>(sp_mat.nnz()); ++i) {
          NumericT x = static_cast<NumericT>(sp_mat_elements[i]);
          unsigned int r = sp_mat_coords[2*i];
          unsigned int c = sp_mat_coords[2*i+1];
          for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {
            NumericT y = d_mat_wrapper( col, c);
            result_wrapper(r, col) += x * y;
          }
        }

      }
      //
      // ELL Matrix
      //
      /** @brief Carries out matrix-vector multiplication with a ell_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
      template<class ScalarType, unsigned int ALIGNMENT>
      void prod_impl(const viennacl::ell_matrix<ScalarType, ALIGNMENT> & mat,
                     const viennacl::vector_base<ScalarType> & vec,
                           viennacl::vector_base<ScalarType> & result)
      {
        ScalarType         * result_buf   = detail::extract_raw_pointer<ScalarType>(result.handle());
        ScalarType   const * vec_buf      = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements     = detail::extract_raw_pointer<ScalarType>(mat.handle());
        unsigned int const * coords       = detail::extract_raw_pointer<unsigned int>(mat.handle2());

        for(vcl_size_t row = 0; row < mat.size1(); ++row)
        {
          ScalarType sum = 0;

          for(unsigned int item_id = 0; item_id < mat.internal_maxnnz(); ++item_id)
          {
            vcl_size_t offset = row + item_id * mat.internal_size1();
            ScalarType val = elements[offset];

            if(val != 0)
            {
              unsigned int col = coords[offset];
              sum += (vec_buf[col * vec.stride() + vec.start()] * val);
            }
          }

          result_buf[row * result.stride() + result.start()] = sum;
        }
      }

      /** @brief Carries out ell_matrix-d_matrix multiplication
      *
      * Implementation of the convenience expression result = prod(sp_mat, d_mat);
      *
      * @param sp_mat     The sparse(ELL) matrix
      * @param d_mat      The dense matrix
      * @param result     The result dense matrix
      */
      template<class ScalarType, typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::ell_matrix<ScalarType, ALIGNMENT> & sp_mat,
                     const viennacl::matrix_base<NumericT, F1> & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result)
      {
        ScalarType   const * sp_mat_elements     = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_coords       = detail::extract_raw_pointer<unsigned int>(sp_mat.handle2());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat);
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat);
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat);
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat);
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat);
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat);
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat);

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
                for (vcl_size_t col = 0; col < d_mat.size2(); ++col)
                  result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
          {
            for (long item_id = 0; item_id < static_cast<long>(sp_mat.maxnnz()); ++item_id)
            {
              vcl_size_t offset = static_cast<vcl_size_t>(row) + item_id * sp_mat.internal_size1();
              NumericT sp_mat_val = static_cast<NumericT>(sp_mat_elements[offset]);
              unsigned int sp_mat_col = sp_mat_coords[offset];

              if( sp_mat_val != 0)
              {
                for (vcl_size_t col = 0; col < d_mat.size2(); ++col)
                  result_wrapper(static_cast<vcl_size_t>(row), col) += sp_mat_val * d_mat_wrapper( sp_mat_col, col);
              }
            }
          }
        }
        else {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col)
            for (long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
                result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col) {

            for(unsigned int item_id = 0; item_id < sp_mat.maxnnz(); ++item_id) {

              for(vcl_size_t row = 0; row < sp_mat.size1(); ++row) {

                vcl_size_t offset = row + item_id * sp_mat.internal_size1();
                NumericT sp_mat_val = static_cast<NumericT>(sp_mat_elements[offset]);
                unsigned int sp_mat_col = sp_mat_coords[offset];

                if( sp_mat_val != 0) {

                  result_wrapper( row, col) += sp_mat_val * d_mat_wrapper( sp_mat_col, col);
                }
              }
            }
          }
        }

      }

      /** @brief Carries out matrix-trans(matrix) multiplication first matrix being sparse ell
      *          and the second dense transposed
      *
      * Implementation of the convenience expression result = prod(sp_mat, trans(d_mat));
      *
      * @param sp_mat             The sparse matrix
      * @param d_mat              The transposed dense matrix
      * @param result             The result matrix
      */
      template<class ScalarType, typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::ell_matrix<ScalarType, ALIGNMENT> & sp_mat,
                     const viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                                        const viennacl::matrix_base<NumericT, F1>,
                                                        viennacl::op_trans > & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result) {

        ScalarType   const * sp_mat_elements     = detail::extract_raw_pointer<ScalarType>(sp_mat.handle());
        unsigned int const * sp_mat_coords       = detail::extract_raw_pointer<unsigned int>(sp_mat.handle2());

        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat.lhs());
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat.lhs());
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat.lhs());
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat.lhs());
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat.lhs());
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat.lhs());
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat.lhs());

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        if ( detail::is_row_major(typename F1::orientation_category()) ) {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for(long row = 0; row < static_cast<long>(sp_mat.size1()); ++row)
            for (vcl_size_t col = 0; col < d_mat.size2(); ++col)
              result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

          for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {

            for(unsigned int item_id = 0; item_id < sp_mat.maxnnz(); ++item_id) {

              for(vcl_size_t row = 0; row < sp_mat.size1(); ++row) {

                vcl_size_t offset = row + item_id * sp_mat.internal_size1();
                NumericT sp_mat_val = static_cast<NumericT>(sp_mat_elements[offset]);
                unsigned int sp_mat_col = sp_mat_coords[offset];

                if( sp_mat_val != 0) {

                  result_wrapper( row, col) += sp_mat_val * d_mat_wrapper( col, sp_mat_col);
                }
              }
            }
          }
        }
        else {
#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for (long col = 0; col < static_cast<long>(d_mat.size2()); ++col)
            for (vcl_size_t row = 0; row < sp_mat.size1(); ++row)
              result_wrapper( row, col) = (NumericT)0; /* filling result with zeros, as the product loops are reordered */

#ifdef VIENNACL_WITH_OPENMP
        #pragma omp parallel for
#endif
          for(long item_id = 0; item_id < static_cast<long>(sp_mat.maxnnz()); ++item_id) {

            for(vcl_size_t row = 0; row < sp_mat.size1(); ++row) {

              vcl_size_t offset = row + item_id * sp_mat.internal_size1();
              NumericT sp_mat_val = static_cast<NumericT>(sp_mat_elements[offset]);
              unsigned int sp_mat_col = sp_mat_coords[offset];

              if( sp_mat_val != 0) {

                for (vcl_size_t col = 0; col < d_mat.size2(); ++col) {

                  result_wrapper( row, col) += sp_mat_val * d_mat_wrapper( col, sp_mat_col);
                }
              }
            }
          }
        }

      }

      //
      // Hybrid Matrix
      //
      /** @brief Carries out matrix-vector multiplication with a hyb_matrix
      *
      * Implementation of the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      * @param result The result vector
      */
      template<class ScalarType, unsigned int ALIGNMENT>
      void prod_impl(const viennacl::hyb_matrix<ScalarType, ALIGNMENT> & mat,
                     const viennacl::vector_base<ScalarType> & vec,
                           viennacl::vector_base<ScalarType> & result)
      {
        ScalarType         * result_buf     = detail::extract_raw_pointer<ScalarType>(result.handle());
        ScalarType   const * vec_buf        = detail::extract_raw_pointer<ScalarType>(vec.handle());
        ScalarType   const * elements       = detail::extract_raw_pointer<ScalarType>(mat.handle());
        unsigned int const * coords         = detail::extract_raw_pointer<unsigned int>(mat.handle2());
        ScalarType   const * csr_elements   = detail::extract_raw_pointer<ScalarType>(mat.handle5());
        unsigned int const * csr_row_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle3());
        unsigned int const * csr_col_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle4());


        for(vcl_size_t row = 0; row < mat.size1(); ++row)
        {
          ScalarType sum = 0;

          //
          // Part 1: Process ELL part
          //
          for(unsigned int item_id = 0; item_id < mat.internal_ellnnz(); ++item_id)
          {
            vcl_size_t offset = row + item_id * mat.internal_size1();
            ScalarType val = elements[offset];

            if(val != 0)
            {
              unsigned int col = coords[offset];
              sum += (vec_buf[col * vec.stride() + vec.start()] * val);
            }
          }

          //
          // Part 2: Process HYB part
          //
          vcl_size_t col_begin = csr_row_buffer[row];
          vcl_size_t col_end   = csr_row_buffer[row + 1];

          for(vcl_size_t item_id = col_begin; item_id < col_end; item_id++)
          {
              sum += (vec_buf[csr_col_buffer[item_id] * vec.stride() + vec.start()] * csr_elements[item_id]);
          }

          result_buf[row * result.stride() + result.start()] = sum;
        }

      }

      //
      // Hybrid Matrix
      //
      /** @brief Carries out sparse-matrix-dense-matrix multiplication with a hyb_matrix
      *
      * Implementation of the convenience expression C = prod(A, B);
      *
      * @param mat    The sparse matrix A
      * @param d_mat  The dense matrix B
      * @param result The dense result matrix C
      */
      template<typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::hyb_matrix<NumericT, ALIGNMENT> & mat,
                     const viennacl::matrix_base<NumericT, F1> & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result)
      {
        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat);
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat);
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat);
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat);
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat);
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat);
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat);

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        NumericT     const * elements       = detail::extract_raw_pointer<NumericT>(mat.handle());
        unsigned int const * coords         = detail::extract_raw_pointer<unsigned int>(mat.handle2());
        NumericT     const * csr_elements   = detail::extract_raw_pointer<NumericT>(mat.handle5());
        unsigned int const * csr_row_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle3());
        unsigned int const * csr_col_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle4());


        for (vcl_size_t result_col = 0; result_col < result.size2(); ++result_col)
        {
          for(vcl_size_t row = 0; row < mat.size1(); ++row)
          {
            NumericT sum = 0;

            //
            // Part 1: Process ELL part
            //
            for(unsigned int item_id = 0; item_id < mat.internal_ellnnz(); ++item_id)
            {
              vcl_size_t offset = row + item_id * mat.internal_size1();
              NumericT val = elements[offset];

              if(val != 0)
              {
                unsigned int col = coords[offset];
                sum += d_mat_wrapper(col, result_col) * val;
              }
            }

            //
            // Part 2: Process HYB/CSR part
            //
            vcl_size_t col_begin = csr_row_buffer[row];
            vcl_size_t col_end   = csr_row_buffer[row + 1];

            for(vcl_size_t item_id = col_begin; item_id < col_end; item_id++)
              sum += d_mat_wrapper(csr_col_buffer[item_id], result_col) * csr_elements[item_id];

            result_wrapper(row, result_col) = sum;
          }
        } // for result_col
      }


      /** @brief Carries out sparse-matrix-transposed-dense-matrix multiplication with a hyb_matrix
      *
      * Implementation of the convenience expression C = prod(A, trans(B));
      *
      * @param mat    The sparse matrix A
      * @param d_mat  The dense matrix B
      * @param result The dense result matrix C
      */
      template<typename NumericT, unsigned int ALIGNMENT, typename F1, typename F2>
      void prod_impl(const viennacl::hyb_matrix<NumericT, ALIGNMENT> & mat,
                     const viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                                        const viennacl::matrix_base<NumericT, F1>,
                                                        viennacl::op_trans > & d_mat,
                           viennacl::matrix_base<NumericT, F2> & result)
      {
        NumericT const * d_mat_data = detail::extract_raw_pointer<NumericT>(d_mat);
        NumericT       * result_data = detail::extract_raw_pointer<NumericT>(result);

        vcl_size_t d_mat_start1 = viennacl::traits::start1(d_mat.lhs());
        vcl_size_t d_mat_start2 = viennacl::traits::start2(d_mat.lhs());
        vcl_size_t d_mat_inc1   = viennacl::traits::stride1(d_mat.lhs());
        vcl_size_t d_mat_inc2   = viennacl::traits::stride2(d_mat.lhs());
        vcl_size_t d_mat_internal_size1  = viennacl::traits::internal_size1(d_mat.lhs());
        vcl_size_t d_mat_internal_size2  = viennacl::traits::internal_size2(d_mat.lhs());

        vcl_size_t result_start1 = viennacl::traits::start1(result);
        vcl_size_t result_start2 = viennacl::traits::start2(result);
        vcl_size_t result_inc1   = viennacl::traits::stride1(result);
        vcl_size_t result_inc2   = viennacl::traits::stride2(result);
        vcl_size_t result_internal_size1  = viennacl::traits::internal_size1(result);
        vcl_size_t result_internal_size2  = viennacl::traits::internal_size2(result);

        detail::matrix_array_wrapper<NumericT const, typename F1::orientation_category, false>
            d_mat_wrapper(d_mat_data, d_mat_start1, d_mat_start2, d_mat_inc1, d_mat_inc2, d_mat_internal_size1, d_mat_internal_size2);
        detail::matrix_array_wrapper<NumericT,       typename F2::orientation_category, false>
            result_wrapper(result_data, result_start1, result_start2, result_inc1, result_inc2, result_internal_size1, result_internal_size2);

        NumericT     const * elements       = detail::extract_raw_pointer<NumericT>(mat.handle());
        unsigned int const * coords         = detail::extract_raw_pointer<unsigned int>(mat.handle2());
        NumericT     const * csr_elements   = detail::extract_raw_pointer<NumericT>(mat.handle5());
        unsigned int const * csr_row_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle3());
        unsigned int const * csr_col_buffer = detail::extract_raw_pointer<unsigned int>(mat.handle4());


        for (vcl_size_t result_col = 0; result_col < result.size2(); ++result_col)
        {
          for(vcl_size_t row = 0; row < mat.size1(); ++row)
          {
            NumericT sum = 0;

            //
            // Part 1: Process ELL part
            //
            for(unsigned int item_id = 0; item_id < mat.internal_ellnnz(); ++item_id)
            {
              vcl_size_t offset = row + item_id * mat.internal_size1();
              NumericT val = elements[offset];

              if(val != 0)
              {
                unsigned int col = coords[offset];
                sum += d_mat_wrapper(result_col, col) * val;
              }
            }

            //
            // Part 2: Process HYB/CSR part
            //
            vcl_size_t col_begin = csr_row_buffer[row];
            vcl_size_t col_end   = csr_row_buffer[row + 1];

            for(vcl_size_t item_id = col_begin; item_id < col_end; item_id++)
              sum += d_mat_wrapper(result_col, csr_col_buffer[item_id]) * csr_elements[item_id];

            result_wrapper(row, result_col) = sum;
          }
        } // for result_col
      }


    } // namespace host_based
  } //namespace linalg
} //namespace viennacl


#endif
