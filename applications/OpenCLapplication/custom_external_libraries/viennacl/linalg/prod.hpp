#ifndef VIENNACL_LINALG_PROD_HPP_
#define VIENNACL_LINALG_PROD_HPP_

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

/** @file viennacl/linalg/prod.hpp
    @brief Generic interface for matrix-vector and matrix-matrix products.
           See viennacl/linalg/vector_operations.hpp, viennacl/linalg/matrix_operations.hpp, and
           viennacl/linalg/sparse_matrix_operations.hpp for implementations.
*/

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"
#include <vector>
#include <map>

namespace viennacl
{
  //
  // generic prod function
  //   uses tag dispatch to identify which algorithm
  //   should be called
  //
  namespace linalg
  {
    #ifdef VIENNACL_WITH_MTL4
    // ----------------------------------------------------
    // mtl4
    //
    template< typename MatrixT, typename VectorT >
    typename viennacl::enable_if< viennacl::is_mtl4< typename viennacl::traits::tag_of< MatrixT >::type >::value,
                                  VectorT>::type
    prod(MatrixT const& matrix, VectorT const& vector)
    {
      return VectorT(matrix * vector);
    }
    #endif

    #ifdef VIENNACL_WITH_EIGEN
    // ----------------------------------------------------
    // Eigen
    //
    template< typename MatrixT, typename VectorT >
    typename viennacl::enable_if< viennacl::is_eigen< typename viennacl::traits::tag_of< MatrixT >::type >::value,
                                  VectorT>::type
    prod(MatrixT const& matrix, VectorT const& vector)
    {
      return matrix * vector;
    }
    #endif

    #ifdef VIENNACL_WITH_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename MatrixT, typename VectorT >
    typename viennacl::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< MatrixT >::type >::value,
                                  VectorT>::type
    prod(MatrixT const& matrix, VectorT const& vector)
    {
      // std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::prod(matrix, vector);
    }
    #endif


    // ----------------------------------------------------
    // STL type
    //

    // dense matrix-vector product:
    template< typename T, typename A1, typename A2, typename VectorT >
    VectorT
    prod(std::vector< std::vector<T, A1>, A2 > const & matrix, VectorT const& vector)
    {
      VectorT result(matrix.size());
      for (typename std::vector<T, A1>::size_type i=0; i<matrix.size(); ++i)
      {
        result[i] = 0; //we will not assume that VectorT is initialized to zero
        for (typename std::vector<T, A1>::size_type j=0; j<matrix[i].size(); ++j)
          result[i] += matrix[i][j] * vector[j];
      }
      return result;
    }

    // sparse matrix-vector product:
    template< typename KEY, typename DATA, typename COMPARE, typename AMAP, typename AVEC, typename VectorT >
    VectorT
    prod(std::vector< std::map<KEY, DATA, COMPARE, AMAP>, AVEC > const& matrix, VectorT const& vector)
    {
      typedef std::vector< std::map<KEY, DATA, COMPARE, AMAP>, AVEC > MatrixType;

      VectorT result(matrix.size());
      for (typename MatrixType::size_type i=0; i<matrix.size(); ++i)
      {
        result[i] = 0; //we will not assume that VectorT is initialized to zero
        for (typename std::map<KEY, DATA, COMPARE, AMAP>::const_iterator row_entries = matrix[i].begin();
             row_entries != matrix[i].end();
             ++row_entries)
          result[i] += row_entries->second * vector[row_entries->first];
      }
      return result;
    }


    /*template< typename MatrixT, typename VectorT >
    VectorT
    prod(MatrixT const& matrix, VectorT const& vector,
         typename viennacl::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                     >::type* dummy = 0)
    {
      // std::cout << "std .. " << std::endl;
      return prod_impl(matrix, vector);
    }*/

    // ----------------------------------------------------
    // VIENNACL
    //

    // standard product:
    template< typename NumericT, typename F1, typename F2>
    viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                 const viennacl::matrix_base<NumericT, F2>,
                                 viennacl::op_mat_mat_prod >
    prod(viennacl::matrix_base<NumericT, F1> const & A,
         viennacl::matrix_base<NumericT, F2> const & B)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                          const viennacl::matrix_base<NumericT, F2>,
                                          viennacl::op_mat_mat_prod >(A, B);
    }

    // right factor is transposed:
    template< typename NumericT, typename F1, typename F2>
    viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                 const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                                                   const viennacl::matrix_base<NumericT, F2>,
                                                                   op_trans>,
                                 viennacl::op_mat_mat_prod >
    prod(viennacl::matrix_base<NumericT, F1> const & A,
         viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                     const viennacl::matrix_base<NumericT, F2>,
                                     op_trans> const & B)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const viennacl::matrix_base<NumericT, F1>,
                                          const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                                                            const viennacl::matrix_base<NumericT, F2>,
                                                                            op_trans>,
                                          viennacl::op_mat_mat_prod >(A, B);
    }

    // left factor transposed:
    template< typename NumericT, typename F1, typename F2>
    viennacl::matrix_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                                                   const viennacl::matrix_base<NumericT, F1>,
                                                                   op_trans>,
                                 const viennacl::matrix_base<NumericT, F2>,
                                 viennacl::op_mat_mat_prod >
    prod(viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                     const viennacl::matrix_base<NumericT, F1>,
                                     op_trans> const & A,
         viennacl::matrix_base<NumericT, F2> const & B)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                                                            const viennacl::matrix_base<NumericT, F1>,
                                                                            op_trans>,
                                          const viennacl::matrix_base<NumericT, F2>,
                                          viennacl::op_mat_mat_prod >(A, B);
    }


    // both factors transposed:
    template< typename NumericT, typename F1, typename F2>
    viennacl::matrix_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                                                   const viennacl::matrix_base<NumericT, F1>,
                                                                   op_trans>,
                                 const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                                                   const viennacl::matrix_base<NumericT, F2>,
                                                                   op_trans>,
                                 viennacl::op_mat_mat_prod >
    prod(viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                     const viennacl::matrix_base<NumericT, F1>,
                                     op_trans> const & A,
         viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                     const viennacl::matrix_base<NumericT, F2>,
                                     op_trans> const & B)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F1>,
                                                                            const viennacl::matrix_base<NumericT, F1>,
                                                                            op_trans>,
                                          const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F2>,
                                                                            const viennacl::matrix_base<NumericT, F2>,
                                                                            op_trans>,
                                          viennacl::op_mat_mat_prod >(A, B);
    }



    // matrix-vector product
    template< typename NumericT, typename F>
    viennacl::vector_expression< const viennacl::matrix_base<NumericT, F>,
                                 const viennacl::vector_base<NumericT>,
                                 viennacl::op_prod >
    prod(viennacl::matrix_base<NumericT, F> const & matrix,
         viennacl::vector_base<NumericT> const & vector)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::vector_expression< const viennacl::matrix_base<NumericT, F>,
                                          const viennacl::vector_base<NumericT>,
                                          viennacl::op_prod >(matrix, vector);
    }

    // transposed matrix-vector product
    template< typename NumericT, typename F>
    viennacl::vector_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F>,
                                                                   const viennacl::matrix_base<NumericT, F>,
                                                                   op_trans>,
                                 const viennacl::vector_base<NumericT>,
                                 viennacl::op_prod >
    prod(viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F>,
                                     const viennacl::matrix_base<NumericT, F>,
                                     op_trans> const & matrix,
         viennacl::vector_base<NumericT> const & vector)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::vector_expression< const viennacl::matrix_expression<const viennacl::matrix_base<NumericT, F>,
                                                                            const viennacl::matrix_base<NumericT, F>,
                                                                            op_trans>,
                                          const viennacl::vector_base<NumericT>,
                                          viennacl::op_prod >(matrix, vector);
    }


    template<typename SparseMatrixType, class SCALARTYPE>
    typename viennacl::enable_if< viennacl::is_any_sparse_matrix<SparseMatrixType>::value,
                                  vector_expression<const SparseMatrixType,
                                                    const vector_base<SCALARTYPE>,
                                                    op_prod >
                                 >::type
    prod(const SparseMatrixType & mat,
         const vector_base<SCALARTYPE> & vec)
    {
      return vector_expression<const SparseMatrixType,
                               const vector_base<SCALARTYPE>,
                               op_prod >(mat, vec);
    }

    template< typename SparseMatrixType, typename SCALARTYPE, typename F1>
    typename viennacl::enable_if< viennacl::is_any_sparse_matrix<SparseMatrixType>::value,
                                  viennacl::matrix_expression<const SparseMatrixType,
                                                              const matrix_base < SCALARTYPE, F1 >,
                                                              op_prod >
                                 >::type
    prod(const SparseMatrixType & sp_mat,
         const viennacl::matrix_base<SCALARTYPE, F1> & d_mat)
    {
      return viennacl::matrix_expression<const SparseMatrixType,
                                         const viennacl::matrix_base < SCALARTYPE, F1 >,
                                         op_prod >(sp_mat, d_mat);
    }

    // right factor is transposed
    template< typename SparseMatrixType, typename SCALARTYPE, typename F1 >
    typename viennacl::enable_if< viennacl::is_any_sparse_matrix<SparseMatrixType>::value,
                                  viennacl::matrix_expression< const SparseMatrixType,
                                                               const viennacl::matrix_expression<const viennacl::matrix_base<SCALARTYPE, F1>,
                                                                                                 const viennacl::matrix_base<SCALARTYPE, F1>,
                                                                                                 op_trans>,
                                                               viennacl::op_prod >
                                  >::type
    prod(const SparseMatrixType & A,
         viennacl::matrix_expression<const viennacl::matrix_base < SCALARTYPE, F1 >,
                                     const viennacl::matrix_base < SCALARTYPE, F1 >,
                                     op_trans> const & B)
    {
      return viennacl::matrix_expression< const SparseMatrixType,
                                          const viennacl::matrix_expression<const viennacl::matrix_base < SCALARTYPE, F1 >,
                                                                            const viennacl::matrix_base < SCALARTYPE, F1 >,
                                                                            op_trans>,
                                          viennacl::op_prod >(A, B);
    }

    template<typename StructuredMatrixType, class SCALARTYPE>
    typename viennacl::enable_if< viennacl::is_any_dense_structured_matrix<StructuredMatrixType>::value,
                                  vector_expression<const StructuredMatrixType,
                                                    const vector_base<SCALARTYPE>,
                                                    op_prod >
                                 >::type
    prod(const StructuredMatrixType & mat,
         const vector_base<SCALARTYPE> & vec)
    {
      return vector_expression<const StructuredMatrixType,
                               const vector_base<SCALARTYPE>,
                               op_prod >(mat, vec);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





