#ifndef VIENNACL_LINALG_DETAIL_ILU_COMMON_HPP_
#define VIENNACL_LINALG_DETAIL_ILU_COMMON_HPP_

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

/** @file viennacl/linalg/detail/ilu/common.hpp
    @brief Common routines used within ILU-type preconditioners
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
    
      /** @brief Increments a row iterator (iteration along increasing row indices) up to a certain row index k.
      * 
      * Generic implementation using the iterator concept from boost::numeric::ublas. Could not find a better way for sparse matrices...
      *
      * @param row_iter   The row iterator
      * @param k      The final row index
      */
      template <typename T>
      void ilu_inc_row_iterator_to_row_index(T & row_iter, unsigned int k)
      {
        while (row_iter.index1() < k)
          ++row_iter;
      }
      
      /** @brief Increments a row iterator (iteration along increasing row indices) up to a certain row index k.
      * 
      * Specialization for the sparse matrix adapter shipped with ViennaCL
      *
      * @param row_iter   The row iterator
      * @param k      The final row index
      */
      template <typename ScalarType>
      void ilu_inc_row_iterator_to_row_index(viennacl::tools::sparse_matrix_adapter<ScalarType> & row_iter, unsigned int k)
      {
        row_iter += k - row_iter.index1();
      }
      
      /** @brief Increments a row iterator (iteration along increasing row indices) up to a certain row index k.
      * 
      * Specialization for the const sparse matrix adapter shipped with ViennaCL
      *
      * @param row_iter   The row iterator
      * @param k      The final row index
      */
      template <typename ScalarType>
      void ilu_inc_row_iterator_to_row_index(viennacl::tools::const_sparse_matrix_adapter<ScalarType> & row_iter, unsigned int k)
      {
        row_iter += k - row_iter.index1();
      }

      /** @brief Generic inplace solution of a unit lower triangular system
      *   
      * @param mat  The system matrix
      * @param vec  The right hand side vector
      */
      template<typename MatrixType, typename VectorType>
      void ilu_inplace_solve(MatrixType const & mat, VectorType & vec, viennacl::linalg::unit_lower_tag)
      {
        typedef typename MatrixType::const_iterator1    InputRowIterator;  //iterate along increasing row index
        typedef typename MatrixType::const_iterator2    InputColIterator;  //iterate along increasing column index
        
        for (InputRowIterator row_iter = mat.begin1(); row_iter != mat.end1(); ++row_iter)
        {
          for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
          {
            if (col_iter.index2() < col_iter.index1())
              vec[col_iter.index1()] -= *col_iter * vec[col_iter.index2()];
          }
        }
      }

      /** @brief Generic inplace solution of a upper triangular system
      *   
      * @param mat  The system matrix
      * @param vec  The right hand side vector
      */
      template<typename MatrixType, typename VectorType>
      void ilu_inplace_solve(MatrixType const & mat, VectorType & vec, viennacl::linalg::upper_tag)
      {
        typedef typename MatrixType::const_reverse_iterator1    InputRowIterator;  //iterate along increasing row index
        typedef typename MatrixType::const_iterator2            InputColIterator;  //iterate along increasing column index
        typedef typename VectorType::value_type                 ScalarType;
        
        ScalarType diagonal_entry = 1.0;
        
        for (InputRowIterator row_iter = mat.rbegin1(); row_iter != mat.rend1(); ++row_iter)
        {
          for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
          {
            if (col_iter.index2() > col_iter.index1())
              vec[col_iter.index1()] -= *col_iter * vec[col_iter.index2()];
            if (col_iter.index2() == col_iter.index1())
              diagonal_entry = *col_iter;
          }
          vec[row_iter.index1()] /= diagonal_entry;
        }
      }

      /** @brief Generic LU substitution
      *   
      * @param mat  The system matrix
      * @param vec  The right hand side vector
      */
      template<typename MatrixType, typename VectorType>
      void ilu_lu_substitute(MatrixType const & mat, VectorType & vec)
      {
        ilu_inplace_solve(mat, vec, unit_lower_tag());
        ilu_inplace_solve(mat, vec, upper_tag());
      }

    } // namespace detail
  } // namespace linalg
} // namespace viennacl




#endif



