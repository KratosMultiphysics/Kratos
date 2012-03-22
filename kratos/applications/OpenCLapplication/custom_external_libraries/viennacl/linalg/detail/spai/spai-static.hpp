#ifndef VIENNACL_LINALG_DETAIL_SPAI_SPAI_STATIC_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_SPAI_STATIC_HPP

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

/** @file viennacl/linalg/detail/spai/spai-static.hpp
    @brief Implementation of a static SPAI. Experimental in 1.2.x.
    
    SPAI code contributed by Nikolay Lukash
*/

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>
//#include "spai-dynamic.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/storage.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/triangular.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"
// ViennaCL includes
#include "viennacl/linalg/prod.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/compressed_matrix_operations.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"

//#include "boost/numeric/ublas/detail/matrix_assign.hpp"

namespace viennacl
{
    namespace linalg
    {
      namespace detail
      {
        namespace spai
        {
        
          /********************************* STATIC SPAI FUNCTIONS******************************************/
          
          /** @brief Projects solution of LS problem onto original column m 
          * @param m_in solution of LS
          * @param J set of non-zero columns 
          * @param m original column of M
          */
          template <typename VectorType, typename SparseVectorType>
          void fanOutVector(const VectorType& m_in, const std::vector<unsigned int>& J, SparseVectorType& m){
              unsigned int  cnt = 0;
              for (size_t i = 0; i < J.size(); ++i) {
                  m[J[i]] = m_in(cnt++);
              }
          }
          /** @brief Solution of linear:R*x=y system by backward substitution
          * @param R uppertriangular matrix 
          * @param y right handside vector
          * @param x solution vector
          */
          template <typename MatrixType, typename VectorType>
          void backwardSolve(const MatrixType& R, const VectorType& y, VectorType& x){
              typedef typename MatrixType::value_type ScalarType;
              for (long i = R.size2()-1; i >= 0 ; i--) {
                  x(i) = y(i);
                  for (size_t j = i+1; j < R.size2(); ++j) {
                      x(i) -= R(i,j)*x(j);
                  }
                  x(i) /= R(i,i);
              }
          }
          /** @brief Perform projection of set I on the unit-vector
          * @param I set of non-zero rows
          * @param y result vector
          * @param ind index of unit vector
          */
          template <typename VectorType, typename ScalarType>
          void projectI(const std::vector<unsigned int>& I, VectorType& y, unsigned int ind){
              for(size_t i = 0; i < I.size(); ++i){
                  //y.resize(y.size()+1);
                  if(I[i] == ind){
                      y(i) = static_cast<ScalarType>(1.0);
                  }
                  else{
                      y(i) = static_cast<ScalarType>(0.0);
                  }
              }
          }
          
          /** @brief Builds index set of projected columns for current column of preconditioner
          * @param v current column of preconditioner
          * @param J output - index set of non-zero columns
          */
          template <typename SparseVectorType>
          void buildColumnIndexSet(const SparseVectorType& v, std::vector<unsigned int>& J){
              //typedef typename VectorType::value_type ScalarType;
              unsigned int tmp_v;
              for(typename SparseVectorType::const_iterator vec_it = v.begin(); vec_it != v.end(); ++vec_it){
                  tmp_v = vec_it->first;
                  J.push_back(vec_it->first);
              }
              std::sort(J.begin(), J.end());
          }
          
          /** @brief Initialize preconditioner with sparcity pattern = p(A)
          * @param A input matrix
          * @param M output matrix - initialized preconditioner
          */
          template <typename SparseMatrixType>
          void initPreconditioner(const SparseMatrixType& A, SparseMatrixType& M){
              typedef typename SparseMatrixType::value_type ScalarType;
              M.resize(A.size1(), A.size2(), false);
              for(typename SparseMatrixType::const_iterator1 row_it = A.begin1(); row_it!= A.end1(); ++row_it){
                  //
                  for(typename SparseMatrixType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it){
                      M(col_it.index1(),col_it.index2()) = static_cast<ScalarType>(1);
                  }
              }
          }
          
          /** @brief Row projection for matrix A(:,J) -> A(I,J), building index set of non-zero rows
          * @param A_v_c input matrix
          * @param J set of non-zero rows
          * @param I output matrix 
          */
          template <typename SparseVectorType>
          void projectRows(const std::vector<SparseVectorType>& A_v_c, const std::vector<unsigned int>& J, std::vector<unsigned int>& I){
              for(size_t i = 0; i < J.size(); ++i){
                  for(typename SparseVectorType::const_iterator col_it = A_v_c[J[i]].begin(); col_it!=A_v_c[J[i]].end(); ++col_it){
                      if(!isInIndexSet(I, col_it->first)){
                          I.push_back(col_it->first);
                      }
                  }
              }
              std::sort(I.begin(), I.end());
          }
        }
      }
    }
}

#endif