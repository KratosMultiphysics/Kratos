#ifndef VIENNACL_TRAITS_SIZE_HPP_
#define VIENNACL_TRAITS_SIZE_HPP_

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

/** @file size.hpp
    @brief Generic size and resize functionality for different vector and matrix types
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/meta/result_of.hpp"

#ifdef VIENNACL_HAVE_UBLAS  
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif

#ifdef VIENNACL_HAVE_EIGEN  
#include <Eigen/Core>
#include <Eigen/Sparse>
#endif

#ifdef VIENNACL_HAVE_MTL4
#include <boost/numeric/mtl/mtl.hpp>
#endif

#include <vector>
#include <map>

namespace viennacl
{

  namespace traits
  {
    //
    // Resize: Change the size of vectors and matrices
    //
    template <typename MatrixType>
    void resize(MatrixType & matrix, size_t rows, size_t cols)
    {
      matrix.resize(rows, cols); 
    }
    
    template <typename VectorType>
    void resize(VectorType & vec, size_t new_size)
    {
      vec.resize(new_size); 
    }
    
    #ifdef VIENNACL_HAVE_UBLAS  
    //ublas needs separate treatment:
    template <typename ScalarType>
    void resize(boost::numeric::ublas::compressed_matrix<ScalarType> & matrix,
                size_t rows,
                size_t cols)
    {
      matrix.resize(rows, cols, false); //Note: omitting third parameter leads to compile time error (not implemented in ublas <= 1.42) 
    }
    #endif  
    
    
    #ifdef VIENNACL_HAVE_MTL4
    template <typename ScalarType>
    void resize(mtl::compressed2D<ScalarType> & matrix,
                size_t rows,
                size_t cols)
    {
      matrix.change_dim(rows, cols);
    }
    
    template <typename ScalarType>
    void resize(mtl::dense_vector<ScalarType> & vec,
                size_t new_size)
    {
      vec.change_dim(new_size);
    }
    #endif      

    #ifdef VIENNACL_HAVE_EIGEN
    inline void resize(Eigen::MatrixXf & m,
                       std::size_t new_rows,
                       std::size_t new_cols)
    {
      m.resize(new_rows, new_cols);
    }
    
    inline void resize(Eigen::MatrixXd & m,
                       std::size_t new_rows,
                       std::size_t new_cols)
    {
      m.resize(new_rows, new_cols);
    }
    
    template <typename T, int options>
    inline void resize(Eigen::SparseMatrix<T, options> & m,
                       std::size_t new_rows,
                       std::size_t new_cols)
    {
      m.resize(new_rows, new_cols);
    }    
    
    inline void resize(Eigen::VectorXf & v,
                       std::size_t new_size)
    {
      v.resize(new_size);
    }
    
    inline void resize(Eigen::VectorXd & v,
                       std::size_t new_size)
    {
      v.resize(new_size);
    }
    #endif


    //
    // size: Returns the length of vectors
    //
    template <typename VectorType>
    typename result_of::size_type<VectorType>::type size(VectorType const & vec)
    {
      return vec.size(); 
    }

    #ifdef VIENNACL_HAVE_MTL4
    template <typename ScalarType>
    typename result_of::size_type< mtl::dense_vector<ScalarType> >::type
    size(mtl::dense_vector<ScalarType> const & vec) { return vec.used_memory(); }
    #endif
    
    #ifdef VIENNACL_HAVE_EIGEN
    inline std::size_t size(Eigen::VectorXf const & v) { return v.rows(); }
    inline std::size_t size(Eigen::VectorXd const & v) { return v.rows(); }
    #endif

    //
    // size1: No. of rows for matrices
    //
    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    size1(MatrixType const & mat) { return mat.size1(); }

    #ifdef VIENNACL_HAVE_EIGEN
    inline std::size_t size1(Eigen::MatrixXf const & m) { return m.rows(); }
    inline std::size_t size1(Eigen::MatrixXd const & m) { return m.rows(); }
    template <typename T, int options>
    inline std::size_t size1(Eigen::SparseMatrix<T, options> & m) { return m.rows(); }    
    #endif

    //
    // size2: No. of columns for matrices
    //
    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    size2(MatrixType const & mat) { return mat.size2(); }
 
    #ifdef VIENNACL_HAVE_EIGEN
    inline std::size_t size2(Eigen::MatrixXf const & m) { return m.cols(); }
    inline std::size_t size2(Eigen::MatrixXd const & m) { return m.cols(); }
    template <typename T, int options>
    inline std::size_t size2(Eigen::SparseMatrix<T, options> & m) { return m.cols(); }    
    #endif
 
    //
    // internal_size: Returns the internal (padded) length of vectors
    //
    template <typename VectorType>
    typename result_of::size_type<VectorType>::type 
    internal_size(VectorType const & vec)
    {
      return vec.internal_size(); 
    }

    template <typename VectorType>
    typename result_of::size_type<VectorType>::type 
    internal_size(viennacl::vector_range<VectorType> const & vec)
    {
      return vec.get().internal_size(); 
    }

    //
    // internal_size1: No. of internal (padded) rows for matrices
    //
    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    internal_size1(MatrixType const & mat) { return mat.internal_size1(); }

    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    internal_size1(viennacl::matrix_range<MatrixType> const & mat) { return mat.get().internal_size1(); }

    //
    // internal_size2: No. of internal (padded) columns for matrices
    //
    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    internal_size2(MatrixType const & mat) { return mat.internal_size2(); }
 
    template <typename MatrixType>
    typename result_of::size_type<MatrixType>::type
    internal_size2(viennacl::matrix_range<MatrixType> const & mat) { return mat.get().internal_size2(); }

 
  } //namespace traits
} //namespace viennacl
    

#endif
