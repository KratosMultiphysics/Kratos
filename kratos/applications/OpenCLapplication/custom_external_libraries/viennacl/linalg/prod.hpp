/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory
======================================================================= */

#ifndef _VIENNACL_PROD_HPP_
#define _VIENNACL_PROD_HPP_

/** @file prod.hpp
    @brief Generic interface for matrix-vector and matrix-matrix products. 
           See viennacl/linalg/vector_operations.hpp, viennacl/linalg/matrix_operations.hpp,
           viennacl/linalg/compressed_matrix_operations.hpp and viennacl/linalg/coordinate_matrix_operations.hpp for implementations.
*/

#include "viennacl/forwards.h"
#include "tag_of.hpp"
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
    #ifdef VIENNACL_HAVE_MTL4
    // ----------------------------------------------------
    // mtl4
    //
    template< typename MatrixT, typename VectorT >
    VectorT 
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_mtl4< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "mtl4 .. " << std::endl;
      return VectorT(matrix * vector);
    }
    #endif
    
    #ifdef VIENNACL_HAVE_EIGEN
    // ----------------------------------------------------
    // Eigen
    //
    template< typename MatrixT, typename VectorT >
    VectorT 
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_eigen< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "ublas .. " << std::endl;
      return matrix * vector;
    }
    #endif
    
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename MatrixT, typename VectorT >
    VectorT 
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
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
    prod_impl(std::vector< std::vector<T, A1>, A2 > const& matrix, VectorT const& vector)
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
    prod_impl(std::vector< std::map<KEY, DATA, COMPARE, AMAP>, AVEC > const& matrix, VectorT const& vector)
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
    
    
    template< typename MatrixT, typename VectorT >
    VectorT 
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "std .. " << std::endl;
      return prod_impl(matrix, vector);
    }

    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename MatrixT, typename NumericT, unsigned int ALIGNMENT >
    viennacl::vector_expression< const MatrixT, 
                                 const viennacl::vector<NumericT, ALIGNMENT>,
                                 viennacl::op_prod >
    prod(MatrixT const& matrix,
         viennacl::vector<NumericT, ALIGNMENT> const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::linalg::prod_impl(matrix, vector);
    }

    template< typename MatrixT, typename NumericT, typename F, unsigned int ALIGNMENT >
    viennacl::matrix_expression< const MatrixT, 
                                 const viennacl::matrix<NumericT, F, ALIGNMENT>,
                                 viennacl::op_prod >
    prod(MatrixT const& matrix_A,
         viennacl::matrix<NumericT, F, ALIGNMENT> const& matrix_B, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const MatrixT, 
                                          const viennacl::matrix<NumericT, F, ALIGNMENT>,
                                          viennacl::op_prod >(matrix_A, matrix_B);
    }

    template< typename MatrixT, typename NumericT, typename F, unsigned int ALIGNMENT >
    viennacl::matrix_expression< const MatrixT, 
                                 const viennacl::matrix_expression< const viennacl::matrix<NumericT, F, ALIGNMENT>, 
                                                                    const viennacl::matrix<NumericT, F, ALIGNMENT>,
                                                                    viennacl::op_trans >,
                                 viennacl::op_prod >
    prod(MatrixT const& matrix_A,
         const viennacl::matrix_expression< const viennacl::matrix<NumericT, F, ALIGNMENT>, 
                                            const viennacl::matrix<NumericT, F, ALIGNMENT>,
                                            viennacl::op_trans > & matrix_B,
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::matrix_expression< const MatrixT, 
                                          const viennacl::matrix_expression< const viennacl::matrix<NumericT, F, ALIGNMENT>, 
                                                                             const viennacl::matrix<NumericT, F, ALIGNMENT>,
                                                                             viennacl::op_trans >,
                                          viennacl::op_prod >(matrix_A, matrix_B);
      //return viennacl::linalg::prod_impl(matrix_A, matrix_B);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





