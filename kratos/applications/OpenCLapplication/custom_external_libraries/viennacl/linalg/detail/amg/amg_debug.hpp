#ifndef VIENNACL_LINALG_DETAIL_AMG_AMG_DEBUG_HPP
#define VIENNACL_LINALG_DETAIL_AMG_AMG_DEBUG_HPP

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

/** @file amg_debug.hpp
    @brief Debug functionality for AMG. To be removed.
    
    AMG code contributed by Markus Wagner
*/

#include <iostream>
#include "viennacl/io/matrix_market.hpp"

#ifdef SMALL_SIZE
#define VIENNACL_AMG_MATRIXTYPE boost::numeric::ublas::matrix<ScalarType>
#else
#define VIENNACL_AMG_MATRIXTYPE MatrixType
#endif

namespace viennacl
{
  namespace linalg
  {    
    namespace detail
    {
      namespace amg
      {

        template <typename MatrixType>
        void printmatrix(MatrixType & mat, int const value=-1)
        {
          typedef typename MatrixType::value_type ScalarType;  
          typedef typename VIENNACL_AMG_MATRIXTYPE::iterator1 InternalRowIterator;
          typedef typename VIENNACL_AMG_MATRIXTYPE::iterator2 InternalColIterator;
          
          #ifdef DEBUG
          VIENNACL_AMG_MATRIXTYPE mat2 = mat;
          
          for (InternalRowIterator row_iter = mat2.begin1(); row_iter != mat2.end1(); ++row_iter)
          {
            for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {     
              std::cout << *col_iter << " ";
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
          #endif
        }

        template <typename VectorType>
        void printvector(VectorType const & vec)
        {
          #ifdef DEBUGBENCH
          for (typename VectorType::const_iterator iter = vec.begin(); iter != vec.end(); ++iter)
          {
            std::cout << *iter << " ";
          }
          std::cout << std::endl;
          #endif
        }

      }
    }
  }
}
#endif
