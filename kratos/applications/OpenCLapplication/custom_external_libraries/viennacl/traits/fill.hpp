#ifndef VIENNACL_TRAITS_FILL_HPP_
#define VIENNACL_TRAITS_FILL_HPP_

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

/** @file fill.hpp
    @brief Generic fill functionality for different matrix types
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/meta/result_of.hpp"

#ifdef VIENNACL_HAVE_EIGEN  
#include <Eigen/Core>
#include <Eigen/Sparse>
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
    template <typename MatrixType, typename SCALARTYPE>
    void fill(MatrixType & matrix, std::size_t row_index, std::size_t col_index, SCALARTYPE value)
    {
      matrix(row_index, col_index) = value; 
    }
    
    #ifdef VIENNACL_HAVE_EIGEN
    template <typename T, int options, typename SCALARTYPE>
    inline void fill(Eigen::SparseMatrix<T, options> & m,
                     std::size_t row_index,
                     std::size_t col_index,
                     SCALARTYPE value
                    )
    {
      m.fill(row_index, col_index) = value;
    }    
    #endif

 
  } //namespace traits
} //namespace viennacl
    

#endif
