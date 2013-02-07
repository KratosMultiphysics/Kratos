#ifndef VIENNACL_LINALG_DETAIL_SPAI_BLOCK_MATRIX_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_BLOCK_MATRIX_HPP

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

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include "viennacl/ocl/backend.hpp"
#include "viennacl/tools/tools.hpp"

/** @file viennacl/linalg/detail/spai/block_matrix.hpp
    @brief Implementation of a bunch of (small) matrices on GPU. Experimental.
    
    SPAI code contributed by Nikolay Lukash
*/

namespace viennacl
{
    namespace linalg
    {
      namespace detail
      {
        namespace spai
        {
          
          /**
          * @brief Represents a contigious matrices on GPU
          */
          
          class block_matrix{
          public:
              block_matrix(){
                  
              }
              /**
              * @brief Returns a handle to the elements
              */
              viennacl::ocl::handle<cl_mem>& handle(){ return _elements; }
              /**
              * @brief Returns a handle to the matrix dimensions
              */
              viennacl::ocl::handle<cl_mem>& handle1() { return _matrix_dimensions; }
              /**
              * @brief Returns a handle to the start indices of matrix
              */
              viennacl::ocl::handle<cl_mem>& handle2() { return _start_block_inds; }
              
              /**
              * @brief Returns a handle to the const elements
              */
              const viennacl::ocl::handle<cl_mem>& handle() const { return _elements; }
              /**
              * @brief Returns a handle to the const matrix dimensions
              */
              const viennacl::ocl::handle<cl_mem>& handle1() const { return _matrix_dimensions; }
              /**
              * @brief Returns a handle to the const start indices of matrix
              */
              const viennacl::ocl::handle<cl_mem>& handle2() const { return _start_block_inds; }
          private:
              //unsigned int _vectorIndex;
              viennacl::ocl::handle<cl_mem> _elements;
              viennacl::ocl::handle<cl_mem> _matrix_dimensions;
              viennacl::ocl::handle<cl_mem> _start_block_inds;
          };
        
        
        }
      }
    }
}
#endif