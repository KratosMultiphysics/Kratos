#ifndef VIENNACL_LINALG_DETAIL_SPAI_BLOCK_VECTOR_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_BLOCK_VECTOR_HPP

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

/** @file viennacl/linalg/detail/spai/block_vector.hpp
    @brief Implementation of a bunch of vectors on GPU. Experimental.
    
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
          * @brief Represents a contigious vector on GPU
          */
          
          class block_vector{
          public:
              block_vector(){
              }
              /**
              * @brief Return handle to the elements
              */
              viennacl::ocl::handle<cl_mem>& handle(){ return _elements; }
              /**
              * @brief Return handle to start indices
              */
              viennacl::ocl::handle<cl_mem>& handle1() { return _start_block_inds; }
              
              /**
              * @brief Return handle to the const elements
              */
              const viennacl::ocl::handle<cl_mem>& handle() const { return _elements; }
              /**
              * @brief Return handle to const start indices
              */
              const viennacl::ocl::handle<cl_mem>& handle1() const { return _start_block_inds; }
          private:
              //unsigned int _vectorIndex;
              viennacl::ocl::handle<cl_mem> _elements;
              viennacl::ocl::handle<cl_mem> _start_block_inds;
          };
        }
      }
    }
}
#endif