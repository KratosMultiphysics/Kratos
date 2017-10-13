#ifndef VIENNACL_LINALG_DETAIL_SPAI_BLOCK_VECTOR_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_BLOCK_VECTOR_HPP

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

          class block_vector
          {
            public:

              /**
              * @brief Return handle to the elements
              */
              viennacl::ocl::handle<cl_mem>& handle(){ return elements_; }
              /**
              * @brief Return handle to start indices
              */
              viennacl::ocl::handle<cl_mem>& handle1() { return start_block_inds_; }

              /**
              * @brief Return handle to the const elements
              */
              const viennacl::ocl::handle<cl_mem>& handle() const { return elements_; }
              /**
              * @brief Return handle to const start indices
              */
              const viennacl::ocl::handle<cl_mem>& handle1() const { return start_block_inds_; }
            private:
              viennacl::ocl::handle<cl_mem> elements_;
              viennacl::ocl::handle<cl_mem> start_block_inds_;
          };
        }
      }
    }
}
#endif
