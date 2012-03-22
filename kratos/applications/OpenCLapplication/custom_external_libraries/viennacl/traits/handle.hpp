#ifndef VIENNACL_TRAITS_HANDLE_HPP_
#define VIENNACL_TRAITS_HANDLE_HPP_

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

/** @file traits/handle.hpp
    @brief Extracts the underlying OpenCL handle from a vector, a matrix, an expression etc.
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

namespace viennacl
{
  namespace traits
  {
    
    // Returns the OpenCL handle of a ViennaCL object
    template <typename T>
    viennacl::ocl::handle<cl_mem> handle(T & obj)
    {
      return obj.handle();
    }

    template <typename T>
    viennacl::ocl::handle<cl_mem> handle(viennacl::vector_range<T> & obj)
    {
      return handle(obj.get());
    }

    template <typename T>
    viennacl::ocl::handle<cl_mem> handle(viennacl::vector_range<T> const & obj)
    {
      return handle(obj.get());
    }

    template <typename T>
    viennacl::ocl::handle<cl_mem> handle(viennacl::matrix_range<T> & obj)
    {
      return handle(obj.get());
    }

    template <typename T>
    viennacl::ocl::handle<cl_mem> handle(viennacl::matrix_range<T> const & obj)
    {
      return handle(obj.get());
    }

  } //namespace traits
} //namespace viennacl
    

#endif
