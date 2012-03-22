#ifndef VIENNACL_META_PREDICATE_HPP_
#define VIENNACL_META_PREDICATE_HPP_

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

/** @file predicate.hpp
    @brief All the predicates used within ViennaCL. Checks for expressions to be vectors, etc.
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"

namespace viennacl
{
    //
    // is_cpu_scalar: checks for float or double
    //
    template <typename T>
    struct is_cpu_scalar
    {
      enum { value = false };
    };
  
    template <>
    struct is_cpu_scalar<float>
    {
      enum { value = true };
    };

    template <>
    struct is_cpu_scalar<double>
    {
      enum { value = true };
    };
    
    //
    // is_scalar: checks for viennacl::scalar
    //
    template <typename T>
    struct is_scalar
    {
      enum { value = false };
    };
  
    template <typename T>
    struct is_scalar<viennacl::scalar<T> >
    {
      enum { value = true };
    };
  
    //
    // is_vector
    //
    template <typename T>
    struct is_vector
    {
      enum { value = false };
    };

    template <typename ScalarType, unsigned int ALIGNMENT>
    struct is_vector<viennacl::vector<ScalarType, ALIGNMENT> >
    {
      enum { value = true };
    };

    template <typename T>
    struct is_vector<viennacl::vector_range<T> >
    {
      enum { value = true };
    };
    
    
    
    //
    // is_matrix
    //
    template <typename T>
    struct is_matrix
    {
      enum { value = false };
    };

    template <typename ScalarType, typename F, unsigned int ALIGNMENT>
    struct is_matrix<viennacl::matrix<ScalarType, F, ALIGNMENT> >
    {
      enum { value = true };
    };

    template <typename T>
    struct is_matrix<viennacl::matrix_range<T> >
    {
      enum { value = true };
    };
    

} //namespace viennacl
    

#endif
