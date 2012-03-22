#ifndef VIENNACL_OCL_UTILS_HPP_
#define VIENNACL_OCL_UTILS_HPP_

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

/** @file viennacl/ocl/utils.hpp
    @brief Provides OpenCL-related utilities.
*/

#include <vector>
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/device.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief Ensures that double precision types are only allocated if it is supported by the device. If double precision is requested for a device not capable of providing that, a double_precision_not_provided_error is thrown.
     */
    template <typename ScalarType>
    struct DOUBLE_PRECISION_CHECKER
    {
      static void apply() {} 
    };
    
    template <>
    struct DOUBLE_PRECISION_CHECKER<double>
    {
      static void apply()
      {
        if (!viennacl::ocl::current_device().double_support())
          throw viennacl::ocl::double_precision_not_provided_error();
      }
    };
    
    

  } //ocl
} //viennacl
#endif
