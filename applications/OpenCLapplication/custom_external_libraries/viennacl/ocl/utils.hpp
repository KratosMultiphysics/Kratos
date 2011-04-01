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

#ifndef _VIENNACL_OCL_UTILS_HPP_
#define _VIENNACL_OCL_UTILS_HPP_

/** @file backend.hpp
    @brief Implementations of the OpenCL backend, where all contexts are stored in.
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
