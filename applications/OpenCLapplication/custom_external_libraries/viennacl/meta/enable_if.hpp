#ifndef VIENNACL_META_ENABLE_IF_HPP_
#define VIENNACL_META_ENABLE_IF_HPP_

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

/** @file enable_if.hpp
    @brief Simple enable-if variant that uses the SFINAE pattern
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"


#include <vector>
#include <map>

namespace viennacl
{
    /** @brief Simple enable-if variant that uses the SFINAE pattern */
    template <bool b, class T = void> 
    struct enable_if
    {
      typedef T   type;
    };

    template <class T> 
    struct enable_if<false, T> {};

} //namespace viennacl
    

#endif
