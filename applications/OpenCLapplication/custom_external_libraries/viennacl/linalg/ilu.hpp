#ifndef VIENNACL_LINALG_ILU_HPP_
#define VIENNACL_LINALG_ILU_HPP_

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

/** @file ilu.hpp
    @brief Implementations of incomplete factorization preconditioners. Convenience header file.
*/

#include "viennacl/linalg/detail/ilu/ilut.hpp"
#include "viennacl/linalg/detail/ilu/ilu0.hpp"
#include "viennacl/linalg/detail/ilu/host_block_ilu.hpp"
//#include "viennacl/linalg/detail/ilu/opencl_block_ilu.hpp" //to be enabled in 1.4.0

#endif



