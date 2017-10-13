#ifndef VIENNACL_RAND_UTILS_HPP_
#define VIENNACL_RAND_UTILS_HPP_

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
/** @file   viennacl/rand/utils.hpp
    @brief  Unused: Helper functionality random number generation. */

/** \cond */

#ifdef VIENNACL_WITH_OPENCL
#include "viennacl/linalg/kernels/rand_kernels.h"

namespace viennacl{

namespace rand{


template<class SCALARTYPE, class DISTRIBUTION>
struct random_matrix_t{
    typedef size_t size_type;
    random_matrix_t(size_type _size1, unsigned int _size2, DISTRIBUTION const & _distribution) : size1(_size1), size2(_size2), distribution(_distribution){
        #ifdef VIENNACL_WITH_OPENCL
        viennacl::linalg::kernels::rand<SCALARTYPE,1>::init();
        #endif
    }
    size_type size1;
    size_type size2;
    DISTRIBUTION distribution;
};


template<class SCALARTYPE, class DISTRIBUTION>
struct random_vector_t{
    typedef size_t size_type;
    random_vector_t(size_type _size, DISTRIBUTION const & _distribution) : size(_size), distribution(_distribution){
        #ifdef VIENNACL_WITH_OPENCL
        viennacl::linalg::kernels::rand<SCALARTYPE,1>::init();
        #endif
    }
    size_type size;
    DISTRIBUTION distribution;
};

template<class ScalarType, class Distribution>
struct buffer_dumper;


}

}

#endif

/** \endcond */

#endif
