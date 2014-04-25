#ifndef VIENNACL_RAND_GAUSSIAN_HPP_
#define VIENNACL_RAND_GAUSSIAN_HPP_

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

#include "viennacl/backend/mem_handle.hpp"
#include "viennacl/rand/utils.hpp"

/** @file   viennacl/rand/gaussian.hpp
    @brief  Unused: Generation of random numbers with a Gaussian normal distribution */

/** \cond */


namespace viennacl{

namespace rand{

struct gaussian_tag{
    gaussian_tag(float _mu = 0, float _sigma = 1) : mu(_mu), sigma(_sigma){ }
    float mu;
    float sigma;
};

template<class ScalarType>
struct buffer_dumper<ScalarType, gaussian_tag>{
    static void dump(viennacl::backend::mem_handle const & buff, gaussian_tag tag, cl_uint start, cl_uint size){
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::rand<ScalarType,1>::program_name(),"dump_gaussian");
      k.global_work_size(0, viennacl::tools::align_to_multiple<unsigned int>(size/2,k.local_work_size(0)));
      viennacl::ocl::enqueue(k(buff.opencl_handle(), start, size, cl_float(tag.mu), cl_float(tag.sigma) , cl_uint(time(0))));
    }
};

}

}

/** \endcond */
#endif
