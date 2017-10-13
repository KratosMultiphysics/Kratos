#ifndef VIENNACL_RAND_UNIFORM_HPP_
#define VIENNACL_RAND_UNIFORM_HPP_

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

/** @file   viennacl/rand/uniform.hpp
    @brief  Unused: Generation of uniformly distributed random numbers. */

/** \cond */

namespace viennacl{

namespace rand{

struct uniform_tag{
    uniform_tag(unsigned int _a = 0, unsigned int _b = 1) : a(_a), b(_b){ }
    float a;
    float b;
};

template<class ScalarType>
struct buffer_dumper<ScalarType, uniform_tag>{
  static void dump(viennacl::backend::mem_handle const & buff, uniform_tag tag, cl_uint start, cl_uint size){
    viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::rand<ScalarType,1>::program_name(),"dump_uniform");
    k.global_work_size(0, viennacl::tools::align_to_multiple<unsigned int>(size,k.local_work_size(0)));
    viennacl::ocl::enqueue(k(buff.opencl_handle(), start, size, cl_float(tag.a), cl_float(tag.b) , cl_uint(time(0))));
  }
};



}

}

/** \endcond */

#endif
