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

/** @file ocl/forwards.h
    @brief This file provides the forward declarations for the OpenCL layer of ViennaCL
*/

#ifndef _VIENNACL_OCL_FORWARDS_H_
#define _VIENNACL_OCL_FORWARDS_H_

#define VIENNACL_OCL_MAX_DEVICE_NUM  8

#include <stddef.h>

namespace viennacl
{
  namespace ocl
  {
    //device type tags (cf. OpenCL standard)
    struct gpu_tag {};
    struct cpu_tag {};
    struct accelerator_tag {};
    struct default_tag {};
    
    
    class kernel;
    class device;
    class command_queue;
    class context;
    class program;

    template<class OCL_TYPE>
    class handle;

    template <typename KernelType>
    void enqueue(KernelType & k, viennacl::ocl::command_queue const & queue);
    
    inline viennacl::ocl::context & current_context();
    inline viennacl::ocl::device const & current_device();
  }
} //namespace viennacl

#endif

/*@}*/
