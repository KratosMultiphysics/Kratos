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

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_ERROR_HPP_
#define _VIENNACL_ERROR_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <string>
#include <iostream>

namespace viennacl
{
  namespace ocl
  {
    void checkError(cl_int err, const std::string & file, const std::string & func, int line);
    std::string getErrorString(cl_int err);

    #define CL_ERR_CHECK(err) viennacl::ocl::checkError(err, __FILE__, __FUNCTION__, __LINE__);

    
    
    void checkError(cl_int err, const std::string & file, const std::string & func, int line)
    {
      if (err != CL_SUCCESS)
      {
        std::cerr << "Error " << err  << " in function " << func << " ( "<< file << ":" << line << " ) " << getErrorString(err) << std::endl;
        //std::cin.ignore(0, '\n');
        //std::cin.get();
        std::cerr << "If you think that this is a bug in ViennaCL, please report it at viennacl-support@lists.sourceforge.net and supply at least the following information:" << std::endl;
        std::cerr << " * Operating System" << std::endl;
        std::cerr << " * Which OpenCL implementation (NVIDIA, ATI, etc.)" << std::endl;
        std::cerr << " * ViennaCL version" << std::endl;
        std::cerr << "Many thanks in advance!" << std::endl;
        std::cerr << "Aborting (forcing seg-fault)..." << std::endl;
        *((int*)(0)) = 0;
        //exit(0);
      }
    }
    
    std::string getErrorString(cl_int err)
    {
      switch (err)
      {
        case CL_SUCCESS:
          return "CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:
          return "CL_DEVICE_NOT_FOUND \n ViennaCL could not find a suitable device. Please check whether an OpenCL implementation is properly installed and a suitable device available.";
        case CL_DEVICE_NOT_AVAILABLE:
          return "CL_DEVICE_NOT_AVAILABLE \n ViennaCL could not use the compute device because it is not available.";
        case CL_COMPILER_NOT_AVAILABLE:
          return "CL_COMPILER_NOT_AVAILABLE \n Your OpenCL framework does not provide an OpenCL compiler. Unfortunately, ViennaCL cannot be used without such a compiler.";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
          return "CL_MEM_OBJECT_ALLOCATION_FAILURE \n ViennaCL could not allocate memory on the device. Most likely the device simply ran out of memory.";
        case CL_OUT_OF_RESOURCES:
          return "CL_OUT_OF_RESOURCES \n ViennaCL tried to launch a compute kernel, but the device does not provide enough resources. Try changing the global and local work item sizes.";
        case CL_OUT_OF_HOST_MEMORY:
          return "CL_OUT_OF_HOST_MEMORY \n The host ran out of memory (usually CPU RAM). Please try again on smaller problems.";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
          return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
          return "CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
          return "CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
          return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
          return "CL_BUILD_PROGRAM_FAILURE \n The OpenCL compiler encountered an error during the compilation of ViennaCL sources. This is most likely a bug in ViennaCL.";
        case CL_MAP_FAILURE:
          return "CL_MAP_FAILURE";

        case CL_INVALID_VALUE:
          return "CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:
          return "CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:
          return "CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:
          return "CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:
          return "CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:
          return "CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:
          return "CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:
          return "CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:
          return "CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
          return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:
          return "CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:
          return "CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:
          return "CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:
          return "CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:
          return "CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:
          return "CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:
          return "CL_INVALID_KERNEL_NAME \n The supplied kernel name is invalid. If you have written your own OpenCL kernel, please check that the correct kernel name is used in the initalization of the kernel object.";
        case CL_INVALID_KERNEL_DEFINITION:
          return "CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
          return "CL_INVALID_KERNEL \n The kernel is invalid. Did you ";
        case CL_INVALID_ARG_INDEX:
          return "CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
          return "CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
          return "CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
          return "CL_INVALID_KERNEL_ARGS \n The supplied kernel arguments do not fit the kernel parameter list. If you have written your own OpenCL kernel, please check that the correct kernel arguments are set in the appropriate order.";
        case CL_INVALID_WORK_DIMENSION:
          return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
          return "CL_INVALID_WORK_GROUP_SIZE \n The supplied work group size is invalid. If you set this value manually, please reconsider your choice.";
        case CL_INVALID_WORK_ITEM_SIZE:
          return "CL_INVALID_WORK_ITEM_SIZE \n The work item size is invalid. If you set this value manually, please reconsider your choice.";
        case CL_INVALID_GLOBAL_OFFSET:
          return "CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:
          return "CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:
          return "CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:
          return "CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:
          return "CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:
          return "CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:
          return "CL_INVALID_MIP_LEVEL";
        //case CL_INVALID_GLOBAL_WORK_SIZE:
        //  return "CL_INVALID_GLOBAL_WORK_SIZE";
      }

      return "NO CL Error \n ViennaCL encountered an unknown OpenCL error. In some cases, this might be due to an invalid global work size, but it can also be due to several compilation errors.";
    }
  } //namespace ocl
} //namespace viennacl

#endif
