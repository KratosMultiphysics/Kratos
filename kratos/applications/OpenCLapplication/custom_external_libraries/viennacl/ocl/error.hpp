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
          return "CL_DEVICE_NOT_FOUND \n ViennaCL could not find a suitable device. Please check whether an OpenCL implementation is properly installed and a suitable device available. \n Aborting...";
        case CL_DEVICE_NOT_AVAILABLE:
          return "CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:
          return "CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
          return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:
          return "CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:
          return "CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
          return "CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
          return "CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
          return "CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
          return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
          return "CL_BUILD_PROGRAM_FAILURE";
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
          return "CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:
          return "CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
          return "CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:
          return "CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
          return "CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
          return "CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
          return "CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:
          return "CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
          return "CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:
          return "CL_INVALID_WORK_ITEM_SIZE";
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

      return "NO CL Error";
    }
  } //namespace ocl
} //namespace viennacl

#endif
