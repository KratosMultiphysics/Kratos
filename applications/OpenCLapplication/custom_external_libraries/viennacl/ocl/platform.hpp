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

#ifndef _VIENNACL_PLATFORM_HPP_
#define _VIENNACL_PLATFORM_HPP_

/** @file platform.hpp
    @brief Implements a OpenCL platform within ViennaCL
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <vector>
#include "viennacl/ocl/forwards.h"
#include "viennacl/ocl/device.hpp"

namespace viennacl
{
  namespace ocl
  {
    class platform
    {
      
      public:
        platform()
        {
          cl_int err;
          cl_uint num_platforms;
          #if defined(VIENNACL_DEBUG_ALL)
          std::cout << "ViennaCL: Getting platform..." << std::endl;
          #endif
          err = clGetPlatformIDs(1, &id_, &num_platforms);
          VIENNACL_ERR_CHECK(err);
          assert(num_platforms > 0 && "ViennaCL: ERROR: No platform found!");          
        }
        
        cl_platform_id id() const
        {
          return id_;
        }
        
        /** @brief Returns an information string */
        std::string info() const
        {
          char buffer[1024];
          cl_int err;
          err = clGetPlatformInfo(id_, CL_PLATFORM_VENDOR, 1024 * sizeof(char), buffer, NULL);
          VIENNACL_ERR_CHECK(err);
          
          std::stringstream ss;
          ss << buffer;
          
          return ss.str();
        }
        
        //////////////////// get device //////////////////
        /** @brief Returns the available devices of the supplied device type */
        std::vector<device> devices(cl_device_type dtype = CL_DEVICE_TYPE_DEFAULT)
        {
          cl_int err;
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_DEVICE)
          std::cout << "ViennaCL: Querying devices available at current platform." << std::endl;
          #endif
          cl_device_id device_ids[VIENNACL_OCL_MAX_DEVICE_NUM];
          cl_uint num_devices;
          err = clGetDeviceIDs(id_, dtype, VIENNACL_OCL_MAX_DEVICE_NUM, device_ids, &num_devices);
          if (err == CL_DEVICE_NOT_FOUND && dtype == CL_DEVICE_TYPE_DEFAULT)
          {
            //workaround for ATI Stream SDK v2.3: No CPUs detected with default device type:
            err = clGetDeviceIDs(id_, CL_DEVICE_TYPE_CPU, VIENNACL_OCL_MAX_DEVICE_NUM, device_ids, &num_devices);
          }
          
          VIENNACL_ERR_CHECK(err);
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_DEVICE)
          std::cout << "ViennaCL: Found " << num_devices << " devices." << std::endl;
          #endif
          
          assert(num_devices > 0 && "Error in viennacl::ocl::platform::devices(): No OpenCL devices available!");
          std::vector<device> devices;
          
          for (cl_uint i=0; i<num_devices; ++i)
            devices.push_back(device(device_ids[i]));

          return devices;
        }
        
      private:
        cl_platform_id id_;
    };
    
  }
}

#endif
