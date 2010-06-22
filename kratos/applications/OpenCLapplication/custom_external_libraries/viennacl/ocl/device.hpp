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

#ifndef _VIENNACL_DEVICE_HPP_
#define _VIENNACL_DEVICE_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <vector>
#include <string>
#include <sstream>
#include <assert.h>
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/error.hpp"
#include "viennacl/tools/tools.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief A class representing a compute device (e.g. a GPU)
    *
    */
    class Device
    {
      friend Device & device();

    public:

      handle<cl_mem> createMemory(cl_mem_flags flags, unsigned int size, void * ptr = NULL)
      {
        if (ptr) flags |= CL_MEM_COPY_HOST_PTR;
        cl_int err;
        handle<cl_mem> mem = clCreateBuffer(_context, flags, size, ptr, &err);
        //assert(err == CL_SUCCESS);
        CL_ERR_CHECK(err);
        return mem;
      }

      template < typename SCALARTYPE, typename A, template <typename, typename> class VectorType >
      handle<cl_mem> createMemory(cl_mem_flags flags, const VectorType<SCALARTYPE, A> & _buffer)
      {
        return createMemory(flags, sizeof(SCALARTYPE) * _buffer.size(), (void*)&_buffer[0]);
      }


      handle<cl_command_queue> & queue() { return default_command_queue; }

      bool double_support() const { return supported_double; }
      
      handle<cl_context> & context() { return _context; }
      
      cl_device_id id() { return device; }
      
      std::string info() const
      {
        std::ostringstream oss;
        char buffer[1024]; buffer[0] = 0;
        cl_int err;
        cl_uint vendor_id;
        cl_ulong local_mem_size;
        cl_ulong global_mem_size;
        
        err = clGetDeviceInfo(device, CL_DEVICE_VENDOR_ID, sizeof(cl_uint), &vendor_id, NULL);
        CL_ERR_CHECK(err);
        oss << "CL Device Vendor ID: " << vendor_id << std::endl;

        err = clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(char)*1024, buffer, NULL);
        CL_ERR_CHECK(err);
        oss << "CL Device Name: " << buffer << std::endl;

        err = clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(char)*1024, buffer, NULL);
        CL_ERR_CHECK(err);
        std::string test = buffer;
        oss << "CL Driver Version: " << test << std::endl;

        oss << "--------------------------------" << std::endl;
        
        oss << "CL Device Max Compute Units: " << _compute_units << std::endl;

//         err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(char)*1024, buffer, NULL);
//         CL_ERR_CHECK(err);
//         oss << "CL Device Max Work Item Dimensions: " << buffer << std::endl;
// 
//         err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(char)*1024, buffer, NULL);
//         CL_ERR_CHECK(err);
//         oss << "CL Device Max Work Item Sizes: " << buffer << std::endl;

        oss << "CL Device Max Work Group Size: " << _max_work_group_size << std::endl;

        err = clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_mem_size, NULL);
        CL_ERR_CHECK(err);
        oss << "CL Device Global Mem Size: " << global_mem_size << std::endl;
        
        err = clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL);
        CL_ERR_CHECK(err);
        oss << "CL Device Local Mem Size: " << local_mem_size << std::endl;
        
        //return info string:
        std::string ret(oss.str());
        return ret;
      }
      
      unsigned int max_work_group_size() const { return _max_work_group_size; }
      cl_uint compute_units() const { return _compute_units; }
      
      unsigned int type() const { return _type; }
      unsigned int work_items_per_group() const { return _work_items_per_group; }
      unsigned int work_groups() const { return _work_groups; }

    private:

      Device()
      {
        cl_int err;
        char buffer[1024];

        err = clGetPlatformIDs(1, &platform, NULL);
        CL_ERR_CHECK(err);
        
        // Try to find a suitable GPU
        err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
        if (err == CL_DEVICE_NOT_FOUND)
        {
          //No suitable GPU found. Try to find a suitable CPU (available via ATI Stream SDK)
          err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_CPU, 1, &device, NULL);
          
          //set default parameters:
          _type = CL_DEVICE_TYPE_CPU;
          _work_items_per_group = 1;
          _work_groups = 2;   //some experiments on a Core 2 Quad 9550 show that two work groups give good results
        }
        else
        {
          _type = CL_DEVICE_TYPE_GPU;
          _work_items_per_group = 128;
          _work_groups = 128;
        }
        CL_ERR_CHECK(err);
        
        //create OpenCL context for device
        _context = clCreateContext(0, 1, &device, NULL, NULL, &err);
        CL_ERR_CHECK(err);
        default_command_queue = clCreateCommandQueue(_context, device, 0, &err);
        CL_ERR_CHECK(err);

        //get extensions and search for double precision
        clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(char)*1024, buffer, NULL);
        std::string extensions(buffer);
        supported_double = false;
        if (extensions.find("cl_khr_fp64") != std::string::npos
        #ifdef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK
          || extensions.find("cl_amd_fp64") != std::string::npos
        #endif
          )
        {
          //std::cout << "Double precision supported" << std::endl;
          supported_double = true;
        }
        
        //query a little bit of info:
        err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &_max_work_group_size, NULL);
        CL_ERR_CHECK(err);
        err = clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &_compute_units, NULL);
        CL_ERR_CHECK(err);
      }

      bool supported_double;

      cl_platform_id platform;
      cl_device_id device;
      handle<cl_context> _context;
      handle<cl_command_queue> default_command_queue;
      
      size_t _max_work_group_size;
      cl_uint _compute_units;
      
      unsigned int _work_items_per_group;
      unsigned int _work_groups;
      unsigned int _type; //device type
    };

    
    /** @brief Singleton pattern. Returns a compute device.
    */
    Device & device()
    {
      static Device dev;
      return dev;
    }
    
    /** @brief Flush the device command queue, i.e. wait until all kernels in the command queue have started.
    */
    void flush()
    {
      clFlush(device().queue());
    }
    
    /** @brief Blocks until all kernels on the device have finished
    */
    void finish()
    {
      clFinish(device().queue());
    }
    
  } //namespace ocl
} //namespace viennacl

#endif
