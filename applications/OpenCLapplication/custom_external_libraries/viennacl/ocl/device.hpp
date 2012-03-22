#ifndef VIENNACL_OCL_DEVICE_HPP_
#define VIENNACL_OCL_DEVICE_HPP_

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

/** @file device.hpp
    @brief Represents an OpenCL device within ViennaCL
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include<stdio.h>

#include <vector>
#include <string>
#include <sstream>
#include <assert.h>
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/error.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief A class representing a compute device (e.g. a GPU)
    *
    */
    class device
    {
      public:
        explicit device() : device_(0) {}
        
        explicit device(cl_device_id dev) : device_(dev)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_DEVICE)
          std::cout << "ViennaCL: Creating device object (CTOR with cl_device_id)" << std::endl;
          #endif
          init(dev);
        }
        
        device(const device & other)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_DEVICE)
          std::cout << "ViennaCL: Creating device object (Copy CTOR)" << std::endl;
          #endif
          device_ = other.device_;
          init(device_);
        }
        
        /** @brief Initializes the class from a given device ID */
        void init(cl_device_id dev)
        {
          cl_int err;

          //query a little bit of info:
          err = clGetDeviceInfo(dev, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &max_work_group_size_, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clGetDeviceInfo(dev, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &compute_units_, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clGetDeviceInfo(dev, CL_DEVICE_TYPE, sizeof(cl_device_type), &type_, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clGetDeviceInfo(dev, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_memory_, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clGetDeviceInfo(dev, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &max_memory_alloc_, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clGetDeviceInfo(dev, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_memory_, NULL);
          VIENNACL_ERR_CHECK(err);
        }

        /** @brief Returns true if the device supports double precision */
        bool double_support() const
        { 
          char buffer[1024];
          bool ret = false;
          
          //get extensions and search for double precision
          clGetDeviceInfo(device_, CL_DEVICE_EXTENSIONS, sizeof(char)*1024, buffer, NULL);
          std::string extensions(buffer);
          if (extensions.find("cl_khr_fp64") != std::string::npos
              || extensions.find("cl_amd_fp64") != std::string::npos)
          {
            ret = true;
          }
          
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_DEVICE)
          std::cout << "ViennaCL: Device extensions: " << std::endl;
          std::cout << extensions << std::endl;
          if (ret)
            std::cout << "ViennaCL: Device " << name() << " supports double precision." << std::endl;
          else
            std::cout << "ViennaCL: No double precision for device " << name() << "." << std::endl;
          #endif
          
          return ret;
        }
        
        std::string double_support_extension() const
        {
          char buffer[1024];
          clGetDeviceInfo(device_, CL_DEVICE_EXTENSIONS, sizeof(char)*1024, buffer, NULL);
          std::string extensions(buffer);
          
          if (extensions.find("cl_amd_fp64") != std::string::npos) //AMD extension
            return "cl_amd_fp64";
          
          if (extensions.find("cl_khr_fp64") != std::string::npos) //Khronos-certified standard extension for double precision
            return "cl_khr_fp64";
          
          return "";
        }
        
        /** @brief Returns the OpenCL device id */
        cl_device_id id() const
        {
          assert(device_ != 0);
          return device_;
        }
        
        /** @brief Returns the device name */
        std::string name() const
        {
          std::ostringstream oss;        
          char buffer[1024]; 
          cl_int err;          
          err = clGetDeviceInfo(device_, CL_DEVICE_NAME, sizeof(char)*1024, &buffer, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << buffer;
          return oss.str();          
        }
        
        /** @brief Returns the driver version */
        std::string driver_version() const
        {
          std::ostringstream oss;
          char buffer[1024]; buffer[0] = 0;
          cl_int err;          
          err = clGetDeviceInfo(device_, CL_DRIVER_VERSION, sizeof(char)*1024, buffer, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << buffer;
          return oss.str();          
        }        
        
        /** @brief Returns the number of compute units on the device */
        cl_uint max_compute_units() const
        {
          return compute_units_;
        }
        
        /** @brief Returns the maximum work group size for the device*/
        size_t max_workgroup_size() const
        {
          return max_work_group_size_;
        }                        

        /** @brief Returns the global memory for the device*/
        cl_ulong global_memory() const
        {
          return global_memory_;
        }           

        /** @brief Returns the local memory for the device*/
        cl_ulong local_memory() const
        {
          return local_memory_;
        }       

        /** @brief Returns the maximum allocable memory for the device*/
        cl_ulong max_allocable_memory() const
        {
          return max_memory_alloc_;
        }           
        
        /** @brief Returns an info string with a few properties of the device */
        std::string info() const
        {
          std::ostringstream oss;
          char buffer[1024]; buffer[0] = 0;
          cl_int err;
          cl_uint vendor_id;
          cl_ulong local_mem_size;
          cl_ulong global_mem_size;
          
          err = clGetDeviceInfo(device_, CL_DEVICE_VENDOR_ID, sizeof(cl_uint), &vendor_id, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << "CL Device Vendor ID: " << vendor_id << std::endl;

          err = clGetDeviceInfo(device_, CL_DEVICE_NAME, sizeof(char)*1024, buffer, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << "CL Device Name: " << buffer << std::endl;

          err = clGetDeviceInfo(device_, CL_DRIVER_VERSION, sizeof(char)*1024, buffer, NULL);
          VIENNACL_ERR_CHECK(err);
          std::string test = buffer;
          oss << "CL Driver Version: " << test << std::endl;

          oss << "--------------------------------" << std::endl;
          
          oss << "CL Device Max Compute Units: " << compute_units_ << std::endl;

  //         err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(char)*1024, buffer, NULL);
  //         CL_ERR_CHECK(err);
  //         oss << "CL Device Max Work Item Dimensions: " << buffer << std::endl;
  // 
  //         err = clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(char)*1024, buffer, NULL);
  //         CL_ERR_CHECK(err);
  //         oss << "CL Device Max Work Item Sizes: " << buffer << std::endl;

          oss << "CL Device Max Work Group Size: " << max_work_group_size_ << std::endl;

          err = clGetDeviceInfo(device_, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_mem_size, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << "CL Device Global Mem Size: " << global_mem_size << std::endl;
          
          err = clGetDeviceInfo(device_, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem_size, NULL);
          VIENNACL_ERR_CHECK(err);
          oss << "CL Device Local Mem Size: " << local_mem_size << std::endl;
          
          //return info string:
          std::string ret(oss.str());
          return ret;
        }
        
        size_t max_work_group_size() const { return max_work_group_size_; }
        cl_uint compute_units() const { return compute_units_; }
        cl_device_type type() const { return type_; }
        
        bool operator==(device const & other) const
        {
          return device_ == other.device_;
        }

        bool operator==(cl_device_id other) const
        {
          return device_ == other;
        }

      private:
        
        cl_device_id    device_;
        size_t          max_work_group_size_;
        cl_uint         compute_units_;
        cl_device_type  type_; //device type
        cl_ulong        max_memory_alloc_;
        cl_ulong        global_memory_;
        cl_ulong        local_memory_;
    };

  } //namespace ocl
} //namespace viennacl

#endif
