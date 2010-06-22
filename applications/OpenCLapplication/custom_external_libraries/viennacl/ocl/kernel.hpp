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

#ifndef _VIENNACL_KERNEL_HPP_
#define _VIENNACL_KERNEL_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/program.hpp"
#include "viennacl/ocl/device.hpp"


namespace viennacl
{
  namespace ocl
  {
    class kernel
    {
    public:
      kernel() : _init_done(false) {}
      //kernel(const handle<cl_kernel> & _kernel_handle) : h(_kernel_handle) {}
      
      void prepareInit(const char * name, program const & prog)
      {
        if (!_init_done)
        {
          _name = name;
          _program = prog;
        }
        else
          std::cout << "WARNING: Kernel " << _name << " already initialized, cannot assign new name '" << name << "' and program." << std::endl;
      }
      
      void setArgument(unsigned int pos, unsigned long val)
      {
        init();
        #ifdef VCL_BUILD_INFO
        //std::cout << "Setting unsigned long kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        setArgument(pos, static_cast<unsigned int>(val));
      }

      template<class TYPE>
      void setArgument(unsigned int pos, TYPE val)
      {
        init();
        #ifdef VCL_BUILD_INFO
        //std::cout << "Setting TYPE kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h, pos, sizeof(TYPE), (void*)&val);
        CL_ERR_CHECK(err);
      }

      void setArgument(unsigned int pos, const handle<cl_mem> & val)
      {
        init();
        #ifdef VCL_BUILD_INFO
        //std::cout << "Setting cl_mem kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h, pos, sizeof(cl_mem), (void*)&val.get());
        CL_ERR_CHECK(err);
      }

      void setLocalBuffer(unsigned int pos, unsigned int size)
      {
        init();
        #ifdef VCL_BUILD_INFO
        //std::cout << "Setting local buffer kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h, pos, size, 0);
        CL_ERR_CHECK(err);
      }
      
      void start1D()
      {
        start1D(device().work_groups() * device().work_items_per_group(),
                device().work_items_per_group());
      }

      void start1D(const unsigned int & global_work_size)
      {
        #ifdef VCL_BUILD_INFO
        std::cout << "Starting kernel '" << _name << "'..." << std::endl;
        #endif
        size_t tmp = global_work_size;
        cl_int err = clEnqueueNDRangeKernel(device().queue(), h, 1, NULL, &tmp, NULL, 0, NULL, NULL);
        CL_ERR_CHECK(err);
        //assert(err == CL_SUCCESS);
      }
      void start1D(const unsigned int & global_work_size, const unsigned int & local_work_size)
      {
        #ifdef VCL_BUILD_INFO
        std::cout << "Starting kernel '" << _name << "'..." << std::endl;
        #endif
        size_t tmp_global = global_work_size;
        size_t tmp_local = local_work_size;
        cl_int err = clEnqueueNDRangeKernel(device().queue(), h, 1, NULL, &tmp_global, &tmp_local, 0, NULL, NULL);
        CL_ERR_CHECK(err);
        /*
        if (err != CL_SUCCESS)
        {
          //std::cout << "Flushing queue, then enqueuing again with half the size..." << std::endl;
          unsigned int new_global_work_size = global_work_size / 2;
          unsigned int new_local_work_size = local_work_size / 2;
          err = clEnqueueNDRangeKernel(device().command_queue(), h, 1, NULL, &new_global_work_size, &new_local_work_size, 0, NULL, NULL);
          if (err != CL_SUCCESS)
          {
            std::cout << "Finishing queue, then enqueuing kernel '" << _name << "' again..." << std::endl;
            clFlush(device().command_queue());
            err = clEnqueueNDRangeKernel(device().command_queue(), h, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
            CL_ERR_CHECK(err);
          }
        } */
        //assert(err == CL_SUCCESS);
      }

      const handle<cl_kernel> & get()
      {
        init();
        return h; 
      }

    private:
      
      //prevent messing around with kernels:
      kernel(const kernel & other) {}
      template <typename T>
      void operator=(T other) {}
      
      void create_kernel()
      {
        cl_int err;
        #ifdef VCL_BUILD_INFO
        std::cout << "Building " << _name << std::endl;
        #endif
        h = clCreateKernel(_program.get(), _name.c_str(), &err);
        CL_ERR_CHECK(err);
      }

      void init()
      {
        if (!_init_done)
        {
          create_kernel();
          _init_done = true;
        }
      }

      handle<cl_kernel> h;
      program _program;
      std::string _name;
      std::string _source;
      bool _init_done;
    };
    
  } //namespace ocl
} //namespace viennacl

#endif
