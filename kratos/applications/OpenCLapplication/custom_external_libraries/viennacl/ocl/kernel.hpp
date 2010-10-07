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
        #ifdef VIENNACL_BUILD_INFO
        //std::cout << "Setting unsigned long kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        setArgument(pos, static_cast<unsigned int>(val));
      }

      template<class TYPE>
      void setArgument(unsigned int pos, TYPE val)
      {
        init();
        #ifdef VIENNACL_BUILD_INFO
        //std::cout << "Setting TYPE kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h.get(), pos, sizeof(TYPE), (void*)&val);
        CL_ERR_CHECK(err);
      }

      void setArgument(unsigned int pos, const handle<cl_mem> & val)
      {
        init();
        #ifdef VIENNACL_BUILD_INFO
        //std::cout << "Setting cl_mem kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h.get(), pos, sizeof(cl_mem), (void*)&val.get());
        CL_ERR_CHECK(err);
      }

      void setLocalBuffer(unsigned int pos, unsigned int size)
      {
        init();
        #ifdef VIENNACL_BUILD_INFO
        //std::cout << "Setting local buffer kernel argument at pos " << pos << " for kernel " << _name << std::endl;
        #endif
        cl_int err = clSetKernelArg(h.get(), pos, size, 0);
        CL_ERR_CHECK(err);
      }
      
      void start1D()
      {
        start1D(_work_groups * _work_items_per_group, _work_items_per_group);
      }

      void start1D(const size_t & global_work_size)
      {
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Starting kernel '" << _name << "'..." << std::endl;
        #endif
        size_t tmp = global_work_size;
        cl_int err = clEnqueueNDRangeKernel(device().queue().get(), h.get(), 1, NULL, &tmp, NULL, 0, NULL, NULL);
        CL_ERR_CHECK(err);
        //assert(err == CL_SUCCESS);
      }
      void start1D(const size_t & global_work_size, const size_t & local_work_size)
      {
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Starting kernel '" << _name << "'..." << std::endl;
        #endif
        size_t tmp_global = global_work_size;
        size_t tmp_local = local_work_size;
        cl_int err = clEnqueueNDRangeKernel(device().queue().get(), h.get(), 1, NULL, &tmp_global, &tmp_local, 0, NULL, NULL);

        if (err != CL_SUCCESS)  //if not successful, try to start with smaller work size
        {
          while (err != CL_SUCCESS && tmp_local > 1)
          {
            //std::cout << "Flushing queue, then enqueuing again with half the size..." << std::endl;
            tmp_global /= 2;
            tmp_local /= 2;

            err = clEnqueueNDRangeKernel(device().queue().get(), h.get(), 1, NULL, &tmp_global, &tmp_local, 0, NULL, NULL);
          }
          
          if (err != CL_SUCCESS)
          {
            //could not start kernel with any parameters
            std::cerr << "Could not start kernel '" << _name << "' with global work size " << tmp_global << " and local work size " << tmp_local << "." << std::endl;
            std::cerr << "Smaller work sizes also failed." << std::endl;
            CL_ERR_CHECK(err);
          }
          else
          {
            //remember parameters:
            _work_items_per_group = tmp_local;
            _work_groups = tmp_global / tmp_local;
            #ifdef VIENNACL_BUILD_INFO
            std::cout << "Kernel '" << _name << "' now uses global work size " << tmp_global << " and local work size " << tmp_local << "."  << std::endl;
            #endif
          }          
        }
        
        CL_ERR_CHECK(err);
        
        //assert(err == CL_SUCCESS);
      }

      const handle<cl_kernel> & get()
      {
        init();
        return h; 
      }

      size_t work_items_per_group() const { return _work_items_per_group; }
      size_t work_groups() const { return _work_groups; }

    private:
      
      //prevent messing around with kernels:
      kernel(const kernel & other) {}
      template <typename T>
      void operator=(T other) {}
      
      void create_kernel()
      {
        cl_int err;
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Building " << _name << std::endl;
        #endif
        h = clCreateKernel(_program.handle().get(), _name.c_str(), &err);
        
        if (err != CL_SUCCESS)
        {
          std::cerr << "Could not build kernel '" << _name << "'." << std::endl;
          #ifdef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK_ON_GPU
          if (_name.find("norm") != std::string::npos)
            std::cerr << "You seem to be using one of the functions norm_1, norm_2, norm_inf and index_norm_inf, which are not available on ATI GPUs! This also affects the GMRES solver." << std::endl;
          #endif
        }
        CL_ERR_CHECK(err);
        
        //set work group sizes and the like:
        if (viennacl::ocl::device().type() == CL_DEVICE_TYPE_CPU)
        {
          _work_items_per_group = 1;
          _work_groups = 2;   //some experiments on a Core 2 Quad 9550 show that two work groups give good results
        }
        else
        {
          _work_items_per_group = 128;
          _work_groups = 128;
        }
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
      size_t _work_items_per_group;
      size_t _work_groups;
      bool _init_done;
    };
    
  } //namespace ocl
} //namespace viennacl

#endif
