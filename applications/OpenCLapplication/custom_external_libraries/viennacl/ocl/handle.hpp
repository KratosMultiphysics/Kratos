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

#ifndef _VIENNACL_HANDLE_HPP_
#define _VIENNACL_HANDLE_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <assert.h>
#include <string>
#include <iostream>
#include "viennacl/ocl/error.hpp"

namespace viennacl
{
  namespace ocl
  {
    /** @brief Helper for OpenCL reference counting used by class handle.
    *   @tparam OCL_TYPE Must be one out of cl_mem, cl_program, cl_kernel, cl_command_queue and cl_context, otherwise a compile time error is thrown.
    */
    template<class OCL_TYPE>
    class handle_inc_dec_helper
    {
      typedef typename OCL_TYPE::ERROR_TEMPLATE_ARGUMENT_FOR_CLASS_INVALID   ErrorType;
    };
    
    
    //cl_mem:
    template <>
    struct handle_inc_dec_helper<cl_mem>
    {
      static void inc(cl_mem & something)
      {
        cl_int err = clRetainMemObject(something);
        VIENNACL_ERR_CHECK(err);
      }
      
      static void dec(cl_mem & something)
      {
        #ifndef __APPLE__
        cl_int err = clReleaseMemObject(something);
        VIENNACL_ERR_CHECK(err);
        #endif
      }
    };
    
    //cl_program:
    template <>
    struct handle_inc_dec_helper<cl_program>
    {
      static void inc(cl_program & something)
      {
        cl_int err = clRetainProgram(something);
        VIENNACL_ERR_CHECK(err);
      }
      
      static void dec(cl_program & something)
      {
        #ifndef __APPLE__
        cl_int err = clReleaseProgram(something);
        VIENNACL_ERR_CHECK(err);
        #endif
      }
    };
    
    //cl_kernel:
    template <>
    struct handle_inc_dec_helper<cl_kernel>
    {
      static void inc(cl_kernel & something)
      {
        cl_int err = clRetainKernel(something);
        VIENNACL_ERR_CHECK(err);
      }
      
      static void dec(cl_kernel & something)
      {
        #ifndef __APPLE__
        cl_int err = clReleaseKernel(something);
        VIENNACL_ERR_CHECK(err);
        #endif
      }
    };

    //cl_command_queue:
    template <>
    struct handle_inc_dec_helper<cl_command_queue>
    {
      static void inc(cl_command_queue & something)
      {
        cl_int err = clRetainCommandQueue(something);
        VIENNACL_ERR_CHECK(err);
      }
      
      static void dec(cl_command_queue & something)
      {
        #ifndef __APPLE__
        cl_int err = clReleaseCommandQueue(something);
        VIENNACL_ERR_CHECK(err);
        #endif
      }
    };
    
    //cl_context:
    template <>
    struct handle_inc_dec_helper<cl_context>
    {
      static void inc(cl_context & something)
      {
        cl_int err = clRetainContext(something);
        VIENNACL_ERR_CHECK(err);
      }
      
      static void dec(cl_context & something)
      {
        #ifndef __APPLE__
        cl_int err = clReleaseContext(something);
        VIENNACL_ERR_CHECK(err);
        #endif
      }
    };
    
    /** @brief Handle class the effectively represents a smart pointer for OpenCL handles */
    template<class OCL_TYPE>
    class handle
    {
    public:
      handle() : something(0) {}
      handle(const OCL_TYPE & _something) : something(_something) {}
      handle(const handle & h) : something(h.something) { if (something != 0) inc(); }
      ~handle() { if (something != 0) dec(); }
      handle & operator=(const handle & h)
      {
        if (something != 0) dec();
        something = h.something;
        inc();
        return *this;
      }
      handle & operator=(const OCL_TYPE & _something)
      {
        if (something != 0) dec();
        something = _something;
        return *this;
      }
      operator OCL_TYPE() const { return something; }
      //const OCL_TYPE & get() const { return something; }
      
      /** @brief Swaps the OpenCL handle of two handle objects */
      handle & swap(handle & other)
      {
        OCL_TYPE tmp = other.something;
        other.something = this->something;
        this->something = tmp;
        return *this;
      }
      
      /** @brief Manually increment the OpenCL reference count. Typically called automatically, but is necessary if user-supplied memory objects are wrapped. */
      void inc() { handle_inc_dec_helper<OCL_TYPE>::inc(something); };
      /** @brief Manually decrement the OpenCL reference count. Typically called automatically, but might be useful with user-supplied memory objects.  */
      void dec() { handle_inc_dec_helper<OCL_TYPE>::dec(something); };
    private:
      OCL_TYPE something;
    };

    
  } //namespace ocl
} //namespace viennacl

#endif
