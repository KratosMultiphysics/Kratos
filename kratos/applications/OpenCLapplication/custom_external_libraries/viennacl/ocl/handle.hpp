#ifndef VIENNACL_OCL_HANDLE_HPP_
#define VIENNACL_OCL_HANDLE_HPP_

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

/** @file viennacl/ocl/handle.hpp
    @brief Implementation of a smart-pointer-like class for handling OpenCL handles.
*/

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
      handle() : h_(0) {}
      handle(const OCL_TYPE & _something) : h_(_something) {}
      handle(const handle & other) : h_(other.h_) { if (h_ != 0) inc(); }
      ~handle() { if (h_ != 0) dec(); }
      handle & operator=(const handle & other)
      {
        if (h_ != 0) 
          dec();
        h_ = other.h_;
        inc();
        return *this;
      }
      handle & operator=(const OCL_TYPE & _something)
      {
        if (h_ != 0) dec();
        h_ = _something;
        return *this;
      }
      
      /** @brief Implicit conversion to the plain OpenCL handle. DEPRECATED and will be removed some time in the future. */
      operator OCL_TYPE() const { return h_; }
      
      const OCL_TYPE & get() const { return h_; }
      
      
      
      /** @brief Swaps the OpenCL handle of two handle objects */
      handle & swap(handle & other)
      {
        OCL_TYPE tmp = other.h_;
        other.h_ = this->h_;
        this->h_ = tmp;
        return *this;
      }
      
      /** @brief Manually increment the OpenCL reference count. Typically called automatically, but is necessary if user-supplied memory objects are wrapped. */
      void inc() { handle_inc_dec_helper<OCL_TYPE>::inc(h_); };
      /** @brief Manually decrement the OpenCL reference count. Typically called automatically, but might be useful with user-supplied memory objects.  */
      void dec() { handle_inc_dec_helper<OCL_TYPE>::dec(h_); };
    private:
      OCL_TYPE h_;
    };

    
  } //namespace ocl
} //namespace viennacl

#endif
