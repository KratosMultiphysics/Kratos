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
    template<class OCL_TYPE>
    class handle
    {
    public:
      handle() : something(0) {}
      handle(const OCL_TYPE & _something) : something(_something) {}
      handle(const handle & h) : something(h.something) { inc(); }
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
      const OCL_TYPE & get() const { return something; }
    private:
      void inc();
      void dec();
      OCL_TYPE something;
    };

    
    ////////////////// specializations ///////////////////
    
    template<>
    void handle<cl_mem>::dec()
    {
      cl_int err = clReleaseMemObject(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_MEM_OBJECT);
      //cout << "DEC cl_mem HANDLE" << endl;
    }

    template<>
    void handle<cl_mem>::inc()
    {
      cl_int err = clRetainMemObject(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_MEM_OBJECT);
      //cout << "INC cl_mem HANDLE" << endl;
    }

    template<>
    void handle<cl_program>::dec()
    {
      cl_int err = clReleaseProgram(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_PROGRAM);
      //cout << "DEC cl_program HANDLE" << endl;
    }

    template<>
    void handle<cl_program>::inc()
    {
      cl_int err = clRetainProgram(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_PROGRAM);
      //cout << "INC cl_program HANDLE" << endl;
    }

    template<>
    void handle<cl_kernel>::dec()
    {
      cl_int err = clReleaseKernel(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_KERNEL);
      //cout << "DEC cl_kernel HANDLE" << endl;
    }

    template<>
    void handle<cl_kernel>::inc()
    {
      cl_int err = clRetainKernel(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_KERNEL);
      //cout << "INC cl_kernel HANDLE" << endl;
    }

    template<>
    void handle<cl_command_queue>::dec()
    {
      cl_int err = clReleaseCommandQueue(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_COMMAND_QUEUE);
      //cout << "DEC cl_command_queue HANDLE" << endl;
    }

    template<>
    void handle<cl_command_queue>::inc()
    {
      cl_int err = clReleaseCommandQueue(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_COMMAND_QUEUE);
      //cout << "INC cl_command_queue HANDLE" << endl;
    }

    template<>
    void handle<cl_context>::dec()
    {
      cl_int err = clReleaseContext(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_CONTEXT);
      //cout << "DEC cl_context HANDLE" << endl;
    }

    template<>
    void handle<cl_context>::inc()
    {
      cl_int err = clReleaseContext(something);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS || err == CL_INVALID_CONTEXT);
      //cout << "INC cl_context HANDLE" << endl;
    }


  /*
    void kernel::start1D(const unsigned int & global_work_size)
    {
      cl_int err = clEnqueueNDRangeKernel(device().cq(), h, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS);
    }

    void kernel::start1D(const unsigned int & global_work_size, const unsigned int & local_work_size)
    {
      cl_int err = clEnqueueNDRangeKernel(device().cq(), h, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);
      CL_ERR_CHECK(err);
      //assert(err == CL_SUCCESS);
    } */
  //  */
  
  } //namespace ocl
} //namespace viennacl

#endif
