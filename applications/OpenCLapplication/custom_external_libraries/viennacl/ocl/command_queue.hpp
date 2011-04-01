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

#ifndef _VIENNACL_COMMAND_QUEUE_HPP_
#define _VIENNACL_COMMAND_QUEUE_HPP_

/** @file command_queue.hpp
    @brief Implementations of command queue representations
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <vector>
#include <string>
#include <sstream>
#include "viennacl/ocl/context.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief A class representing a command queue
    *
    */
    class command_queue
    {
      public:
        command_queue() {};
        command_queue(viennacl::ocl::handle<cl_command_queue> h, cl_device_id dev) : handle_(h) {}
        
        //Copy constructor:
        command_queue(command_queue const & other)
        {
          handle_ = other.handle_;
        }

        //assignment operator:
        command_queue & operator=(command_queue const & other)
        {
          handle_ = other.handle_;
          return *this;
        }
        
        /** @brief Waits until all kernels in the queue have finished their execution */
        void finish() const
        {
          clFinish(handle_);
        }
        
        /** @brief Waits until all kernels in the queue have started their execution */
        void flush() const
        {
          clFlush(handle_);
        }

        viennacl::ocl::handle<cl_command_queue> const & handle() const { return handle_; }

      private:
        
        viennacl::ocl::handle<cl_command_queue> handle_;
    };

 
    
  } //namespace ocl
} //namespace viennacl

#endif
