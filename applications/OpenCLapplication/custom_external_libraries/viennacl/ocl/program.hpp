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

#ifndef _VIENNACL_PROGRAM_HPP_
#define _VIENNACL_PROGRAM_HPP_

/** @file program.hpp
    @brief Implements an OpenCL program class for ViennaCL
*/

#include <string>
#include <vector>
#include "viennacl/ocl/forwards.h"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"

namespace viennacl
{
  namespace ocl
  {
    class program
    {
      friend class kernel;
      
      typedef std::vector<viennacl::ocl::kernel>    KernelContainer;
      
    public:
      program() {}
      program(viennacl::ocl::handle<cl_program> const & h, std::string const & prog_name = std::string()) : handle_(h), name_(prog_name) {}
      
      program(program const & other)
      {
        handle_ = other.handle_;
        name_ = other.name_;
        kernels_ = other.kernels_;
      }
      
      viennacl::ocl::program & operator=(const program & other)
      {
        handle_ = other.handle_;
        name_ = other.name_;
        kernels_ = other.kernels_;
        return *this;
      }

      std::string const & name() const { return name_; }
      
      /** @brief Adds a kernel to the program */
      viennacl::ocl::kernel & add_kernel(std::string const & kernel_name)
      {
        viennacl::ocl::kernel temp(handle_, kernel_name);
        kernels_.push_back(temp);
        return kernels_.back();
      }
      
      /** @brief Returns the kernel with the provided name */
      viennacl::ocl::kernel & get_kernel(std::string const & name)
      {
        //std::cout << "Requiring kernel " << name << " from program " << name_ << std::endl;
        for (KernelContainer::iterator it = kernels_.begin();
              it != kernels_.end();
             ++it)
        {
          if (it->name() == name)
            return *it;
        }
        std::cerr << "ViennaCL: FATAL ERROR: Could not find kernel '" << name << "'" << std::endl;
        assert(!"Kernel not found");
        return kernels_[0];  //return a defined object
      }

    private:
      const viennacl::ocl::handle<cl_program> & handle() const { return handle_; }
      
      viennacl::ocl::handle<cl_program> handle_;
      std::string name_;
      KernelContainer kernels_;
    };
  } //namespace ocl
} //namespace viennacl


#endif
