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

#ifndef _VIENNACL_PROGRAM_HPP_
#define _VIENNACL_PROGRAM_HPP_


#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/device.hpp"

namespace viennacl
{
  namespace ocl
  {
    class program
    {
    public:
      program() {}
      program(const handle<cl_program> & _prog) : h(_prog) {}
      
      program & operator=(const handle<cl_program> & _prog) { h = _prog; return *this; }

      const handle<cl_program> & get() const { return h; }

      void create(const char * source)
      {
        std::string source2(source);
        create(source2);
      }
      
      void create(std::string & source)
      {
        const char * source_text = source.c_str();
        size_t source_size = source.size();
        
        //std::cout << "Building source: " << std::endl;
        //std::cout << _source << std::endl;

        cl_int err;
        h = clCreateProgramWithSource(device().context(), 1, (const char **)&source_text, &source_size, &err);
        //if (err != CL_SUCCESS)
          //std::cerr << _source << std::endl;
        //assert(err == CL_SUCCESS);
        CL_ERR_CHECK(err);
        //Timer watch;
        //double cpu_time = watch.get();
        //err = clBuildProgram(prog, 0, NULL, "-cl-mad-enable", NULL, NULL);
        err = clBuildProgram(h, 0, NULL, NULL, NULL, NULL);
        //std::cout << "Build time: " << watch.get() - cpu_time << std::endl;
        //assert(err == CL_SUCCESS);
        //if (err != CL_SUCCESS)
        //  std::cerr << _source << std::endl;
        
        //additional debug info
        #ifdef VCL_BUILD_INFO
        char buffer[1024];
        cl_build_status status;
        clGetProgramBuildInfo(h, device().id(), CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, NULL);
        clGetProgramBuildInfo(h, device().id(), CL_PROGRAM_BUILD_LOG, sizeof(char)*1024, &buffer, NULL);
        std::cout << "Build Scalar: Err = " << err << " Status = " << status << std::endl;
        std::cout << "Log: " << buffer << std::endl;
        #endif

        CL_ERR_CHECK(err);
      }

    private:
      handle<cl_program> h;
    };
  } //namespace ocl
} //namespace viennacl


#endif
