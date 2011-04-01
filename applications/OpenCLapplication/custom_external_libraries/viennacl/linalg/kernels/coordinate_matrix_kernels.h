#ifndef _VIENNACL_COORDINATE_MATRIX_KERNELS_HPP_
#define _VIENNACL_COORDINATE_MATRIX_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/coordinate_matrix_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct coordinate_matrix;


    /////////////// single precision kernels //////////////// 
   template <>
   struct coordinate_matrix<float, 128>
   {
    static std::string program_name()
    {
      return "f_coordinate_matrix_128";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(coordinate_matrix_align1_vec_mul);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("vec_mul");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct coordinate_matrix<float, 1>
   {
    static std::string program_name()
    {
      return "f_coordinate_matrix_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(coordinate_matrix_align1_vec_mul);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("vec_mul");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct coordinate_matrix<double, 128>
   {
    static std::string program_name()
    {
      return "d_coordinate_matrix_128";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<double>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        viennacl::ocl::platform pf;
        std::string pf_info(pf.info());
        source.append(viennacl::tools::make_double_kernel(coordinate_matrix_align1_vec_mul, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("vec_mul");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct coordinate_matrix<double, 1>
   {
    static std::string program_name()
    {
      return "d_coordinate_matrix_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<double>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        viennacl::ocl::platform pf;
        std::string pf_info(pf.info());
        source.append(viennacl::tools::make_double_kernel(coordinate_matrix_align1_vec_mul, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("vec_mul");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
