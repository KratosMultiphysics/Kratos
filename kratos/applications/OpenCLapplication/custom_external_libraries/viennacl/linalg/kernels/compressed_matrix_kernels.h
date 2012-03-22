#ifndef _VIENNACL_COMPRESSED_MATRIX_KERNELS_HPP_
#define _VIENNACL_COMPRESSED_MATRIX_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/compressed_matrix_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct compressed_matrix;


    /////////////// single precision kernels //////////////// 
   template <>
   struct compressed_matrix<float, 4>
   {
    static std::string program_name()
    {
      return "f_compressed_matrix_4";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(compressed_matrix_align1_row_scaling_1);
        source.append(compressed_matrix_align1_lu_forward);
        source.append(compressed_matrix_align1_bicgstab_kernel2);
        source.append(compressed_matrix_align1_lu_backward);
        source.append(compressed_matrix_align4_vec_mul);
        source.append(compressed_matrix_align1_jacobi);
        source.append(compressed_matrix_align1_jacobi_precond);
        source.append(compressed_matrix_align1_bicgstab_kernel1);
        source.append(compressed_matrix_align1_row_scaling_2);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct compressed_matrix<float, 1>
   {
    static std::string program_name()
    {
      return "f_compressed_matrix_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(compressed_matrix_align1_row_scaling_1);
        source.append(compressed_matrix_align1_lu_forward);
        source.append(compressed_matrix_align1_bicgstab_kernel2);
        source.append(compressed_matrix_align1_lu_backward);
        source.append(compressed_matrix_align1_vec_mul);
        source.append(compressed_matrix_align1_jacobi);
        source.append(compressed_matrix_align1_jacobi_precond);
        source.append(compressed_matrix_align1_bicgstab_kernel1);
        source.append(compressed_matrix_align1_row_scaling_2);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct compressed_matrix<float, 8>
   {
    static std::string program_name()
    {
      return "f_compressed_matrix_8";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(compressed_matrix_align1_row_scaling_1);
        source.append(compressed_matrix_align1_lu_forward);
        source.append(compressed_matrix_align1_bicgstab_kernel2);
        source.append(compressed_matrix_align1_lu_backward);
        source.append(compressed_matrix_align8_vec_mul);
        source.append(compressed_matrix_align1_jacobi);
        source.append(compressed_matrix_align1_jacobi_precond);
        source.append(compressed_matrix_align1_bicgstab_kernel1);
        source.append(compressed_matrix_align1_row_scaling_2);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct compressed_matrix<double, 4>
   {
    static std::string program_name()
    {
      return "d_compressed_matrix_4";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<double>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        std::string fp64_ext = viennacl::ocl::current_device().double_support_extension();
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align4_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi_precond, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_2, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct compressed_matrix<double, 1>
   {
    static std::string program_name()
    {
      return "d_compressed_matrix_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<double>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        std::string fp64_ext = viennacl::ocl::current_device().double_support_extension();
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi_precond, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_2, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct compressed_matrix<double, 8>
   {
    static std::string program_name()
    {
      return "d_compressed_matrix_8";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<double>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        std::string fp64_ext = viennacl::ocl::current_device().double_support_extension();
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align8_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_jacobi_precond, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_bicgstab_kernel1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_row_scaling_2, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("row_scaling_1");
        prog_.add_kernel("lu_forward");
        prog_.add_kernel("bicgstab_kernel2");
        prog_.add_kernel("lu_backward");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("jacobi");
        prog_.add_kernel("jacobi_precond");
        prog_.add_kernel("bicgstab_kernel1");
        prog_.add_kernel("row_scaling_2");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
