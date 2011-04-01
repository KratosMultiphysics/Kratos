#ifndef _VIENNACL_VECTOR_KERNELS_HPP_
#define _VIENNACL_VECTOR_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/vector_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct vector;


    /////////////// single precision kernels //////////////// 
   template <>
   struct vector<float, 16>
   {
    static std::string program_name()
    {
      return "f_vector_16";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(vector_align1_plane_rotation);
        source.append(vector_align4_inplace_mul_add);
        source.append(vector_align1_index_norm_inf);
        source.append(vector_align1_norm_inf);
        source.append(vector_align16_cpu_mult);
        source.append(vector_align16_sub);
        source.append(vector_align4_cpu_inplace_mul_add);
        source.append(vector_align1_sum);
        source.append(vector_align16_inplace_mult);
        source.append(vector_align4_inplace_div_sub);
        source.append(vector_align1_norm_2);
        source.append(vector_align4_cpu_mul_add);
        source.append(vector_align4_inplace_div_add);
        source.append(vector_align1_diag_precond);
        source.append(vector_align4_inplace_mul_sub);
        source.append(vector_align1_cpu_inplace_mult);
        source.append(vector_align16_divide);
        source.append(vector_align4_mul_add);
        source.append(vector_align16_add);
        source.append(vector_align16_inplace_divide);
        source.append(vector_align1_swap);
        source.append(vector_align1_mul_sub);
        source.append(vector_align1_vmax);
        source.append(vector_align1_inner_prod);
        source.append(vector_align16_inplace_sub);
        source.append(vector_align16_inplace_add);
        source.append(vector_align16_mult);
        source.append(vector_align1_sqrt_sum);
        source.append(vector_align1_clear);
        source.append(vector_align1_norm_1);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct vector<float, 1>
   {
    static std::string program_name()
    {
      return "f_vector_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(vector_align1_plane_rotation);
        source.append(vector_align1_inplace_mul_add);
        source.append(vector_align1_index_norm_inf);
        source.append(vector_align1_norm_inf);
        source.append(vector_align1_cpu_mult);
        source.append(vector_align1_sub);
        source.append(vector_align1_cpu_inplace_mul_add);
        source.append(vector_align1_sum);
        source.append(vector_align1_inplace_mult);
        source.append(vector_align1_inplace_div_sub);
        source.append(vector_align1_norm_2);
        source.append(vector_align1_cpu_mul_add);
        source.append(vector_align1_inplace_div_add);
        source.append(vector_align1_diag_precond);
        source.append(vector_align1_inplace_mul_sub);
        source.append(vector_align1_cpu_inplace_mult);
        source.append(vector_align1_divide);
        source.append(vector_align1_mul_add);
        source.append(vector_align1_add);
        source.append(vector_align1_inplace_divide);
        source.append(vector_align1_swap);
        source.append(vector_align1_mul_sub);
        source.append(vector_align1_vmax);
        source.append(vector_align1_inner_prod);
        source.append(vector_align1_inplace_sub);
        source.append(vector_align1_inplace_add);
        source.append(vector_align1_mult);
        source.append(vector_align1_sqrt_sum);
        source.append(vector_align1_clear);
        source.append(vector_align1_norm_1);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct vector<float, 4>
   {
    static std::string program_name()
    {
      return "f_vector_4";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(vector_align1_plane_rotation);
        source.append(vector_align4_inplace_mul_add);
        source.append(vector_align1_index_norm_inf);
        source.append(vector_align1_norm_inf);
        source.append(vector_align1_cpu_mult);
        source.append(vector_align1_sub);
        source.append(vector_align4_cpu_inplace_mul_add);
        source.append(vector_align1_sum);
        source.append(vector_align1_inplace_mult);
        source.append(vector_align4_inplace_div_sub);
        source.append(vector_align1_norm_2);
        source.append(vector_align4_cpu_mul_add);
        source.append(vector_align4_inplace_div_add);
        source.append(vector_align1_diag_precond);
        source.append(vector_align4_inplace_mul_sub);
        source.append(vector_align1_cpu_inplace_mult);
        source.append(vector_align1_divide);
        source.append(vector_align4_mul_add);
        source.append(vector_align1_add);
        source.append(vector_align1_inplace_divide);
        source.append(vector_align1_swap);
        source.append(vector_align1_mul_sub);
        source.append(vector_align1_vmax);
        source.append(vector_align1_inner_prod);
        source.append(vector_align1_inplace_sub);
        source.append(vector_align1_inplace_add);
        source.append(vector_align1_mult);
        source.append(vector_align1_sqrt_sum);
        source.append(vector_align1_clear);
        source.append(vector_align1_norm_1);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct vector<double, 16>
   {
    static std::string program_name()
    {
      return "d_vector_16";
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
        source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_cpu_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_2, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_diag_precond, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_swap, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_vmax, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align16_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sqrt_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_clear, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct vector<double, 1>
   {
    static std::string program_name()
    {
      return "d_vector_1";
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
        source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_div_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_2, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_div_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_diag_precond, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_swap, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_vmax, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sqrt_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_clear, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct vector<double, 4>
   {
    static std::string program_name()
    {
      return "d_vector_4";
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
        source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_inplace_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_2, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_diag_precond, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align4_mul_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_swap, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_vmax, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_mult, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_sqrt_sum, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_clear, pf_info));
        source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("plane_rotation");
        prog_.add_kernel("inplace_mul_add");
        prog_.add_kernel("index_norm_inf");
        prog_.add_kernel("norm_inf");
        prog_.add_kernel("cpu_mult");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_mul_add");
        prog_.add_kernel("sum");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("inplace_div_sub");
        prog_.add_kernel("norm_2");
        prog_.add_kernel("cpu_mul_add");
        prog_.add_kernel("inplace_div_add");
        prog_.add_kernel("diag_precond");
        prog_.add_kernel("inplace_mul_sub");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("divide");
        prog_.add_kernel("mul_add");
        prog_.add_kernel("add");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("swap");
        prog_.add_kernel("mul_sub");
        prog_.add_kernel("vmax");
        prog_.add_kernel("inner_prod");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("mult");
        prog_.add_kernel("sqrt_sum");
        prog_.add_kernel("clear");
        prog_.add_kernel("norm_1");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
