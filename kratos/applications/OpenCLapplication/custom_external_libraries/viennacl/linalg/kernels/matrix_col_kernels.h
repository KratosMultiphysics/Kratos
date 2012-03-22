#ifndef _VIENNACL_MATRIX_COL_KERNELS_HPP_
#define _VIENNACL_MATRIX_COL_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/matrix_col_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct matrix_col;


    /////////////// single precision kernels //////////////// 
   template <>
   struct matrix_col<float, 1>
   {
    static std::string program_name()
    {
      return "f_matrix_col_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(matrix_col_align1_fft_direct);
        source.append(matrix_col_align1_fft_reorder);
        source.append(matrix_col_align1_trans_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_unit_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_rank1_update);
        source.append(matrix_col_align1_inplace_divide);
        source.append(matrix_col_align1_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_lu_factorize);
        source.append(matrix_col_align1_trans_vec_mul);
        source.append(matrix_col_align1_clear);
        source.append(matrix_col_align1_scaled_rank1_update);
        source.append(matrix_col_align1_inplace_mult);
        source.append(matrix_col_align1_cpu_inplace_mult);
        source.append(matrix_col_align1_vec_mul);
        source.append(matrix_col_align1_inplace_sub);
        source.append(matrix_col_align1_fft_radix2_local);
        source.append(matrix_col_align1_inplace_add);
        source.append(matrix_col_align1_fft_radix2);
        source.append(matrix_col_align1_unit_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_unit_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_unit_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_assign);
        source.append(matrix_col_align1_add);
        source.append(matrix_col_align1_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_sub);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_direct");
        prog_.add_kernel("fft_reorder");
        prog_.add_kernel("trans_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("rank1_update");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("lower_triangular_substitute_inplace");
        prog_.add_kernel("lu_factorize");
        prog_.add_kernel("trans_vec_mul");
        prog_.add_kernel("clear");
        prog_.add_kernel("scaled_rank1_update");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("fft_radix2_local");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("fft_radix2");
        prog_.add_kernel("unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_upper_triangular_substitute_inplace");
        prog_.add_kernel("assign");
        prog_.add_kernel("add");
        prog_.add_kernel("upper_triangular_substitute_inplace");
        prog_.add_kernel("sub");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct matrix_col<float, 16>
   {
    static std::string program_name()
    {
      return "f_matrix_col_16";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(matrix_col_align1_fft_direct);
        source.append(matrix_col_align1_fft_reorder);
        source.append(matrix_col_align1_trans_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_unit_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_rank1_update);
        source.append(matrix_col_align1_inplace_divide);
        source.append(matrix_col_align1_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_lu_factorize);
        source.append(matrix_col_align1_trans_vec_mul);
        source.append(matrix_col_align1_clear);
        source.append(matrix_col_align1_scaled_rank1_update);
        source.append(matrix_col_align1_inplace_mult);
        source.append(matrix_col_align1_cpu_inplace_mult);
        source.append(matrix_col_align1_vec_mul);
        source.append(matrix_col_align1_inplace_sub);
        source.append(matrix_col_align1_fft_radix2_local);
        source.append(matrix_col_align1_inplace_add);
        source.append(matrix_col_align1_fft_radix2);
        source.append(matrix_col_align1_unit_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_unit_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_unit_lower_triangular_substitute_inplace);
        source.append(matrix_col_align1_trans_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_assign);
        source.append(matrix_col_align1_add);
        source.append(matrix_col_align1_upper_triangular_substitute_inplace);
        source.append(matrix_col_align1_sub);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_direct");
        prog_.add_kernel("fft_reorder");
        prog_.add_kernel("trans_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("rank1_update");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("lower_triangular_substitute_inplace");
        prog_.add_kernel("lu_factorize");
        prog_.add_kernel("trans_vec_mul");
        prog_.add_kernel("clear");
        prog_.add_kernel("scaled_rank1_update");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("fft_radix2_local");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("fft_radix2");
        prog_.add_kernel("unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_upper_triangular_substitute_inplace");
        prog_.add_kernel("assign");
        prog_.add_kernel("add");
        prog_.add_kernel("upper_triangular_substitute_inplace");
        prog_.add_kernel("sub");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct matrix_col<double, 1>
   {
    static std::string program_name()
    {
      return "d_matrix_col_1";
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
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_direct, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_reorder, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_unit_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_rank1_update, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_divide, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_lu_factorize, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_clear, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_scaled_rank1_update, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_mult, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_cpu_inplace_mult, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_sub, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_radix2_local, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_add, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_radix2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_unit_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_unit_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_unit_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_assign, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_add, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_sub, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_direct");
        prog_.add_kernel("fft_reorder");
        prog_.add_kernel("trans_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("rank1_update");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("lower_triangular_substitute_inplace");
        prog_.add_kernel("lu_factorize");
        prog_.add_kernel("trans_vec_mul");
        prog_.add_kernel("clear");
        prog_.add_kernel("scaled_rank1_update");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("fft_radix2_local");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("fft_radix2");
        prog_.add_kernel("unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_upper_triangular_substitute_inplace");
        prog_.add_kernel("assign");
        prog_.add_kernel("add");
        prog_.add_kernel("upper_triangular_substitute_inplace");
        prog_.add_kernel("sub");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct

   template <>
   struct matrix_col<double, 16>
   {
    static std::string program_name()
    {
      return "d_matrix_col_16";
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
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_direct, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_reorder, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_unit_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_rank1_update, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_divide, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_lu_factorize, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_clear, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_scaled_rank1_update, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_mult, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_cpu_inplace_mult, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_vec_mul, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_sub, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_radix2_local, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_inplace_add, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_fft_radix2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_unit_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_unit_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_unit_lower_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_trans_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_assign, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_add, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_upper_triangular_substitute_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(matrix_col_align1_sub, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_direct");
        prog_.add_kernel("fft_reorder");
        prog_.add_kernel("trans_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("rank1_update");
        prog_.add_kernel("inplace_divide");
        prog_.add_kernel("lower_triangular_substitute_inplace");
        prog_.add_kernel("lu_factorize");
        prog_.add_kernel("trans_vec_mul");
        prog_.add_kernel("clear");
        prog_.add_kernel("scaled_rank1_update");
        prog_.add_kernel("inplace_mult");
        prog_.add_kernel("cpu_inplace_mult");
        prog_.add_kernel("vec_mul");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("fft_radix2_local");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("fft_radix2");
        prog_.add_kernel("unit_upper_triangular_substitute_inplace");
        prog_.add_kernel("trans_unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("unit_lower_triangular_substitute_inplace");
        prog_.add_kernel("trans_upper_triangular_substitute_inplace");
        prog_.add_kernel("assign");
        prog_.add_kernel("add");
        prog_.add_kernel("upper_triangular_substitute_inplace");
        prog_.add_kernel("sub");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
