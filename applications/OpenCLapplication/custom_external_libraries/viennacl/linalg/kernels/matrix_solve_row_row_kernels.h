#ifndef _VIENNACL_MATRIX_SOLVE_ROW_ROW_KERNELS_HPP_
#define _VIENNACL_MATRIX_SOLVE_ROW_ROW_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/matrix_solve_row_row_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct matrix_solve_row_row;


    /////////////// single precision kernels //////////////// 
   template <>
   struct matrix_solve_row_row<float, 1>
   {
    static std::string program_name()
    {
      return "f_matrix_solve_row_row_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(matrix_solve_row_row_align1_upper_trans_solve);
        source.append(matrix_solve_row_row_align1_trans_lower_trans_solve);
        source.append(matrix_solve_row_row_align1_trans_lower_solve);
        source.append(matrix_solve_row_row_align1_lower_solve);
        source.append(matrix_solve_row_row_align1_trans_unit_lower_solve);
        source.append(matrix_solve_row_row_align1_upper_solve);
        source.append(matrix_solve_row_row_align1_trans_upper_trans_solve);
        source.append(matrix_solve_row_row_align1_trans_unit_lower_trans_solve);
        source.append(matrix_solve_row_row_align1_unit_lower_trans_solve);
        source.append(matrix_solve_row_row_align1_unit_upper_solve);
        source.append(matrix_solve_row_row_align1_unit_lower_solve);
        source.append(matrix_solve_row_row_align1_lower_trans_solve);
        source.append(matrix_solve_row_row_align1_trans_unit_upper_solve);
        source.append(matrix_solve_row_row_align1_unit_upper_trans_solve);
        source.append(matrix_solve_row_row_align1_trans_upper_solve);
        source.append(matrix_solve_row_row_align1_trans_unit_upper_trans_solve);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("upper_trans_solve");
        prog_.add_kernel("trans_lower_trans_solve");
        prog_.add_kernel("trans_lower_solve");
        prog_.add_kernel("lower_solve");
        prog_.add_kernel("trans_unit_lower_solve");
        prog_.add_kernel("upper_solve");
        prog_.add_kernel("trans_upper_trans_solve");
        prog_.add_kernel("trans_unit_lower_trans_solve");
        prog_.add_kernel("unit_lower_trans_solve");
        prog_.add_kernel("unit_upper_solve");
        prog_.add_kernel("unit_lower_solve");
        prog_.add_kernel("lower_trans_solve");
        prog_.add_kernel("trans_unit_upper_solve");
        prog_.add_kernel("unit_upper_trans_solve");
        prog_.add_kernel("trans_upper_solve");
        prog_.add_kernel("trans_unit_upper_trans_solve");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct matrix_solve_row_row<double, 1>
   {
    static std::string program_name()
    {
      return "d_matrix_solve_row_row_1";
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
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_upper_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_lower_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_lower_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_lower_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_unit_lower_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_upper_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_upper_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_unit_lower_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_unit_lower_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_unit_upper_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_unit_lower_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_lower_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_unit_upper_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_unit_upper_trans_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_upper_solve, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_solve_row_row_align1_trans_unit_upper_trans_solve, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("upper_trans_solve");
        prog_.add_kernel("trans_lower_trans_solve");
        prog_.add_kernel("trans_lower_solve");
        prog_.add_kernel("lower_solve");
        prog_.add_kernel("trans_unit_lower_solve");
        prog_.add_kernel("upper_solve");
        prog_.add_kernel("trans_upper_trans_solve");
        prog_.add_kernel("trans_unit_lower_trans_solve");
        prog_.add_kernel("unit_lower_trans_solve");
        prog_.add_kernel("unit_upper_solve");
        prog_.add_kernel("unit_lower_solve");
        prog_.add_kernel("lower_trans_solve");
        prog_.add_kernel("trans_unit_upper_solve");
        prog_.add_kernel("unit_upper_trans_solve");
        prog_.add_kernel("trans_upper_solve");
        prog_.add_kernel("trans_unit_upper_trans_solve");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
