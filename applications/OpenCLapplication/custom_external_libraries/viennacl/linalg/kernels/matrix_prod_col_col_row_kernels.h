#ifndef _VIENNACL_MATRIX_PROD_COL_COL_ROW_KERNELS_HPP_
#define _VIENNACL_MATRIX_PROD_COL_COL_ROW_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/matrix_prod_col_col_row_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct matrix_prod_col_col_row;


    /////////////// single precision kernels //////////////// 
   template <>
   struct matrix_prod_col_col_row<float, 1>
   {
    static std::string program_name()
    {
      return "f_matrix_prod_col_col_row_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(matrix_prod_col_col_row_align1_prod_TA);
        source.append(matrix_prod_col_col_row_align1_prod_TT);
        source.append(matrix_prod_col_col_row_align1_prod_AT);
        source.append(matrix_prod_col_col_row_align1_prod_AA);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("prod_TA");
        prog_.add_kernel("prod_TT");
        prog_.add_kernel("prod_AT");
        prog_.add_kernel("prod_AA");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct matrix_prod_col_col_row<double, 1>
   {
    static std::string program_name()
    {
      return "d_matrix_prod_col_col_row_1";
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
        source.append(viennacl::tools::make_double_kernel(matrix_prod_col_col_row_align1_prod_TA, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_prod_col_col_row_align1_prod_TT, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_prod_col_col_row_align1_prod_AT, pf_info));
        source.append(viennacl::tools::make_double_kernel(matrix_prod_col_col_row_align1_prod_AA, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("prod_TA");
        prog_.add_kernel("prod_TT");
        prog_.add_kernel("prod_AT");
        prog_.add_kernel("prod_AA");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
