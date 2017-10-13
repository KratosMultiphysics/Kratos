#ifndef VIENNACL_SVD_KERNELS_HPP_
#define VIENNACL_SVD_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/svd_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct svd;


    /////////////// single precision kernels //////////////// 
   template <>
   struct svd<float, 1>
   {
    static std::string program_name()
    {
      return "f_svd_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(svd_align1_bidiag_pack);
        source.append(svd_align1_house_col);
        source.append(svd_align1_inverse_signs);
        source.append(svd_align1_transpose_inplace);
        source.append(svd_align1_house_row);
        source.append(svd_align1_copy_col);
        source.append(svd_align1_copy_row);
        source.append(svd_align1_givens_prev);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("bidiag_pack");
        prog_.add_kernel("house_col");
        prog_.add_kernel("inverse_signs");
        prog_.add_kernel("transpose_inplace");
        prog_.add_kernel("house_row");
        prog_.add_kernel("copy_col");
        prog_.add_kernel("copy_row");
        prog_.add_kernel("givens_prev");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct svd<double, 1>
   {
    static std::string program_name()
    {
      return "d_svd_1";
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
        source.append(viennacl::tools::make_double_kernel(svd_align1_bidiag_pack, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_house_col, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_inverse_signs, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_transpose_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_house_row, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_copy_col, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_copy_row, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(svd_align1_givens_prev, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("bidiag_pack");
        prog_.add_kernel("house_col");
        prog_.add_kernel("inverse_signs");
        prog_.add_kernel("transpose_inplace");
        prog_.add_kernel("house_row");
        prog_.add_kernel("copy_col");
        prog_.add_kernel("copy_row");
        prog_.add_kernel("givens_prev");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
