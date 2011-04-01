#ifndef _VIENNACL_SCALAR_KERNELS_HPP_
#define _VIENNACL_SCALAR_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/scalar_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct scalar;


    /////////////// single precision kernels //////////////// 
   template <>
   struct scalar<float, 1>
   {
    static std::string program_name()
    {
      return "f_scalar_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle()])
      {
        std::string source;
        source.append(scalar_align1_inplace_mul);
        source.append(scalar_align1_cpu_inplace_mul);
        source.append(scalar_align1_sub);
        source.append(scalar_align1_cpu_inplace_sub);
        source.append(scalar_align1_mul);
        source.append(scalar_align1_inplace_div);
        source.append(scalar_align1_cpu_mul);
        source.append(scalar_align1_divide);
        source.append(scalar_align1_add);
        source.append(scalar_align1_cpu_inplace_add);
        source.append(scalar_align1_cpu_div);
        source.append(scalar_align1_cpu_inplace_div);
        source.append(scalar_align1_cpu_add);
        source.append(scalar_align1_inplace_sub);
        source.append(scalar_align1_inplace_add);
        source.append(scalar_align1_cpu_sub);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("inplace_mul");
        prog_.add_kernel("cpu_inplace_mul");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_sub");
        prog_.add_kernel("mul");
        prog_.add_kernel("inplace_div");
        prog_.add_kernel("cpu_mul");
        prog_.add_kernel("divide");
        prog_.add_kernel("add");
        prog_.add_kernel("cpu_inplace_add");
        prog_.add_kernel("cpu_div");
        prog_.add_kernel("cpu_inplace_div");
        prog_.add_kernel("cpu_add");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("cpu_sub");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct scalar<double, 1>
   {
    static std::string program_name()
    {
      return "d_scalar_1";
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
        source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_mul, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_mul, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_mul, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_div, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_mul, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_divide, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_div, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_div, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_sub, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_add, pf_info));
        source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_sub, pf_info));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("inplace_mul");
        prog_.add_kernel("cpu_inplace_mul");
        prog_.add_kernel("sub");
        prog_.add_kernel("cpu_inplace_sub");
        prog_.add_kernel("mul");
        prog_.add_kernel("inplace_div");
        prog_.add_kernel("cpu_mul");
        prog_.add_kernel("divide");
        prog_.add_kernel("add");
        prog_.add_kernel("cpu_inplace_add");
        prog_.add_kernel("cpu_div");
        prog_.add_kernel("cpu_inplace_div");
        prog_.add_kernel("cpu_add");
        prog_.add_kernel("inplace_sub");
        prog_.add_kernel("inplace_add");
        prog_.add_kernel("cpu_sub");
        init_done[context_.handle()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
