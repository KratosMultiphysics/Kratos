#ifndef _VIENNACL_SPAI_KERNELS_HPP_
#define _VIENNACL_SPAI_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/spai_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct spai;


    /////////////// single precision kernels //////////////// 
   template <>
   struct spai<float, 1>
   {
    static std::string program_name()
    {
      return "f_spai_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(spai_align1_block_qr);
        source.append(spai_align1_block_qr_assembly);
        source.append(spai_align1_block_q_mult);
        source.append(spai_align1_block_qr_assembly_1);
        source.append(spai_align1_block_bv_assembly);
        source.append(spai_align1_block_r_assembly);
        source.append(spai_align1_assemble_blocks);
        source.append(spai_align1_block_least_squares);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("block_qr");
        prog_.add_kernel("block_qr_assembly");
        prog_.add_kernel("block_q_mult");
        prog_.add_kernel("block_qr_assembly_1");
        prog_.add_kernel("block_bv_assembly");
        prog_.add_kernel("block_r_assembly");
        prog_.add_kernel("assemble_blocks");
        prog_.add_kernel("block_least_squares");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct spai<double, 1>
   {
    static std::string program_name()
    {
      return "d_spai_1";
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
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_qr, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_qr_assembly, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_q_mult, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_qr_assembly_1, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_bv_assembly, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_r_assembly, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_assemble_blocks, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(spai_align1_block_least_squares, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("block_qr");
        prog_.add_kernel("block_qr_assembly");
        prog_.add_kernel("block_q_mult");
        prog_.add_kernel("block_qr_assembly_1");
        prog_.add_kernel("block_bv_assembly");
        prog_.add_kernel("block_r_assembly");
        prog_.add_kernel("assemble_blocks");
        prog_.add_kernel("block_least_squares");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
