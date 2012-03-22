#ifndef _VIENNACL_FFT_KERNELS_HPP_
#define _VIENNACL_FFT_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/platform.hpp"
#include "viennacl/ocl/utils.hpp"
#include "viennacl/linalg/kernels/fft_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct fft;


    /////////////// single precision kernels //////////////// 
   template <>
   struct fft<float, 1>
   {
    static std::string program_name()
    {
      return "f_fft_1";
    }
    static void init()
    {
      viennacl::ocl::DOUBLE_PRECISION_CHECKER<float>::apply();
      static std::map<cl_context, bool> init_done;
      viennacl::ocl::context & context_ = viennacl::ocl::current_context();
      if (!init_done[context_.handle().get()])
      {
        std::string source;
        source.append(fft_align1_fft_mult_vec);
        source.append(fft_align1_bluestein_pre);
        source.append(fft_align1_zero2);
        source.append(fft_align1_bluestein_post);
        source.append(fft_align1_complex_to_real);
        source.append(fft_align1_transpose_inplace);
        source.append(fft_align1_transpose);
        source.append(fft_align1_reverse_inplace);
        source.append(fft_align1_real_to_complex);
        source.append(fft_align1_fft_div_vec_scalar);
        source.append(fft_align1_vandermonde_prod);
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_mult_vec");
        prog_.add_kernel("bluestein_pre");
        prog_.add_kernel("zero2");
        prog_.add_kernel("bluestein_post");
        prog_.add_kernel("complex_to_real");
        prog_.add_kernel("transpose_inplace");
        prog_.add_kernel("transpose");
        prog_.add_kernel("reverse_inplace");
        prog_.add_kernel("real_to_complex");
        prog_.add_kernel("fft_div_vec_scalar");
        prog_.add_kernel("vandermonde_prod");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct



    /////////////// double precision kernels //////////////// 
   template <>
   struct fft<double, 1>
   {
    static std::string program_name()
    {
      return "d_fft_1";
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
        source.append(viennacl::tools::make_double_kernel(fft_align1_fft_mult_vec, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_bluestein_pre, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_zero2, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_bluestein_post, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_complex_to_real, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_transpose_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_transpose, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_reverse_inplace, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_real_to_complex, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_fft_div_vec_scalar, fp64_ext));
        source.append(viennacl::tools::make_double_kernel(fft_align1_vandermonde_prod, fp64_ext));
        std::string prog_name = program_name();
        #ifdef VIENNACL_BUILD_INFO
        std::cout << "Creating program " << prog_name << std::endl;
        #endif
        context_.add_program(source, prog_name);
        viennacl::ocl::program & prog_ = context_.get_program(prog_name);
        prog_.add_kernel("fft_mult_vec");
        prog_.add_kernel("bluestein_pre");
        prog_.add_kernel("zero2");
        prog_.add_kernel("bluestein_post");
        prog_.add_kernel("complex_to_real");
        prog_.add_kernel("transpose_inplace");
        prog_.add_kernel("transpose");
        prog_.add_kernel("reverse_inplace");
        prog_.add_kernel("real_to_complex");
        prog_.add_kernel("fft_div_vec_scalar");
        prog_.add_kernel("vandermonde_prod");
        init_done[context_.handle().get()] = true;
       } //if
     } //init
    }; // struct


  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
