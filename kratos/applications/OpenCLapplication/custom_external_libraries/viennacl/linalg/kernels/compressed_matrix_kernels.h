#ifndef _VIENNACL_COMPRESSED_MATRIX_KERNELS_HPP_
#define _VIENNACL_COMPRESSED_MATRIX_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/kernels/compressed_matrix_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct compressed_matrix
   {
    static viennacl::ocl::kernel lu_forward;
    static viennacl::ocl::kernel lu_backward;
    static viennacl::ocl::kernel vec_mul;

    static void init();
   };

   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel compressed_matrix<T, ALIGNMENT>::lu_forward;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel compressed_matrix<T, ALIGNMENT>::lu_backward;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel compressed_matrix<T, ALIGNMENT>::vec_mul;


    /////////////// single precision kernels //////////////// 
    template <> void compressed_matrix<float, 4>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(compressed_matrix_align1_lu_forward);
       source.append(compressed_matrix_align1_lu_backward);
       source.append(compressed_matrix_align4_vec_mul);
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void compressed_matrix<float, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(compressed_matrix_align1_lu_forward);
       source.append(compressed_matrix_align1_lu_backward);
       source.append(compressed_matrix_align1_vec_mul);
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void compressed_matrix<float, 8>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(compressed_matrix_align1_lu_forward);
       source.append(compressed_matrix_align1_lu_backward);
       source.append(compressed_matrix_align8_vec_mul);
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }


    /////////////// double precision kernels //////////////// 
    template <> void compressed_matrix<double, 4>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align4_vec_mul));
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void compressed_matrix<double, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_vec_mul));
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void compressed_matrix<double, 8>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_forward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align1_lu_backward));
       source.append(viennacl::tools::make_double_kernel(compressed_matrix_align8_vec_mul));
       prog.create(source);

       lu_forward.prepareInit("lu_forward", prog);
       lu_backward.prepareInit("lu_backward", prog);
       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
