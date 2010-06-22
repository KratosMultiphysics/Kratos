#ifndef _VIENNACL_COORDINATE_MATRIX_KERNELS_HPP_
#define _VIENNACL_COORDINATE_MATRIX_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/kernels/coordinate_matrix_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct coordinate_matrix
   {
    static viennacl::ocl::kernel vec_mul;

    static void init();
   };

   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel coordinate_matrix<T, ALIGNMENT>::vec_mul;


    /////////////// single precision kernels //////////////// 
    template <> void coordinate_matrix<float, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(coordinate_matrix_align1_vec_mul);
       prog.create(source);

       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void coordinate_matrix<float, 128>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(coordinate_matrix_align1_vec_mul);
       prog.create(source);

       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }


    /////////////// double precision kernels //////////////// 
    template <> void coordinate_matrix<double, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(coordinate_matrix_align1_vec_mul));
       prog.create(source);

       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }
    template <> void coordinate_matrix<double, 128>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(coordinate_matrix_align1_vec_mul));
       prog.create(source);

       vec_mul.prepareInit("vec_mul", prog);
       init_done = true;
     }
    }

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
