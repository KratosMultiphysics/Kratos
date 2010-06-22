#ifndef _VIENNACL_MATRIX_KERNELS_HPP_
#define _VIENNACL_MATRIX_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/kernels/matrix_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct matrix
   {
    static viennacl::ocl::kernel trans_lower_triangular_substitute_inplace;
    static viennacl::ocl::kernel rank1_update;
    static viennacl::ocl::kernel lower_triangular_substitute_inplace;
    static viennacl::ocl::kernel lu_factorize;
    static viennacl::ocl::kernel trans_vec_mul;
    static viennacl::ocl::kernel scaled_rank1_update;
    static viennacl::ocl::kernel vec_mul;
    static viennacl::ocl::kernel unit_upper_triangular_substitute_inplace;
    static viennacl::ocl::kernel unit_lower_triangular_substitute_inplace;
    static viennacl::ocl::kernel trans_upper_triangular_substitute_inplace;
    static viennacl::ocl::kernel upper_triangular_substitute_inplace;

    static void init();
   };

   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::trans_lower_triangular_substitute_inplace;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::rank1_update;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::lower_triangular_substitute_inplace;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::lu_factorize;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::trans_vec_mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::scaled_rank1_update;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::vec_mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::unit_upper_triangular_substitute_inplace;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::unit_lower_triangular_substitute_inplace;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::trans_upper_triangular_substitute_inplace;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel matrix<T, ALIGNMENT>::upper_triangular_substitute_inplace;


    /////////////// single precision kernels //////////////// 
    template <> void matrix<float, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(matrix_align1_trans_lower_triangular_substitute_inplace);
       source.append(matrix_align1_rank1_update);
       source.append(matrix_align1_lower_triangular_substitute_inplace);
       source.append(matrix_align1_lu_factorize);
       source.append(matrix_align1_trans_vec_mul);
       source.append(matrix_align1_scaled_rank1_update);
       source.append(matrix_align1_vec_mul);
       source.append(matrix_align1_unit_upper_triangular_substitute_inplace);
       source.append(matrix_align1_unit_lower_triangular_substitute_inplace);
       source.append(matrix_align1_trans_upper_triangular_substitute_inplace);
       source.append(matrix_align1_upper_triangular_substitute_inplace);
       prog.create(source);

       trans_lower_triangular_substitute_inplace.prepareInit("trans_lower_triangular_substitute_inplace", prog);
       rank1_update.prepareInit("rank1_update", prog);
       lower_triangular_substitute_inplace.prepareInit("lower_triangular_substitute_inplace", prog);
       lu_factorize.prepareInit("lu_factorize", prog);
       trans_vec_mul.prepareInit("trans_vec_mul", prog);
       scaled_rank1_update.prepareInit("scaled_rank1_update", prog);
       vec_mul.prepareInit("vec_mul", prog);
       unit_upper_triangular_substitute_inplace.prepareInit("unit_upper_triangular_substitute_inplace", prog);
       unit_lower_triangular_substitute_inplace.prepareInit("unit_lower_triangular_substitute_inplace", prog);
       trans_upper_triangular_substitute_inplace.prepareInit("trans_upper_triangular_substitute_inplace", prog);
       upper_triangular_substitute_inplace.prepareInit("upper_triangular_substitute_inplace", prog);
       init_done = true;
     }
    }
    template <> void matrix<float, 16>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(matrix_align1_trans_lower_triangular_substitute_inplace);
       source.append(matrix_align1_rank1_update);
       source.append(matrix_align1_lower_triangular_substitute_inplace);
       source.append(matrix_align1_lu_factorize);
       source.append(matrix_align1_trans_vec_mul);
       source.append(matrix_align1_scaled_rank1_update);
       source.append(matrix_align1_vec_mul);
       source.append(matrix_align1_unit_upper_triangular_substitute_inplace);
       source.append(matrix_align1_unit_lower_triangular_substitute_inplace);
       source.append(matrix_align1_trans_upper_triangular_substitute_inplace);
       source.append(matrix_align1_upper_triangular_substitute_inplace);
       prog.create(source);

       trans_lower_triangular_substitute_inplace.prepareInit("trans_lower_triangular_substitute_inplace", prog);
       rank1_update.prepareInit("rank1_update", prog);
       lower_triangular_substitute_inplace.prepareInit("lower_triangular_substitute_inplace", prog);
       lu_factorize.prepareInit("lu_factorize", prog);
       trans_vec_mul.prepareInit("trans_vec_mul", prog);
       scaled_rank1_update.prepareInit("scaled_rank1_update", prog);
       vec_mul.prepareInit("vec_mul", prog);
       unit_upper_triangular_substitute_inplace.prepareInit("unit_upper_triangular_substitute_inplace", prog);
       unit_lower_triangular_substitute_inplace.prepareInit("unit_lower_triangular_substitute_inplace", prog);
       trans_upper_triangular_substitute_inplace.prepareInit("trans_upper_triangular_substitute_inplace", prog);
       upper_triangular_substitute_inplace.prepareInit("upper_triangular_substitute_inplace", prog);
       init_done = true;
     }
    }


    /////////////// double precision kernels //////////////// 
    template <> void matrix<double, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_rank1_update));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_lu_factorize));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_vec_mul));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_scaled_rank1_update));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_vec_mul));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_unit_upper_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_unit_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_upper_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_upper_triangular_substitute_inplace));
       prog.create(source);

       trans_lower_triangular_substitute_inplace.prepareInit("trans_lower_triangular_substitute_inplace", prog);
       rank1_update.prepareInit("rank1_update", prog);
       lower_triangular_substitute_inplace.prepareInit("lower_triangular_substitute_inplace", prog);
       lu_factorize.prepareInit("lu_factorize", prog);
       trans_vec_mul.prepareInit("trans_vec_mul", prog);
       scaled_rank1_update.prepareInit("scaled_rank1_update", prog);
       vec_mul.prepareInit("vec_mul", prog);
       unit_upper_triangular_substitute_inplace.prepareInit("unit_upper_triangular_substitute_inplace", prog);
       unit_lower_triangular_substitute_inplace.prepareInit("unit_lower_triangular_substitute_inplace", prog);
       trans_upper_triangular_substitute_inplace.prepareInit("trans_upper_triangular_substitute_inplace", prog);
       upper_triangular_substitute_inplace.prepareInit("upper_triangular_substitute_inplace", prog);
       init_done = true;
     }
    }
    template <> void matrix<double, 16>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_rank1_update));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_lu_factorize));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_vec_mul));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_scaled_rank1_update));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_vec_mul));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_unit_upper_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_unit_lower_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_trans_upper_triangular_substitute_inplace));
       source.append(viennacl::tools::make_double_kernel(matrix_align1_upper_triangular_substitute_inplace));
       prog.create(source);

       trans_lower_triangular_substitute_inplace.prepareInit("trans_lower_triangular_substitute_inplace", prog);
       rank1_update.prepareInit("rank1_update", prog);
       lower_triangular_substitute_inplace.prepareInit("lower_triangular_substitute_inplace", prog);
       lu_factorize.prepareInit("lu_factorize", prog);
       trans_vec_mul.prepareInit("trans_vec_mul", prog);
       scaled_rank1_update.prepareInit("scaled_rank1_update", prog);
       vec_mul.prepareInit("vec_mul", prog);
       unit_upper_triangular_substitute_inplace.prepareInit("unit_upper_triangular_substitute_inplace", prog);
       unit_lower_triangular_substitute_inplace.prepareInit("unit_lower_triangular_substitute_inplace", prog);
       trans_upper_triangular_substitute_inplace.prepareInit("trans_upper_triangular_substitute_inplace", prog);
       upper_triangular_substitute_inplace.prepareInit("upper_triangular_substitute_inplace", prog);
       init_done = true;
     }
    }

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
