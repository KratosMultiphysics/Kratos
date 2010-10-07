#ifndef _VIENNACL_VECTOR_KERNELS_HPP_
#define _VIENNACL_VECTOR_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/kernels/vector_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct vector
   {
    static viennacl::ocl::kernel inplace_mul_sub;
    static viennacl::ocl::kernel swap;
    static viennacl::ocl::kernel cpu_mul_add;
    static viennacl::ocl::kernel norm_1;
    static viennacl::ocl::kernel norm_2;
    static viennacl::ocl::kernel index_norm_inf;
    static viennacl::ocl::kernel inplace_divide;
    static viennacl::ocl::kernel mul_sub;
    static viennacl::ocl::kernel sum;
    static viennacl::ocl::kernel clear;
    static viennacl::ocl::kernel plane_rotation;
    static viennacl::ocl::kernel cpu_mult;
    static viennacl::ocl::kernel inplace_mult;
    static viennacl::ocl::kernel norm_inf;
    static viennacl::ocl::kernel cpu_inplace_mult;
    static viennacl::ocl::kernel inplace_sub;
    static viennacl::ocl::kernel inplace_div_add;
    static viennacl::ocl::kernel cpu_inplace_mul_add;
    static viennacl::ocl::kernel inplace_add;
    static viennacl::ocl::kernel inplace_div_sub;
    static viennacl::ocl::kernel mul_add;
    static viennacl::ocl::kernel inplace_mul_add;
    static viennacl::ocl::kernel mult;
    static viennacl::ocl::kernel divide;
    static viennacl::ocl::kernel inner_prod;
    static viennacl::ocl::kernel add;
    static viennacl::ocl::kernel sub;

    static void init();
   };

   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_mul_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::swap;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::cpu_mul_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::norm_1;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::norm_2;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::index_norm_inf;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_divide;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::mul_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::sum;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::clear;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::plane_rotation;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::cpu_mult;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_mult;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::norm_inf;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::cpu_inplace_mult;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_div_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::cpu_inplace_mul_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_div_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::mul_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inplace_mul_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::mult;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::divide;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::inner_prod;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel vector<T, ALIGNMENT>::sub;


    /////////////// single precision kernels //////////////// 
    template <> void vector<float, 4>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(vector_align4_inplace_mul_sub);
       source.append(vector_align1_swap);
       source.append(vector_align4_cpu_mul_add);
       source.append(vector_align1_norm_1);
       source.append(vector_align4_norm_2);
       source.append(vector_align1_index_norm_inf);
       source.append(vector_align1_inplace_divide);
       source.append(vector_align1_mul_sub);
       source.append(vector_align1_sum);
       source.append(vector_align1_clear);
       source.append(vector_align1_plane_rotation);
       source.append(vector_align1_cpu_mult);
       source.append(vector_align1_inplace_mult);
       source.append(vector_align1_norm_inf);
       source.append(vector_align1_cpu_inplace_mult);
       source.append(vector_align1_inplace_sub);
       source.append(vector_align4_inplace_div_add);
       source.append(vector_align4_cpu_inplace_mul_add);
       source.append(vector_align1_inplace_add);
       source.append(vector_align4_inplace_div_sub);
       source.append(vector_align4_mul_add);
       source.append(vector_align4_inplace_mul_add);
       source.append(vector_align1_mult);
       source.append(vector_align1_divide);
       source.append(vector_align1_inner_prod);
       source.append(vector_align1_add);
       source.append(vector_align1_sub);
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }
    template <> void vector<float, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(vector_align1_inplace_mul_sub);
       source.append(vector_align1_swap);
       source.append(vector_align1_cpu_mul_add);
       source.append(vector_align1_norm_1);
       source.append(vector_align1_norm_2);
       source.append(vector_align1_index_norm_inf);
       source.append(vector_align1_inplace_divide);
       source.append(vector_align1_mul_sub);
       source.append(vector_align1_sum);
       source.append(vector_align1_clear);
       source.append(vector_align1_plane_rotation);
       source.append(vector_align1_cpu_mult);
       source.append(vector_align1_inplace_mult);
       source.append(vector_align1_norm_inf);
       source.append(vector_align1_cpu_inplace_mult);
       source.append(vector_align1_inplace_sub);
       source.append(vector_align1_inplace_div_add);
       source.append(vector_align1_cpu_inplace_mul_add);
       source.append(vector_align1_inplace_add);
       source.append(vector_align1_inplace_div_sub);
       source.append(vector_align1_mul_add);
       source.append(vector_align1_inplace_mul_add);
       source.append(vector_align1_mult);
       source.append(vector_align1_divide);
       source.append(vector_align1_inner_prod);
       source.append(vector_align1_add);
       source.append(vector_align1_sub);
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }
    template <> void vector<float, 16>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(vector_align4_inplace_mul_sub);
       source.append(vector_align1_swap);
       source.append(vector_align4_cpu_mul_add);
       source.append(vector_align1_norm_1);
       source.append(vector_align4_norm_2);
       source.append(vector_align1_index_norm_inf);
       source.append(vector_align16_inplace_divide);
       source.append(vector_align1_mul_sub);
       source.append(vector_align1_sum);
       source.append(vector_align1_clear);
       source.append(vector_align1_plane_rotation);
       source.append(vector_align16_cpu_mult);
       source.append(vector_align16_inplace_mult);
       source.append(vector_align1_norm_inf);
       source.append(vector_align1_cpu_inplace_mult);
       source.append(vector_align16_inplace_sub);
       source.append(vector_align4_inplace_div_add);
       source.append(vector_align4_cpu_inplace_mul_add);
       source.append(vector_align16_inplace_add);
       source.append(vector_align4_inplace_div_sub);
       source.append(vector_align4_mul_add);
       source.append(vector_align4_inplace_mul_add);
       source.append(vector_align16_mult);
       source.append(vector_align16_divide);
       source.append(vector_align1_inner_prod);
       source.append(vector_align16_add);
       source.append(vector_align16_sub);
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }


    /////////////// double precision kernels //////////////// 
    template <> void vector<double, 4>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_swap));
       source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_mul_add));
       #ifndef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK_ON_GPU
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1));
       source.append(viennacl::tools::make_double_kernel(vector_align4_norm_2));
       source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf));
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf));
       #endif
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_sum));
       source.append(viennacl::tools::make_double_kernel(vector_align1_clear));
       source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align4_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod));
       source.append(viennacl::tools::make_double_kernel(vector_align1_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_sub));
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }
    template <> void vector<double, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_swap));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mul_add));
       #ifndef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK_ON_GPU
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1));
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_2));
       source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf));
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf));
       #endif
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_sum));
       source.append(viennacl::tools::make_double_kernel(vector_align1_clear));
       source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_div_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_div_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod));
       source.append(viennacl::tools::make_double_kernel(vector_align1_add));
       source.append(viennacl::tools::make_double_kernel(vector_align1_sub));
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }
    template <> void vector<double, 16>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_swap));
       source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_mul_add));
       #ifndef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK_ON_GPU
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_1));
       source.append(viennacl::tools::make_double_kernel(vector_align4_norm_2));
       source.append(viennacl::tools::make_double_kernel(vector_align1_index_norm_inf));
       source.append(viennacl::tools::make_double_kernel(vector_align1_norm_inf));
       #endif
       source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_mul_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align1_sum));
       source.append(viennacl::tools::make_double_kernel(vector_align1_clear));
       source.append(viennacl::tools::make_double_kernel(vector_align1_plane_rotation));
       source.append(viennacl::tools::make_double_kernel(vector_align16_cpu_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align1_cpu_inplace_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_cpu_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align16_inplace_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_div_sub));
       source.append(viennacl::tools::make_double_kernel(vector_align4_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align4_inplace_mul_add));
       source.append(viennacl::tools::make_double_kernel(vector_align16_mult));
       source.append(viennacl::tools::make_double_kernel(vector_align16_divide));
       source.append(viennacl::tools::make_double_kernel(vector_align1_inner_prod));
       source.append(viennacl::tools::make_double_kernel(vector_align16_add));
       source.append(viennacl::tools::make_double_kernel(vector_align16_sub));
       prog.create(source);

       inplace_mul_sub.prepareInit("inplace_mul_sub", prog);
       swap.prepareInit("swap", prog);
       cpu_mul_add.prepareInit("cpu_mul_add", prog);
       norm_1.prepareInit("norm_1", prog);
       norm_2.prepareInit("norm_2", prog);
       index_norm_inf.prepareInit("index_norm_inf", prog);
       inplace_divide.prepareInit("inplace_divide", prog);
       mul_sub.prepareInit("mul_sub", prog);
       sum.prepareInit("sum", prog);
       clear.prepareInit("clear", prog);
       plane_rotation.prepareInit("plane_rotation", prog);
       cpu_mult.prepareInit("cpu_mult", prog);
       inplace_mult.prepareInit("inplace_mult", prog);
       norm_inf.prepareInit("norm_inf", prog);
       cpu_inplace_mult.prepareInit("cpu_inplace_mult", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_div_add.prepareInit("inplace_div_add", prog);
       cpu_inplace_mul_add.prepareInit("cpu_inplace_mul_add", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div_sub.prepareInit("inplace_div_sub", prog);
       mul_add.prepareInit("mul_add", prog);
       inplace_mul_add.prepareInit("inplace_mul_add", prog);
       mult.prepareInit("mult", prog);
       divide.prepareInit("divide", prog);
       inner_prod.prepareInit("inner_prod", prog);
       add.prepareInit("add", prog);
       sub.prepareInit("sub", prog);
       init_done = true;
     }
    }

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
