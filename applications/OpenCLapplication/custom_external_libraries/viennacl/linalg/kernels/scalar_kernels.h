#ifndef _VIENNACL_SCALAR_KERNELS_HPP_
#define _VIENNACL_SCALAR_KERNELS_HPP_
#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/linalg/kernels/scalar_source.h"

//Automatically generated file from aux-directory, do not edit manually!
namespace viennacl
{
 namespace linalg
 {
  namespace kernels
  {
   template<class TYPE, unsigned int alignment>
   struct scalar
   {
    static viennacl::ocl::kernel cpu_sub;
    static viennacl::ocl::kernel cpu_inplace_mul;
    static viennacl::ocl::kernel inplace_mul;
    static viennacl::ocl::kernel cpu_inplace_sub;
    static viennacl::ocl::kernel cpu_inplace_add;
    static viennacl::ocl::kernel mul;
    static viennacl::ocl::kernel inplace_sub;
    static viennacl::ocl::kernel inplace_add;
    static viennacl::ocl::kernel inplace_div;
    static viennacl::ocl::kernel cpu_inplace_div;
    static viennacl::ocl::kernel divide;
    static viennacl::ocl::kernel cpu_add;
    static viennacl::ocl::kernel add;
    static viennacl::ocl::kernel cpu_mul;
    static viennacl::ocl::kernel sub;
    static viennacl::ocl::kernel cpu_div;

    static void init();
   };

   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_inplace_mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::inplace_mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_inplace_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_inplace_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::inplace_sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::inplace_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::inplace_div;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_inplace_div;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::divide;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::add;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_mul;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::sub;
   template <typename T, unsigned int ALIGNMENT> viennacl::ocl::kernel scalar<T, ALIGNMENT>::cpu_div;


    /////////////// single precision kernels //////////////// 
    template <> void scalar<float, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(scalar_align1_cpu_sub);
       source.append(scalar_align1_cpu_inplace_mul);
       source.append(scalar_align1_inplace_mul);
       source.append(scalar_align1_cpu_inplace_sub);
       source.append(scalar_align1_cpu_inplace_add);
       source.append(scalar_align1_mul);
       source.append(scalar_align1_inplace_sub);
       source.append(scalar_align1_inplace_add);
       source.append(scalar_align1_inplace_div);
       source.append(scalar_align1_cpu_inplace_div);
       source.append(scalar_align1_divide);
       source.append(scalar_align1_cpu_add);
       source.append(scalar_align1_add);
       source.append(scalar_align1_cpu_mul);
       source.append(scalar_align1_sub);
       source.append(scalar_align1_cpu_div);
       prog.create(source);

       cpu_sub.prepareInit("cpu_sub", prog);
       cpu_inplace_mul.prepareInit("cpu_inplace_mul", prog);
       inplace_mul.prepareInit("inplace_mul", prog);
       cpu_inplace_sub.prepareInit("cpu_inplace_sub", prog);
       cpu_inplace_add.prepareInit("cpu_inplace_add", prog);
       mul.prepareInit("mul", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div.prepareInit("inplace_div", prog);
       cpu_inplace_div.prepareInit("cpu_inplace_div", prog);
       divide.prepareInit("divide", prog);
       cpu_add.prepareInit("cpu_add", prog);
       add.prepareInit("add", prog);
       cpu_mul.prepareInit("cpu_mul", prog);
       sub.prepareInit("sub", prog);
       cpu_div.prepareInit("cpu_div", prog);
       init_done = true;
     }
    }


    /////////////// double precision kernels //////////////// 
    template <> void scalar<double, 1>::init()
    {
     static bool init_done = false;
     if (!init_done)
     {
       std::string source;
       viennacl::ocl::program prog;
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_sub));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_mul));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_mul));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_sub));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_add));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_mul));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_sub));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_add));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_inplace_div));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_inplace_div));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_divide));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_add));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_add));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_mul));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_sub));
       source.append(viennacl::tools::make_double_kernel(scalar_align1_cpu_div));
       prog.create(source);

       cpu_sub.prepareInit("cpu_sub", prog);
       cpu_inplace_mul.prepareInit("cpu_inplace_mul", prog);
       inplace_mul.prepareInit("inplace_mul", prog);
       cpu_inplace_sub.prepareInit("cpu_inplace_sub", prog);
       cpu_inplace_add.prepareInit("cpu_inplace_add", prog);
       mul.prepareInit("mul", prog);
       inplace_sub.prepareInit("inplace_sub", prog);
       inplace_add.prepareInit("inplace_add", prog);
       inplace_div.prepareInit("inplace_div", prog);
       cpu_inplace_div.prepareInit("cpu_inplace_div", prog);
       divide.prepareInit("divide", prog);
       cpu_add.prepareInit("cpu_add", prog);
       add.prepareInit("add", prog);
       cpu_mul.prepareInit("cpu_mul", prog);
       sub.prepareInit("sub", prog);
       cpu_div.prepareInit("cpu_div", prog);
       init_done = true;
     }
    }

  }  //namespace kernels
 }  //namespace linalg
}  //namespace viennacl
#endif
