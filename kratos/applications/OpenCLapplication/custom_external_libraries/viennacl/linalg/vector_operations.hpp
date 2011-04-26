/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory
======================================================================= */

#ifndef _VIENNACL_VECTOR_OPERATIONS_HPP_
#define _VIENNACL_VECTOR_OPERATIONS_HPP_

/** @file vector_operations.hpp
    @brief Implementations of vector operations.
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/kernels/vector_kernels.h"

namespace viennacl
{
  namespace linalg
  {
    /** @brief Addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * @param vec1  The first addend. 
    * @param vec2  The second addend.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
             const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
             viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec1.size() == vec2.size() && vec1.size() == result.size());

      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "add");
      
      viennacl::ocl::enqueue(k(vec1, vec2, result, size));        
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The addend
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                     const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_add");

      viennacl::ocl::enqueue(k(vec1, vec2, size));        
    }

    /** @brief Subtraction of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * result = vec1 - vec2
    *
    * @param vec1  The first operand. 
    * @param vec2  The second operand.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void sub(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
             const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
             viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "sub");

      viennacl::ocl::enqueue(k(vec1, vec2, result, size));        
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The subtracted vector
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                     const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_sub");

      viennacl::ocl::enqueue(k(vec1, vec2, size));        
    }


    //result = vec * scalar
    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec * alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mult(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
              scalar<SCALARTYPE> const & alpha,
              viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      result.resize(vec.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "mult");

      viennacl::ocl::enqueue(k(vec, alpha, result, static_cast<cl_uint>(vec.internal_size())));        
    }

    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec * alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mult(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
              SCALARTYPE alpha,
              viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      result.resize(vec.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "cpu_mult");

      viennacl::ocl::enqueue(k(vec, alpha, result, static_cast<cl_uint>(vec.internal_size())));        
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mult(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                      scalar<SCALARTYPE> const & alpha)
    {
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_mult");

      viennacl::ocl::enqueue(k(vec, alpha, static_cast<cl_uint>(vec.internal_size())));        
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mult(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                      SCALARTYPE alpha)
    {
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "cpu_inplace_mult");

      viennacl::ocl::enqueue(k(vec, alpha, static_cast<cl_uint>(vec.internal_size())));        
    }

    //result = vec / scalar
    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec / alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The (inverse) scaling factor.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void divide(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                scalar<SCALARTYPE> const & alpha,
                viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec.size() == result.size());
      result.resize(vec.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "divide");

      viennacl::ocl::enqueue(k(vec, alpha, result, static_cast<cl_uint>(vec.internal_size())));        
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The (inverse) scaling factor.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_divide(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                        scalar<SCALARTYPE> const & alpha)
    {
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_divide");

      viennacl::ocl::enqueue(k(vec, alpha, static_cast<cl_uint>(vec.internal_size())));        
    }

    //result = factor * vec1 + vec2
    /** @brief Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = alpha * vec1 + vec2, where alpha is a gpu scalar
    *
    * @param vec1    The first added
    * @param alpha  The scaling factor for the first addend.
    * @param vec2    The second added.
    * @param result The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 scalar<SCALARTYPE> const & alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec1.size() == vec2.size() && result.size() == vec1.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "mul_add");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, alpha, vec2, result, size));        
    }

    /** @brief Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = alpha * vec1 + vec2, where alpha is a cpu scalar
    *
    * @param vec1    The first added
    * @param alpha   The scaling factor for the first addend.
    * @param vec2    The second added.
    * @param result  The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 SCALARTYPE alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec1.size() == vec2.size() && result.size() == vec1.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "cpu_mul_add");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, alpha, vec2, result, size));        
    }

    //vec1 += factor * vec2
    /** @brief Inplace Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += alpha * vec2, where alpha is a gpu scalar
    *
    * @param vec1    The first added
    * @param alpha   The scaling factor for the first addend.
    * @param vec2    The second added.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha)
    {
      assert(vec1.size() == vec2.size());
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_mul_add");

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, size));        
    }

    /** @brief Inplace Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += alpha * vec2, where alpha is a cpu scalar
    *
    * @param vec1    The first added
    * @param vec2    The second added.
    * @param alpha   The scaling factor for the first addend.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         SCALARTYPE alpha)
    {
      assert(vec1.size() == vec2.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "cpu_inplace_mul_add");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, size));        
    }

    /** @brief Multiply-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = alpha * vec1 - vec2, where alpha is a gpu scalar
    *
    * @param vec1    The first vector operand
    * @param alpha   The scaling factor for the first vector.
    * @param vec2    The second operand.
    * @param result  The result vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_sub(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 scalar<SCALARTYPE> const & alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result)
    {
      assert(vec1.size() == vec2.size() && result.size() == vec1.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "mul_sub");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, alpha, vec2, result, size));        
    }


    /** @brief Inplace Multiply-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= alpha * vec2, where alpha is a gpu scalar
    *
    * @param vec1    The result vector which is updated
    * @param vec2    The second operand.
    * @param alpha   The scaling factor for the vector update.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha)
    {
      assert(vec1.size() == vec2.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_mul_sub");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, size));        
    }

    /** @brief Inplace divide-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_div_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha)
    {
      assert(vec1.size() == vec2.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_div_add");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, size));        
    }

    /** @brief Inplace divide-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_div_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha)
    {
      assert(vec1.size() == vec2.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inplace_div_sub");
      cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, size));        
    }


    ///////////////////////// Norms and inner product ///////////////////


    //implementation of inner product:
    //namespace {
      /** @brief Computes the inner product of two vectors - implementation. Library users should call inner_prod(vec1, vec2).
      *
      * @param vec1 The first vector
      * @param vec2 The second vector
      * @param result The result scalar (on the gpu)
      */
      template<class SCALARTYPE, unsigned int ALIGNMENT>
      void inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                           const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                           scalar<SCALARTYPE> & result)
      {
        assert(vec1.size() == vec2.size());
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "inner_prod");
        cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));
        unsigned int work_groups = k.global_work_size() / k.local_work_size();
        static viennacl::vector<SCALARTYPE> temp(work_groups);
        
        /*unsigned int pos = 0;
        k.argument(pos++, vec1.handle());
        k.argument(pos++, vec2.handle());
        k.argument(pos++, size);
        k.local_buffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE) * k.local_work_size()));
        k.argument(pos++, temp.handle());*/
        
        //Note: Number of work groups MUST be a power of two!
        //std::cout << work_groups << ", " << k.local_work_size() << ", " << k.global_work_size() << std::endl;
        assert( work_groups * k.local_work_size() == k.global_work_size() );
        assert( (k.global_work_size() / k.local_work_size()) == 1 
               || (k.global_work_size() / k.local_work_size()) == 2 
               || (k.global_work_size() / k.local_work_size()) == 4
               || (k.global_work_size() / k.local_work_size()) == 8
               || (k.global_work_size() / k.local_work_size()) == 16
               || (k.global_work_size() / k.local_work_size()) == 32
               || (k.global_work_size() / k.local_work_size()) == 64
               || (k.global_work_size() / k.local_work_size()) == 128
               || (k.global_work_size() / k.local_work_size()) == 256
               || (k.global_work_size() / k.local_work_size()) == 512 );
               
        viennacl::ocl::enqueue(k(vec1, vec2, size, viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()), temp));        

        viennacl::ocl::kernel & ksum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "sum");
        
        ksum.local_work_size(0, work_groups);
        ksum.global_work_size(0, work_groups);
        viennacl::ocl::enqueue(ksum(temp, result));
      }
    //}

    //public interface of inner product
    /** @brief Computes the inner product of two vectors.
    *
    * @param vec1 The first vector
    * @param vec2 The second vector
    * @return The result
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT1, unsigned int ALIGNMENT2>
    viennacl::scalar_expression< const viennacl::vector<SCALARTYPE, ALIGNMENT1>, 
                                 const viennacl::vector<SCALARTYPE, ALIGNMENT2>,
                                 viennacl::op_inner_prod >
    inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT1> & vec1,
                    const viennacl::vector<SCALARTYPE, ALIGNMENT2> & vec2)
    {
      return viennacl::scalar_expression< const viennacl::vector<SCALARTYPE, ALIGNMENT1>, 
                                          const viennacl::vector<SCALARTYPE, ALIGNMENT2>,
                                          viennacl::op_inner_prod >(vec1, vec2);
    }


    
    /** @brief Computes the l^1-norm of a vector
    *
    * @param vcl_vec The vector
    * @param result The result scalar
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void norm_1_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                     scalar<SCALARTYPE> & result)
    {
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "norm_1");
      cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      
      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }
      
      
      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<SCALARTYPE> temp(work_groups);
        
      //Note: Number of work groups MUST be a power of two!
      //std::cout << work_groups << ", " << k.local_work_size() << ", " << k.global_work_size() << std::endl;
      assert( work_groups * k.local_work_size() == k.global_work_size() );
      assert( (k.global_work_size() / k.local_work_size()) == 1 
             || (k.global_work_size() / k.local_work_size()) == 2 
             || (k.global_work_size() / k.local_work_size()) == 4
             || (k.global_work_size() / k.local_work_size()) == 8
             || (k.global_work_size() / k.local_work_size()) == 16
             || (k.global_work_size() / k.local_work_size()) == 32
             || (k.global_work_size() / k.local_work_size()) == 64
             || (k.global_work_size() / k.local_work_size()) == 128
             || (k.global_work_size() / k.local_work_size()) == 256
             || (k.global_work_size() / k.local_work_size()) == 512 );
               
        viennacl::ocl::enqueue(k(vcl_vec, size, viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()), temp));        
        
        viennacl::ocl::kernel & ksum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "sum");
        
        ksum.local_work_size(0, work_groups);
        ksum.global_work_size(0, work_groups);
        viennacl::ocl::enqueue(ksum(temp, result));
    }

    /** @brief Computes the l^2-norm of a vector - implementation
    *
    * @param vcl_vec The vector
    * @param result The result scalar
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void norm_2_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                     scalar<SCALARTYPE> & result)
    {
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "norm_2");
      cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      
      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }

      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<SCALARTYPE> temp(work_groups);
        
      //Note: Number of work groups MUST be a power of two!
      //std::cout << work_groups << ", " << k.local_work_size() << ", " << k.global_work_size() << std::endl;
      assert( work_groups * k.local_work_size() == k.global_work_size() );
      assert( (k.global_work_size() / k.local_work_size()) == 1 
             || (k.global_work_size() / k.local_work_size()) == 2 
             || (k.global_work_size() / k.local_work_size()) == 4
             || (k.global_work_size() / k.local_work_size()) == 8
             || (k.global_work_size() / k.local_work_size()) == 16
             || (k.global_work_size() / k.local_work_size()) == 32
             || (k.global_work_size() / k.local_work_size()) == 64
             || (k.global_work_size() / k.local_work_size()) == 128
             || (k.global_work_size() / k.local_work_size()) == 256
             || (k.global_work_size() / k.local_work_size()) == 512 );
               
        viennacl::ocl::enqueue(k(vcl_vec, size, viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()), temp));        

        viennacl::ocl::kernel & sqrt_sum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "sqrt_sum");
        
        sqrt_sum.local_work_size(0, work_groups);
        sqrt_sum.global_work_size(0, work_groups);
        viennacl::ocl::enqueue(sqrt_sum(temp, result, work_groups));
    }

    /** @brief Computes the supremum-norm of a vector
    *
    * @param vcl_vec The vector
    * @param result The result scalar
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void norm_inf_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                  scalar<SCALARTYPE> & result)
    {
      cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "norm_inf");

      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }
      
      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<SCALARTYPE> temp(work_groups);
        
      //Note: Number of work groups MUST be a power of two!
      //std::cout << work_groups << ", " << k.local_work_size() << ", " << k.global_work_size() << std::endl;
      assert( work_groups * k.local_work_size() == k.global_work_size() );
      assert( work_groups == 1 
             || work_groups == 2 
             || work_groups == 4
             || work_groups == 8
             || work_groups == 16
             || work_groups == 32
             || work_groups == 64
             || work_groups == 128
             || work_groups == 256
             || work_groups == 512 );
               
        viennacl::ocl::enqueue(k(vcl_vec, size, viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()), temp));
        //viennacl::ocl::get_queue().finish();
        
        //part 2: parallel reduction of reduced kernel:
        viennacl::ocl::kernel & max_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "vmax");
        max_kernel.local_work_size(0, work_groups);
        max_kernel.global_work_size(0, work_groups);
        
        viennacl::ocl::enqueue(max_kernel(temp, result, work_groups));
    }

    //This function should return a CPU scalar, otherwise statements like 
    // vcl_rhs[index_norm_inf(vcl_rhs)] 
    // are ambiguous
    /** @brief Computes the index of the first entry that is equal to the supremum-norm in modulus.
    *
    * @param vcl_vec The vector
    * @return The result. Note that the result must be a CPU scalar (unsigned int), since gpu scalars are floating point types.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    cl_uint index_norm_inf(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec)
    {
      viennacl::ocl::handle<cl_mem> h = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(cl_uint));
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "index_norm_inf");
      cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());

      k.global_work_size(0, k.local_work_size());
      viennacl::ocl::enqueue(k(vcl_vec,
                               size,
                               viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()),
                               viennacl::ocl::local_mem(sizeof(cl_uint) * k.local_work_size()), h));
      
      //read value:
      cl_uint result;
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), h, CL_TRUE, 0, sizeof(cl_uint), &result, 0, NULL, NULL);
      VIENNACL_ERR_CHECK(err);
      return result;
    }
    
    //TODO: Special case vec1 == vec2 allows improvement!!
    /** @brief Computes a plane rotation of two vectors.
    *
    * Computes (x,y) <- (alpha * x + beta * y, -beta * x + alpha * y)
    *
    * @param vec1   The first vector
    * @param vec2   The second vector
    * @param alpha  The first transformation coefficient
    * @param beta   The second transformation coefficient
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void plane_rotation(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                        const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                        SCALARTYPE alpha,
                        SCALARTYPE beta)
    {
      assert(vec1.size() == vec2.size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "plane_rotation");

      viennacl::ocl::enqueue(k(vec1, vec2, alpha, beta, static_cast<cl_uint>(vec1.size())));
    }
    
  } //namespace linalg
} //namespace viennacl


#endif
