#ifndef VIENNACL_VECTOR_OPERATIONS_HPP_
#define VIENNACL_VECTOR_OPERATIONS_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

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
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/handle.hpp"

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
    template <typename V1, typename V2, typename V3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_vector<V3>::value
                                >::type
    add(const V1 & vec1, 
        const V2 & vec2, 
        V3 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && (viennacl::traits::size(vec1) == viennacl::traits::size(result))
             && "Incompatible vector sizes in add()!");

      //unsigned int size = std::min(viennacl::traits::internal_size(vec1),
      //                             viennacl::traits::internal_size(vec2));
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "add");
      
      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)),
                               viennacl::traits::handle(result),  cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)) )
                            );
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The addend
    */
    template <typename V1, typename V2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                >::type
    inplace_add(V1 & vec1,
                const V2 & vec2)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_add()!");
      
      
      //unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_add");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)))
                            );
    }



    /** @brief Subtraction of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * result = vec1 - vec2
    *
    * @param vec1  The first operand. 
    * @param vec2  The second operand.
    * @param result The result vector.
    */
    template <typename V1, typename V2, typename V3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_vector<V3>::value
                                >::type
    sub(const V1 & vec1,
        const V2 & vec2,
        V3 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && (viennacl::traits::size(vec1) == viennacl::traits::size(result))
             && "Incompatible vector sizes in sub()!");
      
      //unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "sub");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)),
                               viennacl::traits::handle(result),  cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)) )
                            );        
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The subtracted vector
    */
    template <typename V1, typename V2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                >::type
    inplace_sub(V1 & vec1,
                const V2 & vec2)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_sub()!");
      
      //unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_sub");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)))
                            );        
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
    template <typename V1, typename S2, typename V3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                  && viennacl::is_vector<V3>::value
                                >::type
    mult(const V1 & vec,
         S2 const & alpha,
         V3 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec) == viennacl::traits::size(result))
             && "Incompatible vector sizes in mult()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "mult");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),
                               alpha,
                               viennacl::traits::handle(result),  cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );        
    }

    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec * alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param result The result vector.
    */
    template <typename V1, typename SCALARTYPE, typename V3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_cpu_scalar<SCALARTYPE>::value
                                  && viennacl::is_vector<V3>::value
                                >::type
    mult(V1 const & vec,
         SCALARTYPE alpha,
         V3 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec) == viennacl::traits::size(result))
             && "Incompatible vector sizes in mult()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "cpu_mult");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),
                               static_cast<value_type>(alpha),
                               viennacl::traits::handle(result),  cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );        
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    */
    template <typename V1, typename S2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                >::type
    inplace_mult(V1 & vec,
                 S2 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_mult");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),
                               alpha)
                            );
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    */
    template <typename V1, typename S2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_cpu_scalar<S2>::value
                                >::type
    inplace_mult(V1 & vec,
                 S2 alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "cpu_inplace_mult");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)), 
                               static_cast<value_type>(alpha))
                            );        
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
    template <typename V1, typename S2, typename V3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                  && viennacl::is_vector<V3>::value
                                >::type
    divide(V1 const & vec,
           S2 const & alpha,
           V3 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec) == viennacl::traits::size(result))
             && "Incompatible vector sizes in divide()!");

      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "divide");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)), 
                               alpha,
                               viennacl::traits::handle(result), cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The (inverse) scaling factor.
    */
    template <typename V1, typename S2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                >::type
    inplace_divide(V1 & vec,
                   S2 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_divide");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)), 
                               alpha) 
                            );
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
    template <typename V1, typename S2, typename V3, typename V4>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                  && viennacl::is_vector<V3>::value
                                  && viennacl::is_vector<V4>::value
                                >::type
    mul_add(V1 const & vec1,
            S2 const & alpha,
            V3 const & vec2,
            V4 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && (viennacl::traits::size(vec1) == viennacl::traits::size(result))
             && "Incompatible vector sizes in mul_add()!");
      
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "mul_add");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               alpha,
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               viennacl::traits::handle(result), cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );        
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
    template <typename V1, typename SCALARTYPE, typename V3, typename V4>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_cpu_scalar<SCALARTYPE>::value
                                  && viennacl::is_vector<V3>::value
                                  && viennacl::is_vector<V4>::value
                                >::type
    mul_add(V1 const & vec1,
            SCALARTYPE alpha,
            V3 const & vec2,
            V4 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && (viennacl::traits::size(vec1) == viennacl::traits::size(result))
             && "Incompatible vector sizes in mul_add()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "cpu_mul_add");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               static_cast<value_type>(alpha),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               viennacl::traits::handle(result), cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );
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
    template <typename V1, typename V2, typename S3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_scalar<S3>::value
                                >::type
    inplace_mul_add(V1 & vec1,
                    V2 const & vec2,
                    S3 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_mul_add()!");
      
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_mul_add");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               alpha));
    }

    /** @brief Inplace Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += alpha * vec2, where alpha is a cpu scalar
    *
    * @param vec1    The first added
    * @param vec2    The second added.
    * @param alpha   The scaling factor for the first addend.
    */
    template <typename V1, typename V2, typename SCALARTYPE>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_cpu_scalar<SCALARTYPE>::value
                                >::type
    inplace_mul_add(V1 & vec1,
                    V2 const & vec2,
                    SCALARTYPE alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_mul_add()!");

      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "cpu_inplace_mul_add");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               value_type(alpha)));
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
    template <typename V1, typename S2, typename V3, typename V4>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_scalar<S2>::value
                                  && viennacl::is_vector<V3>::value
                                  && viennacl::is_vector<V4>::value
                                >::type
    mul_sub(V1 const & vec1,
            S2 const & alpha,
            V3 const & vec2,
            V4 & result)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && (viennacl::traits::size(vec1) == viennacl::traits::size(result))
             && "Incompatible vector sizes in mul_sub()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "mul_sub");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               alpha,
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               viennacl::traits::handle(result), cl_uint(viennacl::traits::start(result)), cl_uint(viennacl::traits::size(result)))
                            );
    }


    /** @brief Inplace Multiply-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= alpha * vec2, where alpha is a gpu scalar
    *
    * @param vec1    The result vector which is updated
    * @param vec2    The second operand.
    * @param alpha   The scaling factor for the vector update.
    */
    template <typename V1, typename V2, typename S3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_scalar<S3>::value
                                >::type
    inplace_mul_sub(V1 & vec1,
                    V2 const & vec2,
                    S3 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_mul_sub()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_mul_sub");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               alpha)
                            );        
    }

    /** @brief Inplace divide-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    */
    template <typename V1, typename V2, typename S3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_scalar<S3>::value
                                >::type
    inplace_div_add(V1 & vec1,
                    V2 const & vec2,
                    S3 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_div_add()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_div_add");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)), 
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)), 
                               alpha)
                            );
    }

    /** @brief Inplace divide-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    */
    template <typename V1, typename V2, typename S3>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_scalar<S3>::value
                                >::type
    inplace_div_sub(V1 & vec1,
                    V2 const & vec2,
                    S3 const & alpha)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inplace_div_sub()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inplace_div_sub");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)),
                               alpha)
                            );        
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
    template <typename V1, typename V2, typename S3>
    void inner_prod_impl(V1 const & vec1,
                         V2 const & vec2,
                         S3 & result,
                         typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                       && viennacl::is_vector<V2>::value
                                                       && viennacl::is_scalar<S3>::value
#ifdef _MSC_VER
                                                     >::type * dummy = 0)
#else
                                                     >::type * dummy)
#endif                                                   
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
    
      assert( (viennacl::traits::size(vec1) == viennacl::traits::size(vec2))
             && "Incompatible vector sizes in inner_prod_impl()!");
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "inner_prod");
      //cl_uint size = static_cast<cl_uint>(std::min(vec1.internal_size(), vec2.internal_size()));
      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      
      static viennacl::vector<value_type> temp(work_groups);
      
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
              
      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)),
                               viennacl::ocl::local_mem(sizeof(value_type) * k.local_work_size()),
                               temp));        

      viennacl::ocl::kernel & ksum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "sum");
      
      ksum.local_work_size(0, work_groups);
      ksum.global_work_size(0, work_groups);
      viennacl::ocl::enqueue(ksum(viennacl::traits::handle(temp), cl_uint(viennacl::traits::start(temp)), cl_uint(viennacl::traits::size(temp)),
                                  result)
                            );
    }

    //public interface of inner product
    /** @brief Computes the inner product of two vectors.
    *
    * @param vec1 The first vector
    * @param vec2 The second vector
    * @return The result
    */
    template <typename V1, typename V2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value,
                                  viennacl::scalar_expression< const V1, 
                                                               const V2,
                                                               viennacl::op_inner_prod >
                                >::type
    inner_prod_impl(V1 const & vec1,
                    V2 const & vec2)
    {
      return viennacl::scalar_expression< const V1, 
                                          const V2,
                                          viennacl::op_inner_prod >(vec1, vec2);
    }


    
    /** @brief Computes the l^1-norm of a vector
    *
    * @param vec The vector
    * @param result The result scalar
    */
    template <typename V1, typename S2>
    void norm_1_impl(V1 const & vec,
                     S2 & result,
                     typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                   && viennacl::is_scalar<S2>::value
#ifdef _MSC_VER
                                                 >::type * dummy = 0)
#else
                                                 >::type * dummy)
#endif                                                   
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "norm_1");
      //cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      
      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }
      
      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<value_type> temp(work_groups);

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
               
      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),                                 
                                viennacl::ocl::local_mem(sizeof(value_type) * k.local_work_size()),
                                temp));        
      
      viennacl::ocl::kernel & ksum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "sum");
      
      ksum.local_work_size(0, work_groups);
      ksum.global_work_size(0, work_groups);
      viennacl::ocl::enqueue(ksum(viennacl::traits::handle(temp), cl_uint(viennacl::traits::start(temp)), cl_uint(viennacl::traits::size(temp)),
                                  result)
                            );
    }

    /** @brief Computes the l^2-norm of a vector - implementation
    *
    * @param vec The vector
    * @param result The result scalar
    */
    template <typename V1, typename S2>
    void norm_2_impl(V1 const & vec,
                     S2 & result,
                     typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                  && viennacl::is_scalar<S2>::value
#ifdef _MSC_VER
                                                 >::type * dummy = 0)
#else
                                                 >::type * dummy)
#endif                                                   
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "norm_2");
      //cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      
      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }

      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<value_type> temp(work_groups);
        
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
               
        viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),                                 
                                 viennacl::ocl::local_mem(sizeof(value_type) * k.local_work_size()),
                                 temp)
                              );

        viennacl::ocl::kernel & sqrt_sum = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "sqrt_sum");
        
        sqrt_sum.local_work_size(0, work_groups);
        sqrt_sum.global_work_size(0, work_groups);
        viennacl::ocl::enqueue(
                        sqrt_sum(viennacl::traits::handle(temp), cl_uint(viennacl::traits::start(temp)), cl_uint(viennacl::traits::size(temp)),
                                 result)
                              );
    }

    /** @brief Computes the supremum-norm of a vector
    *
    * @param vec The vector
    * @param result The result scalar
    */
    template <typename V1, typename S2>
    void norm_inf_impl(V1 const & vec,
                       S2 & result,
                       typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                     && viennacl::is_scalar<S2>::value
#ifdef _MSC_VER
                                                   >::type * dummy = 0)
#else
                                                   >::type * dummy)
#endif                                                   
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      //cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "norm_inf");

      if (k.local_work_size() != k.global_work_size())
      {
        //NOTE: For some reasons the kernel could not be started with several work groups on NVIDIA hardware. This forces us to use as many parallel threads within a single work group as possible
        k.local_work_size(0, viennacl::ocl::current_device().max_work_group_size());
        k.global_work_size(0, viennacl::ocl::current_device().max_work_group_size());
      }
      
      unsigned int work_groups = k.global_work_size() / k.local_work_size();
      viennacl::vector<value_type> temp(work_groups);
        
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
               
      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),                                 
                               viennacl::ocl::local_mem(sizeof(value_type) * k.local_work_size()),
                               temp));
      //viennacl::ocl::get_queue().finish();
      
      //part 2: parallel reduction of reduced kernel:
      viennacl::ocl::kernel & max_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "vmax");
      max_kernel.local_work_size(0, work_groups);
      max_kernel.global_work_size(0, work_groups);
      
      viennacl::ocl::enqueue(
                       max_kernel(viennacl::traits::handle(temp), cl_uint(viennacl::traits::start(temp)), cl_uint(viennacl::traits::size(temp)),
                                  result)
                            );
    }

    //This function should return a CPU scalar, otherwise statements like 
    // vcl_rhs[index_norm_inf(vcl_rhs)] 
    // are ambiguous
    /** @brief Computes the index of the first entry that is equal to the supremum-norm in modulus.
    *
    * @param vec The vector
    * @return The result. Note that the result must be a CPU scalar (unsigned int), since gpu scalars are floating point types.
    */
    template <typename V1>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value,
                                  cl_uint
                                >::type
    index_norm_inf(V1 const & vec)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      viennacl::ocl::handle<cl_mem> h = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(cl_uint));
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<value_type, ALIGNMENT>::program_name(), "index_norm_inf");
      //cl_uint size = static_cast<cl_uint>(vcl_vec.internal_size());

      k.global_work_size(0, k.local_work_size());
      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)),                                 
                               viennacl::ocl::local_mem(sizeof(value_type) * k.local_work_size()),
                               viennacl::ocl::local_mem(sizeof(cl_uint) * k.local_work_size()), h));
      
      //read value:
      cl_uint result;
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), h.get(), CL_TRUE, 0, sizeof(cl_uint), &result, 0, NULL, NULL);
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
    template <typename V1, typename V2, typename SCALARTYPE>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value
                                  && viennacl::is_cpu_scalar<SCALARTYPE>::value
                                >::type
    plane_rotation(V1 & vec1,
                   V2 & vec2,
                   SCALARTYPE alpha,
                   SCALARTYPE beta)
    {
      typedef typename viennacl::result_of::cpu_value_type<V1>::type        value_type;
      
      //TODO: Ensure that correct alignment is chosen for the kernels.
      const unsigned int ALIGNMENT = V1::alignment;
      
      assert(viennacl::traits::size(vec1) == viennacl::traits::size(vec2));
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "plane_rotation");

      viennacl::ocl::enqueue(k(viennacl::traits::handle(vec1), cl_uint(viennacl::traits::start(vec1)), cl_uint(viennacl::traits::size(vec1)),                                 
                               viennacl::traits::handle(vec2), cl_uint(viennacl::traits::start(vec2)), cl_uint(viennacl::traits::size(vec2)),                                 
                               alpha,
                               beta)
                            );
    }
    
  } //namespace linalg
} //namespace viennacl


#endif
