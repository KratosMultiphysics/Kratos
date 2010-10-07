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

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_VECTOR_OPERATIONS_HPP_
#define _VIENNACL_VECTOR_OPERATIONS_HPP_

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
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
             const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
             viennacl::vector<SCALARTYPE, ALIGNMENT> & result, 
             size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The addend
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                     const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
                     size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Subtraction of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * result = vec1 - vec2
    *
    * @param vec1  The first operand. 
    * @param vec2  The second operand.
    * @param result The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void sub(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
             const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
             viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
             size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sub.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Inplace addition of two vectors. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2.
    * 
    * @param vec1  The result. 
    * @param vec2  The subtracted vector
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                     const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                     size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_sub.work_groups() * NUM_THREADS, NUM_THREADS);
    }


    //result = vec * scalar
    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec * alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param result The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mult(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
              scalar<SCALARTYPE> const & alpha,
              viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
              size_t NUM_THREADS = 0)
    {
      result.resize(vec.size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::mult.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec * alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param result The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mult(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
              SCALARTYPE alpha,
              viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
              size_t NUM_THREADS = 0)
    {
      result.resize(vec.size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.setArgument(pos++, alpha);
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_mult.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mult(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                      scalar<SCALARTYPE> const & alpha,
                      size_t NUM_THREADS = 0)
    {
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_mult.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a cpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The scaling factor.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mult(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                      SCALARTYPE alpha,
                      size_t NUM_THREADS = 0)
    {
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.setArgument(pos++, alpha);
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::cpu_inplace_mult.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    //result = vec / scalar
    /** @brief Scales a vector. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = vec / alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The (inverse) scaling factor.
    * @param result The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void divide(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                scalar<SCALARTYPE> const & alpha,
                viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
                size_t NUM_THREADS = 0)
    {
      result.resize(vec.size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::divide.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Scales a vector inplace. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result *= alpha, where alpha is a gpu scalar
    *
    * @param vec    The vector to be scaled.
    * @param alpha  The (inverse) scaling factor.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_divide(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec,
                        scalar<SCALARTYPE> const & alpha,
                        size_t NUM_THREADS = 0)
    {
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.setArgument(pos++, vec.internal_size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.start1D(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::inplace_divide.work_groups() * NUM_THREADS, NUM_THREADS);
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
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 scalar<SCALARTYPE> const & alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
                 size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = alpha * vec1 + vec2, where alpha is a cpu scalar
    *
    * @param vec1    The first added
    * @param alpha   The scaling factor for the first addend.
    * @param vec2    The second added.
    * @param result  The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_add(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 SCALARTYPE alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
                 size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.setArgument(pos++, alpha);
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_mul_add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    //vec1 += factor * vec2
    /** @brief Inplace Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += alpha * vec2, where alpha is a gpu scalar
    *
    * @param vec1    The first added
    * @param alpha   The scaling factor for the first addend.
    * @param vec2    The second added.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha,
                         size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_add.work_groups() * NUM_THREADS, NUM_THREADS);

    }

    /** @brief Inplace Multiply-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += alpha * vec2, where alpha is a cpu scalar
    *
    * @param vec1    The first added
    * @param vec2    The second added.
    * @param alpha   The scaling factor for the first addend.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         SCALARTYPE alpha,
                         size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.setArgument(pos++, alpha);
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::cpu_inplace_mul_add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Multiply-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes result = alpha * vec1 - vec2, where alpha is a gpu scalar
    *
    * @param vec1    The first vector operand
    * @param alpha   The scaling factor for the first vector.
    * @param vec2    The second operand.
    * @param result  The result vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void mul_sub(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                 scalar<SCALARTYPE> const & alpha,
                 const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                 viennacl::vector<SCALARTYPE, ALIGNMENT> & result,
                 size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      result.resize(vec1.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::mul_sub.work_groups() * NUM_THREADS, NUM_THREADS);
    }


    /** @brief Inplace Multiply-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= alpha * vec2, where alpha is a gpu scalar
    *
    * @param vec1    The result vector which is updated
    * @param alpha   The scaling factor for the vector update.
    * @param vec2    The second operand.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_mul_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha,
                         size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_mul_sub.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Inplace divide-add operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 += vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_div_add(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha,
                         size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_add.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Inplace divide-subtract operation. Try to use the overloaded operators for vector instead, unless you want to fine-tune the number of GPU threads involved.
    *
    * Computes vec1 -= vec2 / alpha, where alpha is a gpu scalar
    *
    * @param vec1    The first vector
    * @param vec2    The vector update
    * @param alpha   The scaling factor for the second vector.
    * @param NUM_THREADS The number of threads per work group.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inplace_div_sub(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                         const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                         scalar<SCALARTYPE> const & alpha,
                         size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.setArgument(pos++, alpha.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inplace_div_sub.work_groups() * NUM_THREADS, NUM_THREADS);
    }


    ///////////////////////// Norms and inner product ///////////////////


    //implementation of inner product:
    //namespace {
      /** @brief Computes the inner product of two vectors - implementation. Library users should call inner_prod(vec1, vec2).
      *
      * @param vec1 The first vector
      * @param vec2 The second vector
      * @param result The result scalar (on the gpu)
      * @param NUM_THREADS The number of threads per work group
      */
      template<class SCALARTYPE, unsigned int ALIGNMENT>
      void inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                           const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                           scalar<SCALARTYPE> & result,
                           size_t NUM_THREADS /* see forwards.h */)
      {
        assert(vec1.size() == vec2.size());
        unsigned int size = std::min(vec1.internal_size(), vec2.internal_size());
        unsigned int pos = 0;
        unsigned int group_num = 8;   //Note: group_num MUST be a power of two!
        viennacl::vector<SCALARTYPE> temp(group_num);
        
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setArgument(pos++, vec1.handle());
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setArgument(pos++, vec2.handle());
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setArgument(pos++, size);
        if (NUM_THREADS == 0)
          viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.work_items_per_group()));
        else
          viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*NUM_THREADS));
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.setArgument(pos++, temp.handle());

        if (NUM_THREADS == 0)
          viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.start1D(group_num *
                                                        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.work_items_per_group(),
                                                        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.work_items_per_group());
        else
          viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::inner_prod.start1D(group_num * NUM_THREADS, NUM_THREADS);
        
        pos = 0;
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sum.setArgument(pos++, temp.handle());
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sum.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sum.setArgument(pos++, group_num);
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::sum.start1D(group_num, group_num);
      }
    //}

    //public interface of inner product
    /** @brief Computes the inner product of two vectors.
    *
    * @param vec1 The first vector
    * @param vec2 The second vector
    * @param NUM_THREADS The number of threads per work group
    * @return The result
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT1, unsigned int ALIGNMENT2>
    viennacl::scalar_expression< const viennacl::vector<SCALARTYPE, ALIGNMENT1>, 
                                 const viennacl::vector<SCALARTYPE, ALIGNMENT2>,
                                 viennacl::op_inner_prod >
    inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT1> & vec1,
                    const viennacl::vector<SCALARTYPE, ALIGNMENT2> & vec2,
                    size_t NUM_THREADS /* declared in forwards.h */)
    {
      return viennacl::scalar_expression< const viennacl::vector<SCALARTYPE, ALIGNMENT1>, 
                                          const viennacl::vector<SCALARTYPE, ALIGNMENT2>,
                                          viennacl::op_inner_prod >(vec1, vec2);
    }


    
    /** @brief Computes the l^1-norm of a vector
    *
    * @param vcl_vec The vector
    * @param NUM_THREADS The number of threads per work group
    * @return The result
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    SCALARTYPE norm_1(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                      size_t NUM_THREADS = 0)
    {
      scalar<SCALARTYPE> result;
      unsigned int size = vcl_vec.internal_size();
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.setArgument(pos++, vcl_vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.setArgument(pos++, size);
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.work_items_per_group()));
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*NUM_THREADS));
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.setArgument(pos++, result.handle());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.work_items_per_group(),
                                                                                viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.work_items_per_group());
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_1.start1D(NUM_THREADS,NUM_THREADS);
      
      return result;
    }

    /** @brief Computes the l^2-norm of a vector - implementation
    *
    * @param vcl_vec The vector
    * @param result The result scalar
    * @param NUM_THREADS The number of threads per work group
    * @return The result
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void norm_2_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                     scalar<SCALARTYPE> & result,
                     size_t NUM_THREADS /* declared in forwards.h */)
    {
      unsigned int size = vcl_vec.internal_size();
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.setArgument(pos++, vcl_vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.setArgument(pos++, size);
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.work_items_per_group()));
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*NUM_THREADS));
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.setArgument(pos++, result.handle());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.work_items_per_group(),
                                                                                viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.work_items_per_group());
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_2.start1D(NUM_THREADS, NUM_THREADS);
    }

    /** @brief Computes the supremum-norm of a vector
    *
    * @param vcl_vec The vector
    * @param NUM_THREADS The number of threads per work group
    * @return The result
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    SCALARTYPE norm_inf(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                        size_t NUM_THREADS = 0)
    {
      scalar<SCALARTYPE> result;
      unsigned int size = vcl_vec.internal_size();
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.setArgument(pos++, vcl_vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.setArgument(pos++, size);
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.work_items_per_group()));
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*NUM_THREADS));
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.setArgument(pos++, result.handle());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.work_items_per_group(),
                                                                                  viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.work_items_per_group());
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::norm_inf.start1D(NUM_THREADS, NUM_THREADS);
      
      return result;
    }

    //This function should return a CPU scalar, otherwise statements like 
    // vcl_rhs[index_norm_inf(vcl_rhs)] 
    // are ambiguous
    /** @brief Computes the index of the first entry that is equal to the supremum-norm in modulus.
    *
    * @param vcl_vec The vector
    * @param NUM_THREADS The number of threads per work group
    * @return The result. Note that the result must be a CPU scalar (unsigned int), since gpu scalars are floating point types.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    cl_uint index_norm_inf(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                           size_t NUM_THREADS = 0)
    {
      viennacl::ocl::handle<cl_mem> h = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(cl_uint));
      
      unsigned int size = vcl_vec.internal_size();
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setArgument(pos++, vcl_vec.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setArgument(pos++, size);
      if (NUM_THREADS == 0)
      {
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.work_items_per_group()));
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(cl_uint)*viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.work_items_per_group()));
      }
      else
      {
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE)*NUM_THREADS));
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(cl_uint)*NUM_THREADS));
      }
        
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.setArgument(pos++, h);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.start1D(
                                                  viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.work_items_per_group(),
                                                  viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.work_items_per_group());
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::index_norm_inf.start1D(NUM_THREADS, NUM_THREADS);
      
      //read value:
      cl_uint result;
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), h.get(), CL_TRUE, 0, sizeof(cl_uint), &result, 0, NULL, NULL);
      assert(err == CL_SUCCESS);
      viennacl::ocl::finish();
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
    * @param NUM_THREADS  The number of threads per work group
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void plane_rotation(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                        const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                        SCALARTYPE alpha,
                        SCALARTYPE beta,
                        size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());
      unsigned int size = vec1.size();
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.setArgument(pos++, alpha);
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.setArgument(pos++, beta);
      viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.start1D();
      else
        viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.start1D(viennacl::linalg::kernels::vector<SCALARTYPE,ALIGNMENT>::plane_rotation.work_groups() * NUM_THREADS, NUM_THREADS);
    }
    
  } //namespace linalg
} //namespace viennacl


#endif
