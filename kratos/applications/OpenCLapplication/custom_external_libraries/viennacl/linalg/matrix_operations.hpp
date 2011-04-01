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

#ifndef _VIENNACL_MATRIX_OPERATIONS_HPP_
#define _VIENNACL_MATRIX_OPERATIONS_HPP_

/** @file matrix_operations.hpp
    @brief Implementations of dense matrix related operations. also matrix-vector products.
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/matrix_kernel_class_deducer.hpp"
#include "viennacl/tools/matrix_prod_kernel_class_deducer.hpp"
#include "viennacl/linalg/kernels/vector_kernels.h"
#include "viennacl/linalg/kernels/matrix_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_col_kernels.h"

#include "viennacl/linalg/kernels/matrix_prod_col_col_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_col_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_row_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_col_row_row_kernels.h"

#include "viennacl/linalg/kernels/matrix_prod_row_col_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_col_row_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_row_col_kernels.h"
#include "viennacl/linalg/kernels/matrix_prod_row_row_row_kernels.h"

namespace viennacl
{
  namespace linalg
  {
    
    /** @brief Adds two dense matrices and writes the result to a third matrix
    *
    * This is the implementation of the convenience expression result = mat1 + mat2;
    *
    * @param mat1   The left hand side operand
    * @param mat2   The right hand side operand
    * @param result The resulting matrix
    */
    template<class TYPE, typename F, unsigned int ALIGNMENT>
    void add(const viennacl::matrix<TYPE, F, ALIGNMENT> & mat1, 
             const viennacl::matrix<TYPE, F, ALIGNMENT> & mat2,
             viennacl::matrix<TYPE, F, ALIGNMENT> & result)
    {
      assert(mat1.rows() == mat2.rows());
      assert(mat1.columns() == mat2.columns());
      result.resize(mat1.rows(), mat1.columns());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<TYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "add");
      unsigned int size = std::min(mat1.internal_size(), mat2.internal_size());

      viennacl::ocl::enqueue(k(mat1, mat2, result, size));        
    }


    // A * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication
    *
    * This is used for the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec)
    {
      return viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                         const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                         op_prod >(mat, vec);
    }
    
    /** @brief Carries out matrix-vector multiplication
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    */
    template<class TYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(const viennacl::matrix<TYPE, F, ALIGNMENT> & mat, 
                    const viennacl::vector<TYPE, VECTOR_ALIGNMENT> & vec, 
                          viennacl::vector<TYPE, VECTOR_ALIGNMENT> & result)
    {
      assert(mat.size2() == vec.size());
      // Inplace matrix-vector products like x = prod(A, x) are currently illegal: Introduce a temporary like y = prod(A, x); x = y; instead
      assert(vec.handle() != result.handle() && "No direct inplace matrix-vector product possible. Introduce a temporary!");
      result.resize(mat.size1());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<TYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "vec_mul");
      viennacl::ocl::enqueue(
                             k(mat, mat.size1(), mat.size2(), mat.internal_size1(), mat.internal_size2(), vec, result));    
    }



    // trans(A) * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication with a transposed matrix
    *
    * This is used for the convenience expression result = trans(mat) * vec;
    *
    * @param proxy  The transposed matrix proxy
    * @param vec    The vector
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                   const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                   op_trans>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                                       const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                                       op_trans> & proxy, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec)
    {
      return viennacl::vector_expression<const viennacl::matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                            const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                            op_trans>,
                                         const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                         op_prod >(proxy, vec);
    }

    /** @brief Unwraps the transposed matrix proxy and forwards to trans_prod_impl()
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(const viennacl::matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                      const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                      op_trans> & mat,
                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                          viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & result)
    {
      trans_prod_impl(mat.lhs(), vec, result);
    }
    
    /** @brief Carries out matrix-vector multiplication with a transposed matrix
    *
    * Implementation of the convenience expression result = trans(mat) * vec;
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void trans_prod_impl(const matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                          const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & result)
    {
      assert(mat.size1() == vec.size());  //remember: mat is transposed!
      // Inplace matrix-vector products like x = prod(A, x) are currently illegal: Introduce a temporary like y = prod(A, x); x = y; instead
      assert(vec.handle() != result.handle() && "No direct inplace matrix-vector product possible. Introduce a temporary!");
      result.resize(mat.size2());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "trans_vec_mul");
      
      viennacl::ocl::enqueue(k(mat, mat.size1(), mat.size2(), mat.internal_size1(), mat.internal_size2(), vec, result));        
    }



    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(A, B);
    *
    */
    template<class TYPE, typename F1, typename F2, typename F3, unsigned int ALIGNMENT>
    void prod_impl(const viennacl::matrix<TYPE, F1, ALIGNMENT> & A, 
                    const viennacl::matrix<TYPE, F2, ALIGNMENT> & B, 
                          viennacl::matrix<TYPE, F3, ALIGNMENT> & C, 
                          int block_size = 15) // [JW] added ability to set block size from outside ..
    {
      assert(A.size1() == C.size1());
      assert(A.size2() == B.size1());
      assert(B.size2() == C.size2());
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(C.handle() != A.handle() 
             && C.handle() != B.handle()
             && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      typedef typename viennacl::tools::MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<TYPE, F1, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F2, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F3, ALIGNMENT> >::ResultType    KernelClass;
      KernelClass::init();
      
      //std::cout << "KernelClass::program_name() : " << KernelClass::program_name() << std::endl;
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "prod_AA");
      
      /*k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size1() / 2, block_size / 2));
      k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size2() / 2, block_size / 2));
      k.local_work_size(0, block_size / 2);
      k.local_work_size(1, block_size / 2);*/
      
      k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size1(), block_size));
      k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size2(), block_size));
      k.local_work_size(0, block_size);
      k.local_work_size(1, block_size);
      
      viennacl::ocl::enqueue(
                             k(A, A.size1(), A.size2(), A.internal_size1(), A.internal_size2(),
                               B, B.size1(), B.size2(), B.internal_size1(), B.internal_size2(),
                               C, C.size1(), C.size2(), C.internal_size1(), C.internal_size2(),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size) ));        
    }


    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(trans(A), B);
    *
    */
    template<class TYPE, typename F1, typename F2, typename F3, unsigned int ALIGNMENT>
    void prod_impl(const viennacl::matrix_expression< const matrix<TYPE, F1, ALIGNMENT>,
                                                      const matrix<TYPE, F1, ALIGNMENT>,
                                                      op_trans> & A, 
                    const viennacl::matrix<TYPE, F2, ALIGNMENT> & B, 
                          viennacl::matrix<TYPE, F3, ALIGNMENT> & C)
    {
      assert(A.size2() == C.size1());
      assert(A.size1() == B.size1());
      assert(B.size2() == C.size2());
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(C.handle() != A.lhs().handle() 
             && C.handle() != B.handle()
             && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      int block_size = 15;

      typedef typename viennacl::tools::MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<TYPE, F1, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F2, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F3, ALIGNMENT> >::ResultType    KernelClass;
      KernelClass::init();
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "prod_TA");
      
      k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size1(), block_size));
      k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size2(), block_size));
      
      k.local_work_size(0, block_size);
      k.local_work_size(1, block_size);
      viennacl::ocl::enqueue(
                             k(A.lhs(), A.lhs().size1(), A.lhs().size2(), A.lhs().internal_size1(), A.lhs().internal_size2(),
                               B, B.size1(), B.size2(), B.internal_size1(), B.internal_size2(),
                               C, C.size1(), C.size2(), C.internal_size1(), C.internal_size2(),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size) ));        
    }


    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(A, trans(B));
    *
    */
    template<class TYPE, typename F1, typename F2, typename F3, unsigned int ALIGNMENT>
    void prod_impl(const viennacl::matrix<TYPE, F1, ALIGNMENT> & A, 
                   const viennacl::matrix_expression< const matrix<TYPE, F2, ALIGNMENT>,
                                                      const matrix<TYPE, F2, ALIGNMENT>,
                                                      op_trans> & B,
                   viennacl::matrix<TYPE, F3, ALIGNMENT> & C)
    {
      assert(A.size1() == C.size1());
      assert(A.size2() == B.size2());
      assert(B.size1() == C.size2());
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(C.handle() != A.handle() 
             && C.handle() != B.lhs().handle()
             && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      int block_size = 15;

      typedef typename viennacl::tools::MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<TYPE, F1, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F2, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F3, ALIGNMENT> >::ResultType    KernelClass;
      KernelClass::init();
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "prod_AT");
      
      k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size1(), block_size));
      k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size2(), block_size));
      
      k.local_work_size(0, block_size);
      k.local_work_size(1, block_size);
      viennacl::ocl::enqueue(
                             k(A, A.size1(), A.size2(), A.internal_size1(), A.internal_size2(),
                               B.lhs(), B.lhs().size1(), B.lhs().size2(), B.lhs().internal_size1(), B.lhs().internal_size2(),
                               C, C.size1(), C.size2(), C.internal_size1(), C.internal_size2(),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size) ));        
    }



    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(trans(A), trans(B));
    *
    */
    template<class TYPE, typename F1, typename F2, typename F3, unsigned int ALIGNMENT>
    void prod_impl(const viennacl::matrix_expression< const matrix<TYPE, F1, ALIGNMENT>,
                                                      const matrix<TYPE, F1, ALIGNMENT>,
                                                      op_trans> & A,
                   const viennacl::matrix_expression< const matrix<TYPE, F2, ALIGNMENT>,
                                                      const matrix<TYPE, F2, ALIGNMENT>,
                                                      op_trans> & B,
                   viennacl::matrix<TYPE, F3, ALIGNMENT> & C)
    {
      assert(A.size2() == C.size1());
      assert(A.size1() == B.size2());
      assert(B.size1() == C.size2());
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(C.handle() != A.lhs().handle() 
             && C.handle() != B.lhs().handle()
             && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      int block_size = 15;

      typedef typename viennacl::tools::MATRIX_PROD_KERNEL_CLASS_DEDUCER< viennacl::matrix<TYPE, F1, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F2, ALIGNMENT>,
                                                                          viennacl::matrix<TYPE, F3, ALIGNMENT> >::ResultType    KernelClass;
      KernelClass::init();
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "prod_TT");
      
      k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size1(), block_size));
      k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(C.size2(), block_size));
      
      k.local_work_size(0, block_size);
      k.local_work_size(1, block_size);
      viennacl::ocl::enqueue(
                             k(A.lhs(), A.lhs().size1(), A.lhs().size2(), A.lhs().internal_size1(), A.lhs().internal_size2(),
                               B.lhs(), B.lhs().size1(), B.lhs().size2(), B.lhs().internal_size1(), B.lhs().internal_size2(),
                               C, C.size1(), C.size2(), C.internal_size1(), C.internal_size2(),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size),
                               viennacl::ocl::local_mem(sizeof(TYPE) * block_size * block_size) ));        
    }











    /** @brief Returns a proxy class for the operation mat += vec1 * vec2^T, i.e. a rank 1 update
    *
    * @param vec1    The first vector
    * @param vec2    The second vector
    */
    template<class SCALARTYPE, unsigned int VA1, unsigned int VA2>
    viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VA1>,
                                 const viennacl::vector<SCALARTYPE, VA2>,
                                 op_prod> outer_prod(const viennacl::vector<SCALARTYPE, VA1> & vec1, 
                                                     const viennacl::vector<SCALARTYPE, VA2> & vec2)
    {
      return viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VA1>,
                                          const viennacl::vector<SCALARTYPE, VA2>,
                                          op_prod>(vec1, vec2);
    }
    
    

    /** @brief The implementation of the operation mat += vec1 * vec2^T, i.e. a rank 1 update
    *
    * Implementation of the convenience expression result += outer_prod(vec1, vec2);
    *
    * @param mat1    The matrix to be updated
    * @param vec1    The first vector
    * @param vec2    The second vector
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void rank_1_update(viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat1, 
                       const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                       const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2)
    {
      assert(mat1.size1() == vec1.size());
      assert(mat1.size2() == vec2.size());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "rank1_update");

      viennacl::ocl::enqueue(k(mat1, mat1.size1(), mat1.size2(), mat1.internal_size1(), mat1.internal_size2(), vec1, vec2));        
    }
    
    
    /** @brief The implementation of the operation mat += alpha * vec1 * vec2^T, i.e. a scaled rank 1 update
    *
    * Implementation of the convenience expression result += alpha * outer_prod(vec1, vec2);
    *
    * @param mat1    The matrix to be updated
    * @param val     The scaling factor
    * @param vec1    The first vector
    * @param vec2    The second vector
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void scaled_rank_1_update(viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat1,
                              SCALARTYPE val,
                              const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                              const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2)
    {
      assert(mat1.size1() == vec1.size());
      assert(mat1.size2() == vec2.size());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "scaled_rank1_update");

      viennacl::ocl::enqueue(k(mat1, mat1.size1(), mat1.size2(), mat1.internal_size1(), mat1.internal_size2(), 
                                                           val, vec1, vec2));        
    }
    
  } //namespace linalg


    //v = A * x
    /** @brief Implementation of the operation v1 = A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                          const viennacl::vector<SCALARTYPE, ALIGNMENT>,
                                                                                          viennacl::op_prod> & proxy) 
    {
      // check for the special case x = A * x
      if (proxy.rhs().handle() == this->handle())
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> result(proxy.rhs().size());
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
        *this = result;
        return *this;
      }
      else
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      return *this;
    }

    //v += A * x
    /** @brief Implementation of the operation v1 += A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+=(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.lhs().size1());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      *this += result;
      return *this;
    }

    /** @brief Implementation of the operation v1 -= A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-=(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.lhs().size1());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      *this -= result;
      return *this;
    }
    
    
    //free functions:
    /** @brief Implementation of the operation 'result = v1 + A * v2', where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      result += *this;
      return result;
    }

    /** @brief Implementation of the operation 'result = v1 - A * v2', where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      result = *this - result;
      return result;
    }


    ////////// transposed_matrix_proxy


    //v = trans(A) * x
    /** @brief Implementation of the operation v1 = A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                                   const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                                   op_trans>,
                                                                                          const viennacl::vector<SCALARTYPE, ALIGNMENT>,
                                                                                          viennacl::op_prod> & proxy) 
    {
      // check for the special case x = trans(A) * x
      if (proxy.rhs().handle() == this->handle())
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> result(proxy.rhs().size());
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
        *this = result;
        return *this;
      }
      else
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      return *this;
    }

    //v += A * x
    /** @brief Implementation of the operation v1 += A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+=(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                          const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                          op_trans>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.lhs().size1());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      *this += result;
      return *this;
    }

    /** @brief Implementation of the operation v1 -= A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-=(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                          const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                          op_trans>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.lhs().size1());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      *this -= result;
      return *this;
    }
    
    
    //free functions:
    /** @brief Implementation of the operation 'result = v1 + A * v2', where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                         const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                         op_trans>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      result += *this;
      return result;
    }

    /** @brief Implementation of the operation 'result = v1 - A * v2', where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <typename F, unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                         const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                         op_trans>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      result = *this - result;
      return result;
    }


} //namespace viennacl


#endif
