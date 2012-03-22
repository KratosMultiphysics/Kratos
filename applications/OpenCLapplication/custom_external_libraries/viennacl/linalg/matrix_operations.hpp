#ifndef VIENNACL_MATRIX_OPERATIONS_HPP_
#define VIENNACL_MATRIX_OPERATIONS_HPP_

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
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/predicate.hpp"
#include "viennacl/meta/result_of.hpp"
#include "viennacl/traits/size.hpp"
#include "viennacl/traits/start.hpp"
#include "viennacl/traits/handle.hpp"
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
    //
    ///////////////////////////////////// addition and subtraction///////////////////////////////////////////////
    //
    
    namespace detail
    {
      template<class T1, class T2, class T3>
      typename viennacl::enable_if<   viennacl::is_matrix<T1>::value 
                                   && viennacl::is_matrix<T2>::value 
                                   && viennacl::is_matrix<T3>::value >::type
      add_sub_impl(const T1 & mat1, 
                   const T2 & mat2,
                         T3 & result,
                   std::string kernel_name
                  )
      {
        assert(result.size1() == mat1.size1());
        assert(result.size2() == mat1.size2());
        assert(result.size1() == mat2.size1());
        assert(result.size2() == mat2.size2());

        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< T1 >::ResultType    KernelClass;
        
        std::size_t block_size = 16;
        
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), kernel_name);
        k.global_work_size(0, block_size*block_size);
        k.global_work_size(1, block_size*block_size);
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        viennacl::ocl::enqueue(k(viennacl::traits::handle(mat1), 
                                        cl_uint(viennacl::traits::start1(mat1)),           cl_uint(viennacl::traits::start2(mat1)), 
                                        cl_uint(viennacl::traits::size1(mat1)),            cl_uint(viennacl::traits::size2(mat1)),
                                        cl_uint(viennacl::traits::internal_size1(mat1)),   cl_uint(viennacl::traits::internal_size2(mat1)),
                                viennacl::traits::handle(mat2), 
                                        cl_uint(viennacl::traits::start1(mat2)),           cl_uint(viennacl::traits::start2(mat2)), 
                                        cl_uint(viennacl::traits::size1(mat2)),            cl_uint(viennacl::traits::size2(mat2)),
                                        cl_uint(viennacl::traits::internal_size1(mat2)),   cl_uint(viennacl::traits::internal_size2(mat2)),
                                viennacl::traits::handle(result), 
                                        cl_uint(viennacl::traits::start1(result)),         cl_uint(viennacl::traits::start2(result)), 
                                        cl_uint(viennacl::traits::size1(result)),          cl_uint(viennacl::traits::size2(result)),
                                        cl_uint(viennacl::traits::internal_size1(result)), cl_uint(viennacl::traits::internal_size2(result))
                                )
                              );        
      }
      


      template <typename T1, typename T2>
      typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                    && viennacl::is_matrix<T2>::value
                                  >::type
      inplace_add_sub_impl(T1 & result, T2 const & mat2, std::string kernel_name)
      {
        assert(viennacl::traits::size1(result) == viennacl::traits::size1(mat2));
        assert(viennacl::traits::size2(result) == viennacl::traits::size2(mat2));

        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< T1 >::ResultType    KernelClass;
        
        std::size_t block_size = 16;
        
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), kernel_name);
        k.global_work_size(0, block_size*block_size);
        k.global_work_size(1, block_size*block_size);
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        
        viennacl::ocl::enqueue(k(viennacl::traits::handle(result),
                                        cl_uint(viennacl::traits::start1(result)),         cl_uint(viennacl::traits::start2(result)), 
                                        cl_uint(viennacl::traits::size1(result)),          cl_uint(viennacl::traits::size2(result)),
                                        cl_uint(viennacl::traits::internal_size1(result)), cl_uint(viennacl::traits::internal_size2(result)),
                                viennacl::traits::handle(mat2), 
                                        cl_uint(viennacl::traits::start1(mat2)),            cl_uint(viennacl::traits::start2(mat2)), 
                                        cl_uint(viennacl::traits::size1(mat2)),             cl_uint(viennacl::traits::size2(mat2)),
                                        cl_uint(viennacl::traits::internal_size1(mat2)),    cl_uint(viennacl::traits::internal_size2(mat2))
                                )
                              );
      }
      
    }
    
    /** @brief Adds two dense matrices or submatrices and writes the result to a third matrix or submatrix
    *
    * This is the implementation of the convenience expression result = mat1 + mat2;
    *
    * @param mat1   The left hand side operand
    * @param mat2   The right hand side operand
    * @param result The resulting matrix
    */
    template<class T1, class T2, class T3>
    typename viennacl::enable_if<   viennacl::is_matrix<T1>::value 
                                 && viennacl::is_matrix<T2>::value 
                                 && viennacl::is_matrix<T3>::value >::type
    add(const T1 & mat1, 
        const T2 & mat2,
              T3 & result)
    {
      detail::add_sub_impl(mat1, mat2, result, "add");
    }

    /** @brief Adds a dense matrix or submatrix to another
    *
    * This is the implementation of the convenience expression result += mat1;
    *
    * @param mat2   The addend (either a matrix or a matrix_range)
    * @param result The resulting matrix  (either a matrix or a matrix_range)
    */
    template <typename T1, typename T2>
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                >::type
    inplace_add(T1 & result, T2 const & mat2)
    {
      detail::inplace_add_sub_impl(result, mat2, "inplace_add");
    }



    /** @brief Subtracts two dense matrices or submatrices and writes the result to a third matrix or submatrix
    *
    * This is the implementation of the convenience expression result = mat1 - mat2;
    *
    * @param mat1   The left hand side operand
    * @param mat2   The right hand side operand
    * @param result The resulting matrix
    */
    template<class T1, class T2, class T3>
    typename viennacl::enable_if<   viennacl::is_matrix<T1>::value 
                                 && viennacl::is_matrix<T2>::value 
                                 && viennacl::is_matrix<T3>::value >::type
    sub(const T1 & mat1, 
        const T2 & mat2,
              T3 & result)
    {
      detail::add_sub_impl(mat1, mat2, result, "sub");
    }

    /** @brief Subtracts a dense matrix or submatrix from another
    *
    * This is the implementation of the convenience expression result -= mat1;
    *
    * @param mat2   The addend (either a matrix or a matrix_range)
    * @param result The resulting matrix  (either a matrix or a matrix_range)
    */
    template <typename T1, typename T2>
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                >::type
    inplace_sub(T1 & result, T2 const & mat2)
    {
      detail::inplace_add_sub_impl(result, mat2, "inplace_sub");
    }




    //
    /////////////////////////   inplace multiplication and division /////////////////////////////////
    //

    namespace detail
    {
      template <typename  T1, typename ScalarType>
      typename viennacl::enable_if< viennacl::is_matrix<T1>::value >::type
      inplace_mult_div_impl(T1 & result, 
                            ScalarType val,
                            std::string kernel_name)
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< T1 >::ResultType    KernelClass;
        
        std::size_t block_size = 16;
          
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), kernel_name);
        
        k.global_work_size(0, block_size*block_size);
        k.global_work_size(1, block_size*block_size);
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        
        viennacl::ocl::enqueue(k(viennacl::traits::handle(result),
                                        cl_uint(viennacl::traits::start1(result)),         cl_uint(viennacl::traits::start2(result)), 
                                        cl_uint(viennacl::traits::size1(result)),          cl_uint(viennacl::traits::size2(result)),
                                        cl_uint(viennacl::traits::internal_size1(result)), cl_uint(viennacl::traits::internal_size2(result)),
                                val)
                              );
      }
    }


    /** @brief Multiplies a dense matrix or submatrix by a scalar
    *
    * This is the implementation of the convenience expression matrix *= val;
    *
    * @param result The matrix to be manipulated
    * @param val    The CPU scalar by which all entries of the matrix are multiplied
    */
    template <typename  T1>
    typename viennacl::enable_if< viennacl::is_matrix<T1>::value >::type
    inplace_mult(T1 & result, 
                 typename viennacl::result_of::cpu_value_type< typename T1::value_type >::type val)
    {
      detail::inplace_mult_div_impl(result, val, "cpu_inplace_mult");
    }


    /** @brief Multiplies a dense matrix or submatrix by a scalar
    *
    * This is the implementation of the convenience expression matrix *= val;
    *
    * @param result The matrix to be manipulated
    * @param val    The scalar by which all entries of the matrix are multiplied
    */
    template <typename  T1>
    typename viennacl::enable_if< viennacl::is_matrix<T1>::value >::type
    inplace_mult(T1 & result, 
                 typename T1::value_type val)
    {
      detail::inplace_mult_div_impl(result, val, "inplace_mult");
    }



    /** @brief Divides a dense matrix or submatrix by a scalar
    *
    * This is the implementation of the convenience expression matrix /= val;
    *
    * @param result The matrix to be manipulated
    * @param val    The scalar by which all entries of the matrix are divided
    */
    template <typename  T1>
    typename viennacl::enable_if< viennacl::is_matrix<T1>::value >::type
    inplace_divide(T1 & result, 
                   typename T1::value_type val)
    {
      detail::inplace_mult_div_impl(result, val, "inplace_divide");
    }



    //
    /////////////////////////   matrix-vector products /////////////////////////////////
    //



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
      assert(vec.handle().get() != result.handle().get() && "No direct inplace matrix-vector product possible. Introduce a temporary!");
      result.resize(mat.size1());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<TYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "vec_mul");
      viennacl::ocl::enqueue(
                             k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()),
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2()), vec, result));    
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
      assert(vec.handle().get() != result.handle().get() && "No direct inplace matrix-vector product possible. Introduce a temporary!");
      result.resize(mat.size2());

      typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
      
      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "trans_vec_mul");
      
      viennacl::ocl::enqueue(k(mat, cl_uint(mat.size1()), cl_uint(mat.size2()),
                                    cl_uint(mat.internal_size1()), cl_uint(mat.internal_size2()), vec, result));        
    }





    //
    /////////////////////////   matrix-matrix products /////////////////////////////////
    //
    
    namespace detail
    {
      // C = A * B and possibly transposed variants
      template <typename T1, typename T2, typename T3 >
      void prod(const T1 & A, 
                const T2 & B, 
                T3 & C,
                std::string kernel_name,
                int block_size = 16) // [JW] added ability to set block size from outside ..
      {
        typename viennacl::result_of::cpu_value_type< typename T1::value_type >::type   cpu_value_type;
        
        typedef typename viennacl::tools::MATRIX_PROD_KERNEL_CLASS_DEDUCER< T1, T2, T3 >::ResultType    KernelClass;
        KernelClass::init();
        
        //std::cout << "KernelClass::program_name() : " << KernelClass::program_name() << std::endl;
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), kernel_name);
        
        k.global_work_size(0, viennacl::tools::roundUpToNextMultiple<unsigned int>(viennacl::traits::size1(C), block_size));
        k.global_work_size(1, viennacl::tools::roundUpToNextMultiple<unsigned int>(viennacl::traits::size2(C), block_size));
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        
        viennacl::ocl::enqueue(k(viennacl::traits::handle(A), 
                                        cl_uint(viennacl::traits::start1(A)),           cl_uint(viennacl::traits::start2(A)), 
                                        cl_uint(viennacl::traits::size1(A)),            cl_uint(viennacl::traits::size2(A)),
                                        cl_uint(viennacl::traits::internal_size1(A)),   cl_uint(viennacl::traits::internal_size2(A)),
                                 viennacl::traits::handle(B), 
                                        cl_uint(viennacl::traits::start1(B)),           cl_uint(viennacl::traits::start2(B)), 
                                        cl_uint(viennacl::traits::size1(B)),            cl_uint(viennacl::traits::size2(B)),
                                        cl_uint(viennacl::traits::internal_size1(B)),   cl_uint(viennacl::traits::internal_size2(B)),
                                 viennacl::traits::handle(C), 
                                        cl_uint(viennacl::traits::start1(C)),         cl_uint(viennacl::traits::start2(C)), 
                                        cl_uint(viennacl::traits::size1(C)),          cl_uint(viennacl::traits::size2(C)),
                                        cl_uint(viennacl::traits::internal_size1(C)), cl_uint(viennacl::traits::internal_size2(C)),
                                 viennacl::ocl::local_mem(sizeof(cpu_value_type) * (block_size+1) * block_size),
                                 viennacl::ocl::local_mem(sizeof(cpu_value_type) * (block_size+1) * block_size)
                                )
                              );        
      }
    }


    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(A, B);
    *
    */
    template <typename T1, typename T2, typename T3 >
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                  && viennacl::is_matrix<T3>::value
                                >::type
    prod_impl(const T1 & A, 
              const T2 & B, 
                    T3 & C, 
              int block_size = 16) // [JW] added ability to set block size from outside ..
    {
      assert(viennacl::traits::size1(A) == viennacl::traits::size1(C));
      assert(viennacl::traits::size2(A) == viennacl::traits::size1(B));
      assert(viennacl::traits::size2(B) == viennacl::traits::size2(C));
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(viennacl::traits::handle(C).get() != viennacl::traits::handle(A).get() 
            && viennacl::traits::handle(C).get() != viennacl::traits::handle(B).get()
            && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
        
      detail::prod(A, B, C, "prod_AA", block_size);
    }



    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(trans(A), B);
    *
    */
    template <typename T1, typename T2, typename T3 >
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                  && viennacl::is_matrix<T3>::value
                                >::type
    prod_impl(const viennacl::matrix_expression< const T1,
                                                 const T1,
                                                 op_trans> & A, 
              const T2 & B, 
                    T3 & C, 
              int block_size = 16)
    {
      //std::cout << "size2(A): " << viennacl::traits::size2(A.lhs()) << std::endl;
      //std::cout << "size1(C): " << viennacl::traits::size1(C) << std::endl;
      assert(viennacl::traits::size2(A.lhs()) == viennacl::traits::size1(C));
      assert(viennacl::traits::size1(A.lhs()) == viennacl::traits::size1(B));
      assert(viennacl::traits::size2(B) == viennacl::traits::size2(C));
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(viennacl::traits::handle(C).get() != viennacl::traits::handle(A.lhs()).get() 
            && viennacl::traits::handle(C).get() != viennacl::traits::handle(B).get()
            && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      detail::prod(A.lhs(), B, C, "prod_TA", block_size);
    }




    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(A, trans(B));
    *
    */
    template <typename T1, typename T2, typename T3 >
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                  && viennacl::is_matrix<T3>::value
                                >::type
    prod_impl(const T1 & A, 
              const viennacl::matrix_expression< const T2,
                                                 const T2,
                                                 op_trans> & B,
              T3 & C, 
              int block_size = 16)
    {
      assert(viennacl::traits::size1(A) == viennacl::traits::size1(C));
      assert(viennacl::traits::size2(A) == viennacl::traits::size2(B.lhs()));
      assert(viennacl::traits::size1(B.lhs()) == viennacl::traits::size2(C));
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(viennacl::traits::handle(C).get() != viennacl::traits::handle(A).get() 
            && viennacl::traits::handle(C).get() != viennacl::traits::handle(B.lhs()).get()
            && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      detail::prod(A, B.lhs(), C, "prod_AT", block_size);
    }



    /** @brief Carries out matrix-matrix multiplication
    *
    * Implementation of C = prod(trans(A), trans(B));
    *
    */
    template <typename T1, typename T2, typename T3 >
    typename viennacl::enable_if<    viennacl::is_matrix<T1>::value
                                  && viennacl::is_matrix<T2>::value
                                  && viennacl::is_matrix<T3>::value
                                >::type
    prod_impl(const viennacl::matrix_expression< const T1,
                                                 const T1,
                                                 op_trans> & A,
              const viennacl::matrix_expression< const T2,
                                                 const T2,
                                                 op_trans> & B,
              T3 & C, 
              int block_size = 16)
    {
      assert(viennacl::traits::size2(A.lhs()) == viennacl::traits::size1(C));
      assert(viennacl::traits::size1(A.lhs()) == viennacl::traits::size2(B.lhs()));
      assert(viennacl::traits::size1(B.lhs()) == viennacl::traits::size2(C));
      // Inplace matrix-vector products like B = prod(A, B) are currently illegal: Introduce a temporary like C = prod(A, B); B = C; instead
      assert(viennacl::traits::handle(C).get() != viennacl::traits::handle(A.lhs()).get() 
            && viennacl::traits::handle(C).get() != viennacl::traits::handle(B.lhs()).get()
            && "No direct inplace matrix-matrix product possible. Introduce a temporary!");
      
      detail::prod(A.lhs(), B.lhs(), C, "prod_TT", block_size);
    }




    //
    /////////////////////////   miscellaneous operations /////////////////////////////////
    //






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

      viennacl::ocl::enqueue(k(mat1, cl_uint(mat1.size1()), cl_uint(mat1.size2()),
                                     cl_uint(mat1.internal_size1()), cl_uint(mat1.internal_size2()), vec1, vec2));        
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

      viennacl::ocl::enqueue(k(mat1, cl_uint(mat1.size1()), cl_uint(mat1.size2()),
                                     cl_uint(mat1.internal_size1()), cl_uint(mat1.internal_size2()), 
                                                           val, vec1, vec2));        
    }
    
  } //namespace linalg




  //
  /////////////////////////  Operator overloads /////////////////////////////////
  //





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
    if (proxy.rhs().handle().get() == this->handle().get())
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
    if (proxy.rhs().handle().get() == this->handle().get())
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
