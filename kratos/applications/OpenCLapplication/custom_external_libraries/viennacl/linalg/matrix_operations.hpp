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

#ifndef _VIENNACL_MATRIX_OPERATIONS_HPP_
#define _VIENNACL_MATRIX_OPERATIONS_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/kernels/vector_kernels.h"
#include "viennacl/linalg/kernels/matrix_kernels.h"

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
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class TYPE, typename F, unsigned int ALIGNMENT>
    void add(const viennacl::matrix<TYPE, F, ALIGNMENT> & mat1, 
             const viennacl::matrix<TYPE, F, ALIGNMENT> & mat2,
             viennacl::matrix<TYPE, F, ALIGNMENT> & result,
             unsigned int NUM_THREADS = 0)
    {
      assert(mat1.rows() == mat2.rows());
      assert(mat1.columns() == mat2.columns());
      result.resize(mat1.rows(), mat1.columns());

      unsigned int size = std::min(mat1.internal_size(), mat2.internal_size());
      unsigned int pos = 0;
      viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::init();
      viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.setArgument(pos++, mat1.handle());
      viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.setArgument(pos++, mat2.handle());
      viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.setArgument(pos++, size);

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.start1D();
      else
        viennacl::linalg::kernels::vector<TYPE, ALIGNMENT>::add.start1D(viennacl::ocl::device().work_groups() * NUM_THREADS, NUM_THREADS);
    }


    // A * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication
    *
    * This is used for the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                                    unsigned int NUM_THREADS = 0)
    {
      return viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                         const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                         op_prod >(mat, vec);
    }
    
    /** @brief Carries out matrix-vector multiplication
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param matrix    The matrix
    * @param vec    The vector
    * @param result The result vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class TYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(const viennacl::matrix<TYPE, F, ALIGNMENT> & matrix, 
                    const viennacl::vector<TYPE, VECTOR_ALIGNMENT> & vec, 
                          viennacl::vector<TYPE, VECTOR_ALIGNMENT> & result, 
                    unsigned int NUM_THREADS = 0)
    {
      assert(matrix.size2() == vec.size());
      result.resize(matrix.size1());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, matrix.handle());
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, matrix.internal_size2());
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, vec.size());
      viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.setArgument(pos++, result.size());
      
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.start1D();
      else
        viennacl::linalg::kernels::matrix<TYPE,ALIGNMENT>::vec_mul.start1D(viennacl::ocl::device().work_groups() * NUM_THREADS, NUM_THREADS);
    }
    
    // trans(A) * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication with a transposed matrix
    *
    * This is used for the convenience expression result = trans(mat) * vec;
    *
    * @param proxy  The transposed matrix proxy
    * @param vec    The vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> & proxy, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                                    unsigned int NUM_THREADS = 0)
    {
      return viennacl::vector_expression<const viennacl::transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT>,
                                         const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                         op_prod >(proxy, vec);
    }

    /** @brief Unwraps the transposed matrix proxy and forwards to trans_prod_impl()
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(const viennacl::transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> & mat,
                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                          viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & result,
                    unsigned int NUM_THREADS = 0)
    {
      trans_prod_impl(mat.get_matrix(), vec, result, NUM_THREADS);
    }
    
    /** @brief Carries out matrix-vector multiplication with a transposed matrix
    *
    * Implementation of the convenience expression result = trans(mat) * vec;
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void trans_prod_impl(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat,
                          const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & result,
                          unsigned int NUM_THREADS = 0)
    {
      assert(mat.size1() == vec.size());  //remember: mat is transposed!
      result.resize(mat.size2());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, mat.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, vec.size());
      viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.setArgument(pos++, result.size());
      
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.start1D();
      else
        viennacl::linalg::kernels::matrix<SCALARTYPE,ALIGNMENT>::trans_vec_mul.start1D(viennacl::ocl::device().work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Returns a proxy class for the operation mat += vec1 * vec2^T, i.e. a rank 1 update
    *
    * Implementation of the convenience expression result = trans(mat) * vec;
    *
    * @param vec1    The first vector
    * @param vec2    The second vector
    */
    template<class SCALARTYPE, unsigned int VECTOR_ALIGNMENT>
    viennacl::outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT> outer_prod(const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec1, 
                                                                        const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec2)
    {
      return viennacl::outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT>(vec1, vec2);
    }
    
    

    /** @brief The implementation of the operation mat += vec1 * vec2^T, i.e. a rank 1 update
    *
    * Implementation of the convenience expression result += outer_prod(vec1, vec2);
    *
    * @param mat1    The matrix to be updated
    * @param vec1    The first vector
    * @param vec2    The second vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void rank_1_update(viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat1, 
                       const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                       const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
                       unsigned int NUM_THREADS = 0)
    {
      assert(mat1.size1() == vec1.size());
      assert(mat1.size2() == vec2.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, mat1.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, mat1.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, vec1.size());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.setArgument(pos++, vec2.size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.start1D(viennacl::ocl::device().work_items_per_group(),
                                                                                       viennacl::ocl::device().work_items_per_group());
      else
        viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::rank1_update.start1D(NUM_THREADS, NUM_THREADS);
    }
    
    
    /** @brief The implementation of the operation mat += alpha * vec1 * vec2^T, i.e. a scaled rank 1 update
    *
    * Implementation of the convenience expression result += alpha * outer_prod(vec1, vec2);
    *
    * @param mat1    The matrix to be updated
    * @param val     The scaling factor
    * @param vec1    The first vector
    * @param vec2    The second vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void scaled_rank_1_update(viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat1,
                              SCALARTYPE val,
                              const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1, 
                              const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
                              unsigned int NUM_THREADS = 0)
    {
      assert(mat1.size1() == vec1.size());
      assert(mat1.size2() == vec2.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, mat1.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, val);
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, vec1.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, vec2.handle());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, mat1.internal_size2());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, vec1.size());
      viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.setArgument(pos++, vec2.size());

      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.start1D(viennacl::ocl::device().work_items_per_group(),
                                                                                              viennacl::ocl::device().work_items_per_group());
      else
        viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::scaled_rank1_update.start1D(NUM_THREADS, NUM_THREADS);
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
      if (proxy.get_lhs().handle() == this->handle())
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> result(proxy.get_rhs().size());
        viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
        *this = result;
        return *this;
      }
      else
      {
        viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), *this);
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
      vector<SCALARTYPE, ALIGNMENT> result(proxy.get_lhs().size1());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
      vector<SCALARTYPE, ALIGNMENT> result(proxy.get_lhs().size1());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
      assert(proxy.get_lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
      assert(proxy.get_lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const transposed_matrix_proxy<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                          const viennacl::vector<SCALARTYPE, ALIGNMENT>,
                                                                                          viennacl::op_prod> & proxy) 
    {
      // check for the special case x = trans(A) * x
      if (proxy.get_lhs().get_matrix().handle() == this->handle())
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> result(proxy.get_rhs().size());
        viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
        *this = result;
        return *this;
      }
      else
      {
        viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), *this);
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+=(const vector_expression< const transposed_matrix_proxy<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.get_lhs().size1());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-=(const vector_expression< const transposed_matrix_proxy<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.get_lhs().size1());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+(const vector_expression< const transposed_matrix_proxy<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.get_lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-(const vector_expression< const transposed_matrix_proxy<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.get_lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.get_lhs(), proxy.get_rhs(), result);
      result = *this - result;
      return result;
    }

} //namespace viennacl


#endif
