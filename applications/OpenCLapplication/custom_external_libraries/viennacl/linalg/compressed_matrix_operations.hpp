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

#ifndef _VIENNACL_COMPRESSED_MATRIX_OPERATIONS_HPP_
#define _VIENNACL_COMPRESSED_MATRIX_OPERATIONS_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/kernels/compressed_matrix_kernels.h"

namespace viennacl
{
  namespace linalg
  {
    // A * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication with a compressed_matrix
    *
    * This is used for the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    vector_expression<const compressed_matrix<SCALARTYPE, ALIGNMENT>,
                      const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                      op_prod > prod_impl(const compressed_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                                     const vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                     size_t NUM_THREADS)
    {
      return vector_expression<const compressed_matrix<SCALARTYPE, ALIGNMENT>,
                               const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                               op_prod >(mat, vec);
    }
    
    /** @brief Carries out matrix-vector multiplication with a compressed_matrix
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class TYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(const viennacl::compressed_matrix<TYPE, ALIGNMENT> & mat, 
                   const viennacl::vector<TYPE, VECTOR_ALIGNMENT> & vec,
                         viennacl::vector<TYPE, VECTOR_ALIGNMENT> & result, 
                   size_t NUM_THREADS = 0)
    {
      assert(mat.size1() == result.size());
      assert(mat.size2() == vec.size());

      unsigned int pos = 0;
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, mat.handle1());
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, mat.handle2());
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, mat.handle());
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, result.handle());
      viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.setArgument(pos++, mat.size1());
      
      if (NUM_THREADS == 0)
        viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.start1D();
      else
        viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.start1D(viennacl::linalg::kernels::compressed_matrix<TYPE, ALIGNMENT>::vec_mul.work_groups() * NUM_THREADS, NUM_THREADS);
    }

    /** @brief Inplace solution of a lower triangular compressed_matrix with unit diagonal. Typically used for LU substitutions
    *
    * @param L    The matrix
    * @param vec    The vector
    */
    template<typename SCALARTYPE, unsigned int MAT_ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(compressed_matrix<SCALARTYPE, MAT_ALIGNMENT> const & L, vector<SCALARTYPE, VEC_ALIGNMENT> & vec, viennacl::linalg::unit_lower_tag)
    {
      unsigned int threads = viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.work_items_per_group();
      unsigned int pos = 0;
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setArgument(pos++, L.handle1());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setArgument(pos++, L.handle2());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setArgument(pos++, L.handle());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(int) * (threads+1)));
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE) * threads));
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.setArgument(pos++, L.size1());
      
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_forward.start1D(threads, threads);
    }
    
    /** @brief Convenience functions for result = solve(trans(mat), vec, unit_lower_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param L      The lower triangular sparse matrix
    * @param vec    The load vector, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, unsigned int MAT_ALIGNMENT, unsigned int VEC_ALIGNMENT, typename TAG>
    vector<SCALARTYPE, VEC_ALIGNMENT> solve(compressed_matrix<SCALARTYPE, MAT_ALIGNMENT> const & L,
                                        const vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                                        const viennacl::linalg::unit_lower_tag & tag)
    {
      // do an inplace solve on the result vector:
      vector<SCALARTYPE, VEC_ALIGNMENT> result(vec.size());
      result = vec;

      inplace_solve(L, result, tag);
    
      return result;
    }
    
    
    /** @brief Inplace solution of a upper triangular compressed_matrix. Typically used for LU substitutions
    *
    * @param U      The upper triangular matrix
    * @param vec    The vector
    */
    template<typename SCALARTYPE, unsigned int MAT_ALIGNMENT, unsigned int VEC_ALIGNMENT>
    void inplace_solve(compressed_matrix<SCALARTYPE, MAT_ALIGNMENT> const & U, vector<SCALARTYPE, VEC_ALIGNMENT> & vec, viennacl::linalg::upper_tag)
    {
      unsigned int threads = viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.work_items_per_group();
      
      unsigned int pos = 0;
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setArgument(pos++, U.handle1());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setArgument(pos++, U.handle2());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setArgument(pos++, U.handle());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(int) * (threads+2)));
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setLocalBuffer(pos++, static_cast<unsigned int>(sizeof(SCALARTYPE) * (threads+2)));
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setArgument(pos++, vec.handle());
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.setArgument(pos++, U.size1());
      
      viennacl::linalg::kernels::compressed_matrix<SCALARTYPE,MAT_ALIGNMENT>::lu_backward.start1D(threads, threads);
    }

    /** @brief Convenience functions for result = solve(trans(mat), vec, unit_lower_tag()); Creates a temporary result vector and forwards the request to inplace_solve()
    *
    * @param L      The lower triangular sparse matrix
    * @param vec    The load vector, where the solution is directly written to
    * @param tag    Dispatch tag
    */
    template<typename SCALARTYPE, unsigned int MAT_ALIGNMENT, unsigned int VEC_ALIGNMENT, typename TAG>
    vector<SCALARTYPE, VEC_ALIGNMENT> solve(compressed_matrix<SCALARTYPE, MAT_ALIGNMENT> const & L,
                                        const vector<SCALARTYPE, VEC_ALIGNMENT> & vec,
                                        viennacl::linalg::upper_tag const & tag)
    {
      // do an inplace solve on the result vector:
      vector<SCALARTYPE, VEC_ALIGNMENT> result(vec.size());
      result = vec;
    
      inplace_solve(L, result, tag);
    
      return result;
    }

    
  } //namespace linalg



    //v = A * x
    /** @brief Implementation of the operation v1 = A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                                          const viennacl::vector<SCALARTYPE, ALIGNMENT>,
                                                                                          viennacl::op_prod> & proxy) 
    {
      // check for the special case x = A * x
      if (proxy.get_rhs().handle().get() == this->handle().get())
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+=(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-=(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
