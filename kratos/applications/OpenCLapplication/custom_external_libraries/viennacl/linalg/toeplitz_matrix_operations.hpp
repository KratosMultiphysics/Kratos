#ifndef VIENNACL_LINALG_TOEPLITZ_MATRIX_OPERATIONS_HPP_
#define VIENNACL_LINALG_TOEPLITZ_MATRIX_OPERATIONS_HPP_

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

/** @file toeplitz_matrix_operations.hpp
    @brief Implementations of operations using toeplitz_matrix
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/fft.hpp"

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
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    vector_expression<const toeplitz_matrix<SCALARTYPE, ALIGNMENT>,
                      const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                      op_prod > prod_impl(const toeplitz_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                                     const vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec)
    {
      return vector_expression<const toeplitz_matrix<SCALARTYPE, ALIGNMENT>,
                               const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                               op_prod >(mat, vec);
    }
    
    // A * x
    /** @brief Returns a proxy class that represents matrix-vector multiplication with a toeplitz_matrix
    *
    * This is used for the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param NUM_THREADS Number of threads per work group. Can be used for fine-tuning.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::toeplitz_matrix<SCALARTYPE, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                viennacl::op_prod > prod_impl(const viennacl::toeplitz_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                                                              const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec, 
                                                              size_t NUM_THREADS)
    {
      return viennacl::vector_expression<const viennacl::toeplitz_matrix<SCALARTYPE, ALIGNMENT>,
                               const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                               viennacl::op_prod >(mat, vec);
    }
    
    /** @brief Carries out matrix-vector multiplication with a toeplitz_matrix
    *
    * Implementation of the convenience expression result = prod(mat, vec);
    *
    * @param mat    The matrix
    * @param vec    The vector
    * @param result The result vector
    */
      template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
      void prod_impl(const viennacl::toeplitz_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                     const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec,
                           viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & result)
      {
        assert(mat.size1() == result.size());
        assert(mat.size2() == vec.size());
        
        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tep(mat.elements().size() * 2);
        viennacl::detail::fft::real_to_complex(mat.elements(), tep, mat.elements().size());

        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp(vec.size() * 4);
        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp2(vec.size() * 4);

        tmp.clear();
        copy(vec, tmp);
        viennacl::detail::fft::real_to_complex(tmp, tmp2, vec.size() * 2);
        viennacl::linalg::convolve(tep, tmp2, tmp);
        viennacl::detail::fft::complex_to_real(tmp, tmp2, vec.size() * 2);
        copy(tmp2.begin(), tmp2.begin() + vec.size(), result.begin());
      }

  } //namespace linalg



    /** @brief Implementation of the operation v1 = A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const toeplitz_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+=(const vector_expression< const toeplitz_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-=(const vector_expression< const toeplitz_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                                 const vector<SCALARTYPE, ALIGNMENT>,
                                                                                 op_prod> & proxy) 
    {
      vector<SCALARTYPE, ALIGNMENT> result(proxy.get_lhs().size1());
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator+(const vector_expression< const toeplitz_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.get_lhs().size1() == size());
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
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator-(const vector_expression< const toeplitz_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                                                op_prod> & proxy) 
    {
      assert(proxy.get_lhs().size1() == size());
      vector<SCALARTYPE, ALIGNMENT> result(size());
      viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
      result = *this - result;
      return result;
    }

} //namespace viennacl


#endif
