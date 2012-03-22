#ifndef VIENNACL_FORWARDS_H
#define VIENNACL_FORWARDS_H

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


/** @file viennacl/forwards.h
    @brief This file provides the forward declarations for the main types used within ViennaCL
*/

/**
 @mainpage Source Code Documentation for ViennaCL 1.2.1

 This is the source code documentation of ViennaCL. Detailed information about the functions in ViennaCL can be found here.
 
 For a general overview over the types and functionality provided by ViennaCL, please refer to the file doc/viennacl.pdf

*/

#include <cstddef>
#include "viennacl/ocl/forwards.h"
#include "viennacl/meta/enable_if.hpp"

namespace viennacl
{
  typedef std::size_t                                       vcl_size_t;
  typedef std::ptrdiff_t                                    vcl_ptrdiff_t;
 
  
  /** @brief A tag class representing addition */
  struct op_add;
  /** @brief A tag class representing subtraction */
  struct op_sub;
  /** @brief A tag class representing division */
  struct op_div;
  
  /** @brief A tag class representing inner products of two vectors */
  struct op_inner_prod;

  /** @brief A tag class representing the 1-norm of a vector */
  struct op_norm_1;

  /** @brief A tag class representing the 2-norm of a vector */
  struct op_norm_2;

  /** @brief A tag class representing the 2-norm of a vector */
  struct op_norm_inf;

  /** @brief A tag class representing a matrix-vector product */
  struct op_prod;
  
  /** @brief A tag class representing transposed matrices */
  struct op_trans;
  

  //forward declaration of basic types:
  template<class TYPE>
  class scalar;

  template <typename LHS, typename RHS, typename OP>
  class scalar_expression;

  template <typename SCALARTYPE>
  class entry_proxy;
  
  template <typename LHS, typename RHS, typename OP>
  class vector_expression;

  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class vector_iterator;

  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class const_vector_iterator;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class vector;
  
  //the following forwards are needed for GMRES
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void copy(CPU_ITERATOR const & cpu_begin,
            CPU_ITERATOR const & cpu_end,
            vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin);

  template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
  void copy(const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_begin,
            const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_end,
            vector_iterator<SCALARTYPE, ALIGNMENT_DEST> gpu_dest_begin);
  
  template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
  void copy(const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_begin,
            const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_end,
            const_vector_iterator<SCALARTYPE, ALIGNMENT_DEST> gpu_dest_begin);
  
  
  struct row_major;    
  struct column_major;    
  
  struct row_iteration;
  struct col_iteration;

  template <typename LHS, typename RHS, typename OP>
  class matrix_expression;

  //
  // Matrix types:
  //  
  template <class SCALARTYPE, typename F = row_major, unsigned int ALIGNMENT = 1>
  class matrix;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class compressed_matrix;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 128>
  class coordinate_matrix;    

  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class circulant_matrix;
    
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class hankel_matrix;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class toeplitz_matrix;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class vandermonde_matrix;
  
  //
  // Proxies:
  //
  template <typename SizeType = std::size_t, typename DistanceType = std::ptrdiff_t>
  class basic_range;
  
  typedef basic_range<>  range;
  
  template <typename MatrixType>
  class matrix_range;

  template <typename VectorType>
  class vector_range;
  
  template <typename T>
  struct is_scalar;

  template <typename T>
  struct is_vector;

  template <typename T>
  struct is_matrix;
  
  namespace tools
  {
    //helper for matrix row/col iterators 
    //must be specialized for every viennacl matrix type
    template <typename ROWCOL, typename MATRIXTYPE>
    struct MATRIX_ITERATOR_INCREMENTER
    {
      static void apply(const MATRIXTYPE & mat, unsigned int & row, unsigned int & col)
      {
          typedef typename MATRIXTYPE::ERROR_SPECIALIZATION_FOR_THIS_MATRIX_TYPE_MISSING          ErrorIndicator;
      }
    };
  }
    
  namespace linalg
  {
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void convolve_i(viennacl::vector<SCALARTYPE, ALIGNMENT>& input1,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& input2,
                    viennacl::vector<SCALARTYPE, ALIGNMENT>& output);
    
#ifndef _MSC_VER
    //forward definition of norm_1_impl function
    template <typename V1, typename S2>
    void norm_1_impl(V1 const & vec,
                     S2 & result,
                      typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                    && viennacl::is_scalar<S2>::value
                                                  >::type * dummy = 0);

    //forward definition of norm_2_impl function
    template <typename V1, typename S2>
    void norm_2_impl(V1 const & vec,
                     S2 & result,
                     typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                   && viennacl::is_scalar<S2>::value
                                                 >::type * dummy = 0);

    //forward definition of norm_inf_impl function
    template <typename V1, typename S2>
    void norm_inf_impl(V1 const & vec,
                       S2 & result,
                       typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                     && viennacl::is_scalar<S2>::value
                                                   >::type * dummy = 0);
#endif    
    
    //forward definition of prod_impl functions
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> &, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> &);

    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::compressed_matrix<SCALARTYPE, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::compressed_matrix<SCALARTYPE, ALIGNMENT> & , 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> &);

    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::coordinate_matrix<SCALARTYPE, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::coordinate_matrix<SCALARTYPE, ALIGNMENT> & , 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> &);

                                                    
    //forward definition of inner_prod_impl function
    /*template <typename V1, typename V2>
    typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                  && viennacl::is_vector<V2>::value,
                                  viennacl::scalar_expression< const V1, 
                                                               const V2,
                                                               viennacl::op_inner_prod >
                                >::type
    inner_prod_impl(V1 const & vec1,
                    V2 const & vec2);*/
    
#ifndef _MSC_VER
    template <typename V1, typename V2, typename S3>
    void inner_prod_impl(V1 const & vec1,
                         V2 const & vec2,
                         S3 & result,
                         typename viennacl::enable_if< viennacl::is_vector<V1>::value
                                                       && viennacl::is_vector<V2>::value
                                                       && viennacl::is_scalar<S3>::value
                                                     >::type * dummy = 0);
#endif                                                   
                    
      
    /** @brief A tag class representing a lower triangular matrix */
    struct lower_tag 
    {
      static const char * const name() { return "lower"; }
    };      //lower triangular matrix
    /** @brief A tag class representing an upper triangular matrix */
    struct upper_tag 
    {
      static const char * const name() { return "upper"; }
    };      //upper triangular matrix
    /** @brief A tag class representing a lower triangular matrix with unit diagonal*/
    struct unit_lower_tag
    {
      static const char * const name() { return "unit_lower"; }
    }; //unit lower triangular matrix
    /** @brief A tag class representing an upper triangular matrix with unit diagonal*/
    struct unit_upper_tag
    {
      static const char * const name() { return "unit_upper"; }
    }; //unit upper triangular matrix
    
    //preconditioner tags
    class ilut_tag;
    
    /** @brief A tag class representing the use of no preconditioner */
    class no_precond
    {
      public:
        template <typename VectorType>
        void apply(VectorType & vec) const {}
    };
    
    
  } //namespace linalg
} //namespace viennacl

#endif

/*@}*/
