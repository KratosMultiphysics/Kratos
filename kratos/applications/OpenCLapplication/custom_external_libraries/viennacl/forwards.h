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

/** @file forwards.h
    @brief This file provides the forward declarations for the main types used within ViennaCL
*/

/**
 @mainpage Source Code Documentation for ViennaCL 1.0.5

 This is the source code documentation of ViennaCL. Detailed information about the functions in ViennaCL can be found here.
 
 For a general overview over the types and functionality provided by ViennaCL, please refer to the file doc/viennacl.pdf

*/

#ifndef _VIENNACL_FORWARDS_H_
#define _VIENNACL_FORWARDS_H_

#include <stddef.h>

namespace viennacl
{
  namespace ocl
  {
    class kernel;
    class Device;
    Device & device();
  }
  
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

  /** @brief A tag class representing a matrix-vector product */
  struct op_prod;
  

  //forward declaration of basic types:
  template<class TYPE>
  class scalar;

  template <typename LHS, typename RHS, typename OP>
  class scalar_expression;

  template <typename SCALARTYPE>
  class vector_entry_proxy;
  
  template <typename LHS, typename RHS, typename OP>
  class vector_expression;

  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class vector_iterator;
    
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class vector;
  
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void copy(CPU_ITERATOR const & cpu_begin,
            CPU_ITERATOR const & cpu_end,
            vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_begin);
  
  
  struct row_major;    
  
  struct row_iteration;
  struct col_iteration;
  
  template <typename SCALARTYPE, unsigned int A>
  class outer_prod_proxy;
  
  template <typename SCALARTYPE, unsigned int A>
  class scalar_times_outer_prod_proxy;
  
  template <typename SCALARTYPE, typename F, unsigned int A>
  class transposed_matrix_proxy;
  
  template <class SCALARTYPE, typename F = row_major, unsigned int ALIGNMENT = 1>
  class matrix;

  template <typename SCALARTYPE, typename F, unsigned int A>
  class transposed_matrix_proxy;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 1>
  class compressed_matrix;
  
  template<class SCALARTYPE, unsigned int ALIGNMENT = 128>
  class coordinate_matrix;    
  
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
    //forward definition of norm_2_impl function
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void norm_2_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vcl_vec,
                     scalar<SCALARTYPE> & result,
                     size_t NUM_THREADS = 0);
    
    //forward definition of prod_impl functions
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::matrix<SCALARTYPE, F, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> &, 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> &, 
                                                    size_t NUM_THREADS = 0);

    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::compressed_matrix<SCALARTYPE, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::compressed_matrix<SCALARTYPE, ALIGNMENT> & , 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & , 
                                                    size_t NUM_THREADS = 0);

    template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    viennacl::vector_expression<const viennacl::coordinate_matrix<SCALARTYPE, ALIGNMENT>,
                                const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod > prod_impl(const viennacl::coordinate_matrix<SCALARTYPE, ALIGNMENT> & , 
                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> & , 
                                                    size_t NUM_THREADS = 0);

                                                    
    //forward definition of inner_prod_impl function
    template<class SCALARTYPE, unsigned int ALIGNMENT1, unsigned int ALIGNMENT2>
    viennacl::scalar_expression< const viennacl::vector<SCALARTYPE, ALIGNMENT1>, 
                                 const viennacl::vector<SCALARTYPE, ALIGNMENT2>,
                                 viennacl::op_inner_prod >
    inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT1> &,
                    const viennacl::vector<SCALARTYPE, ALIGNMENT2> &,
                    size_t NUM_THREADS = 0);
                    
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void inner_prod_impl(const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
                          const viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2,
                          scalar<SCALARTYPE> & result,
                          size_t NUM_THREADS = 0);
                    
      
    /** @brief A tag class representing a lower triangular matrix */
    struct lower_tag {};      //lower triangular matrix
    /** @brief A tag class representing an upper triangular matrix */
    struct upper_tag {};      //upper triangular matrix
    /** @brief A tag class representing a lower triangular matrix with unit diagonal*/
    struct unit_lower_tag {}; //unit lower triangular matrix
    /** @brief A tag class representing an upper triangular matrix with unit diagonal*/
    struct unit_upper_tag {}; //unit upper triangular matrix
    
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
