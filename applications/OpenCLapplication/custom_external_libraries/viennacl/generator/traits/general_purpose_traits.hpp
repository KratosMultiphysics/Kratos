#ifndef VIENNACL_GENERATOR_TRAITS_GENERAL_PURPOSE_TRAITS_HPP
#define VIENNACL_GENERATOR_TRAITS_GENERAL_PURPOSE_TRAITS_HPP

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

/** @file viennacl/generator/traits/general_purpose_traits.hpp
 *  @brief Provides a set of metafunctions for the identification of types
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/operation_types.hpp"
#include "viennacl/generator/forwards.h"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/meta_tools/typelist.hpp"

namespace viennacl 
{
  namespace generator
  {
    template <class T>
    struct is_scalar_expression_impl 
    {
      enum { value = 0 };
    };

    template <class T>
    struct is_scalar_expression_impl<result_of::scalar_expression<T> > 
    {
      enum { value = 1};
    };

    template <class T>
    struct is_scalar_expression 
    {
      enum { value = is_scalar_expression_impl<typename result_of::expression_type<T>::Result >::value };
    };

    template <class T>
    struct is_temporary 
    {
      enum { value = 0 } ;
    };

    template <class LHS, class OP, class RHS>
    struct is_temporary<compound_node<LHS,OP,RHS,true> > 
    {
      enum {value = 1};
    };

    template <class REF>
    struct is_temporary<tmp_symbolic_vector<REF> > 
    {
      enum { value = 1};
    };

    template <class T>
    struct is_temporary_kernel_parameter 
    {
      enum { value = is_temporary<T>::value };
    };

    template <class T>
    struct is_temporary_kernel_parameter<inner_prod_impl_t<T> > 
    {
      enum { value = 1 };
    };


    template <class T>
    struct is_regular_kernel_parameter 
    {
      enum { value = 0 };
    };

    template <unsigned int ID,class SCALARTYPE, unsigned int ALIGNMENT>
    struct is_regular_kernel_parameter<symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > 
    {
      enum { value = 1 };
    };

    template <unsigned int ID,class SCALARTYPE, class F, unsigned int ALIGNMENT>
    struct is_regular_kernel_parameter<symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> > 
    {
      enum { value = 1 };
    };

    template <unsigned int ID, class SCALARTYPE>
    struct is_regular_kernel_parameter<cpu_symbolic_scalar<ID, SCALARTYPE> > 
    {
      enum { value = 1 };
    };

    template <unsigned int ID, class SCALARTYPE>
    struct is_regular_kernel_parameter<gpu_symbolic_scalar<ID, SCALARTYPE> > 
    {
      enum { value = 1 };
    };



    template <class T>
    struct is_pure_inner_product_leaf 
    {
      enum { value = 0};
    };

    template <class LHS,class RHS, bool is_temporary>
    struct is_pure_inner_product_leaf<compound_node<LHS,inner_prod_type,RHS, is_temporary> > 
    {
      enum { value = 1};
    };

    template <class T>
    struct is_inner_product_leaf 
    {
      enum { value = is_pure_inner_product_leaf<T>::value };
    };

    template <class LHS,class RHS>
    struct is_inner_product_leaf<compound_node<LHS,scal_mul_type,RHS> > 
    {
      enum { value = ( is_inner_product_leaf<LHS>::value  && is_scalar_expression<RHS>::value && !is_inner_product_leaf<RHS>::value )
                     || ( is_inner_product_leaf<RHS>::value && is_scalar_expression<LHS>::value &&!is_inner_product_leaf<LHS>::value )
           };
    };

    template <class LHS,class RHS>
    struct is_inner_product_leaf<compound_node<LHS,scal_div_type,RHS> >  
    {
      enum { value = ( is_inner_product_leaf<LHS>::value  && is_scalar_expression<RHS>::value && !is_inner_product_leaf<RHS>::value )
                     || ( is_inner_product_leaf<RHS>::value && is_scalar_expression<LHS>::value &&!is_inner_product_leaf<LHS>::value )
           };
    };

    template <class T>
    struct is_pure_product_leaf 
    {
      enum { value = 0};
    };

    template <class LHS,class RHS, bool is_temporary>
    struct is_pure_product_leaf<compound_node<LHS,prod_type,RHS, is_temporary> > 
    {
      enum { value = 1};
    };

    template <class T>
    struct is_product_leaf 
    {
      enum { value = is_pure_product_leaf<T>::value };
    };

    template <class LHS,class RHS>
    struct is_product_leaf<compound_node<LHS,scal_mul_type,RHS> > 
    {
      enum { value = is_product_leaf<LHS>::value
                     ||is_product_leaf<RHS>::value
           };
    };

    template <class T>
    struct is_null_type 
    {
      enum { value = 0 };
    };

    template <>
    struct is_null_type<NullType> 
    {
      enum { value = 1 };
    };

    template <class T>
    struct is_compound 
    {
      enum { value = 0 } ;
    };

    template <class LHS, class OP, class RHS, bool is_temporary>
    struct is_compound<compound_node<LHS,OP,RHS,is_temporary> > 
    {
      enum {value = 1};
    };

    template <class EXPR1, class EXPR2>
    struct is_same_expression_type_impl 
    {
      enum { value = 0 };
    };

    template <class EXPR1, class DESCRIPTOR1, class EXPR2, class DESCRIPTOR2>
    struct is_same_expression_type_impl<result_of::vector_expression<EXPR1,DESCRIPTOR1>,
                                        result_of::vector_expression<EXPR2,DESCRIPTOR2> > 
    {
      private:
        typedef result_of::vector_expression<EXPR1,DESCRIPTOR1> LHS;
        typedef result_of::vector_expression<EXPR2,DESCRIPTOR2> RHS;
      public:
        enum { value = LHS::Alignment == RHS::Alignment };
    };

    template <class EXPR1, class LHS_DESCRIPTOR1, class RHS_DESCRIPTOR1,
              class EXPR2, class LHS_DESCRIPTOR2, class RHS_DESCRIPTOR2>
    struct is_same_expression_type_impl<result_of::matrix_expression<EXPR1,LHS_DESCRIPTOR1,RHS_DESCRIPTOR1>,
                                        result_of::matrix_expression<EXPR2,LHS_DESCRIPTOR2,RHS_DESCRIPTOR2> > 
    {
      private:
        typedef result_of::matrix_expression<EXPR1,LHS_DESCRIPTOR1,RHS_DESCRIPTOR1> LHS;
        typedef result_of::matrix_expression<EXPR2,LHS_DESCRIPTOR2,RHS_DESCRIPTOR2> RHS;
        
      public:
        enum { value = LHS::Alignment == RHS::Alignment };
    };

    template <class EXPR1, class EXPR2>
    struct is_same_expression_type_impl<result_of::scalar_expression<EXPR1>,
                                        result_of::scalar_expression<EXPR2> > 
    {
      enum { value = 1 };
    };

    template<class EXPR1, class EXPR2>
    struct is_same_expression_type 
    {
      enum { value = is_same_expression_type_impl<typename result_of::expression_type<EXPR1>::Result,
                                                  typename result_of::expression_type<EXPR2>::Result>::value
           };
    };

  }
}

#endif


