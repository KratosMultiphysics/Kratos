#ifndef VIENNACL_GENERATOR_TOKENS_MANAGEMENT_HPP
#define VIENNACL_GENERATOR_TOKENS_MANAGEMENT_HPP

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

/** @file viennacl/generator/tokens_management.hpp
 *  @brief Creation and management of the tokens list
 * 
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/operation_types.hpp"
#include "viennacl/generator/tree_operations.hpp"
#include "viennacl/generator/meta_tools/typelist.hpp"

namespace viennacl 
{
  namespace generator
  {

    /////////////////////////////
    ///////// TOKENS ///////////
    ////////////////////////////

    template<class T, bool is_in_temporary_kernel, class TList = NullType, class TokenOp = add_type,  class Enable = void>
    struct extract_tokens 
    {
      typedef TList Result;
    };

    template <class LHS, class OP, class RHS, bool is_temporary, bool is_in_temporary_kernel, class TList, class TokenOp>
    struct extract_tokens<compound_node<LHS, OP, RHS, is_temporary>,
                          is_in_temporary_kernel,
                          TList,
                          TokenOp> 
    {
      private:
        typedef compound_node<LHS,OP,RHS,is_temporary> T;
        typedef typename extract_tokens<LHS, is_in_temporary_kernel, TList, TokenOp>::Result LHS_Result;
        typedef typename extract_tokens<RHS, is_in_temporary_kernel, TList, OP>::Result      RHS_Result;
        typedef typename typelist_utils::fuse<RHS_Result,LHS_Result>::Result                 ResultFalse;
        typedef typename typelist_utils::append<TList, std::pair<T,TokenOp> >::Result        ResulTrue;
        
      public:
        typedef typename get_type_if<ResulTrue,ResultFalse,is_product_leaf<T>::value>::Result Result;
    };

    template <class TList,bool make_operator_inplace>
    struct tokenize_operators 
    {
      private:
        typedef typename TList::Head Head;
        typedef typename get_type_if<typename make_inplace<typename Head::second_type>::Result,
                                     assign_type,
                                     make_operator_inplace>::Result   NewOperator;
        typedef std::pair<typename Head::first_type, NewOperator>     NewHead;
        typedef typename TList::Tail                                  Tail;
        typedef typename tokenize_operators<Tail, true>::Result       NewTail;
        
      public:
        typedef typelist<NewHead, NewTail> Result;
    };

    template <bool make_operator_inplace>
    struct tokenize_operators<NullType,make_operator_inplace> 
    {
      typedef NullType Result;
    };



    template <class T, bool is_in_temporary_kernel>
    struct generate_tokens;

    template <class LHS,class OP, class RHS, bool is_temporary, bool is_in_temporary_kernel>
    struct generate_tokens<compound_node<LHS,OP,RHS,is_temporary>, is_in_temporary_kernel> 
    {
      private:
        typedef typename tree_utils::remove_if<RHS, is_product_leaf>::Result                NewTree;
        typedef std::pair<NewTree,OP>                                                       LinearToken;
        typedef typename extract_tokens<RHS, is_in_temporary_kernel>::Result                Products;
        typedef typename tokenize_operators<Products,
                                            !is_null_type<NewTree>::value>::Result         TokenizedProducts;

      public:
        typedef typelist<LinearToken,TokenizedProducts> Result;
    };

  }
}
#endif
