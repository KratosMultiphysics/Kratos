#ifndef VIENNACL_GENERATOR_TREE_OPERATIONS_HPP
#define VIENNACL_GENERATOR_TREE_OPERATIONS_HPP

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

/** @file viennacl/generator/tree_operations.hpp
 *  @brief Functors for modifying the expression tree.
 * 
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/generator/elementwise_modifier.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"

namespace viennacl 
{
  namespace generator
  {
    namespace tree_utils 
    {

      /*
      * Count if
      */

      template <class T, template<class> class Pred>
      struct count_if 
      {
        enum { value = Pred<T>::value };
      };

      template <class T, std::string (*U)(), template<class> class Pred>
      struct count_if<elementwise_modifier_impl<T,U>, Pred>
      {
        enum { value = Pred<T>::value + count_if<T, Pred>::value };
      };

      template<class T, template<class> class Pred>
      struct count_if<inner_prod_impl_t<T>, Pred> 
      {
        enum { value = Pred<inner_prod_impl_t<T> >::value + count_if<T, Pred>::value };
      };


      template<class LHS, class RHS, class OP, bool is_temporary, template<class> class Pred>
      struct count_if<compound_node<LHS,OP,RHS,is_temporary>,Pred>
      {
        private:
          typedef compound_node<LHS,OP,RHS,is_temporary> T;
          
        public:
          enum { value = Pred<T>::value
                        +  count_if<LHS, Pred>::value
                        +  count_if<RHS, Pred>::value
                };
      };


      /*
      * Count if type
      */

      template<class T, class Searched>
      struct count_if_type 
      {
        enum { value = 0 };
      };

      template<class T>
      struct count_if_type<T,T> 
      {
        enum { value = 1 };
      };

      template<class T, std::string (*U)(), class Searched>
      struct count_if_type<elementwise_modifier_impl<T,U>, Searched> 
      {
        enum { value = count_if_type<T, Searched>::value };
      };

      template <class T, std::string (*U)()>
      struct count_if_type<elementwise_modifier_impl<T,U>, elementwise_modifier_impl<T,U> > 
      {
        enum { value = 1 + count_if_type<T, elementwise_modifier_impl<T,U> >::value };
      };

      template <class LHS, class OP, class RHS, bool is_temporary>
      struct count_if_type<compound_node<LHS, OP, RHS, is_temporary>,
                           compound_node<LHS, OP, RHS, is_temporary> > 
      {
        private:
          typedef compound_node<LHS, OP, RHS, is_temporary> T;
        public:
          enum { value = 1 +  count_if_type<LHS, T>::value
                           +  count_if_type<RHS, T>::value
               };
      };

      template <class LHS, class OP, class RHS, bool is_temporary, class Searched>
      struct count_if_type< compound_node<LHS,OP,RHS,is_temporary>, Searched> 
      {
        enum { value = count_if_type<LHS, Searched>::value
                      +  count_if_type<RHS, Searched>::value
             };
      };


      /*
      * Expand
      */

      template <class LHS, class OP, bool is_temporary, class RHS_LHS, class RHS_OP, class RHS_RHS, bool RHS_is_temporary>
      struct expand_right 
      {
          typedef compound_node< compound_node<LHS, OP, RHS_LHS, RHS_is_temporary>,
                                 RHS_OP,
                                 compound_node<LHS, OP, RHS_RHS, RHS_is_temporary>,
                                 is_temporary>   Result;
      };

      template <class LHS_LHS, class LHS_OP, class LHS_RHS, bool LHS_is_temporary, class OP, class RHS, bool is_temporary>
      struct expand_left 
      {
          typedef compound_node< compound_node<LHS_LHS, OP, RHS, LHS_is_temporary>,
                                 LHS_OP,
                                 compound_node<LHS_RHS, OP, RHS, LHS_is_temporary>,
                                 is_temporary>        Result;
      };

      template <class T>
      struct expand 
      {
        typedef T Result;
      };

      template <class T, std::string (*U)()>
      struct expand< elementwise_modifier_impl<T,U> > 
      {
        private:
          typedef typename expand<T>::Result                 SUB_Result;
        public:
          typedef elementwise_modifier_impl<SUB_Result,U>    Result;
      };

      template<class T>
      struct expand<inner_prod_impl_t<T> > 
      {
        private:
          typedef typename expand<T>::Result      SUB_Result;
        public:
          typedef inner_prod_impl_t<SUB_Result>   Result;
      };


      template<class LHS,class OP,class RHS,bool is_temporary>
      struct expand< compound_node<LHS,OP,RHS,is_temporary> >
      {
        typedef compound_node<typename expand<LHS>::Result, OP, typename expand<RHS>::Result, is_temporary>   Result;
      };

      #define make_right_expandable(__OPERATOR1__ , __OPERATOR2__) \
                      template<class LHS, class RHS_LHS, class RHS_RHS, bool RHS_is_temporary, bool is_temporary>\
                      struct expand< compound_node<LHS, __OPERATOR1__, compound_node<RHS_LHS, __OPERATOR2__, RHS_RHS, RHS_is_temporary>, is_temporary> >\
                      {\
                        typedef typename expand_right<typename expand<LHS>::Result\
                                                    , __OPERATOR1__\
                                                    , is_temporary\
                                                    , typename expand<RHS_LHS>::Result\
                                                    , __OPERATOR2__\
                                                    , typename expand<RHS_RHS>::Result\
                                                    , RHS_is_temporary>::Result Result;\
                      }

      #define make_left_expandable(__OPERATOR1__ , __OPERATOR2__) \
                      template<class LHS_LHS, class LHS_RHS, bool LHS_is_temporary, class RHS, bool is_temporary>\
                      struct expand< compound_node< compound_node<LHS_LHS, __OPERATOR2__ , LHS_RHS , LHS_is_temporary>\
                                                    , __OPERATOR1__\
                                                    , RHS\
                                                    , is_temporary> >\
                      {\
                        typedef typename expand_left< typename expand<LHS_LHS>::Result\
                                                    , __OPERATOR2__\
                                                    , typename expand<LHS_RHS>::Result\
                                                    , LHS_is_temporary\
                                                    , __OPERATOR1__\
                                                    , typename expand<RHS>::Result\
                                                    , is_temporary\
                                                      >	::Result Result;\
                      }

      make_right_expandable ( scal_mul_type,add_type );
      make_right_expandable ( scal_mul_type,sub_type );
      make_left_expandable ( scal_mul_type,add_type );
      make_left_expandable ( scal_mul_type,sub_type );


      #undef make_left_expandable
      #undef make_right_expandable

      ////////////////////////////
      // REGISTER TEMPORARIES  //
      ///////////////////////////

      template <class T>
      struct make_temporary;

      template <unsigned int ID, class SCALARTYPE, unsigned int ALIGNMENT>
      struct make_temporary<symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > 
      {
        typedef tmp_symbolic_vector< symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > Result;
      };

      template <unsigned int ID,typename SCALARTYPE, class F, unsigned int ALIGNMENT>
      struct make_temporary<symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> > {
        typedef tmp_symbolic_matrix< symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> > Result;
      };

      template <class T, bool only_first_order, class Assigned = void, bool is_nested = false>
      struct register_temporaries 
      {
        typedef T Result;
      };

      template <class T, bool only_first_order>
      struct register_temporaries<T, only_first_order, T, true> 
      {
        typedef typename make_temporary<T>::Result     Result;
      };

      template <class T, std::string (*U)(), bool only_first_order, class Assigned, bool is_nested>
      struct register_temporaries<elementwise_modifier_impl<T,U>, only_first_order, Assigned, is_nested> 
      {
        private:
          typedef typename register_temporaries<T, only_first_order, Assigned, is_nested>::Result   SUB_Result;
        public:
          typedef elementwise_modifier_impl<SUB_Result,U>     Result;
      };


      template <class LHS, class OP, class RHS, bool is_temporary, bool only_first_order, class Assigned, bool is_nested>
      struct register_temporaries<compound_node<LHS,OP,RHS,is_temporary>, only_first_order, Assigned, is_nested> 
      {
        private:
          typedef compound_node<LHS,OP,RHS,is_temporary> T;
          static const bool is_non_trivial =  is_pure_product_leaf<T>::value ||is_pure_inner_product_leaf<T>::value;
          typedef typename register_temporaries<LHS, only_first_order, Assigned, is_nested || is_non_trivial>::Result LHS_Result;
          typedef typename register_temporaries<RHS, only_first_order, Assigned, is_nested || is_non_trivial>::Result RHS_Result;

          typedef compound_node<LHS_Result,OP,RHS_Result, is_non_trivial&& ( is_temporary || is_nested )  > RecursiveResult;
          typedef compound_node<LHS,OP,RHS,true> EarlyStoppingResult;
        public:
          typedef typename get_type_if<EarlyStoppingResult, RecursiveResult,is_non_trivial && only_first_order && is_nested>::Result Result;
      };


      ////////////////////////////////
      //////// EXTRACTIF ////////
      ///////////////////////////////


      template <class T, 
                template<class> class Pred,
                template<class, class> class Comp = typelist_utils::true_comp,
                class TList = NullType>
      struct extract_if 
      {
        private:
          typedef typelist<T,TList>    TypeTrue;
          typedef NullType             TypeFalse;
        public:
          typedef typename get_type_if<TypeTrue, TypeFalse, Pred<T>::value>::Result      Result;
      };

      template <class T,
                std::string (*U)(),
                template<class> class Pred,
                template<class,class> class Comp,
                class TList>
      struct extract_if<elementwise_modifier_impl<T,U>, Pred, Comp, TList> 
      {
        private:
          typedef typename extract_if<T, Pred, Comp, TList>::Result         SUB_Result;
        public:
          typedef typename typelist_utils::fuse<TList,SUB_Result>::Result   Result;
      };

      template <class T,
                template<class> class Pred,
                template<class,class> class Comp,
                class TList>
      struct extract_if<inner_prod_impl_t<T>, Pred, Comp, TList > 
      {
        private:
          typedef typename T::LHS LHS;
          typedef typename T::RHS RHS;
          typedef typename extract_if<LHS, Pred, Comp, TList>::Result            LHS_Result;
          typedef typename extract_if<RHS,  Pred, Comp, TList>::Result           RHS_Result;
          typedef typename typelist_utils::fuse<TList, LHS_Result>::Result       TmpResult1;
          typedef typename typelist_utils::fuse<TmpResult1, RHS_Result>::Result  TmpResult2;
          
          typedef TmpResult2                                                                        TypeFalse;
          typedef typename typelist_utils::append<TmpResult2, inner_prod_impl_t<T> >::Result        TypeTrue;
          
        public:
          typedef typename get_type_if<TypeTrue, TypeFalse, Pred< inner_prod_impl_t<T> >::value>::Result   Result;
      };

      template <class LHS, class OP, class RHS, bool is_temporary,
                template<class> class Pred,
                template<class,class> class Comp,
                class TList>
      struct extract_if< compound_node<LHS, OP, RHS, is_temporary>, Pred, Comp, TList>
      {
        private:
          typedef compound_node<LHS,OP,RHS,is_temporary> T;
          typedef typename extract_if<LHS,Pred,Comp,TList>::Result LHS_Result;
          typedef typename extract_if<RHS,Pred,Comp,TList>::Result RHS_Result;

          typedef typename typelist_utils::fuse< typename typelist_utils::fuse<TList, LHS_Result, Comp>::Result,
                                                 RHS_Result, 
                                                 Comp >::Result     TypeFalse;
          typedef typelist<T, TList>                                TypeTrue;
        public:
          typedef typename get_type_if<TypeTrue, TypeFalse, Pred<T>::value>::Result       Result;
      };


      ///////////////////////////////
      //////// FLIP_TREE  ///////////
      ///////////////////////////////

      template <class OP, bool flip>
      struct invert_flip 
      {
        enum { value = flip };
      };

      template <bool flip>
      struct invert_flip<sub_type, flip> 
      {
        enum { value = !flip };
      };

      template <class OP, bool flip>
      struct flip_operator 
      {
          typedef OP Result;
      };

      template <>
      struct flip_operator<sub_type, true> 
      {
          typedef add_type Result;
      };

      template <>
      struct flip_operator<add_type, true> 
      {
          typedef sub_type Result;
      };

      template <class T, bool flip = false>
      struct flip_tree 
      {
          typedef T Result;
      };

      template <class T, std::string (*U)(),  bool flip>
      struct flip_tree <elementwise_modifier_impl<T,U>, flip> 
      {
        private:
          typedef typename flip_tree<T, flip>::Result       SUB_Result;
        public:
          typedef elementwise_modifier_impl<SUB_Result,U>   Result;
      };

      template <class LHS, class OP, class RHS, bool is_temporary, bool flip>
      struct flip_tree< compound_node<LHS, OP, RHS, is_temporary>, flip>
      {
        private:
          typedef typename flip_tree<LHS,flip>::Result LHS_Result;
          typedef typename flip_tree<RHS, invert_flip<OP, flip>::value >::Result RHS_Result;

        public:
          typedef compound_node<LHS_Result, typename flip_operator<OP, flip>::Result , RHS_Result, is_temporary> Result;
      };

      ////////////////////////////////
      //////// REMOVE_IF ////////////
      ///////////////////////////////

      template <class OP, class RHS>
      struct handle_unary_minus
      {
        typedef RHS Result;
      };

      template <class RHS>
      struct handle_unary_minus<sub_type, RHS> 
      {
        typedef compound_node<NullType,sub_type,RHS> Result;
      };

      template <class T>
      struct compound_to_simple 
      {
        typedef T Result;
      };

      template <class LHS, class OP>
      struct compound_to_simple<compound_node<LHS, OP, NullType> > 
      {
        typedef LHS Result;
      };

      template <class OP, class RHS>
      struct compound_to_simple<compound_node<NullType, OP, RHS> > 
      {
        typedef typename handle_unary_minus<OP,RHS>::Result Result;
      };

      template <class OP, class RHS, class Enable=void>
      struct get_new_operator 
      {
        typedef OP Result;
      };

      template <class RHS_OP, class RHS_RHS>
      struct get_new_operator <sub_type, compound_node<NullType, RHS_OP, RHS_RHS> >
      {
        typedef RHS_OP Result;
      };

      template <class T, template<class> class Pred>
      struct remove_if 
      {
        typedef typename get_type_if<NullType,T,Pred<T>::value>::Result    Result;
        typedef typename get_type_if<NullType,T,Pred<T>::value>::Result    TmpTree;
      };

      template <class T, std::string (*U)(), template<class> class Pred>
      struct remove_if<elementwise_modifier_impl<T,U>,Pred > 
      {
        typedef elementwise_modifier_impl<typename remove_if<T,Pred>::Result, U> Result;
      };

      template <class LHS, class OP, class RHS, bool is_temporary, template<class> class Pred>
      struct remove_if<compound_node<LHS,OP,RHS,is_temporary>, Pred> 
      {
        private:
          typedef compound_node<LHS,OP,RHS,is_temporary> T;

          typedef typename remove_if<LHS,Pred>::TmpTree LHS_TmpTree;
          typedef typename remove_if<RHS,Pred>::TmpTree RHS_TmpTree;

          typedef typename compound_to_simple<typename remove_if<LHS,Pred>::Result>::Result LHS_Result;
          typedef typename compound_to_simple<typename remove_if<RHS,Pred>::Result>::Result RHS_Result;

          typedef compound_node<LHS_TmpTree,OP,RHS_TmpTree> TmpTree0;
          typedef typename compound_to_simple<compound_node<LHS_Result,
                                                            typename get_new_operator<OP,RHS_TmpTree>::Result,
                                                            RHS_Result,
                                                            is_temporary> >::Result    Result0;
        public:
          typedef typename get_type_if<NullType, TmpTree0,  Pred<T>::value>::Result    TmpTree;
          typedef typename get_type_if<NullType, Result0,   Pred<T>::value>::Result    Result;
      };
      
    }  // namespace tree_utils
  } // namespace generator
} // namespace viennacl
#endif
