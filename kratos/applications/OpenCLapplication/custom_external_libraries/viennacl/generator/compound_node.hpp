#ifndef VIENNACL_GENERATOR_COMPOUND_NODE_HPP
#define VIENNACL_GENERATOR_COMPOUND_NODE_HPP

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

/** @file compound_node.hpp
 *  @brief Structures corresponding to binary nodes in the expression tree
 *
 *  Generator code contributed by Philippe Tillet
 */

#include <string>
#include <sstream>
#include <set>

#include "viennacl/generator/forwards.h"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/traits/result_of.hpp"

namespace viennacl 
{
  namespace generator
  {

    /**
    * @brief Binary node class for storing expression trees
    * 
    * @tparam LHS_ LHS of the expression
    * @tparam OP_ Operator of the expression
    * @tparam RHS_ RHS of the expression
    * @tparam is_temporary_ Boolean for storing whether the binary node is temporary.
    */
    template<class LHS_, class OP_, class RHS_, bool is_temporary_>
    class compound_node 
    {
      public:
        typedef LHS_  LHS;
        typedef RHS_  RHS;
        typedef OP_   OP;

        static const bool is_temporary = is_temporary_;

        static const std::string name() 
        {
            return LHS::name() + "_" + OP::name() + "_" + RHS::name();
        }
    };

    template<class LHS_, class RHS_, bool is_temporary_>
    class compound_node<LHS_,inner_prod_type,RHS_, is_temporary_> 
    {
      public:
        /**
        * @brief Specialization for the inner product
        */
        typedef LHS_ LHS;
        typedef RHS_ RHS;
        typedef inner_prod_type OP;
        typedef typename result_of::expression_type<RHS>::Result IntermediateType;  //Note: Visual Studio does not allow to combine this line with the next one directly.
        typedef typename IntermediateType::ScalarType ScalarType;

        static const bool is_temporary = is_temporary_;

        enum { id = -2 };

        static const std::string kernel_arguments() 
        {
          return  "__global float * " + name() + '\n';
        }

        static const std::string name() 
        {
          return  LHS::name() + "_inprod_" + RHS::name();
        }

        static const std::string scalar_name() 
        {
          return name() +"_s";
        };

    };

    /**
    * @brief Specialization for the matrix-vector product.
    */
    template<class LHS_, class RHS_, bool is_temporary_>
    class compound_node<LHS_,prod_type,RHS_, is_temporary_> 
    {
      private:
        typedef compound_node<LHS_,prod_type,RHS_, is_temporary_> self_type;

      public:
        typedef LHS_ LHS;
        typedef RHS_ RHS;

        typedef prod_type OP;
        enum { id = LHS::id };

        typedef typename result_of::expression_type<RHS>::Result IntermediateType;    //Note: Visual Studio does not allow to combine this line with the next one directly.
        typedef typename IntermediateType::ScalarType ScalarType;
        static const unsigned int Alignment = result_of::expression_type<RHS>::Result::Alignment;
        static const bool is_temporary = is_temporary_;

        static const std::string name() 
        {
          return LHS::name() + "_prod_" + RHS::name();
        }

        static const std::string size2_name() 
        {
          return "size_"+name();
        }

        static const std::string internal_size2_name() 
        {
          return "internal_size_"+name();
        }
        
        static const std::string name_argument() 
        {
          return " __global " + print_type<ScalarType*,Alignment>::value() + " " + name();
        }

        static const std::string kernel_arguments() 
        {
          return name_argument() + ", unsigned int " + size2_name() + ", unsigned int " + internal_size2_name() + "\n" ;
        }
    };


    /** @brief Addition operator on 2 elements of the same type */
    template<class LHS_TYPE, class RHS_TYPE>
    typename enable_if< is_same_expression_type<LHS_TYPE, RHS_TYPE>,
                        compound_node<LHS_TYPE, add_type, RHS_TYPE> >::type
    operator+ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ) 
    {
      return compound_node<LHS_TYPE, add_type, RHS_TYPE>();
    }

    /** @brief Substraction operator on 2 elements of the same type */
    template<class LHS_TYPE, class RHS_TYPE>
    typename enable_if< is_same_expression_type<LHS_TYPE, RHS_TYPE>,
                        compound_node<LHS_TYPE, sub_type, RHS_TYPE> >::type
    operator- ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ) 
    {
      return compound_node<LHS_TYPE, sub_type, RHS_TYPE>();
    }

    /** @brief Helper for the inner_prod operator */
    template<class LHS, class RHS>
    struct make_inner_prod;

    template<class LHS, class LHS_SIZE_DESCRIPTOR,
             class RHS, class RHS_SIZE_DESCRIPTOR>
    struct make_inner_prod<result_of::vector_expression<LHS, LHS_SIZE_DESCRIPTOR>,
                           result_of::vector_expression<RHS, RHS_SIZE_DESCRIPTOR> > 
    {
      typedef compound_node<LHS,inner_prod_type,RHS,true> Result;
    };


    /** @brief Inner product operator */
    template<class LHS, class RHS>
    compound_node<LHS,inner_prod_type,RHS,true> inner_prod ( LHS vec_expr1,RHS vec_expr2 ) 
    {
      typedef typename result_of::expression_type<LHS>::Result LHS_TYPE;
      typedef typename result_of::expression_type<RHS>::Result RHS_TYPE;
      typename make_inner_prod<LHS_TYPE,RHS_TYPE>::Result result;
      
      return result;;
    }

    /** @brief Product operator */
    template<class LHS, class RHS>
    compound_node<LHS,prod_type,RHS> prod ( LHS vec_expr1,RHS vec_expr2 ) 
    {
      return compound_node<LHS,prod_type,RHS>();
    }

  } // namespace generator
} // namespace viennacl

#endif

