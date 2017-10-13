#ifndef VIENNACL_GENERATOR_MAKE_CODE_EXPRESSION_HPP
#define VIENNACL_GENERATOR_MAKE_CODE_EXPRESSION_HPP

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

/** @file viennacl/generator/make_code/expression.hpp
 *   @brief Directives for generating code for simple expressions.
 *
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/elementwise_modifier.hpp"
#include "viennacl/generator/symbolic_types/symbolic_scalars.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"

namespace viennacl 
{
  namespace generator
  {

    template <class T>
    struct make_expression_code 
    {
      static const std::string value(std::string const & loop_accessor) 
      {
        return  T::name() + '[' + loop_accessor  + ']';
      }
    };

    template <unsigned int ID, class SCALARTYPE>
    struct make_expression_code<cpu_symbolic_scalar<ID,SCALARTYPE> > 
    {
      static const std::string value(std::string const & loop_accessor) 
      {
        return  cpu_symbolic_scalar<ID,SCALARTYPE>::name();
      }
    };

    template <unsigned int ID,class SCALARTYPE>
    struct make_expression_code<gpu_symbolic_scalar<ID,SCALARTYPE> > 
    {
      static const std::string value(std::string const & loop_accessor) 
      {
        return  '*' + gpu_symbolic_scalar<ID,SCALARTYPE>::name();
      }
    };

    template <class LHS, class RHS, bool is_temporary >
    struct make_expression_code<compound_node<LHS,inner_prod_type,RHS, is_temporary> > 
    {
      private:
        typedef compound_node<LHS,inner_prod_type,RHS, is_temporary> T;
        
      public:
        static const std::string value(std::string const & loop_accessor) 
        {
          return T::name() +"_sum";
        }
    };

    template< >
    struct make_expression_code< NullType > 
    {
      static const std::string value(std::string const & loop_accessor) 
      {
          return "0";
      }
    };

    template<class T, std::string (*U)()>
    struct make_expression_code< elementwise_modifier_impl<T, U> > 
    {
      typedef elementwise_modifier_impl<T, U> EW_M;
      static const std::string value ( std::string const & loop_accessor ) 
      {
        return EW_M::modify(make_expression_code<T>::value(loop_accessor));
      }
    };

    template<class LHS, class OP, class RHS >
    struct make_expression_code<compound_node<LHS, OP, RHS, false> > 
    {
      static const std::string value(std::string const & loop_accessor = "k") 
      {
        return make_expression_code<LHS>::value(loop_accessor)
               + OP::expression_string() 
               + make_expression_code<RHS>::value(loop_accessor);
      }
    };

    template<class LHS, class RHS, unsigned int Alignment>
    struct dot_product_impl
    {
      static const std::string value(std::string lhs_loop_id,
                                     std::string rhs_loop_id)
      {
        return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + "," + make_expression_code<RHS>::value(rhs_loop_id) + ")";
      }
    };

    template<class LHS, class RHS>
    struct dot_product_impl<LHS, RHS, 8>
    {
      static const std::string value(std::string lhs_loop_id,
                                     std::string rhs_loop_id)
      {
        return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + ".s0123" + ","
                      + make_expression_code<RHS>::value(rhs_loop_id) + ".s0123 )" 
         +  " + dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s4567" + ","
                      + make_expression_code<RHS>::value(rhs_loop_id) + ".s4567 );"
        ;
      }
    };

    template<class LHS, class RHS>
    struct dot_product_impl<LHS, RHS, 16>
    {
      static const std::string value(std::string lhs_loop_id,std::string rhs_loop_id)
      {
        return "dot(" + make_expression_code<LHS>::value(lhs_loop_id) + ".s0123" + ","
                      + make_expression_code<RHS>::value(rhs_loop_id) + ".s0123)" 
        +"\n	+ dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s4567" + "," 
                      + make_expression_code<RHS>::value(rhs_loop_id) + ".s4567) "
        +"\n	+ dot("	+ make_expression_code<LHS>::value(lhs_loop_id) + ".s89ab" + "," 
                      + make_expression_code<RHS>::value ( rhs_loop_id ) + ".s89ab) "
        +"\n	+ dot("	+ make_expression_code<LHS>::value ( lhs_loop_id ) + ".scdef" + "," 
                      + make_expression_code<RHS>::value ( rhs_loop_id ) + ".scdef)" 
        ;
      }
    };

    template<class LHS, class RHS>
    struct dot_product
    {
      static const std::string value(std::string lhs_loop_id,std::string rhs_loop_id)
      {
        return dot_product_impl<LHS,RHS,LHS::Alignment>::value(lhs_loop_id,rhs_loop_id);
      }
    };

  }

}

#endif


