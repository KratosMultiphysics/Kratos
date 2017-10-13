#ifndef VIENNACL_GENERATOR_MAKE_CODE_MATRIX_VECTOR_PRODUCT_HPP
#define VIENNACL_GENERATOR_MAKE_CODE_MATRIX_VECTOR_PRODUCT_HPP

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

/** @file viennacl/generator/make_code/matrix-vector_product.hpp
 *   @brief Directives for generating code for the matrix-vector product
 *
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/generator/make_code/expression.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/meta_tools/typelist.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/tree_operations.hpp"

namespace viennacl 
{
  namespace generator
  {

    template<class T, class OP, class ASSIGNED>
    struct make_product_code;

    template<class T, class SIZE_DESCRIPTOR, class OP, class ASSIGNED>
    struct make_product_code<result_of::vector_expression<T,SIZE_DESCRIPTOR>, OP, ASSIGNED> 
    {
      private:
        typedef typename tree_utils::remove_if<T, is_pure_product_leaf>::Result                           SCALAR_EXPR;
        typedef typename tree_utils::extract_if<T,is_pure_product_leaf>::Result::Head                     ARG;
        typedef typename generate_tokens<compound_node<NullType,assign_type,SCALAR_EXPR>, false>::Result  Tokens;
        typedef typename ARG::LHS                                   LHS;
        typedef typename ARG::RHS                                   RHS;
        typedef typename result_of::expression_type<LHS>::Result    MatExpr;
        typedef typename MatExpr::ScalarType                        ScalarType;
        typedef typename MatExpr::Layout                            Layout;

        static const unsigned int Alignment = result_of::expression_type<LHS>::Result::Alignment;
        
        static const std::string assign_res(Int2Type<true>) 
        {
          return ASSIGNED::name() + "[ row ]" + OP::expression_string() +  "dot_prod ;";
        }

        static const std::string assign_res(Int2Type<false>)
        {
          return ASSIGNED::name() + "[ row ]" + OP::expression_string() + make_expression_code<SCALAR_EXPR>::value ( "k" ) + "* dot_prod ;";
        }

        static const std::string expression_string() 
        {
          return make_expression_code<LHS>::value() + "*" +  make_expression_code<RHS>::value();
        }

        static const std::string fill_ith_row(viennacl::row_major)
        {
          std::string internal_size_2_expression = MatExpr::internal_size2_expression();
          if(Alignment==1)
            return  " dot_prod +=  " + dot_product<LHS,RHS>::value("row *" + internal_size_2_expression + " + col","col") + ";\n";
          else if (Alignment == 16)
            return " unsigned int scaled_row = row * " + to_string(Alignment) + ";\n"
                  +  "dot_prod.s0 +=  " + dot_product<LHS,RHS>::value("scaled_row *" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s1 +=  " + dot_product<LHS,RHS>::value("(scaled_row+1)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s2 +=  " + dot_product<LHS,RHS>::value("(scaled_row+2)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s3 +=  " + dot_product<LHS,RHS>::value("(scaled_row+3)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s4 +=  " + dot_product<LHS,RHS>::value("(scaled_row+4)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s5 +=  " + dot_product<LHS,RHS>::value("(scaled_row+5)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s6 +=  " + dot_product<LHS,RHS>::value("(scaled_row+6)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s7 +=  " + dot_product<LHS,RHS>::value("(scaled_row+7)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s8 +=  " + dot_product<LHS,RHS>::value("(scaled_row+8)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.s9 +=  " + dot_product<LHS,RHS>::value("(scaled_row+9)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.sa +=  " + dot_product<LHS,RHS>::value("(scaled_row+10)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.sb +=  " + dot_product<LHS,RHS>::value("(scaled_row+11)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.sc +=  " + dot_product<LHS,RHS>::value("(scaled_row+12)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.sd +=  " + dot_product<LHS,RHS>::value("(scaled_row+13)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.se +=  " + dot_product<LHS,RHS>::value("(scaled_row+14)*" + internal_size_2_expression + " + col","col") + ";\n"
                  + "  dot_prod.sf +=  " + dot_product<LHS,RHS>::value("(scaled_row+15)*" + internal_size_2_expression + " + col","col") + ";\n";
          else
            return "ALIGNMENT NOT IMPLEMENTED";
        }

        static const std::string fill_ith_row(viennacl::column_major)
        {
          std::string internal_size_1_expression = MatExpr::internal_size1_expression();
          VIENNACL_STATIC_ASSERT(Alignment==1);
          return   " dot_prod +=  " + dot_product<LHS,RHS>::value("row  + col * " +  internal_size_1_expression, "col") + ";\n";
          
    //       if(Alignment==1)
    //             return   " dot_prod +=  " + dot_product<LHS,RHS>::value("row  + col * " +  internal_size_1_expression, "col") + ";\n";                    ;
    //       else
    //            return "ALIGNMENT NOT IMPLEMENTED";
    //       
        }
        
      public:
        static const std::string value() 
        {
          return
              "for (unsigned int row = get_global_id(0) ; row < " + MatExpr::internal_size1_expression() + " ; row += get_global_size(0))\n"
              "{\n"
              + print_type<ScalarType,Alignment>::value()+" dot_prod = 0;\n"
              "for (unsigned int col = 0; col < " + MatExpr::internal_size2_expression() + "; ++col){\n"
              + fill_ith_row(Layout() )
              + "}\n"
              + assign_res ( Int2Type<is_null_type<SCALAR_EXPR>::value>() ) + "\n"
              + "}\n";
        }

    };
      
    template <class T, class OP, class ASSIGNED>
    struct make_code<T, OP, ASSIGNED, typename enable_if<is_product_leaf<T> >::type> 
    {
      static const std::string value() 
      {
        typedef typename result_of::expression_type<T>::Result U;
        return make_product_code<U,OP,ASSIGNED>::value();
      }
    };

  }
}

#endif


