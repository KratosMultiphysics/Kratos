#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_MATRIX_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_MATRIX_HPP

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

/** @file viennacl/generator/symbolic_types/symbolic_matrix.hpp
 *  @brief Implementation of a symbolic matrix type
 *
 *  Generator code contributed by Philippe Tillet
 */



#include "viennacl/forwards.h"
#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"

namespace viennacl 
{
  namespace generator
  {

    /**
    * @brief Symbolic matrix type
    * 
    * @tparam ID The argument ID of the matrix in the generated code
    * @tparam SCALARTYPE The Scalartype of the matrix in the generated code
    * @tparam F The Layout of the matrix in the generated code
    * @tparam ALIGNMENT The Alignment of the matrix in the generated code
    */
    template<unsigned int ID, typename SCALARTYPE, class F, unsigned int ALIGNMENT>
    class symbolic_matrix 
    {
        typedef symbolic_matrix<ID, SCALARTYPE, F, ALIGNMENT> self_type;

      public:

        enum { id = ID };

        typedef SCALARTYPE ScalarType;

        typedef F Layout;
        
        static const unsigned int Alignment = ALIGNMENT;
        
        typedef viennacl::matrix<ScalarType,F,Alignment> runtime_type;
        
        static const std::string name()
        {
          F layout;
          return "m_a_" + viennacl::generator::to_string(layout) + "_" 
                        + viennacl::generator::to_string(Alignment) + "_"
                        + viennacl::generator::to_string<long>(id);
        }

        static const std::string size1_name() 
        {
          return "size1_" + name();
        }

        static const std::string size2_name() 
        {
          return "size2_" + name();
        }

        static const std::string internal_size1_name() 
        {
          return "internal_size1_" + name();
        }

        static const std::string internal_size2_name() 
        {
          return "internal_size2_" + name();
        }

        static const std::string kernel_arguments() 
        {
          return " __global " + generator::print_type<SCALARTYPE*,Alignment>::value() + " " + name()
                + ", unsigned int " + size1_name()
                + ", unsigned int " + size2_name()
                + ", unsigned int " + internal_size1_name()
                + ", unsigned int " + internal_size2_name()
                + "\n";
        }

        template<typename RHS_TYPE>
        typename enable_if<generator::is_same_expression_type<self_type,RHS_TYPE>,
                           compound_node<self_type, assign_type, RHS_TYPE > >::type
        operator= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,assign_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<generator::is_scalar_expression<RHS_TYPE>,
                           compound_node<self_type, inplace_scal_mul_type, RHS_TYPE > >::type
        operator*= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_scal_mul_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<generator::is_scalar_expression<RHS_TYPE>,
                           compound_node<self_type, inplace_scal_div_type, RHS_TYPE > >::type
        operator/= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_scal_div_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<generator::is_same_expression_type<self_type,RHS_TYPE>,
                           compound_node<self_type, inplace_add_type, RHS_TYPE > >::type
        operator+= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_add_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<generator::is_same_expression_type<self_type,RHS_TYPE>,
                           compound_node<self_type, inplace_sub_type, RHS_TYPE > >::type
        operator-= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_sub_type,RHS_TYPE >();
        }

        operator compound_node<self_type,assign_type,self_type>() 
        {
          return compound_node<self_type,assign_type,self_type>();
        }
    };

    template<unsigned int ID,typename SCALARTYPE, class F, unsigned int ALIGNMENT>
    class tmp_symbolic_matrix<symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> > {};
    
  } // namespace generator
} // namespace viennacl

#endif


