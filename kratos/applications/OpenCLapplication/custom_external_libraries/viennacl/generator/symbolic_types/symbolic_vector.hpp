#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_VECTOR_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_VECTOR_HPP

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

/** @file viennacl/generator/symbolic_types/symbolic_vector.hpp
 *  @brief Implementation of a symbolic vector type
 *
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/vector.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/traits/result_of.hpp"

namespace viennacl 
{
  namespace generator
  {

    /**
     * @brief Symbolic vector type
     * 
     * @tparam ID The argument ID of the vector in the generated code
     * @tparam SCALARTYPE The Scalartype of the vector in the generated code
     * @tparam ALIGNMENT The Alignment of the vector in the generated code
     */
    template <unsigned int ID, typename SCALARTYPE, unsigned int ALIGNMENT>
    class symbolic_vector 
    {
      private:
        typedef symbolic_vector<ID,SCALARTYPE,ALIGNMENT> self_type;

      public:
	
        typedef SCALARTYPE ScalarType;
        
        static const unsigned int Alignment = ALIGNMENT;

        typedef viennacl::vector<ScalarType,Alignment> runtime_type;

        static const unsigned int id = ID;

        static const std::string name() 
        {
          return "v_a" + to_string(Alignment) + "_" + to_string(ID);
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
          return " __global " + print_type<SCALARTYPE*,Alignment>::value() + " " + name();
        }

        static const std::string kernel_arguments() 
        {
          return " __global " + print_type<SCALARTYPE*,Alignment>::value() + " " + name() 
               + ", unsigned int " + size2_name() 
               + ", unsigned int " + internal_size2_name() + "\n" ;
        }

        template<typename RHS_TYPE>
        typename enable_if<is_same_expression_type<self_type,RHS_TYPE>,
                           compound_node<self_type, assign_type, RHS_TYPE > >::type
        operator= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,assign_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<is_scalar_expression<RHS_TYPE>,
                           compound_node<self_type, inplace_scal_mul_type, RHS_TYPE > >::type
        operator*= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_scal_mul_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<is_scalar_expression<RHS_TYPE>,
                           compound_node<self_type, inplace_scal_div_type, RHS_TYPE > >::type
        operator/= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_scal_div_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<is_same_expression_type<self_type,RHS_TYPE>,
                           compound_node<self_type, inplace_add_type, RHS_TYPE > >::type
        operator+= ( RHS_TYPE const & rhs ) const 
        {
          return compound_node<self_type,inplace_add_type,RHS_TYPE >();
        }

        template<typename RHS_TYPE>
        typename enable_if<is_same_expression_type<self_type,RHS_TYPE>,
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

    template< unsigned int ID,class SCALARTYPE,unsigned int ALIGNMENT>
    class tmp_symbolic_vector<symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > 
    {
        typedef symbolic_vector<ID,SCALARTYPE,ALIGNMENT> ARG;

      public:
        typedef SCALARTYPE ScalarType;

        typedef typename symbolic_vector<ID,SCALARTYPE,ALIGNMENT>::runtime_type runtime_type;
        
        static const unsigned int Alignment = ALIGNMENT;

        static const unsigned int id = ID;


        static const std::string name() 
        {
          return "tmp_" + ARG::name();
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
          return " __global " + print_type<SCALARTYPE*,Alignment>::value() + " " + name();
        }

        static const std::string kernel_arguments() 
        {
            return " __global " + print_type<SCALARTYPE*,Alignment>::value() + " " + name() 
                 + ", unsigned int " + size2_name() 
                 + ", unsigned int " + internal_size2_name() + "\n" ;
        }
    };
  }
}

#endif

