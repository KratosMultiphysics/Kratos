#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_SCALARS_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_SYMBOLIC_SCALARS_HPP

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

/** @file viennacl/generator/symbolic_types/symbolic_scalars.hpp
 *  @brief Implementation of the symbolic scalar types.
 *
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/scalar.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/traits/result_of.hpp"

namespace viennacl 
{
  namespace generator
  {

    ///////////////////////////////////////
    /////// REGULAR SYM SCALARS //////////
    //////////////////////////////////////

    /**
    * @brief Symbolic scalar type. Will be passed by value.
    * 
    * @tparam ID The argument ID of the scalar in the generated code
    * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
    */
    template <unsigned int ID, typename SCALARTYPE>
    class cpu_symbolic_scalar
    {
      private:
        typedef cpu_symbolic_scalar<ID,SCALARTYPE> self_type;

      public:

        typedef SCALARTYPE ScalarType;

        typedef ScalarType runtime_type;
        
        enum { id = ID };

        static const std::string name() 
        {
          std::ostringstream oss;
          oss << "c_s" << ID ;
          return oss.str();
        }

        static const std::string kernel_arguments() 
        {
          return print_type<SCALARTYPE,1>::value() + " " + name() + "\n";
        }
    };

    /**
     * @brief Symbolic scalar type. Will be passed by pointer.
     * 
     * @tparam ID The argument ID of the scalar in the generated code
     * @tparam SCALARTYPE The Scalartype of the scalar in the generated code
     */
    template <unsigned int ID, typename SCALARTYPE>
    class gpu_symbolic_scalar 
    {
      private:
        typedef gpu_symbolic_scalar<ID,SCALARTYPE> self_type;

      public:

        typedef SCALARTYPE ScalarType;

        typedef viennacl::scalar<ScalarType> runtime_type;
        
        enum { id = ID };

        static const std::string name() 
        {
          std::ostringstream oss;
          oss << "g_s" << ID ;
          return oss.str();
        }

        static const std::string kernel_arguments() 
        {
          return "__global " + print_type<SCALARTYPE*,1>::value() + " " + name() + "\n" ;
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

    ///////////////////////////////////////
    ///////// SCALAR MULTIPLICATION ///////
    //////////////////////////////////////

    /** @brief Scalar multiplication operator */
    template<class LHS_TYPE, class RHS_TYPE>
    typename enable_if_c<is_scalar_expression<LHS_TYPE>::value || is_scalar_expression<RHS_TYPE>::value,
                         compound_node<LHS_TYPE,scal_mul_type,RHS_TYPE> >::type
    operator* ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ) 
    {
      return compound_node<LHS_TYPE, scal_mul_type,RHS_TYPE> ();
    }

    /** @brief Scalar division operator */
    template<class LHS_TYPE, class RHS_TYPE>
    typename enable_if_c< is_scalar_expression<RHS_TYPE>::value,
                          compound_node<LHS_TYPE,scal_div_type,RHS_TYPE> > ::type
    operator/ ( LHS_TYPE const & lhs, RHS_TYPE const & rhs ) 
    {
      return compound_node<LHS_TYPE,scal_div_type,RHS_TYPE> ();
    }

  }
}
#endif
