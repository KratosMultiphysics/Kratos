#ifndef MAKE_CODE_REGULAR_COMPOUND_NODE_HPP
#define MAKE_CODE_REGULAR_COMPOUND_NODE_HPP

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

/** @file viennacl/generator/make_code/regular_compound_node.hpp
 *  @brief Directives for generating code for the matrix-vector product
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "expression.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/symbolic_types/symbolic_matrix.hpp"

namespace viennacl
{
  namespace generator
  {

    template <class T>
    struct get_loop_bound_impl;

    template <class T, class SIZE_DESCRIPTOR>
    struct get_loop_bound_impl<result_of::vector_expression<T, SIZE_DESCRIPTOR> > 
    {
      static const std::string value() 
      {
        return result_of::vector_expression<T, SIZE_DESCRIPTOR>::internal_size_expression();
      }
    };

    template <class T, class SIZE1_DESCRIPTOR, class SIZE2_DESCRIPTOR>
    struct get_loop_bound_impl<result_of::matrix_expression<T, SIZE1_DESCRIPTOR, SIZE2_DESCRIPTOR> > 
    {
      private:
        typedef result_of::matrix_expression<T, SIZE1_DESCRIPTOR, SIZE2_DESCRIPTOR> Arg;
        
      public:
        static const std::string value() 
        {
          return Arg::internal_size2_expression() + "*" + Arg::internal_size1_expression();
        }
    };

    template <class T>
    struct get_loop_bound 
    {
      static const std::string value() 
      {
        return get_loop_bound_impl<typename result_of::expression_type<T>::Result>::value();
      }
    };

    template <class T, class ASSIGN_OP, class ASSIGNED, class Enable>
    struct make_code 
    {
      static const std::string value() 
      {
        return "for ( unsigned int k = get_global_id(0)"
               " ; k < " + get_loop_bound<ASSIGNED>::value()
               +" ; k += get_global_size(0) ) \n"
               + "{\n"
               + make_expression_code<ASSIGNED>::value("k") + ASSIGN_OP::expression_string() + make_expression_code<T>::value("k") + ";\n"
               + "}\n";
      }
    };

    template<class T, class ASSIGN_OP, unsigned int ASSIGNED_ID, class ASSIGNED_TYPE>
    struct make_code<T, ASSIGN_OP, gpu_symbolic_scalar<ASSIGNED_ID,ASSIGNED_TYPE> > 
    {
      private:
        typedef gpu_symbolic_scalar<ASSIGNED_ID,ASSIGNED_TYPE>  ASSIGNED;
      public:
        static const std::string value() 
        {
          return "if(get_global_id(0) == 0) " 
                 + make_expression_code<ASSIGNED>::value("0") + '\n' 
                 + ASSIGN_OP::expression_string() + make_expression_code<T>::value ( "k" ) + ";\n" ;
        }
    };

  }
}

#endif


