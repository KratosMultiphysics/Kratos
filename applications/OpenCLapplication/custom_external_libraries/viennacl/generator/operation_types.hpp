#ifndef VIENNACL_GENERATOR_OPERATION_TYPES_HPP
#define VIENNACL_GENERATOR_OPERATION_TYPES_HPP

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

/** @file viennacl/generator/operation_types.hpp
 *   @brief Declaration of the types related to the operators
 */

#include <sstream>

namespace viennacl 
{
  namespace generator
  {

    struct assign_type 
    {
      static const std::string expression_string() { return " = "; }
      static const std::string name() { return "eq"; }
    };

    struct add_type 
    {
      static const std::string expression_string() { return " + "; }
      static const std::string name() { return "p"; }
    };

    struct inplace_add_type 
    {
      static const std::string expression_string() { return " += "; }
      static const std::string name() { return "p_eq"; }
    };

    struct sub_type 
    {
      static const std::string expression_string() { return " - "; }
      static const std::string name() { return "m"; }
    };

    struct inplace_sub_type 
    {
      static const std::string expression_string() { return " -= "; }
      static const std::string name() { return "m_eq"; }
    };

    struct scal_mul_type 
    {
      static const std::string expression_string() { return " * "; }
      static const std::string name() { return "mu"; }
    };

    struct inplace_scal_mul_type 
    {
      static const std::string expression_string() { return " *= "; }
      static const std::string name() { return "mu_eq"; }
    };


    struct scal_div_type
    {
      static const std::string expression_string() { return " / "; }
      static const std::string name() { return "d"; }
    };

    struct inplace_scal_div_type 
    {
      static const std::string expression_string() { return " /= "; }
      static const std::string name() { return "d_eq"; }
    };

    struct inner_prod_type 
    {
      static const std::string expression_string() { return "_i_"; }
      static const std::string name() { return "i"; }
    };

    struct prod_type 
    {
      static const std::string expression_string() { return "_p_"; }
      static const std::string name() { return "p"; }
    };

    template<class T>
    struct make_inplace 
    {
      typedef T Result;
    };

    template<>
    struct make_inplace<add_type> 
    {
      typedef inplace_add_type Result;
    };

    template<>
    struct make_inplace<sub_type> 
    {
      typedef inplace_sub_type Result;
    };

    template<>
    struct make_inplace<scal_mul_type> 
    {
      typedef inplace_scal_mul_type Result;
    };

    template<>
    struct make_inplace<scal_div_type> 
    {
      typedef inplace_scal_div_type Result;
    };

  }
}
#endif
