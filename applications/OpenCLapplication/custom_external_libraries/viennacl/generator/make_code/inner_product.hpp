#ifndef VIENNACL_GENERATOR_MAKE_CODE_INNER_PRODUCT_HPP
#define VIENNACL_GENERATOR_MAKE_CODE_INNER_PRODUCT_HPP

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

/** @file viennacl/generator/make_code/inner_product.hpp
 *   @brief Directives for generating code for the inner product.
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/make_code/expression.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/tree_operations.hpp"


namespace viennacl 
{
  namespace generator
  {

    template <class T>
    struct inner_prod_impl_t 
    {
      typedef T PRIOR_TYPE;
      
      static const std::string name() 
      {
        return T::name();
      }
      
      static const std::string kernel_arguments() 
      {
        return T::kernel_arguments();
      }
      enum { id = T::id };
    };

    template <class TOKEN, class OP, class ASSIGNED, class Enable=void>
    struct make_code;

    template <class T, class OP, class ASSIGNED>
    struct make_code<inner_prod_impl_t<T>, OP, ASSIGNED> 
    {
      private:
        typedef typename tree_utils::extract_if<T,is_pure_inner_product_leaf>::Result::Head ARG;
        typedef typename ARG::LHS LHS;
        typedef typename ARG::RHS RHS;

        static const std::string main_size() 
        {
          return result_of::expression_type<LHS>::Result::internal_size_expression();
        }

      public :

        static const std::string value() 
        {
          return  "sum = 0;\n"
                  "for (unsigned int k = (get_group_id(0) * " + main_size() + ")/get_num_groups(0)+ get_local_id(0); k < ((get_group_id(0)+1) * " + main_size() +")/get_num_groups(0); k += get_local_size(0))\n"
                  "  sum += " + dot_product<LHS,RHS>::value("k","k") + ";\n"
                  "shared_memory_ptr[get_local_id(0)] = sum;\n"

                  "for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
                  "  {\n"
                  "    barrier(CLK_LOCAL_MEM_FENCE);\n"
                  "    if (get_local_id(0) < stride)\n"
                  "    shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
                  "  }\n"
                  "barrier(CLK_LOCAL_MEM_FENCE);\n"
                  "if (get_local_id(0) == 0)\n"
                  "  " + ASSIGNED::name() + "[get_group_id(0)] = shared_memory_ptr[0];\n";
        }
        
        viennacl::generator::compound_node< const char*, add_type, const char* > value(const char* arg1);
    };

    template <class T, class OP>
    struct make_code<T, OP, T, typename enable_if<is_inner_product_leaf<T> >::type> 
    {
      private:
        typedef typename tree_utils::extract_if<T,is_pure_inner_product_leaf>::Result::Head ARG;
        typedef typename ARG::LHS LHS;
        typedef typename ARG::RHS RHS;

      public:

        static const std::string value() 
        {
          return  "sum = 0;\n"
                  "local float " + ARG::name() + "_sum;\n"
                  "for (unsigned int i = get_local_id(0) ; i<get_num_groups(0) ; i+=get_local_size(0))\n"
                  "{\n"
                  "   sum+= " +ARG::name() +"[i];\n"
                  "};\n"
                  "shared_memory_ptr[get_local_id(0)]=sum;\n"
                  "for (unsigned int stride = get_local_size(0)/2; stride > 0; stride /= 2)\n"
                  "  {\n"
                  "    barrier(CLK_LOCAL_MEM_FENCE);\n"
                  "    if (get_local_id(0) < stride)\n"
                  "     shared_memory_ptr[get_local_id(0)] += shared_memory_ptr[get_local_id(0)+stride];\n"
                  "  }\n"
      "if(get_local_id(0)==0);\n"
                  +ARG::name() + "_sum = shared_memory_ptr[0];\n"
                  "barrier(CLK_LOCAL_MEM_FENCE);\n";
        }
    };

  }

}

#endif


