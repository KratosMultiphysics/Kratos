#ifndef VIENNACL_GENERATOR_GENERATE_UTILS_HPP
#define VIENNACL_GENERATOR_GENERATE_UTILS_HPP

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */


/** @file viennacl/generator/helpers.hpp
    @brief several code generation helpers
*/

#include <set>

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include "CL/cl.h"
#endif

#include "viennacl/forwards.h"
#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/utils.hpp"
#include "viennacl/generator/forwards.h"

namespace viennacl{

  namespace generator{

    namespace detail{

    /** @brief generate the string for a pointer kernel argument */
      static std::string generate_value_kernel_argument(std::string const & scalartype, std::string const & name){
        return scalartype + ' ' + name + ",";
      }

      /** @brief generate the string for a pointer kernel argument */
      static std::string generate_pointer_kernel_argument(std::string const & address_space, std::string const & scalartype, std::string const & name){
        return address_space +  " " + scalartype + "* " + name + ",";
      }

      /** @brief generate a string from an operation_node_type */
      inline const char * generate(viennacl::scheduler::operation_node_type type){
        // unary expression
        switch(type){
          case viennacl::scheduler::OPERATION_UNARY_ABS_TYPE : return "abs";
          case viennacl::scheduler::OPERATION_UNARY_TRANS_TYPE : return "trans";
          case viennacl::scheduler::OPERATION_BINARY_ASSIGN_TYPE : return "=";
          case viennacl::scheduler::OPERATION_BINARY_INPLACE_ADD_TYPE : return "+=";
          case viennacl::scheduler::OPERATION_BINARY_INPLACE_SUB_TYPE : return "-=";
          case viennacl::scheduler::OPERATION_BINARY_ADD_TYPE : return "+";
          case viennacl::scheduler::OPERATION_BINARY_SUB_TYPE : return "-";
          case viennacl::scheduler::OPERATION_BINARY_MULT_TYPE : return "*";
          case viennacl::scheduler::OPERATION_BINARY_DIV_TYPE : return "/";
          case viennacl::scheduler::OPERATION_BINARY_INNER_PROD_TYPE : return "iprod";
          case viennacl::scheduler::OPERATION_BINARY_MAT_MAT_PROD_TYPE : return "mmprod";
          case viennacl::scheduler::OPERATION_BINARY_MAT_VEC_PROD_TYPE : return "mvprod";
          case viennacl::scheduler::OPERATION_BINARY_ACCESS_TYPE : return "[]";
          default : throw "not implemented";
        }
      }

      /** @brief checks whether an operator is both a binary node and a leaf */
      inline bool is_binary_leaf_operator(viennacl::scheduler::operation_node_type const & op_type) {
        return op_type == viennacl::scheduler::OPERATION_BINARY_INNER_PROD_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_MAT_VEC_PROD_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_MAT_MAT_PROD_TYPE;
      }

      /** @brief checks whether an operator is arithmetic or not */
      inline bool is_arithmetic_operator(viennacl::scheduler::operation_node_type const & op_type) {
        return op_type == viennacl::scheduler::OPERATION_BINARY_ASSIGN_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_ADD_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_DIV_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_ELEMENT_DIV_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_ELEMENT_PROD_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_INPLACE_ADD_TYPE
             ||op_type == viennacl::scheduler::OPERATION_BINARY_INPLACE_SUB_TYPE
//                 ||op_type == viennacl::scheduler::OPERATION_BINARY_INPLACE_DIV_TYPE
//                ||op_type == viennacl::scheduler::OPERATION_BINARY_INPLACE_MULT_TYPE
            ||op_type == viennacl::scheduler::OPERATION_BINARY_MULT_TYPE
            ||op_type == viennacl::scheduler::OPERATION_BINARY_SUB_TYPE;

      }

      /** @brief Recursively execute a functor on a statement */
      template<class Fun>
      static void traverse(viennacl::scheduler::statement const & statement, viennacl::scheduler::statement_node const & root_node, Fun const & fun, bool recurse_binary_leaf /* see forwards.h for default argument */){

        if(root_node.op.type_family==viennacl::scheduler::OPERATION_UNARY_TYPE_FAMILY)
        {
          //Self:
          fun(&statement, &root_node, PARENT_NODE_TYPE);

          //Lhs:
          fun.call_before_expansion();
          if(root_node.lhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.lhs.node_index], fun, recurse_binary_leaf);
          fun(&statement, &root_node, LHS_NODE_TYPE);
          fun.call_after_expansion();
        }
        else if(root_node.op.type_family==viennacl::scheduler::OPERATION_BINARY_TYPE_FAMILY)
        {
          bool deep_recursion = recurse_binary_leaf || !is_binary_leaf_operator(root_node.op.type);

          fun.call_before_expansion();

          //Lhs:
          if(deep_recursion){
            if(root_node.lhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.lhs.node_index], fun, recurse_binary_leaf);
            fun(&statement, &root_node, LHS_NODE_TYPE);
          }

          //Self:
          fun(&statement, &root_node, PARENT_NODE_TYPE);

          //Rhs:
          if(deep_recursion){
            if(root_node.rhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
              traverse(statement, statement.array()[root_node.rhs.node_index], fun, recurse_binary_leaf);
            fun(&statement, &root_node, RHS_NODE_TYPE);
          }

          fun.call_after_expansion();

        }
      }

      /** @brief base functor class for traversing a statement */
      class traversal_functor{
        public:
          void call_before_expansion() const { }
          void call_after_expansion() const { }
      };

      /** @brief functor for generating the prototype of a statement */
      class prototype_generation_traversal : public traversal_functor{
        private:
          std::set<std::string> & already_generated_;
          std::string & str_;
          unsigned int vector_size_;
          mapping_type const & mapping_;
        public:
          prototype_generation_traversal(std::set<std::string> & already_generated, std::string & str, unsigned int vector_size, mapping_type const & mapping) : already_generated_(already_generated), str_(str), vector_size_(vector_size), mapping_(mapping){ }

          void operator()(viennacl::scheduler::statement const *, viennacl::scheduler::statement_node const * root_node, detail::node_type node_type) const {
              if( (node_type==detail::LHS_NODE_TYPE && root_node->lhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
                ||(node_type==detail::RHS_NODE_TYPE && root_node->rhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY) )
                  append_kernel_arguments(already_generated_, str_, vector_size_, *at(mapping_, std::make_pair(root_node,node_type)));
          }
      };

      /** @brief functor for fetching the elements of a statement */
      class fetch_traversal : public traversal_functor{
        private:
          std::set<std::string> & fetched_;
          std::pair<std::string, std::string> index_string_;
          unsigned int vectorization_;
          utils::kernel_generation_stream & stream_;
          mapping_type const & mapping_;
        public:
          fetch_traversal(std::set<std::string> & fetched, std::pair<std::string, std::string> const & index, unsigned int vectorization, utils::kernel_generation_stream & stream, mapping_type const & mapping) : fetched_(fetched), index_string_(index), vectorization_(vectorization), stream_(stream), mapping_(mapping){ }

          void operator()(viennacl::scheduler::statement const *, viennacl::scheduler::statement_node const * root_node, detail::node_type node_type) const {
            if( (node_type==detail::LHS_NODE_TYPE && root_node->lhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
              ||(node_type==detail::RHS_NODE_TYPE && root_node->rhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY) )
              fetch(index_string_, vectorization_, fetched_, stream_, *at(mapping_, std::make_pair(root_node, node_type)));
          }
      };

      /** @brief functor for fetching the LHS of a statement's node
      *
      *   Forwards to fetch_traversal functor if the LHS is not a leaf
      */
      static void fetch_all_lhs(std::set<std::string> & fetched
                                , viennacl::scheduler::statement const & statement
                                , viennacl::scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , vcl_size_t const & vectorization
                                , utils::kernel_generation_stream & stream
                                , detail::mapping_type const & mapping){
        if(root_node.lhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.lhs.node_index], detail::fetch_traversal(fetched, index, static_cast<unsigned int>(vectorization), stream, mapping));
        else
          detail::fetch(index, static_cast<unsigned int>(vectorization),fetched, stream, *at(mapping, std::make_pair(&root_node,detail::LHS_NODE_TYPE)));

      }

      /** @brief functor for fetching the RHS of a statement's node
      *
      *   Forwards to fetch_traversal functor if the RHS is not a leaf
      */
      static void fetch_all_rhs(std::set<std::string> & fetched
                                , viennacl::scheduler::statement const & statement
                                , viennacl::scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , vcl_size_t const & vectorization
                                , utils::kernel_generation_stream & stream
                                , detail::mapping_type const & mapping){
        if(root_node.rhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.rhs.node_index], detail::fetch_traversal(fetched, index, static_cast<unsigned int>(vectorization), stream, mapping));
        else
          detail::fetch(index, static_cast<unsigned int>(vectorization),fetched, stream, *at(mapping, std::make_pair(&root_node,detail::RHS_NODE_TYPE)));

      }


      /** @brief functor for generating the expression string from a statement */
      class expression_generation_traversal : public traversal_functor{
        private:
          std::pair<std::string, std::string> index_string_;
          int vector_element_;
          std::string & str_;
          mapping_type const & mapping_;

        public:
          expression_generation_traversal(std::pair<std::string, std::string> const & index, int vector_element, std::string & str, mapping_type const & mapping) : index_string_(index), vector_element_(vector_element), str_(str), mapping_(mapping){ }

          void call_before_expansion() const { str_+="("; }
          void call_after_expansion() const { str_+=")"; }

          void operator()(viennacl::scheduler::statement const *, viennacl::scheduler::statement_node const * root_node, detail::node_type node_type) const {
            if(node_type==PARENT_NODE_TYPE)
            {
              if(is_binary_leaf_operator(root_node->op.type))
                str_ += generate(index_string_, vector_element_, *at(mapping_, std::make_pair(root_node, node_type)));
              else if(is_arithmetic_operator(root_node->op.type))
                str_ += generate(root_node->op.type);
            }
            else{
              if(node_type==LHS_NODE_TYPE){
                if(root_node->lhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
                  str_ += detail::generate(index_string_,vector_element_, *at(mapping_, std::make_pair(root_node,node_type)));
              }
              else if(node_type==RHS_NODE_TYPE){
                if(root_node->rhs.type_family!=viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
                  str_ += detail::generate(index_string_,vector_element_, *at(mapping_, std::make_pair(root_node,node_type)));
              }
            }
          }
      };

      static void generate_all_lhs(viennacl::scheduler::statement const & statement
                                , viennacl::scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , int vector_element
                                , std::string & str
                                , detail::mapping_type const & mapping){
        if(root_node.lhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.lhs.node_index], detail::expression_generation_traversal(index, vector_element, str, mapping));
        else
          str += detail::generate(index, vector_element,*at(mapping, std::make_pair(&root_node,detail::LHS_NODE_TYPE)));
      }


      static void generate_all_rhs(viennacl::scheduler::statement const & statement
                                , viennacl::scheduler::statement_node const & root_node
                                , std::pair<std::string, std::string> const & index
                                , int vector_element
                                , std::string & str
                                , detail::mapping_type const & mapping){
        if(root_node.rhs.type_family==viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
          detail::traverse(statement, statement.array()[root_node.rhs.node_index], detail::expression_generation_traversal(index, vector_element, str, mapping));
        else
          str += detail::generate(index, vector_element,*at(mapping, std::make_pair(&root_node,detail::RHS_NODE_TYPE)));
      }

    }
  }
}
#endif
