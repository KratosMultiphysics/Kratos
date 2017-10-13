#ifndef VIENNACL_GENERATOR_MAP_GENERATE_PROTOTYPE_HPP
#define VIENNACL_GENERATOR_MAP_GENERATE_PROTOTYPE_HPP

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


/** @file viennacl/generator/map_functor.hpp
    @brief Functor to map the statements to the types defined in mapped_objects.hpp
*/

#include <set>

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/forwards.h"

#include "viennacl/tools/shared_ptr.hpp"

#include "viennacl/generator/helpers.hpp"
#include "viennacl/generator/utils.hpp"
#include "viennacl/generator/mapped_objects.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{

      /** @brief Functor to map the statements to the types defined in mapped_objects.hpp */
      class map_functor : public traversal_functor{
          std::string create_name(unsigned int & current_arg, std::map<void *, vcl_size_t> & memory, void * handle) const{
            if(handle==NULL)
              return "arg" + utils::to_string(current_arg_++);
            if(memory.insert(std::make_pair(handle, current_arg)).second)
              return "arg" + utils::to_string(current_arg_++);
            else
              return "arg" + utils::to_string(memory[handle]);
          }

        public:
          typedef container_ptr_type result_type;

          map_functor(std::map<void *, vcl_size_t> & memory, unsigned int & current_arg, mapping_type & mapping) : memory_(memory), current_arg_(current_arg), mapping_(mapping){ }

          /** @brief Binary leaf */
          template<class T>
          result_type binary_leaf(viennacl::scheduler::statement const * statement, viennacl::scheduler::statement_node const * root_node, mapping_type const * mapping) const {
            T * p = new T("float");

            p->info_.statement = statement;
            p->info_.root_node = root_node;
            p->info_.mapping = mapping;

            return container_ptr_type(p);
          }

          template<class ScalarType>
          result_type operator()(ScalarType const & /*scal*/) const {
            mapped_host_scalar * p = new mapped_host_scalar(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, NULL);
            return container_ptr_type(p);
          }

          /** @brief Scalar mapping */
          template<class ScalarType>
          result_type operator()(scalar<ScalarType> const & scal) const {
            mapped_scalar * p = new mapped_scalar(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, (void*)&scal);
            return container_ptr_type(p);
          }

          /** @brief Vector mapping */
          template<class ScalarType>
          result_type operator()(vector_base<ScalarType> const & vec) const {
            mapped_vector * p = new mapped_vector(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, (void*)&vec);
            if(vec.start() > 0)
              p->start_name_ = p->name_ +"_start";
            if(vec.stride() > 1)
              p->stride_name_ = p->name_ + "_stride";
            return container_ptr_type(p);
          }

          /** @brief Implicit vector mapping */
          template<class ScalarType>
          result_type operator()(implicit_vector_base<ScalarType> const & vec) const {
            mapped_implicit_vector * p = new mapped_implicit_vector(utils::type_to_string<ScalarType>::value());

            if(vec.is_value_static()==false)
              p->value_name_ = create_name(current_arg_, memory_, NULL);
            if(vec.has_index())
              p->value_name_ = create_name(current_arg_, memory_, NULL);
            return container_ptr_type(p);
          }

          /** @brief Matrix mapping */
          template<class ScalarType, class Layout>
          result_type operator()(matrix_base<ScalarType, Layout> const & mat) const {
            mapped_matrix * p = new mapped_matrix(utils::type_to_string<ScalarType>::value());
            p->name_ = create_name(current_arg_, memory_, (void*)&mat);
            p->is_row_major_ = static_cast<bool>(utils::is_same_type<Layout, viennacl::row_major>::value);
            if(mat.start1() > 0)
              p->start1_name_ = p->name_ +"_start1";
            if(mat.stride1() > 1)
              p->stride1_name_ = p->name_ + "_stride1";
            if(mat.start2() > 0)
              p->start2_name_ = p->name_ +"_start2";
            if(mat.stride2() > 1)
              p->stride2_name_ = p->name_ + "_stride2";
            return container_ptr_type(p);
          }

          /** @brief Implicit matrix mapping */
          template<class ScalarType>
          result_type operator()(implicit_matrix_base<ScalarType> const & mat) const {
            mapped_implicit_matrix * p = new mapped_implicit_matrix(utils::type_to_string<ScalarType>::value());

            if(mat.is_value_static()==false)
              p->value_name_ = create_name(current_arg_, memory_, NULL);

            return container_ptr_type(p);
          }

          /** @brief Traversal functor */
          void operator()(viennacl::scheduler::statement const * statement, viennacl::scheduler::statement_node const * root_node, detail::node_type node_type) const {
            const key_type key(root_node, node_type);
            if(node_type == LHS_NODE_TYPE && root_node->lhs.type_family != viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
                 mapping_.insert(mapping_type::value_type(key, utils::call_on_element(root_node->lhs, *this)));
            else if(node_type == RHS_NODE_TYPE && root_node->rhs.type_family != viennacl::scheduler::COMPOSITE_OPERATION_FAMILY)
                 mapping_.insert(mapping_type::value_type(key,  utils::call_on_element(root_node->rhs, *this)));
            else if( node_type== PARENT_NODE_TYPE){
                  viennacl::scheduler::operation_node_type op_type = root_node->op.type;
                if(op_type == viennacl::scheduler::OPERATION_BINARY_INNER_PROD_TYPE)
                  mapping_.insert(mapping_type::value_type(key, binary_leaf<mapped_scalar_reduction>(statement, root_node, &mapping_)));
                else if(op_type == viennacl::scheduler::OPERATION_BINARY_MAT_VEC_PROD_TYPE)
                  mapping_.insert(mapping_type::value_type(key, binary_leaf<mapped_vector_reduction>(statement, root_node, &mapping_)));
                else if(op_type == viennacl::scheduler::OPERATION_BINARY_MAT_MAT_PROD_TYPE)
                  mapping_.insert(mapping_type::value_type(key, binary_leaf<mapped_matrix_product>(statement, root_node, &mapping_)));
            }
          }

        private:
          std::map<void *, vcl_size_t> & memory_;
          unsigned int & current_arg_;
          mapping_type & mapping_;
      };

    }

  }

}
#endif
