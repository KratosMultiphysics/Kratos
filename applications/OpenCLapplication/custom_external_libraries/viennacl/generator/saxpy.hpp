#ifndef VIENNACL_GENERATOR_GENERATE_SAXPY_HPP
#define VIENNACL_GENERATOR_GENERATE_SAXPY_HPP

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


/** @file viennacl/generator/saxpy.hpp
 *
 * @brief Kernel template for the saxpy-like operation
*/

#include <vector>

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/mapped_objects.hpp"
#include "viennacl/generator/helpers.hpp"
#include "viennacl/generator/utils.hpp"

#include "viennacl/generator/profile_base.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl{

  namespace generator{

    /** @brief OpenCL kernel generation class for vector expressions of AXPY type, i.e. x = alpha * y + beta * z, where the number of summands can in principle be arbitrarily large. */
    class vector_saxpy : public profile_base{
      public:
        static std::string csv_format() {
          return "Vec,LSize1,NumGroups1,GlobalDecomposition";
        }

        std::string csv_representation() const{
          std::ostringstream oss;
          oss << vector_size_
              << "," << local_size_1_
              << "," << num_groups_
              << "," << decomposition_;
          return oss.str();
        }

        vector_saxpy(unsigned int v, vcl_size_t gs, vcl_size_t ng, unsigned int d) : profile_base(v, gs, 1, 1), num_groups_(ng), decomposition_(d){ }

        void configure_range_enqueue_arguments(vcl_size_t kernel_id, statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg)  const{
          configure_local_sizes(k, kernel_id);

          k.global_work_size(0,local_size_1_*num_groups_);
          k.global_work_size(1,1);

          scheduler::statement_node const & first_node = statements.front().second;
          viennacl::vcl_size_t N = utils::call_on_vector(first_node.lhs, utils::internal_size_fun());
          k.arg(n_arg++, cl_uint(N/vector_size_));
        }
        void kernel_arguments(statements_type  const & /*statements*/, std::string & arguments_string) const{
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
        }

      private:

        void core(vcl_size_t /*kernel_id*/, utils::kernel_generation_stream& stream, statements_type const & statements, std::vector<detail::mapping_type> const & mapping) const {
          stream << "for(unsigned int i = get_global_id(0) ; i < N ; i += get_global_size(0))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();

          //Fetches entries to registers
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::const_iterator it = mapping.begin() ; it != mapping.end() ; ++it)
            for(detail::mapping_type::const_reverse_iterator iit = it->rbegin() ; iit != it->rend() ; ++iit)
              //Useless to fetch cpu scalars into registers
              if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(iit->second.get()))
                p->fetch( std::make_pair("i","0"), vector_size_, fetched, stream);

          //Generates all the expression, in order
          vcl_size_t i = 0;
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            std::string str;
            detail::traverse(it->first, it->second, detail::expression_generation_traversal(std::make_pair("i","0"), -1, str, mapping[i++]));
            stream << str << ";" << std::endl;
          }

          //Writes back
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it)
             //Gets the mapped object at the LHS of each expression
            if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(at(mapping.at(std::distance(statements.begin(),it)), std::make_pair(&it->second, detail::LHS_NODE_TYPE)).get()))
              p->write_back( std::make_pair("i", "0"), fetched, stream);

          stream.dec_tab();
          stream << "}" << std::endl;
        }

      private:
        vcl_size_t num_groups_;
        unsigned int decomposition_;

    };



    /** @brief OpenCL kernel generation class for matrix expressions of AXPY type, i.e. A = alpha * B + beta * C, where the number of summands can in principle be arbitrarily large. */
    class matrix_saxpy : public profile_base{

        bool invalid_impl(viennacl::ocl::device const & /*dev*/, vcl_size_t /*scalartype_size*/) const{ return false; }
        bool is_slow_impl(viennacl::ocl::device const &) const { return false; }

      public:
        matrix_saxpy(unsigned int v, vcl_size_t gs1, vcl_size_t gs2, vcl_size_t ng1, vcl_size_t ng2, unsigned int d) : profile_base(v, gs1, gs2, 1), num_groups_row_(ng1), num_groups_col_(ng2), decomposition_(d){ }

        static std::string csv_format() {
          return "Vec,LSize1,LSize2,NumGroups1,NumGroups2,GlobalDecomposition";
        }

        std::string csv_representation() const{
          std::ostringstream oss;
          oss << vector_size_
                 << "," << local_size_1_
                 << "," << local_size_2_
                 << "," << num_groups_row_
                 << "," << num_groups_col_
                 << "," << decomposition_;
          return oss.str();
        }

        void configure_range_enqueue_arguments(vcl_size_t kernel_id, statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg)  const{
          configure_local_sizes(k, kernel_id);

          k.global_work_size(0,local_size_1_*num_groups_row_);
          k.global_work_size(1,local_size_2_*num_groups_col_);

          scheduler::statement_node const & first_node = statements.front().second;
          k.arg(n_arg++, cl_uint(utils::call_on_matrix(first_node.lhs, utils::internal_size1_fun())));
          k.arg(n_arg++, cl_uint(utils::call_on_matrix(first_node.lhs, utils::internal_size2_fun())));
        }

        void kernel_arguments(statements_type  const & /*statements*/, std::string & arguments_string) const{
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
        }

      private:
        void core(vcl_size_t /*kernel_id*/, utils::kernel_generation_stream& stream, statements_type const & statements, std::vector<detail::mapping_type> const & mapping) const {

          for(std::vector<detail::mapping_type>::const_iterator it = mapping.begin() ; it != mapping.end() ; ++it){
            for(detail::mapping_type::const_iterator iit = it->begin() ; iit != it->end() ; ++iit){
              if(detail::mapped_matrix * p = dynamic_cast<detail::mapped_matrix*>(iit->second.get()))
                p->bind_sizes("M","N");
            }
          }

          stream << "for(unsigned int i = get_global_id(0) ; i < M ; i += get_global_size(0))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();
          stream << "for(unsigned int j = get_global_id(1) ; j < N ; j += get_global_size(1))" << std::endl;
          stream << "{" << std::endl;
          stream.inc_tab();

          //Fetches entries to registers
          std::set<std::string>  fetched;
          for(std::vector<detail::mapping_type>::const_iterator it = mapping.begin() ; it != mapping.end() ; ++it)
            for(detail::mapping_type::const_reverse_iterator it2 = it->rbegin() ; it2 != it->rend() ; ++it2)
              if(detail::mapped_matrix * p = dynamic_cast<detail::mapped_matrix *>(it2->second.get()))
                p->fetch(std::make_pair("i", "j"), vector_size_, fetched, stream);


          vcl_size_t i = 0;
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            std::string str;
            detail::traverse(it->first, it->second, detail::expression_generation_traversal(std::make_pair("i", "j"), -1, str, mapping[i++]));
            stream << str << ";" << std::endl;
          }

          //Writes back
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            if(detail::mapped_handle * p = dynamic_cast<detail::mapped_handle *>(at(mapping.at(std::distance(statements.begin(),it)), std::make_pair(&it->second,detail::LHS_NODE_TYPE)).get()))
              p->write_back(std::make_pair("i", "j"), fetched, stream);
          }

          stream.dec_tab();
          stream << "}" << std::endl;
          stream.dec_tab();
          stream << "}" << std::endl;
        }

      private:
        vcl_size_t num_groups_row_;
        vcl_size_t num_groups_col_;

        unsigned int decomposition_;
    };
  }

}

#endif
