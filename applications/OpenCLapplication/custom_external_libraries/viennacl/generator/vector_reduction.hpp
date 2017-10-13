#ifndef VIENNACL_GENERATOR_GENERATE_VECTOR_REDUCTION_HPP
#define VIENNACL_GENERATOR_GENERATE_VECTOR_REDUCTION_HPP

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


/** @file viennacl/generator/vector_reduction.hpp
 *
 * @brief Kernel template for the vector reduction operation
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

    /** @brief OpenCL kernel template for reductions resulting in a vector. Example: Computing the row norms of a matrix concurrently. */
    class vector_reduction : public profile_base{

        vcl_size_t lmem_used(vcl_size_t scalartype_size) const {
          return m_*(k_+1)*scalartype_size;
        }

      public:
        /** @brief The user constructor */
        vector_reduction(unsigned int vectorization, unsigned int m, unsigned int k, unsigned int num_groups) : profile_base(vectorization, m, k, 1), m_(m), k_(k), num_groups_(num_groups){ }


        static std::string csv_format() {
          return "Vec,M,K,NumGroups";
        }

        std::string csv_representation() const{
          std::ostringstream oss;
          oss << vector_size_
                 << "," << m_
                 << "," << k_
                 << "," << num_groups_;
          return oss.str();
        }

        unsigned int m() const { return m_; }

        unsigned int k() const { return k_; }

        unsigned int num_groups() const { return num_groups_; }

        void configure_range_enqueue_arguments(vcl_size_t kernel_id, statements_type  const & statements, viennacl::ocl::kernel & kernel, unsigned int & n_arg)  const{

          configure_local_sizes(kernel, kernel_id);
          kernel.global_work_size(0,m_*num_groups_);
          kernel.global_work_size(1,k_);


          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            scheduler::statement::container_type exprs = it->first.array();
            for(scheduler::statement::container_type::iterator iit = exprs.begin() ; iit != exprs.end() ; ++iit){
              if(iit->op.type==scheduler::OPERATION_BINARY_MAT_VEC_PROD_TYPE){
                scheduler::statement_node const * current_node = &(*iit);
                //The LHS of the prod is a matrix
                if(current_node->lhs.type_family==scheduler::MATRIX_TYPE_FAMILY)
                {
                  kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size1_fun())));
                  kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size2_fun())));
                  return;
                }
                else{
                  //The LHS of the prod is a matrix expression
                  current_node = &exprs[current_node->lhs.node_index];
                  if(current_node->lhs.type_family==scheduler::MATRIX_TYPE_FAMILY)
                  {
                    kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size1_fun())));
                    kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size2_fun())));
                    return;
                  }
                  else if(current_node->rhs.type_family==scheduler::MATRIX_TYPE_FAMILY)
                  {
                    kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size1_fun())));
                    kernel.arg(n_arg++, cl_uint(utils::call_on_matrix(current_node->lhs, utils::internal_size2_fun())));
                    return;
                  }
                  else{
                    assert(false && bool("unexpected expression tree"));
                  }
                }
                return;
              }
            }
          }
        }

        void kernel_arguments(statements_type  const & /*statements*/, std::string & arguments_string) const{
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "M");
          arguments_string += detail::generate_value_kernel_argument("unsigned int", "N");
        }

      private:
        void core(vcl_size_t /*kernel_id*/, utils::kernel_generation_stream& stream, statements_type const & statements, std::vector<detail::mapping_type> const & mapping) const {

          std::vector<detail::mapped_vector_reduction*> exprs;
          for(std::vector<detail::mapping_type>::const_iterator it = mapping.begin() ; it != mapping.end() ; ++it){
            for(detail::mapping_type::const_iterator iit = it->begin() ; iit != it->end() ; ++iit){
              if(detail::mapped_vector_reduction * p = dynamic_cast<detail::mapped_vector_reduction*>(iit->second.get()))
                exprs.push_back(p);
              if(detail::mapped_matrix * p = dynamic_cast<detail::mapped_matrix*>(iit->second.get()))
                p->bind_sizes("M","N");
            }
          }

          vcl_size_t lsize1 = m_;
          vcl_size_t lsize2 = k_+1;
          std::string scalartype = "float";
          bool is_lhs_transposed = false;
          if(exprs.front()->root_node().lhs.type_family==scheduler::COMPOSITE_OPERATION_FAMILY)
            if(exprs.front()->statement().array()[exprs.front()->root_node().lhs.node_index].op.type==scheduler::OPERATION_UNARY_TRANS_TYPE)
              is_lhs_transposed = true;

          std::string size1 = "M", size2 = "N";
          if(is_lhs_transposed)
            std::swap(size1, size2);

          for(std::vector<detail::mapped_vector_reduction*>::iterator it = exprs.begin() ; it != exprs.end() ; ++it){
            stream << "__local " <<  (*it)->scalartype() << " buf" << std::distance(exprs.begin(), it) << '[' << lsize1*lsize2 << "];" << std::endl;
          }

          stream << "unsigned int lid0 = get_local_id(0);" << std::endl;
          stream << "unsigned int lid1 = get_local_id(1);" << std::endl;


          stream << "for(unsigned int r = get_global_id(0) ; r < " << size1 << " ; r += get_global_size(0)){" << std::endl;
          stream.inc_tab();

          for(vcl_size_t k = 0 ; k < exprs.size() ; ++k)
            stream << scalartype << " sum" << k << " = 0;" << std::endl;

          stream << "for(unsigned int c = get_local_id(1) ; c < " << size2 << " ; c += get_local_size(1)){" << std::endl;
          stream.inc_tab();

          std::set<std::string>  fetched;

          for(std::vector<detail::mapped_vector_reduction*>::iterator it = exprs.begin() ; it != exprs.end() ; ++it){
            viennacl::scheduler::statement const & statement = (*it)->statement();
            viennacl::scheduler::statement_node const & root_node = (*it)->root_node();
            if(is_lhs_transposed)
              detail::fetch_all_lhs(fetched,statement,root_node, std::make_pair("c", "r"),vector_size_,stream,(*it)->mapping());
            else
              detail::fetch_all_lhs(fetched,statement,root_node, std::make_pair("r", "c"),vector_size_,stream,(*it)->mapping());

            detail::fetch_all_rhs(fetched,statement,root_node, std::make_pair("c", "0"),vector_size_,stream,(*it)->mapping());
          }


          //Update sums;
          for(std::vector<detail::mapped_vector_reduction*>::iterator it = exprs.begin() ; it != exprs.end() ; ++it){
            viennacl::scheduler::statement const & statement = (*it)->statement();
            viennacl::scheduler::statement_node const & root_node = (*it)->root_node();
            std::string str;
            detail::generate_all_lhs(statement,root_node,std::make_pair("i","0"),-1,str,(*it)->mapping());
            str += "*";
            detail::generate_all_rhs(statement,root_node,std::make_pair("i","0"),-1,str,(*it)->mapping());
            stream << " sum" << std::distance(exprs.begin(),it) << " += "  << str << ";" << std::endl;
          }


          stream.dec_tab();
          stream << "}" << std::endl;

          for(vcl_size_t k = 0 ; k < exprs.size() ; ++k){
            stream << "buf" << k << "[lid0*" << lsize2 << "+ lid1] = sum" << k << ";" << std::endl;
          }

          for(unsigned int stride = k_/2 ; stride>1 ; stride /=2){
            stream << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
            stream <<  "if(lid1 < " << stride << ")" ;
            stream << "{" << std::endl;
            stream.inc_tab();

            for(vcl_size_t i = 0 ; i < exprs.size() ; ++i)
              stream << "buf" << i << "[lid0*" << lsize2 << "+ lid1] += buf" << i << "[lid0*" << lsize2 << "+ lid1 + " << stride << "];" << std::endl;

            stream.dec_tab();
            stream << "}" << std::endl;
          }


          stream << "barrier(CLK_LOCAL_MEM_FENCE); " << std::endl;
          stream <<  "if(lid1 == 0)" ;
          stream << "{" << std::endl;
          stream.inc_tab();
          for(vcl_size_t i = 0 ; i < exprs.size() ; ++i){
            stream << "buf" << i << "[lid0*" << lsize2 << "] += buf" << i << "[lid0*" << lsize2 << "+ 1];" << std::endl;
            exprs[i]->access_name("buf"+utils::to_string(i)+"[lid0*"+utils::to_string(lsize2)+"]");
          }
          vcl_size_t i = 0;
          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            std::string str;
            detail::traverse(it->first, it->second, detail::expression_generation_traversal(std::make_pair("r","0"), -1, str, mapping[i++]), false);
            stream << str << ";" << std::endl;
          }
          stream.dec_tab();
          stream << "}" << std::endl;


          stream.dec_tab();
          stream << "}" << std::endl;

        }

      private:
        unsigned int m_;
        unsigned int k_;
        unsigned int num_groups_;
    };
  }
}

#endif
