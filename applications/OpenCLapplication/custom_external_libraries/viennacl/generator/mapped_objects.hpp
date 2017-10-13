#ifndef VIENNACL_GENERATOR_MAPPED_TYPE_HPP
#define VIENNACL_GENERATOR_MAPPED_TYPE_HPP

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


/** @file viennacl/generator/mapped_objects.hpp
    @brief Map ViennaCL objects to generator wrappers
*/

#include <string>

#include "viennacl/scheduler/forwards.h"
#include "viennacl/generator/forwards.h"
#include "viennacl/generator/utils.hpp"

namespace viennacl{

  namespace generator{

    namespace detail{


      /** @brief Base class for mapping viennacl datastructure to generator-friendly structures
       */
      class mapped_object{
        protected:
          /** \cond */
          struct node_info{
              node_info() : mapping(NULL), statement(NULL), root_node(NULL) { }
              mapping_type const * mapping;
              scheduler::statement const * statement;
              scheduler::statement_node const * root_node;
          };
          /** \endcond */
          virtual std::string generate_default(std::pair<std::string, std::string> const & index) const = 0;
          virtual std::string append_vector_size(std::string const & scalartype, unsigned int) const { return scalartype; }

        public:
          mapped_object(std::string const & scalartype) : scalartype_(scalartype){          }
          virtual std::string & append_kernel_arguments(std::set<std::string> &, std::string & str, unsigned int) const{ return str; }
          std::string const & scalartype() const { return scalartype_; }
          void access_name(std::string const & str) { access_name_ = str; }
          std::string const & access_name() const { return access_name_; }
          virtual std::string generate(std::pair<std::string, std::string> const & index, int) const{
            if(!access_name_.empty())
              return access_name_;
            else
              return generate_default(index);
          }
          virtual ~mapped_object(){ }
        protected:
          std::string access_name_;
          std::string scalartype_;
      };

      /** @brief Base class for mapping binary leaves (inner product-based, matrix vector product-base, matrix-matrix product based...)
       */
      class mapped_binary_leaf : public mapped_object{
        public:
          mapped_binary_leaf(std::string const & scalartype) : mapped_object(scalartype){ }
          mapping_type const & mapping() const { return *info_.mapping; }
          scheduler::statement const & statement() const { return *info_.statement; }
          scheduler::statement_node const & root_node() const { return *info_.root_node; }
          std::string generate_default(std::pair<std::string, std::string> const &) const { return "";}
        protected:
          node_info info_;
      };

      /** @brief Mapping of a matrix product */
      class mapped_matrix_product : public mapped_binary_leaf{
          friend class map_functor;
        public:
          mapped_matrix_product(std::string const & scalartype) : mapped_binary_leaf(scalartype){ }
      };

      /** @brief Base class for mapping a reduction */
      class mapped_reduction : public mapped_binary_leaf{
        public:
          mapped_reduction(std::string const & scalartype) : mapped_binary_leaf(scalartype){ }
          viennacl::scheduler::operation_node_type reduction_type() const { return reduction_type_; }
        private:
          viennacl::scheduler::operation_node_type reduction_type_;
      };

      /** @brief Mapping of a scalar reduction (based on inner product) */
      class mapped_scalar_reduction : public mapped_reduction{
          friend class map_functor;
        public:
          mapped_scalar_reduction(std::string const & scalartype) : mapped_reduction(scalartype){ }
      };

      /** @brief Mapping of a vector reduction (based on matrix-vector product) */
      class mapped_vector_reduction : public mapped_reduction{
          friend class map_functor;
        public:
          mapped_vector_reduction(std::string const & scalartype) : mapped_reduction(scalartype){ }
      };

      /** @brief Mapping of a host scalar to a generator class */
      class mapped_host_scalar : public mapped_object{
          friend class map_functor;
          std::string generate_default(std::pair<std::string, std::string> const &) const{ return name_;  }
        public:
          mapped_host_scalar(std::string const & scalartype) : mapped_object(scalartype){ }
          std::string const & name() { return name_; }
          std::string & append_kernel_arguments(std::set<std::string> & already_generated, std::string & str, unsigned int) const{
            if(already_generated.insert(name_).second)
              str += detail::generate_value_kernel_argument(scalartype_, name_);
            return str;
          }

        private:
          std::string name_;
      };

      /** @brief Base class for datastructures passed by pointer */
      class mapped_handle : public mapped_object{
          virtual std::string offset(std::pair<std::string, std::string> const & index) const = 0;
          virtual void append_optional_arguments(std::string &) const{ }
          std::string generate_default(std::pair<std::string, std::string> const & index) const{ return name_  + '[' + offset(index) + ']'; }
        public:
          mapped_handle(std::string const & scalartype) : mapped_object(scalartype){ }

          std::string const & name() const { return name_; }

          void fetch(std::pair<std::string, std::string> const & index, unsigned int vectorization, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string new_access_name = name_ + "_private";
            if(fetched.find(name_)==fetched.end()){
              stream << scalartype_;
              if(vectorization > 1) stream << vectorization;
              stream << " " << new_access_name << " = " << generate_default(index) << ';' << std::endl;
              fetched.insert(name_);
            }
            access_name_ = new_access_name;
          }

          void write_back(std::pair<std::string, std::string> const & index, std::set<std::string> & fetched, utils::kernel_generation_stream & stream) {
            std::string old_access_name = access_name_ ;
            access_name_ = "";
            if(fetched.find(name_)!=fetched.end()){
              stream << generate_default(index) << " = " << old_access_name << ';' << std::endl;
              fetched.erase(name_);
            }
          }

          std::string & append_kernel_arguments(std::set<std::string> & already_generated, std::string & str, unsigned int vector_size) const{
            if(already_generated.insert(name_).second){
              std::string vector_scalartype = append_vector_size(scalartype_, vector_size);
              str += detail::generate_pointer_kernel_argument("__global", vector_scalartype, name_);
              append_optional_arguments(str);
            }
            return str;
          }

        protected:
          std::string name_;
      };

      /** @brief Mapping of a scalar to a generator class */
      class mapped_scalar : public mapped_handle{
          friend class map_functor;
        private:
          std::string offset(std::pair<std::string, std::string> const &)  const { return "0"; }
        public:
          mapped_scalar(std::string const & scalartype) : mapped_handle(scalartype){ }
      };


      /** @brief Base class for mapping buffer-based objects to a generator class */
      class mapped_buffer : public mapped_handle{
        protected:
          std::string append_vector_size(std::string const & scalartype, unsigned int vector_size) const {
            if(vector_size>1)
              return scalartype + utils::to_string(vector_size);
            else
              return scalartype;
          }
        public:
          mapped_buffer(std::string const & scalartype) : mapped_handle(scalartype){ }
          virtual std::string generate(std::pair<std::string, std::string> const & index, int vector_element) const{
            if(vector_element>-1)
              return mapped_object::generate(index, vector_element)+".s"+utils::to_string(vector_element);
            return mapped_object::generate(index, vector_element);
          }

      };

      /** @brief Mapping of a vector to a generator class */
      class mapped_vector : public mapped_buffer{
          friend class map_functor;
          std::string offset(std::pair<std::string, std::string> const & index) const {
            if(info_.statement){
              std::string str;
              detail::generate_all_rhs(*info_.statement, *info_.root_node, index, -1, str, *info_.mapping);
              return str;
            }
            else
              return index.first;
          }

          void append_optional_arguments(std::string & str) const{
            if(!start_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", start_name_);
            if(!stride_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", stride_name_);
            if(!shift_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", shift_name_);
          }
        public:
          mapped_vector(std::string const & scalartype) : mapped_buffer(scalartype){ }
        private:
          node_info info_;

          std::string start_name_;
          std::string stride_name_;
          std::string shift_name_;
      };

      /** @brief Mapping of a matrix to a generator class */
      class mapped_matrix : public mapped_buffer{
          friend class map_functor;
          void append_optional_arguments(std::string & str) const{
            if(!start1_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", start1_name_);
            if(!stride1_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", stride1_name_);
            if(!start2_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", start2_name_);
            if(!stride2_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", stride2_name_);
          }
        public:
          mapped_matrix(std::string const & scalartype) : mapped_buffer(scalartype){ }

          bool is_row_major() const { return is_row_major_; }

          std::string const & size1() const { return size1_; }

          std::string const & size2() const { return size2_; }

          void bind_sizes(std::string const & size1, std::string const & size2) const{
            size1_ = size1;
            size2_ = size2;
          }

          std::string offset(std::pair<std::string, std::string> const & index) const {
            std::string i = index.first;
            std::string j = index.second;
            if(is_row_major_)
              if(j=="0")
                return '(' + i + ')' + '*' + size2_;
              else
                return '(' + i + ')' + '*' + size2_ + "+ (" + j + ')';
            else
              if(i=="0")
                return  "(" + j + ')' + '*' + size1_;
              else
                return  '(' + i + ')' + "+ (" + j + ')' + '*' + size1_;
          }

        private:
          mutable std::string size1_;
          mutable std::string size2_;

          std::string start1_name_;
          std::string stride1_name_;
          std::string shift1_name_;
          std::string start2_name_;
          std::string stride2_name_;
          std::string shift2_name_;
          bool is_row_major_;
      };

      /** @brief Mapping of a implicit vector to a generator class */
      class mapped_implicit_vector : public mapped_object{
          friend class map_functor;
          std::string value_name_;
          std::string index_name_;
        public:
          mapped_implicit_vector(std::string const & scalartype) : mapped_object(scalartype){ }
          std::string generate_default(std::pair<std::string, std::string> const & /*index*/) const{
            return value_name_;
          }
          std::string & append_kernel_arguments(std::set<std::string> & /*already_generated*/, std::string & str, unsigned int /*vector_size*/) const{
            if(!value_name_.empty())
              str += detail::generate_value_kernel_argument(scalartype_, value_name_);
            if(!index_name_.empty())
              str += detail::generate_value_kernel_argument("unsigned int", index_name_);
            return str;
          }
      };

      /** @brief Mapping of a implicit matrix to a generator class */
      class mapped_implicit_matrix : public mapped_object{
          friend class map_functor;
          std::string value_name_;
        public:
          mapped_implicit_matrix(std::string const & scalartype) : mapped_object(scalartype){ }
          std::string generate_default(std::pair<std::string, std::string> const & /* index */) const{
            return value_name_;
          }
          std::string & append_kernel_arguments(std::set<std::string> & /*already generated*/, std::string & str, unsigned int /*vector size*/) const{
            if(!value_name_.empty())
              str += detail::generate_value_kernel_argument(scalartype_, value_name_);
            return str;
          }
      };

      inline std::string generate(std::pair<std::string, std::string> const & index, int vector_element, mapped_object const & s){
        return s.generate(index, vector_element);
      }

      static void fetch(std::pair<std::string, std::string> const & index, unsigned int vectorization, std::set<std::string> & fetched, utils::kernel_generation_stream & stream, mapped_object & s){
        if(mapped_handle * p = dynamic_cast<mapped_handle  *>(&s))
          p->fetch(index, vectorization, fetched, stream);
      }

      static std::string & append_kernel_arguments(std::set<std::string> & already_generated, std::string & str, unsigned int vector_size, mapped_object const & s){
        return s.append_kernel_arguments(already_generated, str, vector_size);
      }

    }

  }

}
#endif
