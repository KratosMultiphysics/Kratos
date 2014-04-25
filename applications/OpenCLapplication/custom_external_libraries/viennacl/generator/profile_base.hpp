#ifndef VIENNACL_GENERATOR_GENERATE_TEMPLATE_BASE_BASE
#define VIENNACL_GENERATOR_GENERATE_TEMPLATE_BASE_BASE

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


/** @file viennacl/generator/profile_base.hpp
 *
 * @brief Base classes for the profiles
*/

#include <list>
#include <set>

#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/device_utils.hpp"
#include "viennacl/ocl/infos.hpp"

#include "viennacl/scheduler/forwards.h"

#include "viennacl/generator/helpers.hpp"
#include "viennacl/generator/map_functor.hpp"

namespace viennacl{

  namespace generator{


    /** @brief Base class for an operation profile */
    class profile_base{
      public:
        typedef std::list< std::pair<scheduler::statement, scheduler::statement_node> > statements_type;

      protected:
        friend std::ostream & operator<<(std::ostream &, profile_base const &);

        virtual bool invalid_impl(viennacl::ocl::device const & /*dev*/, vcl_size_t /*scalartype_size*/) const { return false; }
        virtual bool is_slow_impl(viennacl::ocl::device const &) const { return false; }

        virtual vcl_size_t lmem_used(vcl_size_t /*scalartype_size*/) const { return 0; }

        void configure_local_sizes(viennacl::ocl::kernel & k, vcl_size_t /*kernel_id*/) const {
          k.local_work_size(0,local_size_1_);
          k.local_work_size(1,local_size_2_);
        }

        virtual void print(std::ostream & s) const{
          s << csv_representation();
        }

        /** @brief Generates the body of the associated kernel function
         *
         *  @param kernel_id  If this profile requires multiple kernel, the index for which the core should be generated
         *  @param stream     The output stream the kernel is written to
         *  @param statements the statements for which the code should be generated
         *  @param mapping    the mapping of the statement_nodes to the mapped_objects
         */
        virtual void core(vcl_size_t kernel_id, utils::kernel_generation_stream& stream, statements_type const & statements, std::vector<detail::mapping_type> const & mapping) const = 0;

      public:
        /** @brief The constructor */
        profile_base(unsigned int vectorization, vcl_size_t local_size_1, vcl_size_t local_size_2, vcl_size_t num_kernels) : vector_size_(vectorization), local_size_1_(local_size_1), local_size_2_(local_size_2), num_kernels_(num_kernels){ }

        /** @brief The destructor */
        virtual ~profile_base(){ }

        /** @brief Configures the range and enqueues the arguments associated with the profile */
        virtual void configure_range_enqueue_arguments(vcl_size_t kernel_id, statements_type  const & statements, viennacl::ocl::kernel & k, unsigned int & n_arg) const = 0;

        virtual void kernel_arguments(statements_type  const & statements, std::string & arguments_string) const = 0;

        /** @brief Get the vector size of the kernel */
        unsigned int vector_size() const { return vector_size_; }

        /** @brief csv representation of an operation
         *
         *  Useful when writing to a file */
        virtual std::string csv_representation() const = 0;

        /** @brief returns whether or not the profile is likely to be slow on a particular device
         *  @param dev the given device*/
        bool is_slow(viennacl::ocl::device const & dev) const{
          bool res = false;
          if(dev.type()==CL_DEVICE_TYPE_GPU){
            vcl_size_t warp_size = 32;
            if(dev.vendor_id()==4098)
              warp_size = 64;
            res = static_cast<bool>(((local_size_1_*local_size_2_)%warp_size)>0);
          }
          return res || is_slow_impl(dev);
        }

        /** @brief returns whether or not the profile leads to undefined behavior on particular device
         *  @param dev               the given device
         *  @param scalartype_size   Local memory required to execute the kernel
         */
        bool is_invalid(viennacl::ocl::device const & dev, vcl_size_t scalartype_size) const{
          //Query device informations
          vcl_size_t lmem_available = static_cast<vcl_size_t>(dev.local_mem_size());
          vcl_size_t max_workgroup_size = dev.max_work_group_size();

          std::vector<vcl_size_t> max_work_item_sizes = dev.max_work_item_sizes();
          bool invalid_work_group_sizes = local_size_1_*local_size_2_ > max_workgroup_size
              || local_size_1_ > max_work_item_sizes[0]
              || local_size_2_ > max_work_item_sizes[1]; // uses too much resources

          return  invalid_work_group_sizes
              || lmem_used(scalartype_size)>lmem_available
              || invalid_impl(dev, scalartype_size);
        }

        /** @brief Returns the number of kernels needed by this operation */
        vcl_size_t num_kernels() const{ return num_kernels_; }

        /** @brief Generates the code associated with this profile onto the provided stream
         *  Redirects to the virtual core() method
         *
         *  @param stream Stream onto which the code should be generated
         *  @param device_offset the index of the device in the context (used for the kernel name)
         *  @param statements the statements associated with this profile */
        virtual void operator()(utils::kernel_generation_stream & stream, vcl_size_t device_offset, statements_type const & statements) const {
          std::vector<detail::mapping_type> mapping(statements.size());

          ///Get Prototype, initialize mapping
          std::string prototype;
          std::set<std::string> already_generated;
          kernel_arguments(statements, prototype);

          {
            std::map<void *, vcl_size_t> memory;
            unsigned int current_arg = 0;
            vcl_size_t i = 0;
            for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it)
              detail::traverse(it->first, it->second, detail::map_functor(memory,current_arg,mapping[i++]));
          }

          for(statements_type::const_iterator it = statements.begin() ; it != statements.end() ; ++it){
            detail::traverse(it->first, it->second, detail::prototype_generation_traversal(already_generated, prototype, vector_size(), mapping[std::distance(statements.begin(), it)]));
          }

          prototype.erase(prototype.size()-1); //Last comma pruned

          //Generate
          for(vcl_size_t n = 0 ; n < num_kernels() ; ++n){
            //stream << "__attribute__((vec_type_hint()))" << std::endl;
            stream << " __attribute__((reqd_work_group_size(" << local_size_1_ << "," << local_size_2_ << "," << 1 << ")))" << std::endl;
            stream << "__kernel " << "void " << "kernel_" << device_offset << "_" << n << "(" << std::endl;
            stream << prototype << std::endl;
            stream << ")" << std::endl;

            //core:
            stream << "{" << std::endl;
            stream.inc_tab();
            core(n, stream, statements, mapping);
            stream.dec_tab();
            stream << "}" << std::endl;
          }
        }

      protected:
        unsigned int vector_size_;
        vcl_size_t local_size_1_;
        vcl_size_t local_size_2_;
        vcl_size_t num_kernels_;
    };


    inline std::ostream & operator<<(std::ostream & os, profile_base const & profile){
      profile.print(os);
      return os;
    }

  }

}

#endif
