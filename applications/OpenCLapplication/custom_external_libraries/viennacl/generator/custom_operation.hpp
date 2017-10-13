#ifndef VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP
#define VIENNACL_GENERATOR_CUSTOM_OPERATION_HPP

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

/** @file custom_operation.hpp
 *  @brief User Interface for making custom operations.
 *
 *  Generator code contributed by Philippe Tillet
 */


#include <vector>
#include <set>
#include <algorithm>

#include "viennacl/generator/get_kernels_infos.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/generator/traits/result_of.hpp"
#include "viennacl/generator/meta_tools/utils.hpp"


namespace viennacl 
{
  namespace generator 
  {

    /** @brief A class for making a custom operation */
    class custom_operation 
    { 

      public :

        /** @brief CTor
        * 	
        * @param expression the expression to build the interface for
        * @param program_name_hint the code for this expression will be stored in the program provided by this name
        */
        template<class T>
        custom_operation ( T const & expression, std::string const & program_name_hint="" ) 
        {
          program_name_ = viennacl::generator::program_infos<T>::value (program_name_hint, sources_,runtime_wrappers_);
          create_program ( static_cast<bool> ( viennacl::generator::tree_utils::count_if<T,viennacl::generator::is_inner_product_leaf>::value ) );
        }
        
        /** @brief DTor */
        ~custom_operation()
        {
          for (viennacl::generator::runtime_wrappers_t::iterator it  = runtime_wrappers_.begin();
                                                                 it != runtime_wrappers_.end();
                                                               ++it)
          {
            delete (it->second.second);
          }
        }

        /** @brief Returns the list of the kernels involved in the operation */
        viennacl::generator::KernelsSources const & kernels_sources() const 
        {
          return sources_;
        }

        /** @brief Return the generated sources */
        std::string kernels_source_code() const 
        {
          std::string res;
          for (viennacl::generator::KernelsSources::const_iterator it  = sources_.begin();
                                                                   it != sources_.end();
                                                                 ++it)
          {
            res += it->second + "\n";
          }
          
          return res;
        }

        /** @brief Returns the program name */
        std::string const & program_name() const  { return program_name_; }

        
        /** @brief Convenience for enqueuing the custom operation */
        template<class T0>
        custom_operation & operator() ( T0 const & t0) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1>
        custom_operation & operator() ( T0 const & t0, T1 const & t1 ) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1, class T2>
        custom_operation & operator() ( T0 const & t0, T1 const & t1, T2 const & t2 ) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1, class T2, class T3>
        custom_operation & operator() ( T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3 ) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          user_args_.insert( std::make_pair(3, viennacl::any((T3*)&t3)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1, class T2, class T3, class T4>
        custom_operation & operator() ( T0 & t0, T1 & t1, T2 & t2, T3 & t3, T4 & t4 ) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          user_args_.insert( std::make_pair(3, viennacl::any((T3*)&t3)) );
          user_args_.insert( std::make_pair(4, viennacl::any((T4*)&t4)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1, class T2, class T3, class T4, class T5>
        custom_operation & operator() ( T0 & t0, T1 & t1, T2 & t2, T3 & t3, T4 & t4, T5 & t5 ) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          user_args_.insert( std::make_pair(3, viennacl::any((T3*)&t3)) );
          user_args_.insert( std::make_pair(4, viennacl::any((T4*)&t4)) );
          user_args_.insert( std::make_pair(5, viennacl::any((T5*)&t5)) );
          add_operation_arguments();
          return *this;
        }

        /** @brief Convenience for enqueuing the custom operation */
        template<class T0, class T1, class T2, class T3, class T4, class T5, class T6>
        custom_operation & operator() ( T0 & t0, T1 & t1, T2 & t2, T3 & t3, T4 & t4, T5 & t5, T6 & t6) 
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          user_args_.insert( std::make_pair(3, viennacl::any((T3*)&t3)) );
          user_args_.insert( std::make_pair(4, viennacl::any((T4*)&t4)) );
          user_args_.insert( std::make_pair(5, viennacl::any((T5*)&t5)) );
          user_args_.insert( std::make_pair(6, viennacl::any((T6*)&t6)) );
          add_operation_arguments();
          return *this;
        }
        
        /** @brief Convenience for enqueuing the custom operation */
        template <class T0, class T1, class T2, class T3, class T4, class T5, class T6, class T7>
        custom_operation & operator() ( T0 & t0, T1 & t1, T2 & t2, T3 & t3, T4 & t4, T5 & t5, T6 & t6, T7 & t7 )
        {
          user_args_.insert( std::make_pair(0, viennacl::any((T0*)&t0)) );
          user_args_.insert( std::make_pair(1, viennacl::any((T1*)&t1)) );
          user_args_.insert( std::make_pair(2, viennacl::any((T2*)&t2)) );
          user_args_.insert( std::make_pair(3, viennacl::any((T3*)&t3)) );
          user_args_.insert( std::make_pair(4, viennacl::any((T4*)&t4)) );
          user_args_.insert( std::make_pair(5, viennacl::any((T5*)&t5)) );
          user_args_.insert( std::make_pair(6, viennacl::any((T6*)&t6)) );
          user_args_.insert( std::make_pair(7, viennacl::any((T7*)&t7)) );
          add_operation_arguments();
          return *this;
        }

      private:

        void create_program ( bool include_sum_kernel )
        {
          std::string kernels_string;
          for (viennacl::generator::KernelsSources::iterator it  = sources_.begin(); 
                                                             it != sources_.end(); 
                                                           ++it ) 
          {
            kernels_string += it->second + "\n";
          }
          
          viennacl::ocl::program& program = viennacl::ocl::current_context().add_program(kernels_string, program_name_);
          
          for (viennacl::generator::KernelsSources::iterator it  = sources_.begin(); 
                                                             it != sources_.end(); 
                                                           ++it) 
          {
            program.add_kernel(it->first);
          }
        }

        
        void add_operation_arguments() 
        {
          for (generator::runtime_wrappers_t::iterator it  = runtime_wrappers_.begin();
                                                       it != runtime_wrappers_.end();
                                                     ++it) 
          {
            std::string const & kernel_name = it->first;
            viennacl::ocl::kernel& current_kernel = viennacl::ocl::current_context().get_program(program_name_).get_kernel(kernel_name);
            const unsigned int arg_pos = it->second.first;
            generator::result_of::runtime_wrapper * current_arg = it->second.second;
        #ifdef VIENNACL_DEBUG_CUSTOM_OPERATION
            std::cout << "Enqueuing : Kernel " << kernel_name << " Argument : " << current_arg->name() << " | Pos : " << arg_pos << std::endl;
        #endif
            current_arg->enqueue(arg_pos,current_kernel,user_args_,temporaries_);
          }
        }

    private :
        typedef std::map<std::string, unsigned int> CurrentArgsContainer;
        CurrentArgsContainer current_args_pos_;
            
        std::map<unsigned int, viennacl::any> user_args_;
        
        std::string program_name_;
        
        std::vector<viennacl::ocl::local_mem> lmem_;
        
        viennacl::generator::KernelsSources sources_;
        
        viennacl::generator::runtime_wrappers_t runtime_wrappers_;
            
        std::map<std::string, viennacl::ocl::handle<cl_mem> > temporaries_;
    };


    inline void enqueue_custom_op(viennacl::generator::custom_operation & op, viennacl::ocl::command_queue const & queue)
    {
      for(std::map<std::string,std::string>::const_iterator it = op.kernels_sources().begin(); it != op.kernels_sources().end() ; ++it)
      {
        std::string current_kernel_name = it->first;
        #ifdef VIENNACL_DEBUG_CUSTOM_OPERATION
        std::cout << "Enqueueing " << current_kernel_name << std::endl;
        #endif
        enqueue(viennacl::ocl::current_context().get_program(op.program_name()).get_kernel(current_kernel_name));
      }
    }
    
  }
}

#endif
