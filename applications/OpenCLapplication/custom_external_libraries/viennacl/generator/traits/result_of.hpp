#ifndef VIENNACL_GENERATOR_TRAITS_RESULT_OF_HPP
#define VIENNACL_GENERATOR_TRAITS_RESULT_OF_HPP

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

/** @file viennacl/generator/traits/result_of.hpp
 *  @brief Provides a set of metafunctions for type deductions within the kernel generator framework.
 *
 *  Generator code contributed by Philippe Tillet
 */

#include <string>

#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/forwards.h"
#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/elementwise_modifier.hpp"
#include "viennacl/ocl/local_mem.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/forwards.h"
#include "CL/cl.h"

namespace viennacl
{
  namespace generator
  {
    namespace result_of 
    {

      class runtime_wrapper
      {
        protected:
          bool is_temporary_;
          std::string name_;
          int arg_id_;
          
        public:
          runtime_wrapper(bool _is_temporary, std::string const & _name, int _arg_id) 
            : is_temporary_(_is_temporary), name_(_name), arg_id_(_arg_id) {}
          virtual ~runtime_wrapper() {}
            
          bool is_temporary() const { return is_temporary_; }
          int arg_id() const { return arg_id_; }
          std::string name() const { return name_; }

          virtual void enqueue(unsigned int arg_pos, 
                               viennacl::ocl::kernel & k,
                               std::map<unsigned int, viennacl::any> & runtime_args,
                               std::map<std::string, viennacl::ocl::handle<cl_mem> > & temporaries) = 0;
      };

      class shared_memory_wrapper : public runtime_wrapper
      {
        public:
          shared_memory_wrapper() : runtime_wrapper(true, "shared_memory_ptr", -1 ){ }
  
          void enqueue(unsigned int arg_pos,
                       viennacl::ocl::kernel & k,
                       std::map<unsigned int, viennacl::any> & runtime_args,
                       std::map<std::string, viennacl::ocl::handle<cl_mem> > & temporaries)
          {
            unsigned int lmem_size = k.local_work_size();
            k.arg(arg_pos, viennacl::ocl::local_mem(lmem_size*sizeof(float)));
          }
      
      };

      template <class T, class SIZE_T>
      struct vector_runtime_wrapper : public runtime_wrapper 
      {
        private:
          unsigned int size_id_;
          
          template<typename ScalarType, unsigned int Alignment>
          typename SIZE_T::size_type size(viennacl::vector<ScalarType,Alignment> * size_arg) { return size_arg->size(); }

          template<typename ScalarType, class F, unsigned int Alignment>
          typename SIZE_T::size_type size(viennacl::matrix<ScalarType,F,Alignment> * size_arg) { return size_arg->size2(); }
          
          template<typename ScalarType, unsigned int Alignment>
          typename SIZE_T::size_type internal_size(viennacl::vector<ScalarType,Alignment> * size_arg) { return size_arg->internal_size(); }

          template<typename ScalarType, class F, unsigned int Alignment>
          typename SIZE_T::size_type internal_size(viennacl::matrix<ScalarType,F,Alignment> * size_arg) { return size_arg->internal_size2(); }
          
        public:
          vector_runtime_wrapper(bool _is_temporary, std::string const & _name, int _arg_id, unsigned int _size_id) 
            : runtime_wrapper(_is_temporary,_name,_arg_id),size_id_(_size_id) {}
            
          void enqueue(unsigned int arg_pos,
                       viennacl::ocl::kernel & k,
                       std::map<unsigned int, viennacl::any> & runtime_args,
                       std::map<std::string, 
                       viennacl::ocl::handle<cl_mem> > & temporaries)
          { 
            SIZE_T * size_arg = viennacl::any_cast<SIZE_T * >(runtime_args[size_id_]);
            viennacl::ocl::handle<cl_mem> handle = NULL;
            if(is_temporary_)
            {
              if(temporaries.find(name_)==temporaries.end())
              {
                temporaries.insert(
                  std::make_pair(name_,
                                 viennacl::ocl::handle<cl_mem>(
                                 viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                size_arg->internal_size()*sizeof(typename T::value_type))
                                                                                )
                                )
                                  );
              }
              handle = temporaries[name_];
            }
            else
            {
              T * current_arg = viennacl::any_cast<T * >(runtime_args[arg_id_]);
              handle = current_arg->handle();
            }
            k.arg(arg_pos, handle );
            k.arg(arg_pos+1,cl_uint(size(size_arg)));
            k.arg(arg_pos+2,cl_uint(internal_size(size_arg)));
          }
      };

      template <class T, class SIZE_DESCRIPTOR>
      struct vector_expression : public runtime_wrapper
      {
        typedef T type;
        typedef typename SIZE_DESCRIPTOR::ScalarType ScalarType;
        static const unsigned int Alignment = SIZE_DESCRIPTOR::Alignment;

        static runtime_wrapper * runtime_descriptor()
        {
          return new vector_runtime_wrapper<viennacl::vector<ScalarType,Alignment>,
                                            typename SIZE_DESCRIPTOR::runtime_type>(viennacl::generator::is_temporary<T>::value,
                                                                                    T::name(),
                                                                                    T::id,SIZE_DESCRIPTOR::id);
        }
        
        static const std::string size_expression() 
        {
          return SIZE_DESCRIPTOR::size2_name();
        }
        
        static const std::string internal_size_expression() 
        {
          return SIZE_DESCRIPTOR::internal_size2_name() + "/" + to_string(Alignment);
        }
      };

      template <class T, class SIZE1_T, class SIZE2_T>
      struct matrix_runtime_wrapper : public runtime_wrapper
      {
        private:
          unsigned int size1_id_;
          unsigned int size2_id_;
        public:
          matrix_runtime_wrapper(bool _is_temporary, 
                                 std::string const & _name,
                                 int _arg_id,
                                 unsigned int _size1_id,
                                 unsigned int _size2_id) 
                                : runtime_wrapper(_is_temporary,_name,_arg_id), size1_id_(_size1_id), size2_id_(_size2_id) {}
                                
          unsigned int n_elements(){ return size1_id_*size2_id_; }
          
          void enqueue(unsigned int arg_pos,
                       viennacl::ocl::kernel & k,
                       std::map<unsigned int, viennacl::any> & runtime_args,
                       std::map<std::string,
                       viennacl::ocl::handle<cl_mem> > & temporaries)
          { 
            if (is_temporary_) {}
            
            T * current_arg = any_cast<T * >(runtime_args[arg_id_]);
            SIZE1_T * size1_arg = any_cast<SIZE1_T * >(runtime_args[size1_id_]);
            SIZE2_T * size2_arg = any_cast<SIZE2_T * >(runtime_args[size2_id_]);
            k.arg(arg_pos, current_arg->handle());
            k.arg(arg_pos+1,cl_uint(size1_arg->size1()));
            k.arg(arg_pos+2,cl_uint(size2_arg->size2()));
            k.arg(arg_pos+3,cl_uint(size1_arg->internal_size1()));
            k.arg(arg_pos+4,cl_uint(size2_arg->internal_size2()));
          }
      };
          
      template <class T, class SIZE1_DESCRIPTOR, class SIZE2_DESCRIPTOR>
      struct matrix_expression 
      {
        typedef typename SIZE1_DESCRIPTOR::ScalarType ScalarType;
        typedef typename SIZE1_DESCRIPTOR::Layout Layout;
        static const unsigned int Alignment = SIZE1_DESCRIPTOR::Alignment;
        
        static runtime_wrapper * runtime_descriptor()
        {
          return new matrix_runtime_wrapper<viennacl::matrix<ScalarType,Layout,Alignment>,
                                            typename SIZE1_DESCRIPTOR::runtime_type,
                                            typename SIZE2_DESCRIPTOR::runtime_type>(is_temporary<T>::value,T::name(),
                                                                                     T::id,SIZE1_DESCRIPTOR::id,
                                                                                     SIZE2_DESCRIPTOR::id);
        }
        
        static const std::string size1_expression() 
        {
          return SIZE1_DESCRIPTOR::size1_name();
        }

        static const std::string size2_expression() 
        {
          return SIZE2_DESCRIPTOR::size2_name();
        }

        static const std::string internal_size1_expression() 
        {
          return SIZE1_DESCRIPTOR::internal_size1_name() + "/" + to_string(Alignment);
        }

        static const std::string internal_size2_expression() 
        {
          return SIZE2_DESCRIPTOR::internal_size2_name() + "/" + to_string(Alignment);
        }

        typedef T type;
      };

      template <class T>
      struct scalar_size_descriptor
      {
	      static unsigned int size(viennacl::ocl::kernel & k) { return 1; }
      };

      template <class LHS, class RHS, bool is_temporary>
      struct scalar_size_descriptor<compound_node<LHS,inner_prod_type,RHS,is_temporary> >
      {
        static unsigned int size(viennacl::ocl::kernel & k)
        {
          return k.global_work_size(0)/k.local_work_size(0);
        }
      };

      template <class T>
      struct scalar_runtime_wrapper: public runtime_wrapper
      {
        typedef typename T::ScalarType ScalarType;
        
        scalar_runtime_wrapper(bool _is_temporary, std::string const & _name, int _arg_id) : runtime_wrapper(_is_temporary,_name,_arg_id){}
        
        void enqueue(unsigned int arg_pos,
                     viennacl::ocl::kernel & k,
                     std::map<unsigned int,
                     viennacl::any> & runtime_args, 
                     std::map<std::string, 
                     viennacl::ocl::handle<cl_mem> > & temporaries)
        {
          if(is_temporary_)
          {
            if(temporaries.find(name_)==temporaries.end()) 
            {
              temporaries.insert(
                        std::make_pair(name_,
                viennacl::ocl::handle<cl_mem>(
                viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                scalar_size_descriptor<T>::size(k)*sizeof(ScalarType))
                                              )
                )
              );
            }
            k.arg(arg_pos, temporaries[name_]);
          }
          
          if(arg_id_==-2)
                  k.arg(arg_pos, temporaries[name_]);
          else
          {
            viennacl::scalar<ScalarType>* current_arg = any_cast<viennacl::scalar<ScalarType> * >(runtime_args[arg_id_]);
            k.arg(arg_pos, current_arg->handle());
          }
    
        }
      };

      template <unsigned int ID, class ScalarType>
      struct scalar_runtime_wrapper<viennacl::generator::cpu_symbolic_scalar<ID, ScalarType> >: public runtime_wrapper
      {
        scalar_runtime_wrapper(bool _is_temporary, std::string const & _name, int _arg_id) : runtime_wrapper(_is_temporary,_name,_arg_id){ }
        
        void enqueue(unsigned int arg_pos,
                     viennacl::ocl::kernel & k,
                     std::map<unsigned int, viennacl::any> & runtime_args,
                     std::map<std::string, viennacl::ocl::handle<cl_mem> > & temporaries)
        {
          ScalarType* current_arg = any_cast<ScalarType * >(runtime_args[arg_id_]);
          k.arg(arg_pos, cl_float(*current_arg));
        }
      };
          
      template <class T>
      struct scalar_expression 
      {
        typedef typename T::ScalarType ScalarType;
        
        static runtime_wrapper * runtime_descriptor()
        {
          return new scalar_runtime_wrapper<T>(is_temporary<T>::value,T::name(),T::id);
        }
      };

      /*
       * Compound Nodes - General case
       */
      template <class T>
      struct expression_type 
      {
        typedef NullType Result;
      };

      template <class LHS, class OP, class RHS, bool is_temporary>
      struct expression_type<compound_node<LHS,OP,RHS,is_temporary> > 
      {
        private:
          typedef typename expression_type<LHS>::Result LHS_Result;
          typedef typename expression_type<RHS>::Result RHS_Result;
          
        public:
          typedef typename expression_type<compound_node<LHS_Result, OP, RHS_Result,is_temporary> >::Result Result;
      };

      /*
       * Compound Nodes - usual operators
       */
      template <class LHS, class LHS_SIZE_DESCRIPTOR ,class OP ,class RHS, class RHS_SIZE_DESCRIPTOR ,bool is_temporary>
      struct expression_type<compound_node<vector_expression<LHS,LHS_SIZE_DESCRIPTOR>,
                                           OP,
                                           vector_expression<RHS,RHS_SIZE_DESCRIPTOR>,
                                           is_temporary>
                            >
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
          
        public:
          typedef vector_expression<T, LHS_SIZE_DESCRIPTOR> Result;
      };


      template <class LHS, class LHS_SIZE1_DESCRIPTOR, class LHS_SIZE2_DESCRIPTOR,
                class OP,
                class RHS, class RHS_SIZE1_DESCRIPTOR, class RHS_SIZE2_DESCRIPTOR,
                bool is_temporary>
      struct expression_type<compound_node<matrix_expression<LHS, LHS_SIZE1_DESCRIPTOR, LHS_SIZE2_DESCRIPTOR>,
                                           OP,
                                           matrix_expression<RHS, RHS_SIZE1_DESCRIPTOR, RHS_SIZE2_DESCRIPTOR>,
                                           is_temporary> 
                             > 
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
          
        public:
          typedef matrix_expression<T, LHS_SIZE1_DESCRIPTOR, LHS_SIZE2_DESCRIPTOR> Result;
      };

      template <class LHS, class OP, class RHS, bool is_temporary>
      struct expression_type<compound_node<scalar_expression<LHS>, 
                                           OP,
                                           scalar_expression<RHS>,
                                           is_temporary> 
                            > 
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
          
        public:
          typedef scalar_expression<T> Result;
      };

      /*
       * Scalar Operators
       */
      template <class LHS, class LHS_SIZE_DESCRIPTOR,
                class OP,
                class RHS,
                bool is_temporary>
      struct  expression_type<compound_node<vector_expression<LHS,LHS_SIZE_DESCRIPTOR>,
                                            OP,
                                            scalar_expression<RHS>,
                                            is_temporary> > 
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
        public:
          typedef vector_expression<T, LHS_SIZE_DESCRIPTOR> Result;
      };

      template <class LHS,
                class OP,
                class RHS, class RHS_SIZE_DESCRIPTOR,
                bool is_temporary>
      struct expression_type<compound_node<scalar_expression<LHS>,
                                           OP,
                                           vector_expression<RHS,RHS_SIZE_DESCRIPTOR>,
                                           is_temporary>
                            > 
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
        public:
          typedef vector_expression<T, RHS_SIZE_DESCRIPTOR> Result;
      };


      template <class LHS, class LHS_SIZE1_DESCRIPTOR, class LHS_SIZE2_DESCRIPTOR,
                class OP,
                class RHS, bool is_temporary>
      struct expression_type<compound_node<matrix_expression<LHS,LHS_SIZE1_DESCRIPTOR,LHS_SIZE2_DESCRIPTOR>,
                                           OP,
                                           scalar_expression<RHS>,
                                           is_temporary>
                            > 
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
        public:
          typedef matrix_expression<T, LHS_SIZE1_DESCRIPTOR, LHS_SIZE2_DESCRIPTOR> Result;
      };

      template <class LHS, 
                class OP,
                class RHS, class RHS_SIZE1_DESCRIPTOR, class RHS_SIZE2_DESCRIPTOR,
                bool is_temporary>
      struct expression_type<compound_node<scalar_expression<LHS>,
                                           OP,
                                           matrix_expression<RHS,RHS_SIZE1_DESCRIPTOR, RHS_SIZE2_DESCRIPTOR>,
                                           is_temporary>
                            >
      {
        private:
          typedef compound_node<LHS ,OP, RHS, is_temporary> T;
        public:
          typedef matrix_expression<T, RHS_SIZE1_DESCRIPTOR, RHS_SIZE2_DESCRIPTOR> Result;
      };


      /*
       * Compound Nodes - Non Trivial Operators
       */

      //Matrix-Vector product
      template <class LHS, class LHS_SIZE1_DESCRIPTOR, class LHS_SIZE2_DESCRIPTOR,
                class RHS, class RHS_SIZE_DESCRIPTOR,
                bool is_temporary>
      struct expression_type<compound_node<matrix_expression<LHS, LHS_SIZE1_DESCRIPTOR, LHS_SIZE2_DESCRIPTOR>,
                                           prod_type,
                                           vector_expression<RHS,RHS_SIZE_DESCRIPTOR>,
                                           is_temporary> 
                            >
      {
        typedef vector_expression<compound_node<LHS,prod_type,RHS,is_temporary>, LHS_SIZE1_DESCRIPTOR > Result;
      };

      template <class T>
      struct expression_type<inner_prod_impl_t<T> >
      {
	      typedef scalar_expression<T> Result;
      };

      //Matrix-Matrix product
      template <class LHS, class LHS_SIZE1_DESCRIPTOR, class LHS_SIZE2_DESCRIPTOR,
                class RHS, class RHS_SIZE1_DESCRIPTOR, class RHS_SIZE2_DESCRIPTOR,
                bool is_temporary>
      struct expression_type<compound_node<matrix_expression<LHS, LHS_SIZE1_DESCRIPTOR, LHS_SIZE2_DESCRIPTOR>,
                                           prod_type,
                                           matrix_expression<RHS,RHS_SIZE1_DESCRIPTOR,RHS_SIZE2_DESCRIPTOR>,
                                           is_temporary> 
                            >
      {
        typedef matrix_expression<compound_node<LHS,prod_type,RHS,is_temporary>, LHS_SIZE1_DESCRIPTOR, RHS_SIZE2_DESCRIPTOR > Result;
      };

      //Inner product
      template <class LHS, class LHS_SIZE_DESCRIPTOR,
                class RHS, class RHS_SIZE_DESCRIPTOR,
                bool is_temporary>
      struct expression_type< compound_node<vector_expression<LHS,LHS_SIZE_DESCRIPTOR>,
                                            inner_prod_type,
                                            vector_expression<RHS,RHS_SIZE_DESCRIPTOR>,
                                            is_temporary>
                            >
      {
        typedef scalar_expression<compound_node<LHS,inner_prod_type,RHS,is_temporary> > Result;
      };


      /*
       * Elementwise Modifiers
       */
      template <class T, std::string (*U)()>
      struct expression_type< elementwise_modifier_impl<T,U> > 
      {
        typedef typename expression_type<T>::Result Result;
      };

      template <class T, class SIZE_DESCRIPTOR>
      struct expression_type< vector_expression<T,SIZE_DESCRIPTOR> > 
      {
        typedef typename expression_type<T>::Result Result;
      };

      template <class T, class SIZE1_DESCRIPTOR, class SIZE2_DESCRIPTOR>
      struct expression_type< matrix_expression<T,SIZE1_DESCRIPTOR,SIZE2_DESCRIPTOR> > 
      {
        typedef typename expression_type<T>::Result Result;
      };

      template <class T>
      struct expression_type< scalar_expression<T> > 
      {
        typedef typename expression_type<T>::Result Result;
      };

      /*
       * Symbolic Vectors
       */

      template <unsigned int ID,typename SCALARTYPE, unsigned int ALIGNMENT>
      struct expression_type< symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > 
      {
        typedef vector_expression<symbolic_vector<ID,SCALARTYPE,ALIGNMENT>,
                                  symbolic_vector<ID,SCALARTYPE,ALIGNMENT> > Result;
      };

      template <class Ref>
      struct  expression_type<tmp_symbolic_vector<Ref> > 
      {
        typedef vector_expression<tmp_symbolic_vector<Ref>, Ref> Result;
      };

      /*
       * Symbolic Matrices
       */

      template <unsigned int ID,typename SCALARTYPE, class F, unsigned int ALIGNMENT>
      struct expression_type<symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> > 
      {
        private:
          typedef symbolic_matrix<ID,SCALARTYPE,F,ALIGNMENT> T;
        public:
          typedef matrix_expression<T, T, T> Result;
      };

      template <class Ref>
      struct expression_type<tmp_symbolic_matrix<Ref> > 
      {
        typedef matrix_expression<tmp_symbolic_matrix<Ref>, Ref, Ref > Result;
      };

      /*
       * Symbolic Scalars
       */

      template <unsigned int ID, typename SCALARTYPE>
      struct expression_type<cpu_symbolic_scalar<ID, SCALARTYPE> > 
      {
        typedef scalar_expression<cpu_symbolic_scalar<ID, SCALARTYPE> > Result;
      };

      template <unsigned int ID, typename SCALARTYPE>
      struct expression_type<gpu_symbolic_scalar<ID, SCALARTYPE> > 
      {
        typedef scalar_expression< gpu_symbolic_scalar<ID, SCALARTYPE> > Result;
      };


    }
  }
}

#endif


