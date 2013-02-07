#ifndef VIENNACL_GENERATOR_CREATE_KERNEL_HPP
#define VIENNACL_GENERATOR_CREATE_KERNEL_HPP

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

/** @file viennacl/generator/get_kernels_infos.hpp
 *  @brief Provides information about kernels
 * 
 *  Generator code contributed by Philippe Tillet
 */

// #include "kernel_utils.hpp"

#include <map>

#include "viennacl/generator/compound_node.hpp"
#include "viennacl/generator/operation_types.hpp"
#include "viennacl/generator/symbolic_types/symbolic_matrix.hpp"
#include "viennacl/generator/symbolic_types/symbolic_vector.hpp"
#include "viennacl/generator/symbolic_types/symbolic_scalars.hpp"
#include "viennacl/generator/tree_operations.hpp"
#include "viennacl/generator/tokens_management.hpp"
#include "viennacl/generator/make_code/make_code.hpp"
#include "viennacl/generator/meta_tools/typelist.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"
#include "viennacl/generator/traits/result_of.hpp"

namespace viennacl 
{
  namespace generator
  {

    template <class T, bool is_first = true>
    struct arguments_list;

    template <bool is_first>
    struct arguments_list<NullType, is_first > 
    {
      static const std::string  string_value() { return ""; }
    };


    template <class Head, class Tail, bool is_first >
    struct arguments_list<typelist<Head,Tail>, is_first > 
    {
      private:
        static const std::string add_comma ( Int2Type<false> ) { return ", "; }
        static const std::string add_comma ( Int2Type<true> )  { return ""; }

      public:
        static const std::string string_value() 
        {
            return add_comma ( Int2Type<is_first>() )
                  + Head::kernel_arguments()
                  + arguments_list<Tail,false>::string_value();
        }
    };

    template<class T>
    struct requires_local_buffer 
    {
      enum { value = is_inner_product_leaf<T>::value };
    };

    template<class T>
    struct requires_local_buffer<inner_prod_impl_t<T> > 
    {
      enum { value = 1 };
    };

    template<class T>
    struct requires_local_buffer_list;

    template<class Head, class Tail>
    struct requires_local_buffer_list<typelist<Head, Tail> > 
    {
      enum { value = static_cast<bool> ( tree_utils::count_if<Head, requires_local_buffer>::value
                                        || requires_local_buffer_list<Tail>::value )
           };
    };

    template<class Head>
    struct requires_local_buffer_list<typelist<Head, NullType> > 
    {
      enum { value = static_cast<bool> ( tree_utils::count_if<Head, requires_local_buffer >::value ) };
    };

    template<class TLIST, class ASSIGNED>
    struct calculate_tokens 
    {
      typedef typename TLIST::Head::first_type current_token;
      typedef typename TLIST::Head::second_type token_op;

      static const std::string value() 
      {
        return make_code<current_token, token_op, ASSIGNED>::value()
              + calculate_tokens<typename TLIST::Tail, ASSIGNED>::value();
      }
    };

    template<class ASSIGNED>
    struct calculate_tokens<NullType, ASSIGNED> 
    {
      static const std::string value() { return ""; }
    };

    template<class T>
    struct get_temporary_dependancies
    {
      typedef typename get_temporary_dependancies<typename result_of::expression_type<T>::Result>::Result   Result;
    };

    template<class T, class SIZE_DESCRIPTOR>
    struct get_temporary_dependancies<result_of::vector_expression<T,SIZE_DESCRIPTOR> >
    {
      typedef SIZE_DESCRIPTOR Result;
    };

    template<class T>
    struct get_temporary_dependancies<result_of::scalar_expression<T> >
    {
      typedef NullType Result;
    };

    template<>
    struct get_temporary_dependancies<NullType>
    {
      typedef NullType Result;
    };

    template<class Head, class Tail>
    struct get_temporary_dependancies<typelist<Head,Tail> >
    {
      typedef typename typelist_utils::append<typename get_temporary_dependancies<Tail>::Result,
                                              typename get_temporary_dependancies<Head>::Result>::Result   Result;
    };


    template<class T>
    struct get_kernel_arguments;

    template<class LHS, class OP, class RHS, bool _is_temporary>
    struct get_kernel_arguments<compound_node<LHS,OP,RHS,_is_temporary> > 
    {
      typedef compound_node<LHS,OP,RHS,_is_temporary> Arg;
      typedef typename tree_utils::extract_if<Arg, is_regular_kernel_parameter, typelist_utils::compare1>::Result RegularLeafs;
      typedef typename tree_utils::extract_if<Arg, is_temporary_kernel_parameter>::Result                         TemporaryLeafs;
      typedef typename get_temporary_dependancies<TemporaryLeafs>::Result                                         TemporaryDependancies;
      typedef typename typelist_utils::fuse<RegularLeafs, TemporaryLeafs, typelist_utils::compare1>::Result       TmpResult0;
      typedef typename typelist_utils::fuse<TmpResult0, TemporaryDependancies, typelist_utils::compare1>::Result  TmpResult1;
      typedef typename typelist_utils::no_duplicates<TmpResult1>::Result Result;
    };

    template<class Head, class Tail>
    struct get_kernel_arguments<typelist<Head,Tail> > 
    {
      typedef typename typelist_utils::fuse<typename get_kernel_arguments<Head>::Result, 
                                            typename get_kernel_arguments<Tail>::Result,
                                            typelist_utils::compare1>::Result                TmpResult;
      typedef typename typelist_utils::no_duplicates<TmpResult>::Result                      Result;
    };

    template<>
    struct get_kernel_arguments<NullType> 
    {
      typedef NullType Result;
    };

    template<class TreeList>
    struct kernel_header;

    template<class Head, class Tail>
    struct kernel_header<typelist<Head, Tail> > 
    {
      private:
        typedef typelist<Head, Tail> Arg;
        typedef typename tree_utils::expand<Head>::Result             ExpandedHead;
        typedef typename tree_utils::flip_tree<ExpandedHead>::Result  NewHead;
        typedef typename get_kernel_arguments<Arg>::Result            Arguments;
        
        static const std::string shared_memory ( Int2Type<false> ) { return ""; }
        static const std::string shared_memory ( Int2Type<true> )  { return ",__local float* shared_memory_ptr\n"; }

      public:
        static const std::string value ( std::string const & name )
        {
          return "__kernel void " + name + "(\n"
                + arguments_list<Arguments>::string_value()
                + shared_memory ( Int2Type<requires_local_buffer_list<Arg>::value>() )
                + ")\n";
        }
    };


    template<class TreeList, bool is_in_temporary_kernel, bool is_first = true>
    struct kernel_core;

    template<class TList>
    struct finalize_inner_products 
    {
      static const std::string value() { return ""; }
    };

    template<class Head, class Tail>
    struct finalize_inner_products<typelist<Head,Tail> > 
    {
      static const std::string value() 
      {
          return make_code<Head,assign_type,Head>::value() + finalize_inner_products<Tail>::value();
      }
    };

    template<class Head, class Tail, bool is_in_temporary_kernel, bool is_first>
    struct kernel_core<typelist<Head, Tail>, is_in_temporary_kernel, is_first > 
    {
      private:
        typedef typelist<Head, Tail> Arg;
        typedef typename tree_utils::expand<Head>::Result ExpandedHead;
        typedef typename tree_utils::flip_tree<ExpandedHead>::Result NewHead;
        typedef typename generate_tokens<NewHead, is_in_temporary_kernel>::Result Tokens;
        typedef typename tree_utils::extract_if<typename NewHead::RHS,is_inner_product_leaf>::Result InProdsT;
        typedef typename typelist_utils::no_duplicates<InProdsT>::Result InProds;
        typedef typename Head::LHS LHS;

        static const std::string additional_declarations ( Int2Type<true> ) {
            return  "float sum;\n";
        }

        static const std::string additional_declarations ( Int2Type<false> ) { return  "" ; }

        
        static const std::string head ( Int2Type<true> ) {
            return  "{\n"
                    + additional_declarations ( Int2Type<requires_local_buffer_list<Arg>::value>() );
        }

        static const std::string head ( Int2Type<false> ) { return  "\n" ; }

      public:
        static const std::string value() 
        {
            return head ( Int2Type<is_first>() )
                  + finalize_inner_products<InProds>::value()
                  + calculate_tokens<Tokens, LHS>::value()
                  + kernel_core<Tail, is_in_temporary_kernel, false>::value();
        }
    };

    template<bool is_in_temporary_kernel>
    struct kernel_core<NullType, is_in_temporary_kernel, false> 
    {
      static const std::string value() { return "}"; }
    };

    template<class T>
    struct remove_temporary 
    {
      typedef T Result;
    };

    template<class Ref>
    struct remove_temporary<tmp_symbolic_vector<Ref> > 
    {
      typedef Ref Result;
    };

    template<class LHS, class OP, class RHS>
    struct remove_temporary<compound_node<LHS,OP,RHS,true> > 
    {
      typedef compound_node<LHS,OP,RHS> Result;
    };

    template<class TreeList, class Assigned>
    struct get_all_temporaries 
    {
      private:
        typedef typename TreeList::Head Head;
        typedef typename TreeList::Tail Tail;
        typedef typename remove_temporary<Head>::Result                                    NewHead;
        typedef typename tree_utils::register_temporaries<NewHead, true, Assigned>::Result Registered;
        typedef typename tree_utils::extract_if<Registered,is_temporary>::Result           Temporaries;
        typedef typename get_all_temporaries<Tail, Assigned>::Result                       NewList;
        
      public:
        typedef typename typelist_utils::fuse<Temporaries, NewList>::Result      Result;
    };

    template<class Assigned>
    struct get_all_temporaries<NullType, Assigned> 
    {
      typedef NullType Result;
    };

    template<class T>
    struct Unroll
    {
      typedef NullType Result;
    };

    template<class Head, class Tail>
    struct Unroll<typelist<Head, Tail> > 
    {
      typedef typename typelist_utils::fuse<Head,
                                            typename Unroll<Tail>::Result >::Result    Result;
    };

    template<class HeadHead, class HeadTail, class Tail>
    struct Unroll<typelist<typelist<HeadHead, HeadTail>, Tail> > 
    {
      typedef typename typelist_utils::fuse<typelist<HeadHead, HeadTail>,
                                            typename Unroll<Tail>::Result >::Result    Result;
    };

    template<>
    struct Unroll<typelist<NullType, NullType> > 
    {
      typedef NullType Result;
    };

    template<class T>
    struct find_prior_implementations 
    {
      typedef T Result;
    };

    template<class LHS, class RHS>
    struct find_prior_implementations<compound_node<LHS, inner_prod_type, RHS,true> > 
    {
      typedef inner_prod_impl_t<compound_node<LHS, inner_prod_type, RHS,true> > Result;
    };

    template<class Head, class Tail>
    struct find_prior_implementations<typelist<Head,Tail> > 
    {
      private:
        typedef typename find_prior_implementations<Head>::Result NewHead;
        typedef typename find_prior_implementations<Tail>::Result NewTail;
        
      public:
        typedef typelist<NewHead,NewTail> Result;
    };

    template<class TreeList, class Assigned>
    struct register_kernels 
    {
      private:
        typedef typename get_all_temporaries<TreeList, Assigned>::Result Temporaries;
        typedef typename register_kernels<Temporaries, Assigned>::Result CurrentList;
        typedef typename typelist_utils::erase<Temporaries,
                                               typename Unroll<CurrentList>::Result>::Result  NextTemporaries;
        typedef typename find_prior_implementations<NextTemporaries>::Result                  NextList;

      public:
        typedef typename typelist_utils::append<CurrentList, NextList>::Result Result;
    };

    template<class Assigned>
    struct register_kernels<NullType, Assigned> 
    {
      typedef typelist<NullType,NullType> Result;
    };

    template<class TreeList, bool is_in_temporary_kernel, class Enable = void>
    struct kernel_string 
    {
      static const std::string value(std::string name) 
      {
        return std::string ( kernel_header<TreeList>::value(name)
                             + kernel_core<TreeList, is_in_temporary_kernel>::value() );
      }
    };

    template<class T>
    struct make_impl;

    template<class LHS_, class RHS_, bool is_temporary_>
    struct make_impl<compound_node< LHS_, inner_prod_type, RHS_, is_temporary_ > > 
    {
      typedef inner_prod_impl_t<compound_node< LHS_, inner_prod_type, RHS_, is_temporary_ > >   Result;
    };

    template<class Temporary>
    struct format_temporaries;

    template<class Head, class Tail>
    struct format_temporaries<typelist<Head,Tail> > 
    {
        typedef compound_node<Head, assign_type, typename remove_temporary<Head>::Result>                    NewHead;
        typedef typename typelist_utils::append<typename format_temporaries<Tail>::Result, NewHead>::Result  Result;
    };

    template<>
    struct format_temporaries<NullType> 
    {
      typedef NullType Result;
    };

    typedef std::map<std::string,std::string> KernelsSources;

    template<class TemporaryKernelsList, class MainOperation, int Start, int End>
    struct fill_sources
    {
      typedef typename format_temporaries<typename typelist_utils::type_at< TemporaryKernelsList,Start>::Result >::Result CurrentList;
      
      static void execute(KernelsSources & sources, std::string const & operation_name)
      {
        std::string current_kernel_name("__" + operation_name + "_kernel"  
                                        + to_string(typelist_utils::length<TemporaryKernelsList>::value - 1 - Start));
        sources.insert( std::make_pair( current_kernel_name,
                                        kernel_string< CurrentList, true>::value(current_kernel_name) 
                                      )
                      );
        fill_sources<TemporaryKernelsList, MainOperation, Start+1, End>::execute(sources, operation_name);
      }
    };

    template<class TemporaryKernelsList, class MainOperation, int End>
    struct fill_sources<TemporaryKernelsList, MainOperation, End, End>
    {
      static void execute(KernelsSources & sources, std::string const & operation_name)
      {
        sources.insert(std::make_pair(operation_name,
                                      kernel_string<MainOperation, false>::value(operation_name) 
                                     )
                      );
      }
    };

    typedef std::multimap<std::string, std::pair<unsigned int, result_of::runtime_wrapper*> > runtime_wrappers_t;

    template<class U>
    struct foreach_functor
    {
      static void execute(unsigned int & arg_pos,  runtime_wrappers_t & runtime_wrappers, std::string const & name) 
      {
        foreach_functor<typename result_of::expression_type<U>::Result >::execute(arg_pos, runtime_wrappers, name);
      }
    };

    template<>
    struct foreach_functor<NullType>;

    template<class T>
    struct foreach_functor<result_of::scalar_expression<T> >
    {
      static void execute(unsigned int & arg_pos, runtime_wrappers_t & runtime_wrappers, std::string const & name) 
      {
        runtime_wrappers.insert(runtime_wrappers_t::value_type(name,
                                                               std::make_pair(arg_pos,
                                                                              result_of::scalar_expression<T>::runtime_descriptor())
                                                              )
                               );
        arg_pos += 1;
      }
    };

    template<class T, class SIZE_DESCRIPTOR>
    struct foreach_functor<result_of::vector_expression<T,SIZE_DESCRIPTOR> >
    {
      static void execute(unsigned int & arg_pos, runtime_wrappers_t & runtime_wrappers, std::string const & name) 
      {
        runtime_wrappers.insert(runtime_wrappers_t::value_type(name,
                                                               std::make_pair(arg_pos,
                                                                              result_of::vector_expression<T,SIZE_DESCRIPTOR>::runtime_descriptor())
                                                              )
                               );
        arg_pos += 3;
      }
    };

    template<class T, class SIZE1_DESCRIPTOR, class SIZE2_DESCRIPTOR>
    struct foreach_functor<result_of::matrix_expression<T,SIZE1_DESCRIPTOR, SIZE2_DESCRIPTOR> >
    {
      static void execute(unsigned int & arg_pos, runtime_wrappers_t & runtime_wrappers, std::string const & name ) 
      {
        runtime_wrappers.insert(runtime_wrappers_t::value_type(name,
                                                               std::make_pair(arg_pos,
                                                                              result_of::matrix_expression<T,SIZE1_DESCRIPTOR,SIZE2_DESCRIPTOR>::runtime_descriptor())
                                                              )
                               );
        arg_pos += 5;
      }
    };

    template<class TemporaryKernelsList, class MainOperation, int Start, int End>
    struct fill_args
    {
      typedef typename format_temporaries<typename typelist_utils::type_at<TemporaryKernelsList,Start>::Result >::Result CurrentList;
      typedef typename get_kernel_arguments<CurrentList>::Result Arguments;
      
      static void execute(runtime_wrappers_t & runtime_wrappers, std::string const & operation_name)
      {
        unsigned int arg_pos = 0;
        std::string current_kernel_name("__"+operation_name+"_kernel"
                                        + to_string(typelist_utils::length<TemporaryKernelsList>::value - 1 - Start));
        
        typelist_utils::ForEach<Arguments, foreach_functor>::execute(arg_pos,runtime_wrappers,current_kernel_name);
        
        if(requires_local_buffer_list<CurrentList>::value)
        {
          runtime_wrappers.insert(runtime_wrappers_t::value_type(current_kernel_name,
                                                                 std::make_pair(arg_pos,
                                                                                new result_of::shared_memory_wrapper())
                                                                )
                                 );
        }
        
        fill_args<TemporaryKernelsList,MainOperation,Start+1,End>::execute(runtime_wrappers,operation_name);
      }
    };

    template<class TemporaryKernelsList, class MainOperation, int End>
    struct fill_args<TemporaryKernelsList, MainOperation, End, End>
    {
      private:
        typedef MainOperation CurrentList;
        typedef typename get_kernel_arguments<CurrentList>::Result Arguments;
            
      public:
        static void execute(runtime_wrappers_t & runtime_wrappers, std::string const & operation_name)
        {
          unsigned int arg_pos = 0;  
          typelist_utils::ForEach<Arguments, foreach_functor>::execute(arg_pos, runtime_wrappers, operation_name);
          if(requires_local_buffer_list<CurrentList>::value)
          {
            runtime_wrappers.insert(runtime_wrappers_t::value_type(operation_name,
                                                                   std::make_pair(arg_pos,
                                                                                  new result_of::shared_memory_wrapper())
                                                                  )
                                   );
          }
        }
    };

    template<class ARG>
    struct program_infos
    {
      typedef typename tree_utils::register_temporaries<ARG,false, typename ARG::LHS>::Result   NewARG;
      typedef typelist<NewARG,NullType>                                                         MainOperation_Init;
      typedef typename register_kernels<MainOperation_Init, typename ARG::LHS>::Result          KernelsList;
      
      static std::string value(std::string const & name_hint, KernelsSources & sources, runtime_wrappers_t & runtime_wrappers)
      {
        std::string operation_name = ARG::name();
        std::string program_name( (!name_hint.empty()) ? name_hint : operation_name );
        
        fill_sources<KernelsList,
                     MainOperation_Init,
                     0,
                     typelist_utils::length<KernelsList>::value - 1>::execute(sources,operation_name);
                     
        fill_args<KernelsList,
                  MainOperation_Init,
                  0,
                  typelist_utils::length<KernelsList>::value - 1>::execute(runtime_wrappers,operation_name);
                  
        return program_name;
      }
    };



  } // namespace generator
} // namespace viennacl
#endif
