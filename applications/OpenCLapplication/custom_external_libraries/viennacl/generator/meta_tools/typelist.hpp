#ifndef VIENNACL_GENERATOR_META_TOOLS_TYPELIST_HPP
#define VIENNACL_GENERATOR_META_TOOLS_TYPELIST_HPP

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

/** @file typelist.hpp
 *  @brief Generic implementation of the typelist
 *
 *  Generator code contributed by Philippe Tillet
 */


#include "viennacl/generator/meta_tools/utils.hpp"
#include "viennacl/generator/traits/general_purpose_traits.hpp"


namespace viennacl 
{
  namespace generator
  {
    template <class T,class U>
    struct typelist
    {
      typedef T Head;
      typedef U Tail;

      static const std::string name()
      {
        return Head::name() + "  ;  " + Tail::name();
      }
    };

    namespace typelist_utils 
    {

      /*
       * Is Empty
       */

      template
      <
          typename T1  = NullType, typename T2  = NullType, typename T3  = NullType,
          typename T4  = NullType, typename T5  = NullType, typename T6  = NullType,
          typename T7  = NullType, typename T8  = NullType, typename T9  = NullType,
          typename T10 = NullType, typename T11 = NullType, typename T12 = NullType,
          typename T13 = NullType, typename T14 = NullType, typename T15 = NullType,
          typename T16 = NullType, typename T17 = NullType, typename T18 = NullType
      > 
      struct make_typelist
      {
        private:
          typedef typename make_typelist
          <
              T2 , T3 , T4 , 
              T5 , T6 , T7 , 
              T8 , T9 , T10, 
              T11, T12, T13,
              T14, T15, T16, 
              T17, T18
          >
          ::Result TailResult;

        public:
          typedef typelist<T1, TailResult> Result;
      };

      template <>
      struct make_typelist<>
      {
          typedef NullType Result;
      };

      template <class TList>
      struct is_empty 
      {
          enum { value = 0 };
      };

      template <>
      struct is_empty<NullType> 
      {
          enum { value = 1 };
      };


      /*
       * FOREACH
       */


      template <class TList,template<class> class Functor>
      struct ForEach;

      template <template<class> class Functor>
      struct ForEach<NullType,Functor> 
      {
        static void execute() {}
        
        template <class T1>
        static void execute(T1  & t1) {}
        
        template <class T1, class T2>
        static void execute(T1  & t1, T2  & t2) {}
        
        template <class T1, class T2, class T3>
        static void execute(T1  & t1, T2  & t2, T3 & t3) {}
        
        template <class T1, class T2, class T3, class T4>
        static void execute(T1  & t1, T2  & t2, T3 & t3, T4 & t4) {}
        
        template <class T1, class T2, class T3, class T4, class T5>
        static void execute(T1  & t1, T2  & t2, T3 & t3, T4 & t4, T5 & t5) {}
          
      };

      template <class T, class U,template<class> class Functor>
      struct ForEach< typelist<T, U>, Functor >
      {
        static void execute()
        {
          Functor<T>::execute();
          ForEach<U, Functor>::execute();
        }
        
        template <class T1>
        static void execute(T1  & t1)
        {
          Functor<T>::execute(t1);
          ForEach<U,Functor>::execute(t1);
        }
        
        template <class T1, class T2>
        static void execute(T1  & t1, T2  &t2)
        {
          Functor<T>::execute(t1,t2);
          ForEach<U,Functor>::execute(t1,t2);
        }
        
        template <class T1, class T2, class T3>
        static void execute(T1  & t1, T2  & t2, T3 & t3)
        {
          Functor<T>::execute(t1,t2,t3);
          ForEach<U,Functor>::execute(t1,t2,t3);
        }
        
        template <class T1, class T2, class T3, class T4>
        static void execute(T1  & t1, T2  & t2, T3 & t3, T4 & t4)
        {
          Functor<T>::execute(t1,t2,t3,t4);
          ForEach<U,Functor>::execute(t1,t2,t3,t4);
        }
        
        template <class T1, class T2, class T3, class T4, class T5>
        static void execute(T1  & t1, T2  & t2, T3 & t3, T4 & t4, T5 & t5)
        {
          Functor<T>::execute(t1,t2,t3,t4,t5);
          ForEach<U,Functor>::execute(t1,t2,t3,t4,t5);
        }
      };


      /*
       * length
       */


      template <class TList>
      struct length;

      template <>
      struct length<NullType> 
      {
        enum { value = 0 };
      };

      template <class T, class U>
      struct length< typelist<T, U> > 
      {
        enum { value = 1 + length<U>::value };
      };

      /*
       * type_at
       */

      template <class TList, unsigned int i>
      struct type_at;

      template <class Head, class Tail>
      struct type_at<typelist<Head, Tail>, 0> 
      {
        typedef Head Result;
      };

      template <class Head, class Tail, unsigned int i>
      struct type_at<typelist<Head, Tail>, i> 
      {
        typedef typename type_at<Tail, i - 1>::Result Result;
      };

      /*
       * index_of
       */

      template <class TList, class T>
      struct index_of;

      template <class T>
      struct index_of<NullType, T> 
      {
        enum { value = -1 };
      };

      template <class T, class Tail>
      struct index_of<typelist<T, Tail>, T> 
      {
        enum { value = 0 };
      };

      template <class Head, class Tail, class T>
      struct index_of<typelist<Head, Tail>, T> 
      {
        private:
          enum { temp = index_of<Tail, T>::value };
          
        public:
          enum { value = temp == -1 ? -1 : 1 + temp };
      };

      /*
       * append
       */

      template <class T1, class T2>
      struct compare1 
      {
        enum { value = static_cast<int> ( T1::id ) < static_cast<int> ( T2::id ) };
      };

      template <class T>
      struct compare1<NullType, T> 
      {
        enum { value = 0 };
      };


      template <class T1, class T2>
      struct true_comp 
      {
        enum { value = 1 };
      };

      template <class TList, class T, template<class,class> class Compare = true_comp>
      struct append;

      template <template<class,class> class Compare>
      struct append<NullType, NullType, Compare> 
      {
        typedef NullType Result;
      };

      template <class T, template<class,class> class Compare>
      struct append<NullType, T, Compare> 
      {
        typedef typelist<T,NullType> Result;
      };

      template <class Head, class Tail, template<class,class> class Compare>
      struct append<NullType, typelist<Head, Tail>, Compare > 
      {
        typedef typelist<Head, Tail> Result;
      };

      template <class Head, class Tail, template<class,class> class Compare>
      struct append<typelist<Head, Tail>, NullType, Compare > 
      {
        typedef typelist<Head, Tail> Result;
      };

      template <class Head, class Tail, class T, template<class,class> class Compare>
      struct append<typelist<Head,Tail>, T, Compare> 
      {
        private:
          typedef typelist<Head, typename append<Tail, T, Compare>::Result> TypeCompareFalse;
          typedef typelist<T, typelist<Head,Tail> > TypeCompareTrue;
          
        public:
          typedef typename get_type_if<TypeCompareTrue,TypeCompareFalse,Compare<T,Head>::value >::Result Result;
      };

      /*
       * fuse
       */

      template <class TList, class T, template<class,class> class Compare = true_comp>
      struct fuse 
      {
        typedef typename append<TList, T, Compare>::Result Result;
      };

      template <class Head1, class Tail1, class Head2, class Tail2, template<class,class> class Compare >
      struct fuse<typelist<Head1, Tail1>, typelist<Head2,Tail2>, Compare > 
      {
        private:
          typedef typename append< typelist<Head1,Tail1> , Head2, Compare>::Result NewResult;

        public:
          typedef typename fuse< NewResult, Tail2, Compare >::Result Result;
      };

      /*
       * erase
       */


      template<class TList, class T>
      struct erase;

      template <class T>
      struct erase<NullType, T> 
      {
        typedef NullType Result;
      };

      template <class T, class Tail>
      struct erase<typelist<T, Tail>, T>
      {
        typedef Tail Result;
      };

      template <class Head, class Tail, class T>
      struct erase<typelist<Head, Tail>, T> 
      {
        typedef typelist<Head,
                         typename erase<Tail, T>::Result>        Result;
      };

      template <class Head, class Tail, class Head2, class Tail2>
      struct erase<typelist<Head, Tail>, typelist<Head2, Tail2> > 
      {
        typedef typename erase< typename erase<typelist<Head,Tail>, Head2>::Result, Tail2 >::Result Result;
      };

      /*
       * No duplicate
       */

      template<class TList>
      struct no_duplicates;

      template <>
      struct no_duplicates<NullType> 
      {
        typedef NullType Result;
      };

      template <class Head, class Tail>
      struct no_duplicates< typelist<Head, Tail> > 
      {
        private:
          typedef typename no_duplicates<Tail>::Result L1;
          typedef typename erase<L1, Head>::Result L2;
          
        public:
          typedef typelist<Head, L2> Result;
      };

    }
  }
}

#endif
