#ifndef VIENNACL_GENERATOR_META_TOOLS_UTILS_HPP
#define VIENNACL_GENERATOR_META_TOOLS_UTILS_HPP

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

/** @file viennacl/generator/meta_tools/utils.hpp
 *  @brief Various metaprogramming utilities
 *
 *  Generator code contributed by Philippe Tillet
 */


#include <algorithm>
#include <typeinfo>
#include <iostream>
#include "viennacl/matrix.hpp"

#define VIENNACL_STATIC_ASSERT( x ) typedef char __STATIC_ASSERT__[( x )?1:-1]

namespace viennacl 
{

  class any;

  template<class T>
  T any_cast(any& a);

  class value_base
  {
    public:
      virtual ~value_base() { }
      virtual value_base* clone() const = 0;
      virtual std::type_info const & type() const = 0;
  };

  template <class T>
  class value : public value_base
  {
      friend T any_cast<>(any& a);
      
      T t;
      
    public:
      value(const T& t_) : t(t_) { }
      value_base* clone() const
      {
          return new value(t);
      }
      
      std::type_info const &type() const {
          return typeid(T);
      }
  };  

  class any
  {
      template<class T>
      friend T any_cast(any & a);
      
      value_base* v;
      
    public:
      any() : v(0) { }
      
      template <class value_type>
      any(const value_type& v_) : v(new value<value_type>(v_)) { }

      any(any const & other) : v(other.v ? other.v->clone() : 0) {}

      any& operator=(const any& other)
      {
          if(&other != this)
          {
              any copy(other);
              swap(copy);
          }
          return *this;
      }

      void swap(any& other)
      {
          std::swap(v, other.v);
      }

      std::type_info const & type()
      {
        return v->type();
      }
      
      ~any() { delete v; }
  };

  class bad_any_cast : public std::bad_cast
  {
    public:
      virtual const char * what() const throw()
      {
          return "viennacl::bad_any_cast: "
                 "failed conversion using viennacl::any_cast";
      }
  };

  template <class T>
  T any_cast(any& a)
  {
    value<T>* v = dynamic_cast<value<T>*>(a.v);
    
    if(v == 0)
      throw bad_any_cast();
    else
      return v->t;
  }


  namespace generator
  {
    struct NullType 
    {
      static const std::string name() 
      {
          return "Null\n" ;
      }
    };

    template <class T>
    inline std::string to_string ( T const t ) 
    {
      std::stringstream ss;
      ss << t;
      return ss.str();
    }

    inline std::string to_string(viennacl::row_major    const) { return "rowmajor"; }
    inline std::string to_string(viennacl::column_major const) { return "columnmajor"; }


    template <int v>
    struct Int2Type 
    {
      enum { value = v };
    };

    template<class TypeTrue, class TypeFalse, bool cond>
    struct get_type_if 
    {
      typedef TypeTrue Result;
    };

    template<class TypeTrue, class TypeFalse>
    struct get_type_if<TypeTrue, TypeFalse, false> 
    {
      typedef TypeFalse Result;
    };

    template<class T, class U>
    struct are_same_type
    {
      enum { value = 0 };
    };

    template<class T>
    struct are_same_type<T,T>
    {
      enum { value = 1 };
    };


        
    template <bool B, class T = void>
    struct enable_if_c 
    {
      typedef T type;
    };

    template <class T>
    struct enable_if_c<false, T> {};

    template <class Cond, class T = void>
    struct enable_if : public enable_if_c<Cond::value, T> {};


    template <bool B, class T = void>
    struct disable_if_c 
    {
      typedef T type;
    };

    template <class T>
    struct disable_if_c<true, T> {};

    template <class Cond, class T = void>
    struct disable_if : public disable_if_c<Cond::value, T> {};



    template<class T>
    struct print_align1_type;

    template<>
    struct print_align1_type<int> 
    {
      static const std::string value() { return "int"; }
    };

    template<>
    struct print_align1_type<unsigned int>
    {
      static const std::string value() { return "unsigned int"; }
    };

    template<>
    struct print_align1_type<long> 
    {
      static const std::string value() { return "long"; }
    };

    template<>
    struct print_align1_type<unsigned long> 
    {
      static const std::string value() { return "long"; }
    };

    template<>
    struct print_align1_type<float> 
    {
      static const std::string value() { return "float"; }
    };

    template<>
    struct print_align1_type<double> 
    {
      static const std::string value() { return "double"; }
    };

    template<typename T, unsigned int ALIGNMENT>
    struct print_aligned_type 
    {
	    static const std::string value() 
	    {
        return print_align1_type<T>::value() + to_string ( ALIGNMENT );
      }
    };

    template<typename T>
    struct print_aligned_type<T, 1>
    {
	    static const std::string value() 
	    {
        return print_align1_type<T>::value();
      }
    };

    template<typename T, unsigned int ALIGNMENT>
    struct print_type 
    {
      static const std::string value() 
      {
        return print_aligned_type<T,ALIGNMENT>::value();
      }
    };

    template<typename T, unsigned int ALIGNMENT>
    struct print_type<T*, ALIGNMENT> 
    {
      static const std::string value() 
      {
        return print_type<T,ALIGNMENT>::value() + "*" ;
      }
    };

  }
}

#endif


