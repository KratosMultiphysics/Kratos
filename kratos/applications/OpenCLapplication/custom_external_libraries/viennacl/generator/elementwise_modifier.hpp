#ifndef VIENNACL_GENERATOR_ELEMENTWISE_MODIFIER_HPP
#define VIENNACL_GENERATOR_ELEMENTWISE_MODIFIER_HPP

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

/** @file viennacl/generator/elementwise_modifier.hpp
 *   @brief Contains the stuffs related to the elementwise_modifier
 * 
 *  Generator code contributed by Philippe Tillet
 */

#include <typeinfo>
#include <string>
#include <algorithm>

#include "viennacl/generator/forwards.h"

namespace viennacl 
{
  namespace generator
  {

    /**
    * @brief Implementation of the elementwise_modifier
    * 
    * @tparam T the underlying expression to modify
    * @tparam U the function returning the modifier's expression
    */
    template<class T, std::string (*U)()>
    struct elementwise_modifier_impl
    {
      private:
        static std::string expr_name()
        {
          std::string res = U();
          std::replace(res.begin(),res.end(),'/','o');
          std::replace(res.begin(),res.end(),'*','x');
          std::replace(res.begin(),res.end(),'+','a');
          std::replace(res.begin(),res.end(),'-','s');
          std::replace(res.begin(),res.end(),' ','_');
          std::replace(res.begin(),res.end(),'(','p');
          std::replace(res.begin(),res.end(),')','p');
          return res;
        }
        
      public:
        typedef T PRIOR_TYPE;

        enum { id = -2 };

        static std::string name()
        {
          return expr_name() + '_' + T::name();
        }
        
        static std::string modify(std::string const & replacer) 
        {
          std::string result(U());
          int pos;
          while( (pos = result.find('X')) != std::string::npos )
          {
            result.replace(pos, 1, '(' + replacer + ')' );
          }
          
          return result;
        }
    };

    /** @brief Operator for creating an elementwise_modifier from an expression */
    template<std::string (*U)(),class T>
    elementwise_modifier_impl<T,U> elementwise_modifier( T const & t ) 
    {
      return elementwise_modifier_impl<T,U>();
    }

  }
}

#endif
