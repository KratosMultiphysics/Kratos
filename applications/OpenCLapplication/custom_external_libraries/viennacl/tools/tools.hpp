/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_TOOLS_HPP_
#define _VIENNACL_TOOLS_HPP_

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/tools/adapter.hpp"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {
    
    //supply suitable increment functions for the iterators:
    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    struct MATRIX_ITERATOR_INCREMENTER<viennacl::row_iteration, viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
    {
      static void apply(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, unsigned int & row, unsigned int & col) { ++row; }
    };

    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    struct MATRIX_ITERATOR_INCREMENTER<viennacl::col_iteration, viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
    {
      static void apply(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, unsigned int & row, unsigned int & col) { ++col; }
    };

    
    
    template <bool b, class T = void> 
    struct enable_if
    {
      typedef T   type;
    };

    template <class T> 
    struct enable_if<false, T> {};
    
    
    /** @brief A guard that checks whether the floating point type of GPU types is either float or double */
    template <typename T>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT
    {
        typedef typename T::ERROR_SCALAR_MUST_HAVE_TEMPLATE_ARGUMENT_FLOAT_OR_DOUBLE  ResultType;
    };
    
    template <>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT<float>
    {
        typedef float  ResultType;
    };
    
    template <>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT<double>
    {
        typedef double  ResultType;
    };

    
    
    /** @brief Reads a text from a file into a std::string
    *
    * @param filename   The filename
    * @return The text read from the file
    */
    std::string readTextFromFile(const std::string & filename)
    {
      std::ifstream f(filename.c_str());
      if (!f) return std::string();

      std::stringstream result;
      std::string tmp;
      while (getline(f, tmp))
        result << tmp << std::endl;

      return result.str();
    }

    /** @brief Replaces all occurances of a substring by another stringstream
    *
    * @param text   The string to search in
    * @param to_search  The substring to search for
    * @param to_replace The replacement for found substrings
    * @return The resulting string
    */
    std::string strReplace(const std::string & text, std::string to_search, std::string to_replace)
    {
      std::string::size_type pos=0;
        std::string result;
        std::string::size_type found;
        while( (found = text.find(to_search, pos)) != std::string::npos )
        {
        result.append(text.substr(pos,found-pos));
            result.append(to_replace);
            pos = found + to_search.length();
        }
        if (pos < text.length())
            result.append(text.substr(pos));
        return result;
    }

    /** @brief Rounds an integer to the next multiple of another integer
    *
    * @tparam INT_TYPE  The integer type
    * @param to_reach   The integer to be rounded up (ceil operation)
    * @param base       The base
    * @return The smallest multiple of 'base' such that to_reach <= base
    */
    template <class INT_TYPE>
    INT_TYPE roundUpToNextMultiple(INT_TYPE to_reach, INT_TYPE base)
    {
      if (to_reach % base == 0) return to_reach;
      return ((to_reach / base) + 1) * base;
    }
    
    
    /** @brief Create a double precision kernel out of a single precision kernel
    *
    * @param source   The source string
    * @return   The double precision kernel
    */
    std::string make_double_kernel(std::string const & source)
    {
      #ifdef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK
        std::string result = "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n\n";
      #else
        std::string result = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
      #endif
      result.append(strReplace(source, "float", "double"));
      return result;
    }

    /** @brief Removes the const qualifier from a type */
    template <typename T>
    struct CONST_REMOVER
    {
      typedef T   ResultType;
    };

    template <typename T>
    struct CONST_REMOVER<const T>
    {
      typedef T   ResultType;
    };


    /** @brief Extracts the vector type from one of the two arguments. Used for the vector_expression type.
    *
    * @tparam LHS   The left hand side operand of the vector_expression
    * @tparam RHS   The right hand side operand of the vector_expression
    */
    template <typename LHS, typename RHS>
    struct VECTOR_EXTRACTOR_IMPL
    {
      typedef typename LHS::ERROR_COULD_NOT_EXTRACT_VECTOR_INFORMATION_FROM_VECTOR_EXPRESSION  ResultType;
    };
    
    template <typename LHS, typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<LHS, viennacl::vector<ScalarType, A> >
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    template <typename RHS, typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<viennacl::vector<ScalarType, A>, RHS>
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    //resolve ambiguities for previous cases:
    template <typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<viennacl::vector<ScalarType, A>, viennacl::vector<ScalarType, A> >
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    template <typename LHS, typename RHS>
    struct VECTOR_EXTRACTOR
    {
      typedef typename VECTOR_EXTRACTOR_IMPL<typename CONST_REMOVER<LHS>::ResultType,
                                              typename CONST_REMOVER<RHS>::ResultType>::ResultType      ResultType;
    };

    /** @brief Deduces the size of the resulting vector represented by a vector_expression from the operands
    *
    * @tparam LHS   The left hand side operand
    * @tparam RHS   The right hand side operand
    * @tparam OP    The operation tag
    */
    template <typename LHS, typename RHS, typename OP>
    struct VECTOR_SIZE_DEDUCER
    {
      //Standard case: LHS is the vector type and carries the correct size
      static unsigned int size(LHS & lhs, RHS & rhs) { return lhs.size(); }
    };
    
    //special case: matrix-vector product: Return the number of rows of the matrix
    template <typename MatrixType, typename ScalarType, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<MatrixType, viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static unsigned int size(MatrixType & lhs, viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    //special case: matrix-vector product: Return the number of rows of the matrix
    template <typename MatrixType, typename ScalarType, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<MatrixType, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static unsigned int size(MatrixType & lhs, const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    /** @brief Obtain the cpu scalar type from a type, including a GPU type like viennacl::scalar<T>
    *
    * @tparam T   Either a CPU scalar type or a GPU scalar type
    */
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER
    {
      typedef T       ResultType;
    };
    
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::scalar<T> >
    {
      typedef T       ResultType;
    };
    
  } //namespace tools
} //namespace viennacl
    

#endif
