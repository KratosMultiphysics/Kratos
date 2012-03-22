#ifndef VIENNACL_TOOLS_TOOLS_HPP_
#define VIENNACL_TOOLS_TOOLS_HPP_

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

/** @file tools.hpp
    @brief Various little tools used here and there in ViennaCL.
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/tools/adapter.hpp"


#ifdef VIENNACL_HAVE_UBLAS  
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#endif

#ifdef VIENNACL_HAVE_EIGEN  
#include <Eigen/Core>
#include <Eigen/Sparse>
#endif

#ifdef VIENNACL_HAVE_MTL4
#include <boost/numeric/mtl/mtl.hpp>
#endif

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {
    
    /** @brief Supply suitable increment functions for the iterators: */
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
    inline std::string readTextFromFile(const std::string & filename)
    {
      std::ifstream f(filename.c_str());
      if (!f) return std::string();

      std::stringstream result;
      std::string tmp;
      while (std::getline(f, tmp))
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
    inline std::string strReplace(const std::string & text, std::string to_search, std::string to_replace)
    {
      std::string::size_type pos = 0;
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
    * @param source          The source string
    * @param platform_info   An info string that contains the OpenCL platform vendor
    * @return   The double precision kernel
    */
    inline std::string make_double_kernel(std::string const & source, std::string const & fp_extension)
    {
      std::stringstream ss;
      ss << "#pragma OPENCL EXTENSION " << fp_extension << " : enable\n\n";
      
      std::string result = ss.str();
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
      //take care: using a plain, naive .size() on the left hand side type can cause subtle side-effects!
    };

    //Standard case: LHS is the vector type and carries the correct size
    template <typename ScalarType, unsigned int A, typename RHS>
    struct VECTOR_SIZE_DEDUCER<const viennacl::vector<ScalarType, A>, RHS, viennacl::op_prod>
    {
      static size_t size(const viennacl::vector<ScalarType, A> & lhs,
                         const RHS & rhs) { return lhs.size(); }
    };

    template <typename ScalarType, unsigned int A, typename RHS>
    struct VECTOR_SIZE_DEDUCER<const viennacl::vector<ScalarType, A>, RHS, viennacl::op_div>
    {
      static size_t size(const viennacl::vector<ScalarType, A> & lhs,
                         const RHS & rhs) { return lhs.size(); }
    };
    
    //special case: matrix-vector product: Return the number of rows of the matrix
    template <typename ScalarType, typename F, unsigned int Amat, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<const viennacl::matrix<ScalarType, F, Amat>, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static size_t size(const viennacl::matrix<ScalarType, F, Amat> & lhs,
                         const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    template <typename ScalarType, unsigned int Amat, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<const viennacl::circulant_matrix<ScalarType, Amat>, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static size_t size(const viennacl::circulant_matrix<ScalarType, Amat> & lhs,
                         const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };
    
    template <typename ScalarType, unsigned int Amat, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<const viennacl::compressed_matrix<ScalarType, Amat>, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static size_t size(const viennacl::compressed_matrix<ScalarType, Amat> & lhs,
                         const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    template <typename ScalarType, unsigned int Amat, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<const viennacl::coordinate_matrix<ScalarType, Amat>, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static size_t size(const viennacl::coordinate_matrix<ScalarType, Amat> & lhs,
                         const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };
    
    //special case: transposed matrix-vector product: Return the number of cols(!) of the matrix
    template <typename ScalarType, typename F, unsigned int Amat, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<const viennacl::matrix_expression< const viennacl::matrix<ScalarType, F, Amat>,
                                                                  const viennacl::matrix<ScalarType, F, Amat>,
                                                                  op_trans>,
                               const viennacl::vector<ScalarType, A>,
                               viennacl::op_prod>
    {
      static size_t size(const viennacl::matrix_expression< const viennacl::matrix<ScalarType, F, Amat>,
                                                            const viennacl::matrix<ScalarType, F, Amat>,
                                                            op_trans> & lhs,
                         const viennacl::vector<ScalarType, A> & rhs) { return lhs.lhs().size2(); }
    };

    
    
    
    
    /** @brief Obtain the cpu scalar type from a type, including a GPU type like viennacl::scalar<T>
    *
    * @tparam T   Either a CPU scalar type or a GPU scalar type
    */
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER
    {
      //force compiler error if type cannot be deduced
      //typedef T       ResultType;
    };

    template <>
    struct CPU_SCALAR_TYPE_DEDUCER< float >
    {
      typedef float       ResultType;
    };

    template <>
    struct CPU_SCALAR_TYPE_DEDUCER< double >
    {
      typedef double       ResultType;
    };
    
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::scalar<T> >
    {
      typedef T       ResultType;
    };

    template <typename T, unsigned int A>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::vector<T, A> >
    {
      typedef T       ResultType;
    };

    template <typename T, typename F, unsigned int A>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::matrix<T, F, A> >
    {
      typedef T       ResultType;
    };

    
    template <typename T, typename F, unsigned int A>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::matrix_expression<const matrix<T, F, A>, const matrix<T, F, A>, op_trans> >
    {
      typedef T       ResultType;
    };

        
  } //namespace tools
} //namespace viennacl
    

#endif
