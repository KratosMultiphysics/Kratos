#ifndef VIENNACL_META_RESULT_OF_HPP_
#define VIENNACL_META_RESULT_OF_HPP_

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

/** @file result_of.hpp
    @brief A collection of compile time type deductions
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"


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
    namespace result_of
    {
      //
      // Retrieve size_type 
      //
      template <typename T>
      struct size_type
      {
        typedef typename T::size_type   type;
      };

      #ifdef VIENNACL_HAVE_EIGEN
      template <class T, int a, int b, int c, int d, int e>
      struct size_type< Eigen::Matrix<T, a, b, c, d, e> >
      {
        typedef std::size_t   type;
      };
      
      template <>
      struct size_type<Eigen::VectorXf>
      {
        typedef std::size_t   type;
      };
      
      template <>
      struct size_type<Eigen::VectorXd>
      {
        typedef std::size_t   type;
      };

      template <typename T, int options>
      struct size_type<Eigen::SparseMatrix<T, options> >
      {
        typedef std::size_t   type;
      };
      #endif
      
      //
      // Retrieve value_type:
      //
      template <typename T>
      struct value_type
      {
        typedef typename T::value_type    type; 
      };

      //
      // Retrieve cpu value_type:
      //
      template <typename T>
      struct cpu_value_type
      {
        typedef typename T::ERROR_CANNOT_DEDUCE_CPU_SCALAR_TYPE_FOR_T    type; 
      };

      template <>
      struct cpu_value_type<float>
      {
        typedef float    type; 
      };
      
      template <>
      struct cpu_value_type<double>
      {
        typedef double    type; 
      };
      
      template <typename T>
      struct cpu_value_type<viennacl::scalar<T> >
      {
        typedef T    type; 
      };

      template <typename T, unsigned int ALIGNMENT>
      struct cpu_value_type<viennacl::vector<T, ALIGNMENT> >
      {
        typedef T    type; 
      };

      template <typename T>
      struct cpu_value_type<viennacl::vector_range<T> >
      {
        typedef typename cpu_value_type<T>::type    type; 
      };
      
      template <typename T1, typename T2, typename OP>
      struct cpu_value_type<viennacl::vector_expression<T1, T2, OP> >
      {
        typedef typename cpu_value_type<T1>::type    type; 
      };
      
      
      
      template <typename T, typename F, unsigned int ALIGNMENT>
      struct cpu_value_type<viennacl::matrix<T, F, ALIGNMENT> >
      {
        typedef T    type; 
      };
      
      template <typename T>
      struct cpu_value_type<viennacl::matrix_range<T> >
      {
        typedef typename cpu_value_type<T>::type    type; 
      };

      template <typename T1, typename T2, typename OP>
      struct cpu_value_type<viennacl::matrix_expression<T1, T2, OP> >
      {
        typedef typename cpu_value_type<T1>::type    type; 
      };
      
      
    #ifdef VIENNACL_HAVE_EIGEN  
      template <>
      struct value_type<Eigen::MatrixXf>
      {
        typedef Eigen::MatrixXf::RealScalar    type; 
      };
      
      template <>
      struct value_type<Eigen::MatrixXd>
      {
        typedef Eigen::MatrixXd::RealScalar    type; 
      };

      template <typename ScalarType, int option>
      struct value_type<Eigen::SparseMatrix<ScalarType, option> >
      {
        typedef ScalarType    type; 
      };

      template <>
      struct value_type<Eigen::VectorXf>
      {
        typedef Eigen::VectorXf::RealScalar    type; 
      };

      template <>
      struct value_type<Eigen::VectorXd>
      {
        typedef Eigen::VectorXd::RealScalar    type; 
      };
      
    #endif
      
      
      
      template <typename T>
      struct matrix_expression_internal_storage
      {
        typedef T &     type;
      };
     
      template <>
      struct matrix_expression_internal_storage<const float>
      {
        typedef float type;
      };
      
      template <>
      struct matrix_expression_internal_storage<const double>
      {
        typedef double type;
      };
      
      
    } //namespace result_of
} //namespace viennacl
    

#endif
