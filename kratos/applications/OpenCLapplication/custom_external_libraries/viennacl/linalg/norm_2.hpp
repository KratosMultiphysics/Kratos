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
======================================================================= */

#ifndef _VIENNACL_NORM_2_HPP_
#define _VIENNACL_NORM_2_HPP_

/** @file norm_2.hpp
    @brief Generic interface for the l^2-norm. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include <math.h>    //for sqrt()
#include "viennacl/forwards.h"
#include "tag_of.hpp"

namespace viennacl
{
  //
  // generic norm_2 function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {
    #ifdef VIENNACL_HAVE_MTL4
    // ----------------------------------------------------
    // MTL4
    //
      #if defined(_MSC_VER) && _MSC_VER < 1500        //Visual Studio 2005 needs special treatment
      template <typename ScalarType>
      ScalarType norm_2(mtl::dense_vector<ScalarType> const & v)
      {
        // std::cout << "mtl4 .. " << std::endl;
        return mtl::two_norm(v);
      }
      
      #else
      template< typename VectorT >
      typename VectorT::value_type
      norm_2(VectorT const& v, 
          typename viennacl::tools::enable_if< viennacl::is_mtl4< typename viennacl::traits::tag_of< VectorT >::type >::value
                                              >::type* dummy = 0)
      {
        // std::cout << "mtl4 .. " << std::endl;
        return mtl::two_norm(v);
      }
      #endif
    #endif
    
    
    #ifdef VIENNACL_HAVE_EIGEN
    // ----------------------------------------------------
    // EIGEN
    //
      #if defined(_MSC_VER) && _MSC_VER < 1500        //Visual Studio 2005 needs special treatment
      float norm_2(Eigen::VectorXf const & v)
      {
        // std::cout << "eigen .. " << std::endl;
        return v.norm();
      }
      
      double norm_2(Eigen::VectorXd const & v)
      {
        // std::cout << "eigen .. " << std::endl;
        return v.norm();
      }
      
      #else
      template< typename VectorT >
      typename VectorT::RealScalar
      norm_2(VectorT const& v, 
          typename viennacl::tools::enable_if< viennacl::is_eigen< typename viennacl::traits::tag_of< VectorT >::type >::value
                                              >::type* dummy = 0)
      {
        // std::cout << "ublas .. " << std::endl;
        return v.norm();
      }
      #endif
    #endif
    
    
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
      #if defined(_MSC_VER) && _MSC_VER < 1500        //Visual Studio 2005 needs special treatment
      template< typename ScalarType >
      ScalarType
      norm_2(boost::numeric::ublas::vector<ScalarType> const & v)
      {
        // std::cout << "ublas .. " << std::endl;
        return boost::numeric::ublas::norm_2(v);
      }
      #else
      template< typename VectorT >
      typename VectorT::value_type
      norm_2(VectorT const& v, 
          typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT >::type >::value
                                              >::type* dummy = 0)
      {
        // std::cout << "ublas .. " << std::endl;
        return boost::numeric::ublas::norm_2(v);
      }
      #endif
    #endif
    
    
    // ----------------------------------------------------
    // STL
    //
    template< typename VectorT>
    typename VectorT::value_type
    norm_2(VectorT const& v1,
         typename viennacl::tools::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< VectorT >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "stl .. " << std::endl;
      typename VectorT::value_type result = 0;
      for (typename VectorT::size_type i=0; i<v1.size(); ++i)
        result += v1[i] * v1[i];
      
      return sqrt(result);
    }
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                 const viennacl::vector<ScalarType, alignment>,
                                 viennacl::op_norm_2 >
    norm_2(viennacl::vector<ScalarType, alignment> const & v, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment> >::type >::value
                                            >::type* dummy = 0)
    {
       //std::cout << "viennacl .. " << std::endl;
      return viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                          const viennacl::vector<ScalarType, alignment>,
                                          viennacl::op_norm_2 >(v, v);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





