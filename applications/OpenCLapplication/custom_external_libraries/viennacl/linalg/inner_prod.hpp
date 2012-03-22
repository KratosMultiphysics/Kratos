#ifndef VIENNACL_LINALG_INNER_PROD_HPP_
#define VIENNACL_LINALG_INNER_PROD_HPP_

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

/** @file inner_prod.hpp
    @brief Generic interface for the computation of inner products. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"

namespace viennacl
{
  //
  // generic inner_prod function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {
    
    #ifdef VIENNACL_HAVE_EIGEN
    // ----------------------------------------------------
    // EIGEN
    //
      #if defined(_MSC_VER) && _MSC_VER < 1500        //Visual Studio 2005 needs special treatment
      float
      inner_prod(Eigen::VectorXf const & v1,
                 Eigen::VectorXf const & v2)
      {
        return v1 * v2;
      }
      
      double
      inner_prod(Eigen::VectorXd const & v1,
                 Eigen::VectorXd const & v2)
      {
        return v1 * v2;
      }
      
      #else    
      template< typename VectorT1, typename VectorT2 >
      typename VectorT1::RealScalar
      inner_prod(VectorT1 const& v1, VectorT2 const& v2, 
          typename viennacl::enable_if< viennacl::is_eigen< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                                              >::type* dummy = 0)
      {
        //std::cout << "eigen .. " << std::endl;
        return v1.dot(v2);
      }
      #endif
    #endif
    
    #ifdef VIENNACL_HAVE_MTL4
    // ----------------------------------------------------
    // MTL4
    //
      #if defined(_MSC_VER) && _MSC_VER < 1500        //Visual Studio 2005 needs special treatment
      template <typename ScalarType>
      ScalarType inner_prod(mtl::dense_vector<ScalarType> const & v1,
                            mtl::dense_vector<ScalarType> const & v2)
      {
        return mtl::dot(v1, v2);
      }
      #else    
      template< typename VectorT1, typename VectorT2 >
      typename VectorT1::value_type
      inner_prod(VectorT1 const& v1, VectorT2 const& v2, 
          typename viennacl::enable_if< viennacl::is_mtl4< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                                              >::type* dummy = 0)
      {
        //std::cout << "mtl4 .. " << std::endl;
        return mtl::dot(v1, v2);
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
      inner_prod(boost::numeric::ublas::vector<ScalarType> const & v1,
                 boost::numeric::ublas::vector<ScalarType> const & v2)
      {
        // std::cout << "ublas .. " << std::endl;
        return boost::numeric::ublas::inner_prod(v1, v2);
      }
      #else    
      template< typename VectorT1, typename VectorT2 >
      typename VectorT1::value_type
      inner_prod(VectorT1 const& v1, VectorT2 const& v2, 
          typename viennacl::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                                              >::type* dummy = 0)
      {
        //std::cout << "ublas .. " << std::endl;
        return boost::numeric::ublas::inner_prod(v1, v2);
      }
      #endif
    #endif

    // ----------------------------------------------------
    // STL
    //
    template< typename VectorT1, typename VectorT2 >
    typename VectorT1::value_type
    inner_prod(VectorT1 const& v1, VectorT2 const& v2, 
         typename viennacl::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                                     >::type* dummy = 0)
    {
      assert(v1.size() == v2.size());
      //std::cout << "stl .. " << std::endl;
      typename VectorT1::value_type result = 0;
      for (typename VectorT1::size_type i=0; i<v1.size(); ++i)
        result += v1[i] * v2[i];
      
      return result;
    }

    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment1, unsigned int alignment2 >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment1>, 
                                 const viennacl::vector<ScalarType, alignment2>,
                                 viennacl::op_inner_prod >
    inner_prod(viennacl::vector<ScalarType, alignment1> const & vector1, viennacl::vector<ScalarType, alignment2> const & vector2, 
         typename viennacl::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment1> >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "viennacl .. " << std::endl;
      return viennacl::linalg::inner_prod_impl(vector1, vector2);
    }
  } // end namespace linalg
} // end namespace viennacl
#endif


