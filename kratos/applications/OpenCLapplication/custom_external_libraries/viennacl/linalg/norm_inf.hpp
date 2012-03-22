#ifndef VIENNACL_LINALG_NORM_INF_HPP_
#define VIENNACL_LINALG_NORM_INF_HPP_

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

/** @file norm_inf.hpp
    @brief Generic interface for the l^infty-norm. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include <math.h>    //for sqrt()
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/meta/enable_if.hpp"
#include "viennacl/meta/tag_of.hpp"

namespace viennacl
{
  //
  // generic norm_inf function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {
    
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename VectorT >
    typename viennacl::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT >::type >::value,
                                  typename VectorT::value_type      
                                >::type    
    norm_inf(VectorT const& v1)
    {
      // std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::norm_inf(v1);
    }
    #endif
    
    
    // ----------------------------------------------------
    // STL
    //
    template< typename VectorT>
    typename VectorT::value_type
    norm_inf(VectorT const& v1,
         typename viennacl::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< VectorT >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "stl .. " << std::endl;
      typename VectorT::value_type result = 0;
      for (typename VectorT::size_type i=0; i<v1.size(); ++i)
      {
        if (fabs(v1[i]) > result)
          result = fabs(v1[i]);
      }
      
      return result;
    }
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                 const viennacl::vector<ScalarType, alignment>,
                                 viennacl::op_norm_inf >
    norm_inf(viennacl::vector<ScalarType, alignment> const & v1, 
         typename viennacl::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment> >::type >::value
                                            >::type* dummy = 0)
    {
       //std::cout << "viennacl .. " << std::endl;
      return viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                          const viennacl::vector<ScalarType, alignment>,
                                          viennacl::op_norm_inf >(v1, v1);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





