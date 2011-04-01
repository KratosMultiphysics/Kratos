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

#ifndef _VIENNACL_NORM_1_HPP_
#define _VIENNACL_NORM_1_HPP_

/** @file norm_1.hpp
    @brief Generic interface for the l^1-norm. See viennacl/linalg/vector_operations.hpp for implementations.
*/

#include <math.h>    //for sqrt()
#include "viennacl/forwards.h"
#include "tag_of.hpp"

namespace viennacl
{
  //
  // generic norm_1 function
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
    typename VectorT::value_type
    norm_1(VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::norm_1(vector);
    }
    #endif
    
    
    // ----------------------------------------------------
    // STL
    //
    template< typename VectorT>
    typename VectorT::value_type
    norm_1(VectorT const& v1,
         typename viennacl::tools::enable_if< viennacl::is_stl< typename viennacl::traits::tag_of< VectorT >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "stl .. " << std::endl;
      typename VectorT::value_type result = 0;
      for (typename VectorT::size_type i=0; i<v1.size(); ++i)
        result += fabs(v1[i]);
      
      return result;
    }
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                 const viennacl::vector<ScalarType, alignment>,
                                 viennacl::op_norm_1 >
    norm_1(viennacl::vector<ScalarType, alignment> const & vector, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment> >::type >::value
                                            >::type* dummy = 0)
    {
      return viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                          const viennacl::vector<ScalarType, alignment>,
                                          viennacl::op_norm_1 >(vector, vector);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





