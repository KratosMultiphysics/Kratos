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

   file changelog: - Jun 3, 2010   New from scratch for v1.0.1
======================================================================= */

#ifndef _VIENNACL_NORM_2_HPP_
#define _VIENNACL_NORM_2_HPP_

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
    
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename VectorT >
    typename VectorT::value_type
    norm_2(VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::norm_2(vector);
    }
    #endif
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                 const viennacl::vector<ScalarType, alignment>,
                                 viennacl::op_norm_2 >
    norm_2(viennacl::vector<ScalarType, alignment> const & vector, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment> >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment>, 
                                          const viennacl::vector<ScalarType, alignment>,
                                          viennacl::op_norm_2 >(vector, vector);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





