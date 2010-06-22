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

#ifndef _VIENNACL_INNERPROD_HPP_
#define _VIENNACL_INNERPROD_HPP_

#include "viennacl/forwards.h"
#include "tag_of.hpp"

namespace viennacl
{
  //
  // generic inner_prod function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename VectorT1, typename VectorT2 >
    typename VectorT1::value_type
    inner_prod(VectorT1 const& matrix, VectorT2 const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< VectorT1 >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::inner_prod(matrix, vector);
    }
    #endif
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename ScalarType, unsigned int alignment1, unsigned int alignment2 >
    viennacl::scalar_expression< const viennacl::vector<ScalarType, alignment1>, 
                                 const viennacl::vector<ScalarType, alignment2>,
                                 viennacl::op_inner_prod >
    inner_prod(viennacl::vector<ScalarType, alignment1> const & vector1, viennacl::vector<ScalarType, alignment2> const & vector2, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< viennacl::vector<ScalarType, alignment1> >::type >::value
                                            >::type* dummy = 0)
    {
      //std::cout << "viennacl .. " << std::endl;
      return viennacl::linalg::inner_prod_impl(vector1, vector2);
    }
  } // end namespace linalg
} // end namespace viennacl
#endif


