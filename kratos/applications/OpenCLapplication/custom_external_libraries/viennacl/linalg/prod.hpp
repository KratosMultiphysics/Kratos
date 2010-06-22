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

#ifndef _VIENNACL_PROD_HPP_
#define _VIENNACL_PROD_HPP_

#include "viennacl/forwards.h"
#include "tag_of.hpp"


namespace viennacl
{
  //
  // generic prod function
  //   uses tag dispatch to identify which algorithm
  //   should be called 
  //
  namespace linalg 
  {
    
    #ifdef VIENNACL_HAVE_UBLAS
    // ----------------------------------------------------
    // UBLAS
    //
    template< typename MatrixT, typename VectorT >
    VectorT 
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_ublas< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "ublas .. " << std::endl;
      return boost::numeric::ublas::prod(matrix, vector);
    }
    #endif
    
    // ----------------------------------------------------
    // VIENNACL
    //
    template< typename MatrixT, typename VectorT >
    viennacl::vector_expression< const MatrixT, 
      const VectorT, viennacl::op_prod >
    prod(MatrixT const& matrix, VectorT const& vector, 
         typename viennacl::tools::enable_if< viennacl::is_viennacl< typename viennacl::traits::tag_of< MatrixT >::type >::value
                                            >::type* dummy = 0)
    {
      // std::cout << "viennacl .. " << std::endl;
      return viennacl::linalg::prod_impl(matrix, vector);
    }

  } // end namespace linalg
} // end namespace viennacl
#endif





