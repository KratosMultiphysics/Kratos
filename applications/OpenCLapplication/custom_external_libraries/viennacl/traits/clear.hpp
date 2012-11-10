#ifndef VIENNACL_TRAITS_CLEAR_HPP_
#define VIENNACL_TRAITS_CLEAR_HPP_

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

/** @file viennacl/traits/clear.hpp
    @brief Generic clear functionality for different vector and matrix types
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

#include "viennacl/traits/size.hpp"

#include <vector>
#include <map>

namespace viennacl
{
  namespace traits
  {
    
    //clear:
    template <typename VectorType>
    void clear(VectorType & vec)
    {
      typedef typename viennacl::result_of::size_type<VectorType>::type  size_type;
      
      for (size_type i=0; i<viennacl::traits::size(vec); ++i)
        vec[i] = 0;  //TODO: Quantity access can also be wrapped...
    }

    template <typename ScalarType, unsigned int ALIGNMENT>
    void clear(viennacl::vector<ScalarType, ALIGNMENT> & vec)
    {
      vec.clear();
    }
  } //namespace traits
} //namespace viennacl
    

#endif
