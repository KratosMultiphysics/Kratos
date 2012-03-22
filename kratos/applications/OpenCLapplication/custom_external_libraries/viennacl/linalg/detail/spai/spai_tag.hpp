#ifndef VIENNACL_LINALG_DETAIL_SPAI_SPAI_TAG_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_SPAI_TAG_HPP

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


/** @file viennacl/linalg/detail/spai/spai_tag.hpp
    @brief Implementation of the spai tag holding SPAI configuration parameters. Experimental in 1.2.x.
    
    SPAI code contributed by Nikolay Lukash
*/



//#define __NO_STD_VECTOR
#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include <sstream>
#include "viennacl/ocl/backend.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/storage.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"
#include "boost/numeric/ublas/detail/matrix_assign.hpp"
//#include "boost/thread/thread.hpp"
#include "viennacl/linalg/detail/spai/block_matrix.hpp"
#include "viennacl/linalg/detail/spai/block_vector.hpp"

namespace viennacl
{
    namespace linalg
    {
      namespace detail
      {
        namespace spai
        {
        
          /** @brief A tag for SPAI
          * Contains values for the algorithm.
          * Must be passed to spai_precond constructor
          */
          class spai_tag
          {
              /** @brief Constructor
              * @param residual_norm_threshold Calculate until the norm of the residual falls below this threshold
              * @param iteration_limit maximum number of iterations
              * @param residual_threshold determines starting threshold in residual vector for including new indices into set J
              * @param is_static determines if static version of SPAI should be used
              * @param is_right determines if left or right preconditioner should be used
              */
            public:
              spai_tag(
                      double residual_norm_threshold = 1e-3,
                      unsigned int iteration_limit = 5,
                      double residual_threshold = 1e-2,
                      bool is_static = false,
                      bool is_right = false) :
              _residual_norm_threshold(residual_norm_threshold),
              _iteration_limit(iteration_limit),
              _residual_threshold(residual_threshold),
              _is_static(is_static),
              _is_right(is_right) {};
              
              double getResidualNormThreshold() const
              { return _residual_norm_threshold; }
              
              double getResidualThreshold() const
              { return _residual_threshold; }
              
              unsigned int getIterationLimit () const
              { return _iteration_limit; }
              
              bool getIsStatic() const
              { return _is_static; }
              
              bool getIsRight() const
              { return _is_right; }
              
              long getBegInd() const
              { return _beg_ind; }
              
              long getEndInd() const 
              { return _end_ind; }
              
              
              
              void setResidualNormThreshold(double residual_norm_threshold)
              {
                  if(residual_norm_threshold > 0)
                      _residual_norm_threshold = residual_norm_threshold;
              }
              
              void setResidualThreshold(double residual_threshold)
              {
                  if(residual_threshold > 0)
                      _residual_threshold = residual_threshold;
              }
              
              void setIterationLimit(unsigned int iteration_limit)
              {
                  if(iteration_limit > 0)
                      _iteration_limit = iteration_limit;
              }
              
              void setIsRight(bool is_right) { _is_right = is_right; }
              
              void setIsStatic(bool is_static) { _is_static = is_static; }
              
              void setBegInd(long beg_ind) { _beg_ind = beg_ind; }
              
              void setEndInd(long end_ind){ _end_ind = end_ind; }
              
              
            private:
              double _residual_norm_threshold;
              unsigned int _iteration_limit;
              long _beg_ind, _end_ind;
              double _residual_threshold;
              bool _is_static;
              bool _is_right;
          };
        
        }
      }
    }
}
#endif
