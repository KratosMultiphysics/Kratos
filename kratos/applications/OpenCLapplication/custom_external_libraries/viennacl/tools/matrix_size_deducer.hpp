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

#ifndef _VIENNACL_TOOLS_MATRIX_SIZE_DEDUCER_HPP_
#define _VIENNACL_TOOLS_MATRIX_SIZE_DEDUCER_HPP_

/** @file matrix_size_deducer.hpp
    @brief Helper implementations that deduce the dimensions of the supplied matrix-valued expressions.
*/

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"
#include "viennacl/tools/adapter.hpp"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {

    /** @brief Deduces the size of the resulting vector represented by a vector_expression from the operands
    *
    * @tparam LHS   The left hand side operand
    * @tparam RHS   The right hand side operand
    * @tparam OP    The operation tag
    */
    template <typename LHS, typename RHS, typename OP>
    struct MATRIX_SIZE_DEDUCER
    {
      //Standard case: size1 from lhs, size2 from rhs (fits most cases)
      static unsigned int size1(LHS & lhs, RHS & rhs) { return lhs.size1(); }
      static unsigned int size2(LHS & lhs, RHS & rhs) { return rhs.size2(); }
    };
    
    //special case: outer vector product:
    template <typename ScalarType, unsigned int A1, unsigned int A2>
    struct MATRIX_SIZE_DEDUCER<viennacl::vector<ScalarType, A1>,
                               viennacl::vector<ScalarType, A2>,
                               viennacl::op_prod>
    {
      static unsigned int size1(viennacl::vector<ScalarType, A1> & lhs,
                                viennacl::vector<ScalarType, A2> & rhs) { return lhs.size1(); }

      static unsigned int size2(viennacl::vector<ScalarType, A1> & lhs,
                                viennacl::vector<ScalarType, A2> & rhs) { return rhs.size2(); }
    };

    //special case: transposed matrix-matrix product: Return the number of rows of the matrix
    template <typename MatrixType, typename ScalarType, unsigned int A>
    struct MATRIX_SIZE_DEDUCER<MatrixType, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static unsigned int size(MatrixType & lhs, const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    
  }
}

#endif

