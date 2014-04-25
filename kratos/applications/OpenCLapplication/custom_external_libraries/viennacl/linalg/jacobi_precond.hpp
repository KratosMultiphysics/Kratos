#ifndef VIENNACL_LINALG_JACOBI_PRECOND_HPP_
#define VIENNACL_LINALG_JACOBI_PRECOND_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/jacobi_precond.hpp
    @brief Implementation of a simple Jacobi preconditioner
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/sparse_matrix_operations.hpp"
#include "viennacl/linalg/row_scaling.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for a jacobi preconditioner
    */
    class jacobi_tag {};


    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines. Generic version for non-ViennaCL matrices.
    */
    template <typename MatrixType,
              bool is_viennacl = detail::row_scaling_for_viennacl<MatrixType>::value >
    class jacobi_precond
    {
      typedef typename MatrixType::value_type      ScalarType;

      public:
        jacobi_precond(MatrixType const & mat, jacobi_tag const &) : diag_A(viennacl::traits::size1(mat))
        {
          init(mat);
        }

        void init(MatrixType const & mat)
        {
          diag_A.resize(viennacl::traits::size1(mat));  //resize without preserving values

          for (typename MatrixType::const_iterator1 row_it = mat.begin1();
                row_it != mat.end1();
                ++row_it)
          {
            bool diag_found = false;
            for (typename MatrixType::const_iterator2 col_it = row_it.begin();
                  col_it != row_it.end();
                  ++col_it)
            {
              if (col_it.index1() == col_it.index2())
              {
                diag_A[col_it.index1()] = *col_it;
                diag_found = true;
              }
            }
            if (!diag_found)
              throw "ViennaCL: Zero in diagonal encountered while setting up Jacobi preconditioner!";
          }
        }


        /** @brief Apply to res = b - Ax, i.e. jacobi applied vec (right hand side),  */
        template <typename VectorType>
        void apply(VectorType & vec) const
        {
          assert(viennacl::traits::size(diag_A) == viennacl::traits::size(vec) && bool("Size mismatch"));
          for (vcl_size_t i=0; i<diag_A.size(); ++i)
            vec[i] /= diag_A[i];
        }

      private:
        std::vector<ScalarType> diag_A;
    };


    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename MatrixType>
    class jacobi_precond< MatrixType, true>
    {
        typedef typename viennacl::result_of::cpu_value_type<typename MatrixType::value_type>::type  ScalarType;

      public:
        jacobi_precond(MatrixType const & mat, jacobi_tag const &) : diag_A(mat.size1(), viennacl::traits::context(mat))
        {
          init(mat);
        }


        void init(MatrixType const & mat)
        {
          detail::row_info(mat, diag_A, detail::SPARSE_ROW_DIAGONAL);
        }


        template <unsigned int ALIGNMENT>
        void apply(viennacl::vector<ScalarType, ALIGNMENT> & vec) const
        {
          assert(viennacl::traits::size(diag_A) == viennacl::traits::size(vec) && bool("Size mismatch"));
          vec = element_div(vec, diag_A);
        }

      private:
        viennacl::vector<ScalarType> diag_A;
    };

  }
}




#endif



