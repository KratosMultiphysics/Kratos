
#ifndef VIENNACL_LINALG_DETAIL_ILU0_HPP_
#define VIENNACL_LINALG_DETAIL_ILU0_HPP_

/* =========================================================================
   Copyright (c) 2010-2011, Institute for Microelectronics,
   Institute for Analysis and Scientific Computing,
   TU Wien.

   -----------------
   ViennaCL - The Vienna Computing Library
   -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/linalg/detail/ilu/ilu0.hpp
  @brief Implementations of incomplete factorization preconditioners with static nonzero pattern. Contributed by Evan Bollig.

  ILU0 (Incomplete LU with zero fill-in) 
  - All preconditioner nonzeros exist at locations that were nonzero in the input matrix. 
  - The number of nonzeros in the output preconditioner are exactly the same number as the input matrix

 Evan Bollig 3/30/12
 
 Adapted from viennacl/linalg/detail/ilut.hpp

*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/detail/ilu/common.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for incomplete LU factorization with threshold (ILUT)
    */
    class ilu0_tag
    {
      public:
        /** @brief The constructor.
          *
          * @param row_start     The starting row for the block to which we apply ILU
          * @param row_end       The end column of the block to which we apply ILU
          */
        ilu0_tag(unsigned int row_start = 0, unsigned int row_end = -1)
            : _row_start(row_start),  
            _row_end(row_end) {}
              
      public: 
        unsigned int _row_start, _row_end;
    };


    /** @brief Implementation of a ILU-preconditioner with static pattern
      *
      * refer to the Algorithm in Saad's book (1996 edition)
      *
      *  @param input   The input matrix. Type requirements: const_iterator1 for iteration along rows, const_iterator2 for iteration along columns
      *  @param output  The output matrix. Type requirements: const_iterator1 for iteration along rows, const_iterator2 for iteration along columns and write access via operator()
      *  @param tag     An ilu0_tag in order to dispatch among several other preconditioners.
      */
    template<typename MatrixType, typename LUType>
    void precondition(MatrixType const & input, LUType & output, ilu0_tag const & tag)
    {
      typedef std::map<unsigned int, double>          SparseVector;
      typedef typename SparseVector::iterator         SparseVectorIterator;
      typedef typename MatrixType::const_iterator1    InputRowIterator;  //iterate along increasing row index
      typedef typename MatrixType::const_iterator2    InputColIterator;  //iterate along increasing column index
      typedef typename LUType::iterator1              OutputRowIterator;  //iterate along increasing row index
      typedef typename LUType::iterator2              OutputColIterator;  //iterate along increasing column index

      output.clear();
      assert(input.size1() == output.size1());
      assert(input.size2() == output.size2());
      output.resize(static_cast<unsigned int>(input.size1()), static_cast<unsigned int>(input.size2()), false);
      SparseVector w;


      std::map<double, unsigned int> temp_map;

      // For i = 2, ... , N, DO
      for (InputRowIterator row_iter = input.begin1(); row_iter != input.end1(); ++row_iter)
      {
        w.clear();
        for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
        {
          // Only work on the block described by (row_start:row_end, row_start:row_end)
          if ((static_cast<unsigned int>(row_iter.index1()) >= tag._row_start) && (static_cast<unsigned int>(row_iter.index1()) < tag._row_end))
          {
              if ((static_cast<unsigned int>(col_iter.index2()) >= tag._row_start) && (static_cast<unsigned int>(col_iter.index2()) < tag._row_end))
              {
                  w[static_cast<unsigned int>(col_iter.index2())] = *col_iter;
              }
          } 
          else 
          {
              // Put identity on the excluded diagonal
              w[static_cast<unsigned int>(row_iter.index1())] = 1.; 
          }
        }

        //line 3:
        OutputRowIterator row_iter_out = output.begin1();
        for (SparseVectorIterator k = w.begin(); k != w.end(); ++k)
        {
          unsigned int index_k = k->first;
          // Enforce i = 2 and 
          if (index_k >= static_cast<unsigned int>(row_iter.index1()))
              break;

          detail::ilu_inc_row_iterator_to_row_index(row_iter_out, index_k);

          //line 3: temp = a_ik = a_ik / a_kk
          double temp = k->second / output(index_k, index_k);
          if (output(index_k, index_k) == 0.0)
          {
              std::cerr << "ViennaCL: FATAL ERROR in ILUT(): Diagonal entry is zero in row " << index_k << "!" << std::endl;

          }

          for (OutputColIterator j = row_iter_out.begin(); j != row_iter_out.end(); ++j)
          {
              // Only fill if it a nonzero element of the input matrix
              if (input(row_iter.index1(), j.index2())) {
                  // Follow standard ILU algorithm (i.e., for j = k+1, ... , N)
                  if (j.index2() > index_k) 
                  {
                      // set a_ij
                      w[j.index2()] -= temp * *j;
                  }
              }
          }
          // Set a_ik
          w[index_k] = temp;
          
        } //for k

        // Write rows back to LU factor output
        unsigned int k_count = 0; 
        for (SparseVectorIterator k = w.begin(); k != w.end(); ++k )
        {
          output(static_cast<unsigned int>(row_iter.index1()), k->first) = static_cast<typename LUType::value_type>(w[k->first]);
          k_count ++; 
        }
      } //for i
    }


    /** @brief ILUT preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class ilu0_precond
    {
        typedef typename MatrixType::value_type      ScalarType;

      public:
        ilu0_precond(MatrixType const & mat, ilu0_tag const & tag) : _tag(tag), LU(mat.size1())
        {
            //initialize preconditioner:
            //std::cout << "Start CPU precond" << std::endl;
            init(mat);          
            //std::cout << "End CPU precond" << std::endl;
        }

        template <typename VectorType>
            void apply(VectorType & vec) const
            {
                viennacl::tools::const_sparse_matrix_adapter<ScalarType> LU_const_adapter(LU);
                viennacl::linalg::detail::ilu_lu_substitute(LU_const_adapter, vec);
            }

      private:
        void init(MatrixType const & mat)
        {
            viennacl::tools::sparse_matrix_adapter<ScalarType>       LU_adapter(LU);
            viennacl::linalg::precondition(mat, LU_adapter, _tag);
        }

        ilu0_tag const & _tag;
        
        public: std::vector< std::map<unsigned int, ScalarType> > LU;
    };


    /** @brief ILUT preconditioner class, can be supplied to solve()-routines.
      *
      *  Specialization for compressed_matrix
      */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT>
    class ilu0_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT> >
    {
        typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>   MatrixType;

        public:
        ilu0_precond(MatrixType const & mat, ilu0_tag const & tag) : _tag(tag), LU(mat.size1())
        {
            //initialize preconditioner:
            //std::cout << "Start GPU precond" << std::endl;
            init(mat);          
            //std::cout << "End GPU precond" << std::endl;
        }

        void apply(vector<ScalarType> & vec) const
        {
            copy(vec, temp_vec);
            //lu_substitute(LU, vec);
            viennacl::tools::const_sparse_matrix_adapter<ScalarType> LU_const_adapter(LU);
            viennacl::linalg::detail::ilu_lu_substitute(LU_const_adapter, temp_vec);

            copy(temp_vec, vec);
        }

        private:
        void init(MatrixType const & mat)
        {
            std::vector< std::map<unsigned int, ScalarType> > temp(mat.size1());
            //std::vector< std::map<unsigned int, ScalarType> > LU_cpu(mat.size1());

            //copy to cpu:
            copy(mat, temp);

            viennacl::tools::const_sparse_matrix_adapter<ScalarType>       temp_adapter(temp);
            viennacl::tools::sparse_matrix_adapter<ScalarType>       LU_adapter(LU);
            viennacl::linalg::precondition(temp_adapter, LU_adapter, _tag);

            temp_vec.resize(mat.size1());

        }

        ilu0_tag const & _tag;
        //MatrixType LU;
        public: std::vector< std::map<unsigned int, ScalarType> > LU;
        private: mutable std::vector<ScalarType> temp_vec;
    };

  }
}




#endif



