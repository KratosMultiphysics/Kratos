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

#ifndef _VIENNACL_ILU_HPP_
#define _VIENNACL_ILU_HPP_

/** @file ilu.hpp
    @brief Implementations of incomplete factorization preconditioners
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {
    
    /** @brief A tag for incomplete LU factorization with threshold (ILUT)
    */
    class ilut_tag
    {
      public:
        /** @brief The constructor.
        *
        * @param entries_per_row  Number of nonzero entries per row in L and U. Note that L and U are stored in a single matrix, thus there are 2*entries_per_row in total.
        * @param drop_tolerance   The drop tolerance for ILUT
        */
        ilut_tag(unsigned int entries_per_row = 10,
                 double drop_tolerance = 1e-3) : _entries_per_row(entries_per_row), _drop_tolerance(drop_tolerance) {}; 

        void set_drop_tolerance(double tol)
        {
          if (tol > 0)
            _drop_tolerance = tol;
        }
        double get_drop_tolerance() const { return _drop_tolerance; }
        
        void set_entries_per_row(unsigned int e)
        {
          if (e > 0)
            _entries_per_row = e;
        }

        unsigned int get_entries_per_row() const { return _entries_per_row; }

      private:
        unsigned int _entries_per_row;
        double _drop_tolerance;
    };
    
    
    /** @brief Implementation of a ILU-preconditioner with threshold
    *
    * refer to Algorithm 10.6 by Saad's book (1996 edition)
    *
    *  @param input   The input matrix. Type requirements: const_iterator1 for iteration along rows, const_iterator2 for iteration along columns
    *  @param output  The output matrix. Type requirements: const_iterator1 for iteration along rows, const_iterator2 for iteration along columns and write access via operator()
    *  @param tag     An ilut_tag in order to dispatch among several other preconditioners.
    */
    template<typename MatrixType, typename LUType>
    void precondition(MatrixType const & input, LUType & output, ilut_tag const & tag)
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
      
      for (InputRowIterator row_iter = input.begin1(); row_iter != input.end1(); ++row_iter)
      {
    /*    if (i%10 == 0)
      std::cout << i << std::endl;*/
        
        //line 2:
        w.clear();
        for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
          w[static_cast<unsigned int>(col_iter.index2())] = *col_iter;

        //line 3:
        OutputRowIterator row_iter_out = output.begin1();
        for (SparseVectorIterator k = w.begin(); k != w.end();)
        {
          unsigned int index_k = k->first;
          if (index_k >= static_cast<unsigned int>(row_iter.index1()))
            break;
          
          while (row_iter_out.index1() < index_k)
            ++row_iter_out;
          
          //line 4:
          double temp = k->second / output(index_k, index_k);
          if (output(index_k, index_k) == 0.0)
          {
            std::cerr << "ViennaCL: FATAL ERROR in ILUT(): Diagonal entry is zero in row " << index_k << "!" << std::endl;
          }
          
          //line 5: (dropping rule to w_k)
          if ( fabs(temp) > tag.get_drop_tolerance())
          {
            //line 7:
            for (OutputColIterator j = row_iter_out.begin(); j != row_iter_out.end(); ++j)
            {
              if (j.index2() > index_k) //attention: manipulation of w(k->first) would invalidate iterator!
              {
                w[j.index2()] -= temp * *j;
              }
            }
            ++k;  //attention: manipulation of w(k->first) would invalidate iterator!
            w[index_k] = temp;// - temp * A(index_k, index_k);
          }
          else
            ++k;
        } //for k
        
        //Line 10: Apply a dropping rule to w
        //Step 1: Sort all entries:
        temp_map.clear();
        for (SparseVectorIterator k = w.begin(); k != w.end(); )
        {
          if (fabs(k->second) < tag.get_drop_tolerance())
          { 
            long index = k->first;
            ++k;
            w.erase(index);
          }
          else
          {
            double temp = fabs(k->second);
            while (temp_map.find(temp) != temp_map.end())
              temp *= 1.00000001; //make entry slightly larger to maintain uniqueness of the entry
            temp_map[temp] = k->first;
            ++k;
          }
        }

        //Lines 10-12: write the largest p values to L and U
        unsigned int written_L = 0;
        unsigned int written_U = 0;
        for (typename std::map<double, unsigned int>::reverse_iterator iter = temp_map.rbegin(); iter != temp_map.rend(); ++iter)
        {
          if (iter->second > static_cast<unsigned int>(row_iter.index1())) //entry for U
          {
            if (written_U < tag.get_entries_per_row())
            {
              output(static_cast<unsigned int>(row_iter.index1()), iter->second) = static_cast<typename LUType::value_type>(w[iter->second]);
              ++written_U;
            }
          }
          else if (iter->second == static_cast<unsigned int>(row_iter.index1()))
          {
            output(iter->second, iter->second) = static_cast<typename LUType::value_type>(w[static_cast<unsigned int>(row_iter.index1())]);
          }
          else //entry for L
          {
            if (written_L < tag.get_entries_per_row())
            {
              output(static_cast<unsigned int>(row_iter.index1()), iter->second) = static_cast<typename LUType::value_type>(w[iter->second]);
              ++written_L;
            }
          }
        }
      } //for i
    }


    /** @brief Generic inplace solution of a unit lower triangular system
    *   
    * @param mat  The system matrix
    * @param vec  The right hand side vector
    */
    template<typename MatrixType, typename VectorType>
    void inplace_solve(MatrixType const & mat, VectorType & vec, viennacl::linalg::unit_lower_tag)
    {
      typedef typename MatrixType::const_iterator1    InputRowIterator;  //iterate along increasing row index
      typedef typename MatrixType::const_iterator2    InputColIterator;  //iterate along increasing column index
      
      for (InputRowIterator row_iter = mat.begin1(); row_iter != mat.end1(); ++row_iter)
      {
        for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
        {
          if (col_iter.index2() < col_iter.index1())
            vec[col_iter.index1()] -= *col_iter * vec[col_iter.index2()];
        }
      }
    }

    /** @brief Generic inplace solution of a upper triangular system
    *   
    * @param mat  The system matrix
    * @param vec  The right hand side vector
    */
    template<typename MatrixType, typename VectorType>
    void inplace_solve(MatrixType const & mat, VectorType & vec, viennacl::linalg::upper_tag)
    {
      typedef typename MatrixType::const_reverse_iterator1    InputRowIterator;  //iterate along increasing row index
      typedef typename MatrixType::const_iterator2            InputColIterator;  //iterate along increasing column index
      typedef typename VectorType::value_type                 ScalarType;
      
      ScalarType diagonal_entry = 1.0;
      
      for (InputRowIterator row_iter = mat.rbegin1(); row_iter != mat.rend1(); ++row_iter)
      {
        for (InputColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
        {
          if (col_iter.index2() > col_iter.index1())
            vec[col_iter.index1()] -= *col_iter * vec[col_iter.index2()];
          if (col_iter.index2() == col_iter.index1())
            diagonal_entry = *col_iter;
        }
        vec[row_iter.index1()] /= diagonal_entry;
      }
    }

    /** @brief Generic LU substitution
    *   
    * @param mat  The system matrix
    * @param vec  The right hand side vector
    */
    template<typename MatrixType, typename VectorType>
    void lu_substitute(MatrixType const & mat, VectorType & vec)
    {
      inplace_solve(mat, vec, unit_lower_tag());
      inplace_solve(mat, vec, upper_tag());
    }


    /** @brief ILUT preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class ilut_precond
    {
      typedef typename MatrixType::value_type      ScalarType;
      
      public:
        ilut_precond(MatrixType const & mat, ilut_tag const & tag) : _tag(tag), LU(mat.size1())
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
          viennacl::linalg::lu_substitute(LU_const_adapter, vec);
        }
        
      private:
        void init(MatrixType const & mat)
        {
          viennacl::tools::sparse_matrix_adapter<ScalarType>       LU_adapter(LU);
          viennacl::linalg::precondition(mat, LU_adapter, _tag);
        }
        
        ilut_tag const & _tag;
        std::vector< std::map<unsigned int, ScalarType> > LU;
    };

    
    /** @brief ILUT preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT>
    class ilut_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT> >
    {
      typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>   MatrixType;
      
      public:
        ilut_precond(MatrixType const & mat, ilut_tag const & tag) : _tag(tag), LU(mat.size1())
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
          viennacl::linalg::lu_substitute(LU_const_adapter, temp_vec);
          
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
          
          //copy resulting preconditioner back to gpu:
          //copy(LU_cpu, LU);
        }
        
        ilut_tag const & _tag;
        //MatrixType LU;
        std::vector< std::map<unsigned int, ScalarType> > LU;
        mutable std::vector<ScalarType> temp_vec;
    };

  }
}




#endif



