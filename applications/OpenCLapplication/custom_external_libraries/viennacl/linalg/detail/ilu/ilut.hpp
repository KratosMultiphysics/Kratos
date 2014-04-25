#ifndef VIENNACL_LINALG_DETAIL_ILUT_HPP_
#define VIENNACL_LINALG_DETAIL_ILUT_HPP_

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

/** @file viennacl/linalg/detail/ilu/ilut.hpp
    @brief Implementations of an incomplete factorization preconditioner with threshold (ILUT)
*/

#include <vector>
#include <cmath>
#include <iostream>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"

#include "viennacl/linalg/detail/ilu/common.hpp"
#include "viennacl/compressed_matrix.hpp"

#include "viennacl/linalg/host_based/common.hpp"

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
        * @param entries_per_row        Number of nonzero entries per row in L and U. Note that L and U are stored in a single matrix, thus there are 2*entries_per_row in total.
        * @param drop_tolerance         The drop tolerance for ILUT
        * @param with_level_scheduling  Flag for enabling level scheduling on GPUs.
        */
        ilut_tag(unsigned int entries_per_row = 20,
                 double drop_tolerance = 1e-4,
                 bool with_level_scheduling = false) : entries_per_row_(entries_per_row), drop_tolerance_(drop_tolerance), use_level_scheduling_(with_level_scheduling) {}

        void set_drop_tolerance(double tol)
        {
          if (tol > 0)
            drop_tolerance_ = tol;
        }
        double get_drop_tolerance() const { return drop_tolerance_; }

        void set_entries_per_row(unsigned int e)
        {
          if (e > 0)
            entries_per_row_ = e;
        }

        unsigned int get_entries_per_row() const { return entries_per_row_; }

        bool use_level_scheduling() const { return use_level_scheduling_; }
        void use_level_scheduling(bool b) { use_level_scheduling_ = b; }

      private:
        unsigned int entries_per_row_;
        double drop_tolerance_;
        bool use_level_scheduling_;
    };


    /** @brief Dispatcher overload for extracting the row of nonzeros of a compressed matrix */
    template <typename ScalarType, typename SizeType, typename SparseVector>
    ScalarType setup_w(viennacl::compressed_matrix<ScalarType> const & A,
                       SizeType row,
                       SparseVector & w)
    {
      assert( (A.handle1().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );
      assert( (A.handle2().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );
      assert( (A.handle().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );

      ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(A.handle());
      unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(A.handle1());
      unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(A.handle2());

      SizeType row_i_begin = static_cast<SizeType>(row_buffer[row]);
      SizeType row_i_end   = static_cast<SizeType>(row_buffer[row+1]);
      ScalarType row_norm = 0;
      for (SizeType buf_index_i = row_i_begin; buf_index_i < row_i_end; ++buf_index_i) //Note: We do not assume that the column indices within a row are sorted
      {
        ScalarType entry = elements[buf_index_i];
        w[col_buffer[buf_index_i]] = entry;
        row_norm += entry * entry;
      }
      return std::sqrt(row_norm);
    }

    /** @brief Dispatcher overload for extracting the row of nonzeros of a STL-grown sparse matrix */
    template <typename ScalarType, typename SizeType, typename SparseVector>
    ScalarType setup_w(std::vector< std::map<SizeType, ScalarType> > const & A,
                       SizeType row,
                       SparseVector & w)
    {
      ScalarType row_norm = 0;
      w = A[row];
      for (typename std::map<SizeType, ScalarType>::const_iterator iter_w  = w.begin(); iter_w != w.end(); ++iter_w)
        row_norm += iter_w->second * iter_w->second;

      return std::sqrt(row_norm);
    }


    /** @brief Implementation of a ILU-preconditioner with threshold. Optimized implementation for compressed_matrix.
    *
    * refer to Algorithm 10.6 by Saad's book (1996 edition)
    *
    *  @param A       The input matrix. Either a compressed_matrix or of type std::vector< std::map<T, U> >
    *  @param output  The output matrix. Type requirements: const_iterator1 for iteration along rows, const_iterator2 for iteration along columns and write access via operator()
    *  @param tag     An ilut_tag in order to dispatch among several other preconditioners.
    */
    template<typename SparseMatrixType, typename ScalarType, typename SizeType>
    void precondition(SparseMatrixType const & A,
                      std::vector< std::map<SizeType, ScalarType> > & output,
                      ilut_tag const & tag)
    {
      typedef std::map<SizeType, ScalarType>          SparseVector;
      typedef typename SparseVector::iterator         SparseVectorIterator;
      typedef typename std::map<SizeType, ScalarType>::const_iterator   OutputRowConstIterator;
      typedef std::multimap<ScalarType, std::pair<SizeType, ScalarType> >  TemporarySortMap;

      assert(viennacl::traits::size1(A) == output.size() && bool("Output matrix size mismatch") );

      SparseVector w;
      TemporarySortMap temp_map;

      for (SizeType i=0; i<viennacl::traits::size1(A); ++i)  // Line 1
      {
    /*    if (i%10 == 0)
      std::cout << i << std::endl;*/

        //line 2: set up w
        ScalarType row_norm = setup_w(A, i, w);
        ScalarType tau_i = static_cast<ScalarType>(tag.get_drop_tolerance()) * row_norm;

        //line 3:
        for (SparseVectorIterator w_k = w.begin(); w_k != w.end(); ++w_k)
        {
          SizeType k = w_k->first;
          if (k >= i)
            break;

          //line 4:
          ScalarType a_kk = output[k][k];
          if (a_kk == 0)
          {
            std::cerr << "ViennaCL: FATAL ERROR in ILUT(): Diagonal entry is zero in row " << k
                      << " while processing line " << i << "!" << std::endl;
            throw "ILUT zero diagonal!";
          }

          ScalarType w_k_entry = w_k->second / a_kk;
          w_k->second = w_k_entry;

          //line 5: (dropping rule to w_k)
          if ( std::fabs(w_k_entry) > tau_i)
          {
            //line 7:
            for (OutputRowConstIterator u_k = output[k].begin(); u_k != output[k].end(); ++u_k)
            {
              if (u_k->first > k)
                w[u_k->first] -= w_k_entry * u_k->second;
            }
          }
          //else
          //  w.erase(k);

        } //for w_k

        //Line 10: Apply a dropping rule to w
        //Sort entries which are kept
        temp_map.clear();
        for (SparseVectorIterator w_k = w.begin(); w_k != w.end(); ++w_k)
        {
          SizeType k = w_k->first;
          ScalarType w_k_entry = w_k->second;

          ScalarType abs_w_k = std::fabs(w_k_entry);
          if ( (abs_w_k > tau_i) || (k == i) )//do not drop diagonal element!
          {

            if (abs_w_k == 0) // this can only happen for diagonal entry
              throw "Triangular factor in ILUT singular!";

            temp_map.insert(std::make_pair(abs_w_k, std::make_pair(k, w_k_entry)));
          }
        }

        //Lines 10-12: write the largest p values to L and U
        SizeType written_L = 0;
        SizeType written_U = 0;
        for (typename TemporarySortMap::reverse_iterator iter = temp_map.rbegin(); iter != temp_map.rend(); ++iter)
        {
          std::map<SizeType, ScalarType> & row_i = output[i];
          SizeType j = (iter->second).first;
          ScalarType w_j_entry = (iter->second).second;

          if (j < i) // Line 11: entry for L
          {
            if (written_L < tag.get_entries_per_row())
            {
              row_i[j] = w_j_entry;
              ++written_L;
            }
          }
          else if (j == i)  // Diagonal entry is always kept
          {
            row_i[j] = w_j_entry;
          }
          else //Line 12: entry for U
          {
            if (written_U < tag.get_entries_per_row())
            {
              row_i[j] = w_j_entry;
              ++written_U;
            }
          }
        }

        w.clear(); //Line 13

      } //for i
    }


    /** @brief ILUT preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class ilut_precond
    {
      typedef typename MatrixType::value_type      ScalarType;

      public:
        ilut_precond(MatrixType const & mat, ilut_tag const & tag) : tag_(tag), LU(mat.size1(), mat.size2())
        {
          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);
          //std::cout << "End CPU precond" << std::endl;
        }

        template <typename VectorType>
        void apply(VectorType & vec) const
        {
          //Note: Since vec can be a rather arbitrary vector type, we call the more generic version in the backend manually:
          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(LU.handle1());
          unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(LU.handle2());
          ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(LU.handle());

          viennacl::linalg::host_based::detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec, LU.size2(), unit_lower_tag());
          viennacl::linalg::host_based::detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec, LU.size2(), upper_tag());
        }

      private:
        void init(MatrixType const & mat)
        {
          viennacl::context host_context(viennacl::MAIN_MEMORY);
          viennacl::compressed_matrix<ScalarType> temp;
          viennacl::switch_memory_context(temp, host_context);

          viennacl::copy(mat, temp);

          std::vector< std::map<unsigned int, ScalarType> > LU_temp(mat.size1());

          viennacl::linalg::precondition(temp, LU_temp, tag_);

          viennacl::switch_memory_context(LU, host_context);
          viennacl::copy(LU_temp, LU);
        }

        ilut_tag const & tag_;
        viennacl::compressed_matrix<ScalarType> LU;
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
        ilut_precond(MatrixType const & mat, ilut_tag const & tag) : tag_(tag), LU(mat.size1(), mat.size2())
        {
          //initialize preconditioner:
          //std::cout << "Start GPU precond" << std::endl;
          init(mat);
          //std::cout << "End GPU precond" << std::endl;
        }

        void apply(vector<ScalarType> & vec) const
        {
          if (vec.handle().get_active_handle_id() != viennacl::MAIN_MEMORY)
          {
            if (tag_.use_level_scheduling())
            {
              //std::cout << "Using multifrontal on GPU..." << std::endl;
              detail::level_scheduling_substitute(vec,
                                                  multifrontal_L_row_index_arrays_,
                                                  multifrontal_L_row_buffers_,
                                                  multifrontal_L_col_buffers_,
                                                  multifrontal_L_element_buffers_,
                                                  multifrontal_L_row_elimination_num_list_);

              vec = viennacl::linalg::element_div(vec, multifrontal_U_diagonal_);

              detail::level_scheduling_substitute(vec,
                                                  multifrontal_U_row_index_arrays_,
                                                  multifrontal_U_row_buffers_,
                                                  multifrontal_U_col_buffers_,
                                                  multifrontal_U_element_buffers_,
                                                  multifrontal_U_row_elimination_num_list_);
            }
            else
            {
              viennacl::context host_context(viennacl::MAIN_MEMORY);
              viennacl::context old_context = viennacl::traits::context(vec);
              viennacl::switch_memory_context(vec, host_context);
              viennacl::linalg::inplace_solve(LU, vec, unit_lower_tag());
              viennacl::linalg::inplace_solve(LU, vec, upper_tag());
              viennacl::switch_memory_context(vec, old_context);
            }
          }
          else //apply ILUT directly:
          {
            viennacl::linalg::inplace_solve(LU, vec, unit_lower_tag());
            viennacl::linalg::inplace_solve(LU, vec, upper_tag());
          }
        }

      private:
        void init(MatrixType const & mat)
        {
          viennacl::context host_context(viennacl::MAIN_MEMORY);
          viennacl::switch_memory_context(LU, host_context);

          std::vector< std::map<unsigned int, ScalarType> > LU_temp(mat.size1());

          if (viennacl::traits::context(mat).memory_type() == viennacl::MAIN_MEMORY)
          {
            viennacl::linalg::precondition(mat, LU_temp, tag_);
          }
          else //we need to copy to CPU
          {
            viennacl::compressed_matrix<ScalarType> cpu_mat(mat.size1(), mat.size2());
            viennacl::switch_memory_context(cpu_mat, host_context);

            cpu_mat = mat;

            viennacl::linalg::precondition(cpu_mat, LU_temp, tag_);
          }

          viennacl::copy(LU_temp, LU);

          if (!tag_.use_level_scheduling())
            return;

          //
          // multifrontal part:
          //

          viennacl::switch_memory_context(multifrontal_U_diagonal_, host_context);
          multifrontal_U_diagonal_.resize(LU.size1(), false);
          host_based::detail::row_info(LU, multifrontal_U_diagonal_, viennacl::linalg::detail::SPARSE_ROW_DIAGONAL);

          detail::level_scheduling_setup_L(LU,
                                           multifrontal_U_diagonal_, //dummy
                                           multifrontal_L_row_index_arrays_,
                                           multifrontal_L_row_buffers_,
                                           multifrontal_L_col_buffers_,
                                           multifrontal_L_element_buffers_,
                                           multifrontal_L_row_elimination_num_list_);


          detail::level_scheduling_setup_U(LU,
                                           multifrontal_U_diagonal_,
                                           multifrontal_U_row_index_arrays_,
                                           multifrontal_U_row_buffers_,
                                           multifrontal_U_col_buffers_,
                                           multifrontal_U_element_buffers_,
                                           multifrontal_U_row_elimination_num_list_);

          //
          // Bring to device if necessary:
          //

          // L:

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_L_row_index_arrays_.begin();
                                                                             it != multifrontal_L_row_index_arrays_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_L_row_buffers_.begin();
                                                                             it != multifrontal_L_row_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_L_col_buffers_.begin();
                                                                             it != multifrontal_L_col_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_L_element_buffers_.begin();
                                                                             it != multifrontal_L_element_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<ScalarType>(*it, viennacl::traits::context(mat));


          // U:

          viennacl::switch_memory_context(multifrontal_U_diagonal_, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_U_row_index_arrays_.begin();
                                                                             it != multifrontal_U_row_index_arrays_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_U_row_buffers_.begin();
                                                                             it != multifrontal_U_row_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_U_col_buffers_.begin();
                                                                             it != multifrontal_U_col_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<unsigned int>(*it, viennacl::traits::context(mat));

          for (typename std::list< viennacl::backend::mem_handle >::iterator it  = multifrontal_U_element_buffers_.begin();
                                                                             it != multifrontal_U_element_buffers_.end();
                                                                           ++it)
            viennacl::backend::switch_memory_context<ScalarType>(*it, viennacl::traits::context(mat));


        }

        ilut_tag const & tag_;
        viennacl::compressed_matrix<ScalarType> LU;

        std::list< viennacl::backend::mem_handle > multifrontal_L_row_index_arrays_;
        std::list< viennacl::backend::mem_handle > multifrontal_L_row_buffers_;
        std::list< viennacl::backend::mem_handle > multifrontal_L_col_buffers_;
        std::list< viennacl::backend::mem_handle > multifrontal_L_element_buffers_;
        std::list< vcl_size_t > multifrontal_L_row_elimination_num_list_;

        viennacl::vector<ScalarType> multifrontal_U_diagonal_;
        std::list< viennacl::backend::mem_handle > multifrontal_U_row_index_arrays_;
        std::list< viennacl::backend::mem_handle > multifrontal_U_row_buffers_;
        std::list< viennacl::backend::mem_handle > multifrontal_U_col_buffers_;
        std::list< viennacl::backend::mem_handle > multifrontal_U_element_buffers_;
        std::list< vcl_size_t > multifrontal_U_row_elimination_num_list_;
    };

  }
}




#endif



