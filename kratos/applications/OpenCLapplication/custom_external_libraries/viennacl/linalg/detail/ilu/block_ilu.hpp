#ifndef VIENNACL_LINALG_DETAIL_BLOCK_ILU_HPP_
#define VIENNACL_LINALG_DETAIL_BLOCK_ILU_HPP_

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

/** @file viennacl/linalg/detail/ilu/block_ilu.hpp
    @brief Implementations of incomplete block factorization preconditioners
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/detail/ilu/common.hpp"
#include "viennacl/linalg/detail/ilu/ilu0.hpp"
#include "viennacl/linalg/detail/ilu/ilut.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      /** @brief Helper range class for representing a subvector of a larger buffer. */
      template <typename VectorType, typename ValueType, typename SizeType = vcl_size_t>
      class ilu_vector_range
      {
        public:
          //typedef typename VectorType::value_type      value_type;
          //typedef typename VectorType::size_type       size_type;

          ilu_vector_range(VectorType & v,
                           SizeType start_index,
                           SizeType vec_size
                          ) : vec_(v), start_(start_index), size_(vec_size) {}

          ValueType & operator()(SizeType index)
          {
            assert(index < size_ && bool("Index out of bounds!"));
            return vec_[start_ + index];
          }

          ValueType & operator[](SizeType index)
          {
            assert(index < size_ && bool("Index out of bounds!"));
            return vec_[start_ + index];
          }

          SizeType size() const { return size_; }

        private:
          VectorType & vec_;
          SizeType start_;
          SizeType size_;
      };

      /** @brief Extracts a diagonal block from a larger system matrix
        *
        * @param A                   The full matrix
        * @param diagonal_block_A    The output matrix, to which the extracted block is written to
        * @param start_index         First row- and column-index of the block
        * @param stop_index          First row- and column-index beyond the block
        */
      template <typename ScalarType>
      void extract_block_matrix(viennacl::compressed_matrix<ScalarType> const & A,
                                viennacl::compressed_matrix<ScalarType> & diagonal_block_A,
                                vcl_size_t start_index,
                                vcl_size_t stop_index
                                )
      {

        assert( (A.handle1().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );
        assert( (A.handle2().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );
        assert( (A.handle().get_active_handle_id() == viennacl::MAIN_MEMORY) && bool("System matrix must reside in main memory for ILU0") );

        ScalarType   const * A_elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(A.handle());
        unsigned int const * A_row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(A.handle1());
        unsigned int const * A_col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(A.handle2());

        ScalarType   * output_elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(diagonal_block_A.handle());
        unsigned int * output_row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(diagonal_block_A.handle1());
        unsigned int * output_col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(diagonal_block_A.handle2());

        vcl_size_t output_counter = 0;
        for (vcl_size_t row = start_index; row < stop_index; ++row)
        {
          unsigned int buffer_col_start = A_row_buffer[row];
          unsigned int buffer_col_end   = A_row_buffer[row+1];

          output_row_buffer[row - start_index] = static_cast<unsigned int>(output_counter);

          for (unsigned int buf_index = buffer_col_start; buf_index < buffer_col_end; ++buf_index)
          {
            unsigned int col = A_col_buffer[buf_index];
            if (col < start_index)
              continue;

            if (col >= static_cast<unsigned int>(stop_index))
              continue;

            output_col_buffer[output_counter] = static_cast<unsigned int>(col - start_index);
            output_elements[output_counter] = A_elements[buf_index];
            ++output_counter;
          }
          output_row_buffer[row - start_index + 1] = static_cast<unsigned int>(output_counter);
        }
      }


    }

    /** @brief A block ILU preconditioner class, can be supplied to solve()-routines
     *
     * @tparam MatrixType   Type of the system matrix
     * @tparam ILUTag       Type of the tag identifiying the ILU preconditioner to be used on each block.
    */
    template <typename MatrixType, typename ILUTag>
    class block_ilu_precond
    {
      typedef typename MatrixType::value_type      ScalarType;

      public:
        typedef std::vector<std::pair<vcl_size_t, vcl_size_t> >    index_vector_type;   //the pair refers to index range [a, b) of each block


        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          vcl_size_t num_blocks = 8
                         ) : tag_(tag), LU_blocks(num_blocks)
        {

          // Set up vector of block indices:
          block_indices_.resize(num_blocks);
          for (vcl_size_t i=0; i<num_blocks; ++i)
          {
            vcl_size_t start_index = (   i  * mat.size1()) / num_blocks;
            vcl_size_t stop_index  = ((i+1) * mat.size1()) / num_blocks;

            block_indices_[i] = std::pair<vcl_size_t, vcl_size_t>(start_index, stop_index);
          }

          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);
          //std::cout << "End CPU precond" << std::endl;
        }

        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          index_vector_type const & block_boundaries
                         ) : tag_(tag), block_indices_(block_boundaries), LU_blocks(block_boundaries.size())
        {
          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);
          //std::cout << "End CPU precond" << std::endl;
        }


        template <typename VectorType>
        void apply(VectorType & vec) const
        {
          for (vcl_size_t i=0; i<block_indices_.size(); ++i)
          {
            detail::ilu_vector_range<VectorType, ScalarType>  vec_range(vec, block_indices_[i].first, LU_blocks[i].size2());

            unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(LU_blocks[i].handle1());
            unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(LU_blocks[i].handle2());
            ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(LU_blocks[i].handle());

            viennacl::linalg::host_based::detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_range, LU_blocks[i].size2(), unit_lower_tag());
            viennacl::linalg::host_based::detail::csr_inplace_solve<ScalarType>(row_buffer, col_buffer, elements, vec_range, LU_blocks[i].size2(), upper_tag());

          }
        }

      private:
        void init(MatrixType const & A)
        {
          viennacl::context host_context(viennacl::MAIN_MEMORY);
          viennacl::compressed_matrix<ScalarType> mat(host_context);

          viennacl::copy(A, mat);

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(mat.handle1());

#ifdef VIENNACL_WITH_OPENMP
          #pragma omp parallel for
#endif
          for (long i=0; i<static_cast<long>(block_indices_.size()); ++i)
          {
            // Step 1: Extract blocks
            vcl_size_t block_size = block_indices_[i].second - block_indices_[i].first;
            vcl_size_t block_nnz  = row_buffer[block_indices_[i].second] - row_buffer[block_indices_[i].first];
            viennacl::compressed_matrix<ScalarType> mat_block(block_size, block_size, block_nnz, host_context);

            detail::extract_block_matrix(mat, mat_block, block_indices_[i].first, block_indices_[i].second);

            // Step 2: Precondition blocks:
            viennacl::switch_memory_context(LU_blocks[i], host_context);
            preconditioner_dispatch(mat_block, LU_blocks[i], tag_);
          }

        }

        void preconditioner_dispatch(viennacl::compressed_matrix<ScalarType> const & mat_block,
                                     viennacl::compressed_matrix<ScalarType> & LU,
                                     viennacl::linalg::ilu0_tag)
        {
          LU = mat_block;
          viennacl::linalg::precondition(LU, tag_);
        }

        void preconditioner_dispatch(viennacl::compressed_matrix<ScalarType> const & mat_block,
                                     viennacl::compressed_matrix<ScalarType> & LU,
                                     viennacl::linalg::ilut_tag)
        {
          std::vector< std::map<unsigned int, ScalarType> > temp(mat_block.size1());

          viennacl::linalg::precondition(mat_block, temp, tag_);

          viennacl::copy(temp, LU);
        }

        ILUTag const & tag_;
        index_vector_type block_indices_;
        std::vector< viennacl::compressed_matrix<ScalarType> > LU_blocks;
    };





    /** @brief ILUT preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT, typename ILUTag>
    class block_ilu_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT>, ILUTag >
    {
        typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>        MatrixType;
        //typedef std::vector<ScalarType>                             STLVectorType;

      public:
        typedef std::vector<std::pair<vcl_size_t, vcl_size_t> >    index_vector_type;   //the pair refers to index range [a, b) of each block


        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          vcl_size_t num_blocks = 8
                         ) : tag_(tag),
                             block_indices_(num_blocks),
                             gpu_block_indices(),
                             gpu_L_trans(0,0, viennacl::traits::context(mat)),
                             gpu_U_trans(0,0, viennacl::traits::context(mat)),
                             gpu_D(mat.size1(), viennacl::traits::context(mat)),
                             LU_blocks(num_blocks)
        {
          // Set up vector of block indices:
          block_indices_.resize(num_blocks);
          for (vcl_size_t i=0; i<num_blocks; ++i)
          {
            vcl_size_t start_index = (   i  * mat.size1()) / num_blocks;
            vcl_size_t stop_index  = ((i+1) * mat.size1()) / num_blocks;

            block_indices_[i] = std::pair<vcl_size_t, vcl_size_t>(start_index, stop_index);
          }

          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);
          //std::cout << "End CPU precond" << std::endl;
        }

        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          index_vector_type const & block_boundaries
                         ) : tag_(tag),
                             block_indices_(block_boundaries),
                             gpu_block_indices(viennacl::traits::context(mat)),
                             gpu_L_trans(0,0,viennacl::traits::context(mat)),
                             gpu_U_trans(0,0,viennacl::traits::context(mat)),
                             gpu_D(0,viennacl::traits::context(mat)),
                             LU_blocks(block_boundaries.size())
        {
          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);
          //std::cout << "End CPU precond" << std::endl;
        }


        void apply(vector<ScalarType> & vec) const
        {
          viennacl::linalg::detail::block_inplace_solve(trans(gpu_L_trans), gpu_block_indices, block_indices_.size(), gpu_D,
                                                        vec,
                                                        viennacl::linalg::unit_lower_tag());

          viennacl::linalg::detail::block_inplace_solve(trans(gpu_U_trans), gpu_block_indices, block_indices_.size(), gpu_D,
                                                        vec,
                                                        viennacl::linalg::upper_tag());

          //apply_cpu(vec);
        }


      private:

        void init(MatrixType const & A)
        {
          viennacl::context host_context(viennacl::MAIN_MEMORY);
          viennacl::compressed_matrix<ScalarType> mat(host_context);

          mat = A;

          unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(mat.handle1());

#ifdef VIENNACL_WITH_OPENMP
          #pragma omp parallel for
#endif
          for (long i=0; i<static_cast<long>(block_indices_.size()); ++i)
          {
            // Step 1: Extract blocks
            vcl_size_t block_size = block_indices_[i].second - block_indices_[i].first;
            vcl_size_t block_nnz  = row_buffer[block_indices_[i].second] - row_buffer[block_indices_[i].first];
            viennacl::compressed_matrix<ScalarType> mat_block(block_size, block_size, block_nnz, host_context);

            detail::extract_block_matrix(mat, mat_block, block_indices_[i].first, block_indices_[i].second);

            // Step 2: Precondition blocks:
            viennacl::switch_memory_context(LU_blocks[i], host_context);
            preconditioner_dispatch(mat_block, LU_blocks[i], tag_);
          }

          /*
           * copy resulting preconditioner back to GPU:
           */

          viennacl::switch_memory_context(gpu_L_trans, viennacl::traits::context(A));
          viennacl::switch_memory_context(gpu_U_trans, viennacl::traits::context(A));
          viennacl::switch_memory_context(gpu_D, viennacl::traits::context(A));

          viennacl::backend::typesafe_host_array<unsigned int> block_indices_uint(gpu_block_indices, 2 * block_indices_.size());
          for (vcl_size_t i=0; i<block_indices_.size(); ++i)
          {
            block_indices_uint.set(2*i, block_indices_[i].first);
            block_indices_uint.set(2*i + 1, block_indices_[i].second);
          }

          viennacl::backend::memory_create(gpu_block_indices, block_indices_uint.raw_size(), viennacl::traits::context(A), block_indices_uint.get());

          blocks_to_device(mat.size1());

        }

        // Copy computed preconditioned blocks to OpenCL device
        void blocks_to_device(vcl_size_t matrix_size)
        {
          std::vector< std::map<unsigned int, ScalarType> > L_transposed(matrix_size);
          std::vector< std::map<unsigned int, ScalarType> > U_transposed(matrix_size);
          std::vector<ScalarType> entries_D(matrix_size);

          //
          // Transpose individual blocks into a single large matrix:
          //
          for (vcl_size_t block_index = 0; block_index < LU_blocks.size(); ++block_index)
          {
            MatrixType const & current_block = LU_blocks[block_index];

            unsigned int const * row_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(current_block.handle1());
            unsigned int const * col_buffer = viennacl::linalg::host_based::detail::extract_raw_pointer<unsigned int>(current_block.handle2());
            ScalarType   const * elements   = viennacl::linalg::host_based::detail::extract_raw_pointer<ScalarType>(current_block.handle());

            vcl_size_t block_start = block_indices_[block_index].first;

            //transpose L and U:
            for (vcl_size_t row = 0; row < current_block.size1(); ++row)
            {
              unsigned int buffer_col_start = row_buffer[row];
              unsigned int buffer_col_end   = row_buffer[row+1];

              for (unsigned int buf_index = buffer_col_start; buf_index < buffer_col_end; ++buf_index)
              {
                unsigned int col = col_buffer[buf_index];

                if (row > col) //entry for L
                  L_transposed[col + block_start][static_cast<unsigned int>(row + block_start)] = elements[buf_index];
                else if (row == col)
                  entries_D[row + block_start] = elements[buf_index];
                else //entry for U
                  U_transposed[col + block_start][static_cast<unsigned int>(row + block_start)] = elements[buf_index];
              }
            }
          }

          //
          // Move data to GPU:
          //
          tools::const_sparse_matrix_adapter<ScalarType, unsigned int> adapted_L_transposed(L_transposed, matrix_size, matrix_size);
          tools::const_sparse_matrix_adapter<ScalarType, unsigned int> adapted_U_transposed(U_transposed, matrix_size, matrix_size);
          viennacl::copy(adapted_L_transposed, gpu_L_trans);
          viennacl::copy(adapted_U_transposed, gpu_U_trans);
          viennacl::copy(entries_D, gpu_D);
        }

        void preconditioner_dispatch(viennacl::compressed_matrix<ScalarType> const & mat_block,
                                     viennacl::compressed_matrix<ScalarType> & LU,
                                     viennacl::linalg::ilu0_tag)
        {
          LU = mat_block;
          viennacl::linalg::precondition(LU, tag_);
        }

        void preconditioner_dispatch(viennacl::compressed_matrix<ScalarType> const & mat_block,
                                     viennacl::compressed_matrix<ScalarType> & LU,
                                     viennacl::linalg::ilut_tag)
        {
          std::vector< std::map<unsigned int, ScalarType> > temp(mat_block.size1());

          viennacl::linalg::precondition(mat_block, temp, tag_);

          viennacl::copy(temp, LU);
        }


        ILUTag const & tag_;
        index_vector_type block_indices_;
        viennacl::backend::mem_handle gpu_block_indices;
        viennacl::compressed_matrix<ScalarType> gpu_L_trans;
        viennacl::compressed_matrix<ScalarType> gpu_U_trans;
        viennacl::vector<ScalarType> gpu_D;

        std::vector< MatrixType > LU_blocks;
    };


  }
}




#endif



