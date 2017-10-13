#ifndef VIENNACL_LINALG_DETAIL_OPENCL_BLOCK_ILU_HPP_
#define VIENNACL_LINALG_DETAIL_OPENCL_BLOCK_ILU_HPP_

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

/** @file viennacl/linalg/detail/ilu/opencl_block_ilu.hpp
    @brief Implementations of incomplete block factorization preconditioners using OpenCL
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/detail/ilu/common.hpp"
#include "viennacl/linalg/detail/ilu/ilu0.hpp"
#include "viennacl/linalg/detail/ilu/ilut.hpp"
#include "viennacl/linalg/detail/ilu/host_block_ilu.hpp"

#include "viennacl/linalg/kernels/ilu_kernels.h"

#include <map>

namespace viennacl
{
  namespace linalg
  {
    
    /** @brief ILUT preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT, typename ILUTag>
    class block_ilu_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT>, ILUTag >
    {
        typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>        MatrixType;
        typedef std::vector< std::map<unsigned int, ScalarType> >   InternalMatrixType;
        typedef std::vector<ScalarType>                             STLVectorType;
      
      public:
        typedef std::vector<std::pair<std::size_t, std::size_t> >    index_vector_type;   //the pair refers to index range [a, b) of each block
          
        
        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          std::size_t num_blocks = 4
                         ) : tag_(tag), block_indices_(num_blocks),
                             gpu_block_indices(0),
                             gpu_L_trans(0,0),
                             gpu_U_trans(0,0),
                             gpu_D(mat.size1()),
                             LU_blocks(num_blocks),
                             temp_vec_(mat.size1())
        {
          viennacl::linalg::kernels::ilu<ScalarType, 1>::init();
          
          // Set up vector of block indices:
          block_indices_.resize(num_blocks);
          for (std::size_t i=0; i<num_blocks; ++i)
          {
            std::size_t start_index = (   i  * mat.size1()) / num_blocks;
            std::size_t stop_index  = ((i+1) * mat.size1()) / num_blocks;
            
            block_indices_[i] = std::pair<std::size_t, std::size_t>(start_index, stop_index);
          }
          
          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);          
          //std::cout << "End CPU precond" << std::endl;
        }

        block_ilu_precond(MatrixType const & mat,
                          ILUTag const & tag,
                          index_vector_type const & block_boundaries
                         ) : tag_(tag), block_indices_(block_boundaries),
                             gpu_block_indices(0),
                             gpu_L_trans(0,0),
                             gpu_U_trans(0,0),
                             gpu_D(mat.size1()),
                             LU_blocks(block_boundaries.size()),
                             temp_vec_(mat.size1())
        {
          viennacl::linalg::kernels::ilu<ScalarType, 1>::init();
          
          //initialize preconditioner:
          //std::cout << "Start CPU precond" << std::endl;
          init(mat);          
          //std::cout << "End CPU precond" << std::endl;
        }
        
        
        void apply(vector<ScalarType> & vec) const
        {
          apply_cpu(vec);
        }

        // GPU application (default)
        void apply_gpu(vector<ScalarType> & vec) const
        {
          viennacl::ocl::kernel & block_ilut_kernel =
               viennacl::ocl::get_kernel(viennacl::linalg::kernels::ilu<ScalarType, 1>::program_name(), "block_ilu_substitute");

          viennacl::ocl::enqueue(block_ilut_kernel(gpu_L_trans.handle1(),  //L
                                                   gpu_L_trans.handle2(),
                                                   gpu_L_trans.handle(),
                                                   gpu_U_trans.handle1(),  //U
                                                   gpu_U_trans.handle2(),
                                                   gpu_U_trans.handle(),
                                                   gpu_D,                  //D
                                                   gpu_block_indices,
                                                   vec,
                                                   static_cast<cl_uint>(vec.size())));
        }

        // CPU fallback:
        void apply_cpu(vector<ScalarType> & vec) const
        {
          viennacl::copy(vec, temp_vec_);
          
          for (std::size_t i=0; i<block_indices_.size(); ++i)
          {
            viennacl::tools::const_sparse_matrix_adapter<ScalarType> LU_const_adapter(LU_blocks[i],
                                                                                      LU_blocks[i].size(),
                                                                                      LU_blocks[i].size());
            detail::ilu_vector_range<STLVectorType>  vec_range(temp_vec_,
                                                            block_indices_[i].first,
                                                            LU_blocks[i].size());
            viennacl::linalg::detail::ilu_lu_substitute(LU_const_adapter, vec_range);
          }
                    
          viennacl::copy(temp_vec_, vec);
        }
        
      private:
        void init(MatrixType const & mat)
        {
          InternalMatrixType temp(mat.size1());
          //std::vector< std::map<unsigned int, ScalarType> > LU_cpu(mat.size1());

          //copy to cpu:
          viennacl::copy(mat, temp);
          
          for (std::size_t i=0; i<block_indices_.size(); ++i)
          {
            // Step 1: Extract blocks
            std::size_t block_size = block_indices_[i].second - block_indices_[i].first;
            InternalMatrixType mat_block(block_size);
            viennacl::tools::const_sparse_matrix_adapter<ScalarType>  temp_adapter(temp, temp.size(), temp.size());
            detail::extract_block_matrix(temp_adapter, mat_block, block_indices_[i].first, block_indices_[i].second);
            
            
            // Step 2: Precondition blocks:
            viennacl::tools::const_sparse_matrix_adapter<ScalarType>  mat_block_adapter(mat_block, block_size, block_size);
            viennacl::tools::sparse_matrix_adapter<ScalarType>        LU_adapter(LU_blocks[i], block_size, block_size);
            viennacl::linalg::precondition(mat_block_adapter, LU_adapter, tag_);
          }
          
          /*
           * copy resulting preconditioner back to GPU:
           */
          
          std::vector<cl_uint> block_indices_uint(2 * block_indices_.size());
          for (std::size_t i=0; i<block_indices_.size(); ++i)
          {
            block_indices_uint[2*i]     = static_cast<cl_uint>(block_indices_[i].first);
            block_indices_uint[2*i + 1] = static_cast<cl_uint>(block_indices_[i].second);
          }

          gpu_block_indices = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                             sizeof(cl_uint) * block_indices_uint.size(),
                                                                             &(block_indices_uint[0]) );
          
          blocks_to_GPU(mat.size1());

          //set kernel parameters:
          viennacl::ocl::kernel & block_ilut_kernel =
               viennacl::ocl::get_kernel(viennacl::linalg::kernels::ilu<ScalarType, 1>::program_name(), "block_ilu_substitute");
          
          block_ilut_kernel.global_work_size(0, 128 * block_indices_.size() );
          block_ilut_kernel.local_work_size(0, 128);
        }
        
        // Copy computed preconditioned blocks to OpenCL device
        void blocks_to_GPU(std::size_t matrix_size)
        {
          std::vector< std::map<cl_uint, ScalarType> > L_transposed(matrix_size);
          std::vector< std::map<cl_uint, ScalarType> > U_transposed(matrix_size);
          std::vector<ScalarType> entries_D(matrix_size);
          
          //
          // Transpose individual blocks into a single large matrix:
          //
          for (std::size_t block_index = 0; block_index < LU_blocks.size(); ++block_index)
          {
            InternalMatrixType const & current_block = LU_blocks[block_index];
            
            std::size_t block_start = block_indices_[block_index].first;
            
            //transpose L and U:
            for (std::size_t row = 0; row < current_block.size(); ++row)
            {
              for (typename InternalMatrixType::value_type::const_iterator col_iter  = current_block[row].begin();
                                                                           col_iter != current_block[row].end();
                                                                         ++col_iter)
              {
                if (row > col_iter->first) //entry for L
                  L_transposed[col_iter->first + block_start][row + block_start] = col_iter->second;
                else if (row == col_iter->first)
                  entries_D[row + block_start] = col_iter->second;
                else //entry for U
                  U_transposed[col_iter->first + block_start][row + block_start] = col_iter->second;
              }
            }
          }
          
          //
          // Move data to GPU:
          //
          viennacl::copy(L_transposed, gpu_L_trans);
          viennacl::copy(U_transposed, gpu_U_trans);
          viennacl::copy(entries_D, gpu_D);
        }
        
        
        ILUTag const & tag_;
        index_vector_type block_indices_;
        viennacl::ocl::handle<cl_mem> gpu_block_indices;
        viennacl::compressed_matrix<ScalarType> gpu_L_trans;
        viennacl::compressed_matrix<ScalarType> gpu_U_trans;
        viennacl::vector<ScalarType> gpu_D;
        
        std::vector< InternalMatrixType > LU_blocks;
        mutable STLVectorType temp_vec_;
    };

  }
}




#endif



