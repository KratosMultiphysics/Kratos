#ifndef VIENNACL_LINALG_DETAIL_SPAI_SPAI_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_SPAI_HPP

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

/** @file viennacl/linalg/detail/spai/spai.hpp
    @brief Main implementation of SPAI (not FSPAI)
*/

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>

//local includes
#include "viennacl/linalg/detail/spai/spai_tag.hpp"
#include "viennacl/linalg/qr.hpp"
#include "viennacl/linalg/detail/spai/spai-dynamic.hpp"
#include "viennacl/linalg/detail/spai/spai-static.hpp"
#include "viennacl/linalg/detail/spai/sparse_vector.hpp"
#include "viennacl/linalg/detail/spai/block_matrix.hpp"
#include "viennacl/linalg/detail/spai/block_vector.hpp"

//boost includes
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/storage.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/lu.hpp"
#include "boost/numeric/ublas/triangular.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"

// ViennaCL includes
#include "viennacl/linalg/prod.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/compressed_matrix_operations.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/kernels/spai_source.h"
#include "viennacl/linalg/kernels/spai_kernels.h"



#define VIENNACL_SPAI_K_b 20

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      namespace spai
      {
        
        //debug function for print
        template<typename SparseVectorType>
        void print_sparse_vector(const SparseVectorType& v){
            for(typename SparseVectorType::const_iterator vec_it = v.begin(); vec_it!= v.end(); ++vec_it){
                std::cout<<"[ "<<vec_it->first<<" ]:"<<vec_it->second<<std::endl;
            }
        }
        template<typename DenseMatrixType>
        void print_matrix( DenseMatrixType & m){
            for(int i = 0; i < m.size2(); ++i){
                for(int j = 0; j < m.size1(); ++j){
                    std::cout<<m(j, i)<<" ";
                }
                std::cout<<std::endl;
            }
        }
        
        /** @brief Add two sparse vectors res_v = b*v
         * @param v initial sparse vector
         * @param b scalar
         * @param res_v output vector 
         */ 
        template<typename SparseVectorType, typename ScalarType>
        void add_sparse_vectors(const SparseVectorType& v, const ScalarType b,  SparseVectorType& res_v){
            for(typename SparseVectorType::const_iterator v_it = v.begin(); v_it != v.end(); ++v_it){
                res_v[v_it->first] += b*v_it->second;
            }
        }
        //sparse-matrix - vector product
        /** @brief Computation of residual res = A*v - e
         * @param A_v_c column major vectorized input sparse matrix
         * @param v sparse vector, in this case new column of preconditioner matrix 
         * @param ind index for current column
         * @param res residual
         */ 
        template<typename SparseVectorType, typename ScalarType>
        void compute_spai_residual(const std::vector<SparseVectorType>& A_v_c, const SparseVectorType& v,
                                   const unsigned int ind, SparseVectorType& res){
            for(typename SparseVectorType::const_iterator v_it = v.begin(); v_it != v.end(); ++v_it){
                add_sparse_vectors(A_v_c[v_it->first], v_it->second, res);
            }
            res[ind] -= static_cast<ScalarType>(1);
        }
        
        /** @brief Setting up index set of columns and rows for certain column
         * @param A_v_c column major vectorized initial sparse matrix  
         * @param v current column of preconditioner matrix
         * @param J set of column indices
         * @param I set of row indices
         */ 
        template<typename SparseVectorType>
        void build_index_set(const std::vector<SparseVectorType>& A_v_c, const SparseVectorType& v, std::vector<unsigned int>& J, 
                             std::vector<unsigned int>& I){
            buildColumnIndexSet(v, J);
            projectRows(A_v_c, J, I);
        }
        
        /************************************************** CPU BLOCK SET UP ***************************************/
        /** @brief Setting up blocks and QR factorizing them on CPU
         * @param A initial sparse matrix 
         * @param A_v_c column major vectorized initial sparse matrix  
         * @param M_v initialized preconditioner
         * @param g_I container of row indices
         * @param g_J container of column indices
         * @param g_A_I_J container of dense matrices -> R matrices after QR factorization
         * @param g_b_v container of vectors beta, necessary for Q recovery
         */ 
        template<typename SparseMatrixType, typename DenseMatrixType, typename SparseVectorType, typename VectorType>
        void block_set_up(const SparseMatrixType& A,
                          const std::vector<SparseVectorType>& A_v_c,
                          const std::vector<SparseVectorType>& M_v, 
                          std::vector<std::vector<unsigned int> >& g_I, 
                          std::vector<std::vector<unsigned int> >& g_J,
                          std::vector<DenseMatrixType>& g_A_I_J,
                          std::vector<VectorType>& g_b_v){
#ifdef _OPENMP
            #pragma omp parallel for
#endif            
            for(std::size_t i = 0; i < M_v.size(); ++i){
                build_index_set(A_v_c, M_v[i], g_J[i], g_I[i]);
                initProjectSubMatrix(A, g_J[i], g_I[i], g_A_I_J[i]);
                //print_matrix(g_A_I_J[i]);
                single_qr(g_A_I_J[i], g_b_v[i]);
                //print_matrix(g_A_I_J[i]);
            }            
        }
        
        /** @brief Setting up index set of columns and rows for all columns
         * @param A_v_c column major vectorized initial sparse matrix  
         * @param M_v initialized preconditioner
         * @param g_J container of column indices
         * @param g_I container of row indices
         */ 
        template<typename SparseVectorType>
        void index_set_up(const std::vector<SparseVectorType>& A_v_c,
                          const std::vector<SparseVectorType>& M_v, 
                          std::vector<std::vector<unsigned int> >& g_J, 
                          std::vector<std::vector<unsigned int> >& g_I)
        {
#ifdef _OPENMP
            #pragma omp parallel for
#endif            
            for(std::size_t i = 0; i < M_v.size(); ++i){
                build_index_set(A_v_c, M_v[i], g_J[i], g_I[i]);
            }
        }
        
        /************************************************** GPU BLOCK SET UP ***************************************/
        /** @brief Setting up blocks and QR factorizing them on GPU
         * @param A initial sparse matrix 
         * @param A_v_c column major vectorized initial sparse matrix  
         * @param M_v initialized preconditioner
         * @param g_is_update container that indicates which blocks are active
         * @param g_I container of row indices
         * @param g_J container of column indices
         * @param g_A_I_J container of dense matrices -> R matrices after QR factorization
         * @param g_bv container of vectors beta, necessary for Q recovery
         * @param cur_iter number of current iteration
         */ 
        template<typename ScalarType, unsigned int MAT_ALIGNMENT, typename SparseVectorType>
        void block_set_up(const viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>& A,
                          const std::vector<SparseVectorType>& A_v_c, 
                          const std::vector<SparseVectorType>& M_v,
                          std::vector<cl_uint> g_is_update,
                          std::vector<std::vector<unsigned int> >& g_I, 
                          std::vector<std::vector<unsigned int> >& g_J, 
                          block_matrix & g_A_I_J, 
                          block_vector & g_bv,
                          const unsigned int cur_iter)
        {
            bool is_empty_block;
            //build index set
            index_set_up(A_v_c, M_v, g_J, g_I);
            block_assembly(A, g_J, g_I, g_A_I_J, g_is_update, is_empty_block, cur_iter);
            block_qr<ScalarType>(g_I, g_J, g_A_I_J, g_bv, g_is_update, cur_iter);
            
        }
        
        
        /***************************************************************************************************/
        /******************************** SOLVING LS PROBLEMS ON GPU ***************************************/
        /***************************************************************************************************/
        /** @brief Elicitation of sparse vector m for particular column from m_in - contigious vector for all columns  
         * @param m_in contigious sparse vector for all columns
         * @param start_m_ind start index of particular vector   
         * @param J column index set
         * @param m sparse vector for particular column
         */ 
        template<typename ScalarType, typename SparseVectorType>
        void custom_fan_out(const std::vector<ScalarType> & m_in,
                            unsigned int start_m_ind,
                            const std::vector<unsigned int> & J,
                            SparseVectorType & m)
        {
            unsigned int  cnt = 0;
            for (std::size_t i = 0; i < J.size(); ++i) {
                m[J[i]] = m_in[start_m_ind + cnt++];
            }
        }
        
        

        //GPU based least square problem
        /** @brief Solution of Least square problem on GPU 
         * @param A_v_c column-major vectorized initial sparse matrix 
         * @param M_v column-major vectorized sparse preconditioner matrix
         * @param g_I container of row set indices
         * @param g_J container of column set indices
         * @param g_A_I_J_vcl contigious matrix that consists of blocks A(I_k, J_k)
         * @param g_bv_vcl contigious vector that consists of betas, necessary for Q recovery
         * @param g_res container of residuals 
         * @param g_is_update container with indicators which blocks are active
         * @param tag spai tag
         * @param cur_iter current iteration number
         */ 
        template<typename SparseVectorType, typename ScalarType>
        void least_square_solve(std::vector<SparseVectorType> & A_v_c, 
                                std::vector<SparseVectorType> & M_v,
                                std::vector<std::vector<unsigned int> >& g_I, 
                                std::vector<std::vector<unsigned int> > & g_J,
                                block_matrix & g_A_I_J_vcl,
                                block_vector & g_bv_vcl,
                                std::vector<SparseVectorType> & g_res,
                                std::vector<cl_uint> & g_is_update,
                                const spai_tag & tag,
                                const unsigned int cur_iter){
            unsigned int y_sz, m_sz;
            std::vector<cl_uint> y_inds(M_v.size() + 1, static_cast<cl_uint>(0));
            std::vector<cl_uint> m_inds(M_v.size() + 1, static_cast<cl_uint>(0));
            get_size(g_I, y_sz);
            init_start_inds(g_I, y_inds);
            init_start_inds(g_J, m_inds);
            //create y_v
            std::vector<ScalarType> y_v(y_sz, static_cast<ScalarType>(0));
            for(std::size_t i = 0; i < M_v.size(); ++i){
                for(std::size_t j = 0; j < g_I[i].size(); ++j){
                    if(g_I[i][j] == i)
                        y_v[y_inds[i] + j] = static_cast<ScalarType>(1.0);
                }
            }
            //compute m_v
            get_size(g_J, m_sz);
            std::vector<ScalarType> m_v(m_sz, static_cast<cl_uint>(0));
            
            //acquire kernel
            /*if(cur_iter == 0){
                std::string ls_kernel_file_name = "kernels/spai/ls_g.cl";
                std::string ls_kernel_source;
                read_kernel_from_file(ls_kernel_file_name, ls_kernel_source);
                //compilation of a kernel
                viennacl::ocl::program & ls_prog = viennacl::ocl::current_context().add_program(ls_kernel_source.c_str(), "ls_kernel_source");
                //least square kernel
                ls_prog.add_kernel("block_least_squares");
            }*/
            block_vector y_v_vcl;
            block_vector m_v_vcl;
            //prepearing memory for least square problem on GPU
            y_v_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                              static_cast<unsigned int>(sizeof(ScalarType)*y_v.size()), 
                                                                              &(y_v[0]));
            m_v_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                              static_cast<unsigned int>(sizeof(ScalarType)*m_v.size()),
                                                                              &(m_v[0]));
            y_v_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                               static_cast<unsigned int>(sizeof(cl_uint)*(g_I.size() + 1)), 
                                                                               &(y_inds[0]));
            viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                    static_cast<unsigned int>(sizeof(cl_uint)*(g_is_update.size())),
                                                                                        &(g_is_update[0]));
            viennacl::ocl::kernel& ls_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_least_squares");
            ls_kernel.local_work_size(0, 1);
            ls_kernel.global_work_size(0, 256);
            viennacl::ocl::enqueue(ls_kernel(g_A_I_J_vcl.handle(), g_A_I_J_vcl.handle2(), g_bv_vcl.handle(), g_bv_vcl.handle1(), m_v_vcl.handle(), 
                                             y_v_vcl.handle(), y_v_vcl.handle1(), 
                                             g_A_I_J_vcl.handle1(), g_is_update_vcl,
                                             //viennacl::ocl::local_mem(static_cast<unsigned int>(sizeof(ScalarType)*(local_r_n*local_c_n))),
                                             static_cast<unsigned int>(M_v.size())));
            //copy vector m_v back from GPU to CPU
            cl_int vcl_err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                                 m_v_vcl.handle().get(), CL_TRUE, 0, 
                                                 sizeof(ScalarType)*(m_v.size()),
                                                 &(m_v[0]), 0, NULL, NULL);
            VIENNACL_ERR_CHECK(vcl_err);
            //fan out vector in parallel
            //#pragma omp parallel for
            for(std::size_t i = 0; i < M_v.size(); ++i){
                if(g_is_update[i]){
                    //faned out onto sparse vector
                    custom_fan_out(m_v, m_inds[i], g_J[i], M_v[i]);
                    g_res[i].clear();
                    compute_spai_residual<SparseVectorType, ScalarType>(A_v_c,  M_v[i], static_cast<unsigned int>(i), g_res[i]);
                    ScalarType res_norm = 0;
                    //compute norm of res - just to make sure that this implementatino works correct
                    sparse_norm_2(g_res[i], res_norm);
                    //std::cout<<"Residual norm of column #: "<<i<<std::endl;
                    //std::cout<<res_norm<<std::endl;
                    //std::cout<<"************************"<<std::endl;
                    g_is_update[i] = (res_norm > tag.getResidualNormThreshold())&& (!tag.getIsStatic())?(1):(0);
                    
                }
            }
        }
        
        //CPU based least square problems
        /** @brief Solution of Least square problem on CPU 
         * @param A_v_c column-major vectorized initial sparse matrix 
         * @param g_R blocks for least square solution
         * @param g_b_v vectors beta, necessary for Q recovery
         * @param g_I container of row index set for all columns of matrix M
         * @param g_J container of column index set for all columns of matrix M
         * @param g_res container of residuals
         * @param g_is_update container with indicators which blocks are active
         * @param M_v column-major vectorized sparse matrix, final preconditioner
         * @param tag spai tag
         */ 
        template<typename SparseVectorType, typename DenseMatrixType, typename VectorType>
        void least_square_solve(const std::vector<SparseVectorType>& A_v_c,
                                std::vector<DenseMatrixType>& g_R, 
                                std::vector<VectorType>& g_b_v,
                                std::vector<std::vector<unsigned int> >& g_I, 
                                std::vector<std::vector<unsigned int> >& g_J,
                                std::vector<SparseVectorType>& g_res, 
                                std::vector<bool>& g_is_update, 
                                std::vector<SparseVectorType>& M_v,
                                const spai_tag& tag){
            typedef typename DenseMatrixType::value_type ScalarType;
            //VectorType m_new, y;
#ifdef _OPENMP
            #pragma omp parallel for
#endif            
            for(std::size_t i = 0; i < M_v.size(); ++i){
                if(g_is_update[i]){
                    VectorType y = boost::numeric::ublas::zero_vector<ScalarType>(g_I[i].size());
                    //std::cout<<y<<std::endl;
                    projectI<VectorType, ScalarType>(g_I[i], y, static_cast<unsigned int>(tag.getBegInd() + i));
                    apply_q_trans_vec(g_R[i], g_b_v[i], y);
                    VectorType m_new =  boost::numeric::ublas::zero_vector<ScalarType>(g_R[i].size2());
                    backwardSolve(g_R[i], y, m_new);
                    fanOutVector(m_new, g_J[i], M_v[i]);
                    g_res[i].clear();
                    compute_spai_residual<SparseVectorType, ScalarType>(A_v_c,  M_v[i], static_cast<unsigned int>(tag.getBegInd() + i), g_res[i]); 
                    ScalarType res_norm = 0;
                    sparse_norm_2(g_res[i], res_norm);
//                    std::cout<<"Residual norm of column #: "<<i<<std::endl;
//                    std::cout<<res_norm<<std::endl;
//                    std::cout<<"************************"<<std::endl;
                    g_is_update[i] = (res_norm > tag.getResidualNormThreshold())&& (!tag.getIsStatic());
                }
            }
        } 
        
        //************************************ UPDATE CHECK ***************************************************//
        template<typename VectorType>
        bool is_all_update(VectorType& parallel_is_update){
            
            for(unsigned int i = 0; i < parallel_is_update.size(); ++i){
                if(parallel_is_update[i])
                    return true;
            }
            return false;
        }
        
        //********************************** MATRIX VECTORIZATION ***********************************************//
        //Matrix vectorization, column based approach
        /** @brief Solution of Least square problem on CPU 
         * @param M_in input sparse, boost::numeric::ublas::compressed_matrix 
         * @param M_v array of sparse vectors 
         */ 
        template<typename SparseMatrixType, typename SparseVectorType>
        void vectorize_column_matrix(const SparseMatrixType& M_in, std::vector<SparseVectorType>& M_v){
            for(typename SparseMatrixType::const_iterator1 row_it = M_in.begin1(); row_it!= M_in.end1(); ++row_it){
                //
                for(typename SparseMatrixType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it){
                    M_v[static_cast<unsigned int>(col_it.index2())][static_cast<unsigned int>(col_it.index1())] = *col_it; 
                }
                //std::cout<<std::endl;
            }
        }
        
        //Matrix vectorization row based approach
        template<typename SparseMatrixType, typename SparseVectorType>
        void vectorize_row_matrix(const SparseMatrixType& M_in, std::vector<SparseVectorType>& M_v){
            for(typename SparseMatrixType::const_iterator1 row_it = M_in.begin1(); row_it!= M_in.end1(); ++row_it){
                for(typename SparseMatrixType::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it){
                    M_v[static_cast<unsigned int>(col_it.index1())][static_cast<unsigned int>(col_it.index2())] = *col_it; 
                }
            }
        }
        
        //************************************* BLOCK ASSEMBLY CODE *********************************************//
        
        
        
        void write_set_to_array(const std::vector<std::vector<unsigned int> >& ind_set, std::vector<cl_uint>& a){
            unsigned int cnt = 0;
            //unsigned int tmp;
            for(size_t i = 0; i < ind_set.size(); ++i){
                for(size_t j = 0; j < ind_set[i].size(); ++j){
                    a[cnt++] = static_cast<cl_uint>(ind_set[i][j]); 
                }
            }
        }
        
        
        
        //assembling blocks on GPU
        /** @brief Assembly of blocks on GPU by a gived set of row indices: g_I and column indices: g_J 
         * @param A intial sparse matrix  
         * @param g_J container of column index set
         * @param g_I container of row index set
         * @param g_A_I_J_vcl contigious blocks A(I, J) using GPU memory
         * @param g_is_update container with indicators which blocks are active
         * @param is_empty_block parameter that indicates if no block were assembled 
         * @param cur_iter current iteration number 
         */ 
        template<typename ScalarType, unsigned int MAT_ALIGNMENT>
        void block_assembly(const viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>& A, const std::vector<std::vector<unsigned int> >& g_J, 
                            const std::vector<std::vector<unsigned int> >& g_I,
                            block_matrix& g_A_I_J_vcl, 
                            std::vector<cl_uint>& g_is_update,
                            bool& is_empty_block,
                            const unsigned int cur_iter){
            //computing start indices for index sets and start indices for block matrices
            unsigned int sz_I, sz_J, sz_blocks;
            std::vector<cl_uint> matrix_dims(g_I.size()*2, static_cast<cl_uint>(0));
            std::vector<cl_uint> i_ind(g_I.size() + 1, static_cast<cl_uint>(0));
            std::vector<cl_uint> j_ind(g_I.size() + 1, static_cast<cl_uint>(0));
            std::vector<cl_uint> blocks_ind(g_I.size() + 1, static_cast<cl_uint>(0));
            // 
            init_start_inds(g_J, j_ind);
            init_start_inds(g_I, i_ind);
            //
            get_size(g_J, sz_J);
            get_size(g_I, sz_I);
            std::vector<cl_uint> I_set(sz_I, static_cast<cl_uint>(0));
            //
            std::vector<cl_uint> J_set(sz_J, static_cast<cl_uint>(0));
            // computing size for blocks
            // writing set to arrays
            write_set_to_array(g_I, I_set);
            write_set_to_array(g_J, J_set);
            // if block for assembly does exist
            if(I_set.size() > 0 && J_set.size() > 0){
                compute_blocks_size(g_I, g_J, sz_blocks, blocks_ind, matrix_dims);
                std::vector<ScalarType> con_A_I_J(sz_blocks, static_cast<ScalarType>(0));
                
                /*if(cur_iter == 0){
                    std::string block_asm_file_name = "kernels/spai/block_assembly_g.cl";
                    std::string block_asm_source;
                    read_kernel_from_file(block_asm_file_name, block_asm_source);
                    viennacl::ocl::program & block_asm_prog = viennacl::ocl::current_context().add_program(block_asm_source.c_str(), 
                                                                                                           "block_assembly_kernel_source");
                    
                    block_asm_prog.add_kernel("assemble_blocks");
                }*/
                block_vector set_I_vcl, set_J_vcl;
                //init memory on GPU
                //contigious g_A_I_J
                g_A_I_J_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                      static_cast<unsigned int>(sizeof(ScalarType)*(sz_blocks)), 
                                                                                      &(con_A_I_J[0]));
                //matrix_dimensions
                g_A_I_J_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                   static_cast<unsigned int>(sizeof(cl_uint)*2*static_cast<cl_uint>(g_I.size())), 
                                                                                   &(matrix_dims[0]));
                //start_block inds
                g_A_I_J_vcl.handle2() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                       static_cast<unsigned int>(sizeof(cl_uint)*(g_I.size() + 1)), 
                                                                                       &(blocks_ind[0]));
                //set_I
                set_I_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                    static_cast<unsigned int>(sizeof(cl_uint)*sz_I), 
                                                                                    &(I_set[0]));
                //set_J
                set_J_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                    static_cast<unsigned int>(sizeof(cl_uint)*sz_J), 
                                                                                    &(J_set[0]));
                //i_ind
                set_I_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                     static_cast<unsigned int>(sizeof(cl_uint)*(g_I.size() + 1)), 
                                                                                     &(i_ind[0]));
                //j_ind
                set_J_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                     static_cast<unsigned int>(sizeof(cl_uint)*(g_I.size() + 1)), 
                                                                                     &(j_ind[0]));
                
                viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                                               static_cast<unsigned int>(sizeof(cl_uint)*g_is_update.size()),
                                                                                                               &(g_is_update[0]));
                viennacl::ocl::kernel& assembly_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "assemble_blocks");
                assembly_kernel.local_work_size(0, 1);
                assembly_kernel.global_work_size(0, 256);
                viennacl::ocl::enqueue(assembly_kernel(A.handle1(), A.handle2(), A.handle(), 
                                                       set_I_vcl.handle(), set_J_vcl.handle(), set_I_vcl.handle1(), 
                                                       set_J_vcl.handle1(), 
                                                       g_A_I_J_vcl.handle2(), g_A_I_J_vcl.handle1(), g_A_I_J_vcl.handle(),
                                                       g_is_update_vcl,
                                                       static_cast<unsigned int>(g_I.size())));
                is_empty_block = false;
            }else{ 
                is_empty_block = true;
            }
        }
        
        /************************************************************************************************************************/
        
        /** @brief Insertion of vectorized matrix column into original sparse matrix 
         * @param M_v column-major vectorized matrix 
         * @param M original sparse matrix
         * @param is_right indicates if matrix should be transposed in the output
         */ 
        template<typename SparseMatrixType, typename SparseVectorType>
        void insert_sparse_columns(const std::vector<SparseVectorType>& M_v,
                                   SparseMatrixType& M,
                                   bool is_right){
            if (is_right)
            {
              for(unsigned int i = 0; i < M_v.size(); ++i){
                  for(typename SparseVectorType::const_iterator vec_it = M_v[i].begin(); vec_it!=M_v[i].end(); ++vec_it){
                      M(vec_it->first, i) = vec_it->second;
                  }
              }
            }
            else  //transposed fill of M
            {
              for(unsigned int i = 0; i < M_v.size(); ++i){
                  for(typename SparseVectorType::const_iterator vec_it = M_v[i].begin(); vec_it!=M_v[i].end(); ++vec_it){
                      M(i, vec_it->first) = vec_it->second;
                  }
              }
            }
        }
        
        /** @brief Transposition of sparse matrix
         * @param A_in intial sparse matrix  
         * @param A output transposed matrix 
         */ 
        template<typename MatrixType>
        void sparse_transpose(const MatrixType& A_in, MatrixType& A){
            typedef typename MatrixType::value_type ScalarType;
            std::vector<std::map<size_t, ScalarType> >   temp_A(A_in.size2());
            A.resize(A_in.size2(), A_in.size1(), false);
            
            for (typename MatrixType::const_iterator1 row_it = A_in.begin1();
                 row_it != A_in.end1();
                 ++row_it)
            {
                for (typename MatrixType::const_iterator2 col_it = row_it.begin();
                     col_it != row_it.end();
                     ++col_it)
                {
                    temp_A[col_it.index2()][col_it.index1()] = *col_it;
                }
            }
            
            for (size_t i=0; i<temp_A.size(); ++i)
            {
                for (typename std::map<size_t, ScalarType>::const_iterator it = temp_A[i].begin();
                     it != temp_A[i].end();
                     ++it)
                    A(i, it->first) = it->second;
            }
        }
        
        
        
        
//        template<typename SparseVectorType>
//        void custom_copy(std::vector<SparseVectorType> & M_v, std::vector<SparseVectorType> & l_M_v, const unsigned int beg_ind){
//            for(int i = 0; i < l_M_v.size(); ++i){
//                l_M_v[i] = M_v[i + beg_ind];
//            }
//        }
        
        //CPU version
        /** @brief Construction of SPAI preconditioner on CPU 
         * @param A initial sparse matrix   
         * @param M output preconditioner
         * @param tag spai tag
         */ 
        template <typename MatrixType>
        void computeSPAI(const MatrixType & A, MatrixType & M, spai_tag & tag){
            typedef typename MatrixType::value_type ScalarType;
            typedef typename boost::numeric::ublas::vector<ScalarType> VectorType;
            typedef typename viennacl::linalg::detail::spai::sparse_vector<ScalarType> SparseVectorType;
            typedef typename boost::numeric::ublas::matrix<ScalarType> DenseMatrixType;
            //sparse matrix transpose...
            unsigned int cur_iter = 0;
            tag.setBegInd(0); tag.setEndInd(VIENNACL_SPAI_K_b);
            bool go_on = true;
            std::vector<SparseVectorType> A_v_c(M.size2());
            std::vector<SparseVectorType> M_v(M.size2());
            vectorize_column_matrix(A, A_v_c);
            vectorize_column_matrix(M, M_v);
            
            
            while(go_on){
                go_on = (tag.getEndInd() < static_cast<long>(M.size2()));
                cur_iter = 0;
                unsigned int l_sz = tag.getEndInd() - tag.getBegInd();
                //std::vector<bool> g_is_update(M.size2(), true);
                std::vector<bool> g_is_update(l_sz, true);
                //init is update
                //init_parallel_is_update(g_is_update);
                //std::vector<SparseVectorType> A_v_c(K);
                //std::vector<SparseVectorType> M_v(K);
                //vectorization of marices
                //print_matrix(M_v);
                std::vector<SparseVectorType> l_M_v(l_sz);
                //custom_copy(M_v, l_M_v, beg_ind);
                std::copy(M_v.begin() + tag.getBegInd(), M_v.begin() + tag.getEndInd(), l_M_v.begin()); 
                //print_matrix(l_M_v);
                //std::vector<SparseVectorType> l_A_v_c(K);
                //custom_copy(A_v_c, l_A_v_c, beg_ind);
                //std::copy(A_v_c.begin() + beg_ind, A_v_c.begin() + end_ind, l_A_v_c.begin());
                //print_matrix(l_A_v_c);
                //vectorize_row_matrix(A, A_v_r);
                //working blocks
                //std::vector<DenseMatrixType> g_A_I_J(M.size2())                
                std::vector<DenseMatrixType> g_A_I_J(l_sz);
                //std::vector<VectorType> g_b_v(M.size2());
                std::vector<VectorType> g_b_v(l_sz);
                //std::vector<SparseVectorType> g_res(M.size2());
                std::vector<SparseVectorType> g_res(l_sz);
                //std::vector<std::vector<unsigned int> > g_I(M.size2());
                std::vector<std::vector<unsigned int> > g_I(l_sz);
                //std::vector<std::vector<unsigned int> > g_J(M.size2());
                std::vector<std::vector<unsigned int> > g_J(l_sz);
                while((cur_iter < tag.getIterationLimit())&&is_all_update(g_is_update)){
                    // SET UP THE BLOCKS..
                    // PHASE ONE
                    if(cur_iter == 0) block_set_up(A, A_v_c, l_M_v,  g_I, g_J, g_A_I_J, g_b_v);
                    else block_update(A, A_v_c, g_res, g_is_update, g_I, g_J, g_b_v, g_A_I_J, tag);
                    
                    //PHASE TWO, LEAST SQUARE SOLUTION
                    least_square_solve(A_v_c, g_A_I_J, g_b_v, g_I, g_J, g_res, g_is_update, l_M_v, tag);
                    
                    if(tag.getIsStatic()) break;
                    cur_iter++;
                    
                    
                }
                std::copy(l_M_v.begin(), l_M_v.end(), M_v.begin() + tag.getBegInd());
                tag.setBegInd(tag.getEndInd());//beg_ind = end_ind;
                tag.setEndInd(std::min(static_cast<long>(tag.getBegInd() + VIENNACL_SPAI_K_b), static_cast<long>(M.size2())));
                //std::copy(l_M_v.begin(), l_M_v.end(), M_v.begin() + tag.getBegInd()); 
                
            }
            M.resize(M.size1(), M.size2(), false);
            insert_sparse_columns(M_v, M, tag.getIsRight());
        }

        
        //GPU - based version
        /** @brief Construction of SPAI preconditioner on GPU
         * @param A initial sparse matrix   
         * @param cpu_A copy of initial matrix on CPU
         * @param cpu_M output preconditioner on CPU
         * @param M output preconditioner
         * @param tag SPAI tag class with parameters
         */ 
        template <typename ScalarType, unsigned int MAT_ALIGNMENT>
        void computeSPAI(const viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>& A, //input
                         const boost::numeric::ublas::compressed_matrix<ScalarType>& cpu_A,
                         boost::numeric::ublas::compressed_matrix<ScalarType>& cpu_M, //output
                         viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>& M,
                         const spai_tag& tag){
            typedef typename boost::numeric::ublas::vector<ScalarType> VectorType;
            typedef typename viennacl::linalg::detail::spai::sparse_vector<ScalarType> SparseVectorType;
            typedef typename boost::numeric::ublas::matrix<ScalarType> DenseMatrixType;
            typedef typename boost::numeric::ublas::compressed_matrix<ScalarType> CPUMatrixType;
            //typedef typename viennacl::compressed_matrix<ScalarType> GPUSparseMatrixType;
            //sparse matrix transpose...
            unsigned int cur_iter = 0;
            std::vector<cl_uint> g_is_update(cpu_M.size2(), static_cast<cl_uint>(1));
            //init is update
            //init_parallel_is_update(g_is_update);
            std::vector<SparseVectorType> A_v_c(cpu_M.size2());
            std::vector<SparseVectorType> M_v(cpu_M.size2());
            vectorize_column_matrix(cpu_A, A_v_c);
            vectorize_column_matrix(cpu_M, M_v);
            std::vector<SparseVectorType> g_res(cpu_M.size2());
            std::vector<std::vector<unsigned int> > g_I(cpu_M.size2());
            std::vector<std::vector<unsigned int> > g_J(cpu_M.size2());
            
            //OpenCL variables
            block_matrix g_A_I_J_vcl;
            block_vector g_bv_vcl;
            while((cur_iter < tag.getIterationLimit())&&is_all_update(g_is_update)){
                // SET UP THE BLOCKS..
                // PHASE ONE..
                //timer.start();
                //index set up on CPU
                if(cur_iter == 0) block_set_up(A, A_v_c, M_v, g_is_update, g_I, g_J, g_A_I_J_vcl, g_bv_vcl, cur_iter);
                else block_update(A, A_v_c, g_is_update, g_res, g_J, g_I, g_A_I_J_vcl, g_bv_vcl, tag, cur_iter);
                //std::cout<<"Phase 2 timing: "<<timer.get()<<std::endl;
                //PERFORM LEAST SQUARE problems solution
                //PHASE TWO
                //timer.start();
                least_square_solve<SparseVectorType, ScalarType>(A_v_c, M_v, g_I, g_J, g_A_I_J_vcl, g_bv_vcl, g_res, g_is_update, tag, cur_iter);
                //std::cout<<"Phase 3 timing: "<<timer.get()<<std::endl;
                if(tag.getIsStatic()) break;
                cur_iter++;
            }
            cpu_M.resize(cpu_M.size1(), cpu_M.size2(), false);
            insert_sparse_columns(M_v, cpu_M, tag.getIsRight());
            //copy back to GPU
            M.resize(static_cast<unsigned int>(cpu_M.size1()), static_cast<unsigned int>(cpu_M.size2()));
            viennacl::copy(cpu_M, M);
        }
        
      }        
    }
  }
}
#endif