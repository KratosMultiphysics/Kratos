#ifndef VIENNACL_LINALG_DETAIL_SPAI_SPAI_DYNAMIC_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_SPAI_DYNAMIC_HPP

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

/** @file viennacl/linalg/detail/spai/spai-dynamic.hpp
    @brief Implementation of a dynamic SPAI. Provides the routines for automatic pattern updates Experimental in 1.2.x.
    
    SPAI code contributed by Nikolay Lukash
*/

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>
//#include "block_matrix.hpp"
//#include "block_vector.hpp"
//#include "benchmark-utils.hpp"
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
#include "viennacl/linalg/cg.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/ilu.hpp"
#include "viennacl/ocl/backend.hpp"

#include "viennacl/linalg/detail/spai/block_matrix.hpp"
#include "viennacl/linalg/detail/spai/block_vector.hpp"
#include "viennacl/linalg/detail/spai/qr.hpp"
#include "viennacl/linalg/detail/spai/spai_tag.hpp"
#include "viennacl/linalg/kernels/spai_source.h"
#include "viennacl/linalg/kernels/spai_kernels.h"

namespace viennacl
{
    namespace linalg
    {
      namespace detail
      {
        namespace spai
        {
        
          typedef std::pair<unsigned int, double> PairT;
          struct CompareSecond{
              bool operator()(const PairT& left, const PairT& right)
              {
                  return static_cast<double>(left.second) > static_cast<double>(right.second);
              }
          };
          
          
          /** @brief Initializes Dense matrix from sparse one
          * @param A_in Riginal sparse matrix
          * @param J Set of column indices
          * @param I Set of row indices 
          * @param A_out dense matrix output
          */
          template<typename SparseMatrixType, typename DenseMatrixType>
          void initProjectSubMatrix(const SparseMatrixType& A_in, const std::vector<unsigned int>& J, std::vector<unsigned int>& I,
                                    DenseMatrixType& A_out){
              typedef typename DenseMatrixType::value_type ScalarType;
              A_out.resize(I.size(), J.size(), false);
              for(size_t j = 0; j < J.size(); ++j){
                  for(size_t i = 0; i < I.size(); ++i){
                      A_out(i,j) = A_in(I[i],J[j]);
                  }
              }
          }
          
          /** @brief Determines if element ind is in set {J}
          * @param J current set
          * @param ind current element
          */
          bool isInIndexSet(const std::vector<unsigned int>& J, const unsigned int& ind){
              return (std::find(J.begin(), J.end(), ind) != J.end());
          }
          
          /** @brief Composition of new matrix R, that is going to be used in Least Square problem solving
          * @param A matrix Q'*A(I, \\tilde J), where \\tilde J - set of new column indices
          * @param R_n matrix A_Iu_J_u after QR factorization
          * @param R previously composed matrix R 
          */
          template<typename MatrixType>
          void composeNewR(const MatrixType& A, const MatrixType& R_n, MatrixType& R){
              typedef typename MatrixType::value_type ScalarType;
              size_t row_n = R_n.size1() - (A.size1() - R.size2()); 
              MatrixType C = boost::numeric::ublas::zero_matrix<ScalarType>(R.size1() + row_n, R.size2() + A.size2());
              //write original R to new Composite R
              boost::numeric::ublas::project(C, boost::numeric::ublas::range(0,R.size1()), boost::numeric::ublas::range(0, R.size2())) += R;
              //write upper part of Q'*A_I_\hatJ, all columns and number of rows that equals to R.size2()
              boost::numeric::ublas::project(C, boost::numeric::ublas::range(0, R.size2()), boost::numeric::ublas::range(R.size2(), 
                                                                                                                        R.size2() + A.size2())) += 
              boost::numeric::ublas::project(A, boost::numeric::ublas::range(0, R.size2()), boost::numeric::ublas::range(0, A.size2()));
              //adding decomposed(QR) block to Composite R
              if(R_n.size1() > 0 && R_n.size2() > 0)
                  boost::numeric::ublas::project(C, boost::numeric::ublas::range(R.size2(), R.size1() + row_n),
                                                boost::numeric::ublas::range(R.size2(), R.size2() + A.size2())) += R_n;
              R = C;
          }
          
          /** @brief Composition of new vector of coefficients beta from QR factorizations(necessary for Q recovery)
          * @param v_n new vector from last QR factorization
          * @param v composition of previous vectors from QR factorizations
          */
          template<typename VectorType>
          void composeNewVector(const VectorType& v_n, VectorType& v){
              typedef typename VectorType::value_type ScalarType;
              VectorType w  = boost::numeric::ublas::zero_vector<ScalarType>(v.size() + v_n.size());
              boost::numeric::ublas::project(w, boost::numeric::ublas::range(0, v.size())) += v;
              boost::numeric::ublas::project(w, boost::numeric::ublas::range(v.size(), v.size() + v_n.size())) += v_n;
              v = w;
          }
          
          /** @brief Computation of Euclidean norm for sparse vector
          * @param v initial sparse vector
          * @param norm scalar that represents Euclidean norm
          */
          template<typename SparseVectorType, typename ScalarType>
          void sparse_norm_2(const SparseVectorType& v, ScalarType& norm){
              for(typename SparseVectorType::const_iterator vec_it  = v.begin(); vec_it != v.end(); ++vec_it){
                  norm += (vec_it->second)*(vec_it->second);
              }
              norm = std::sqrt(norm);
          }
          
          /** @brief Dot product of two sparse vectors
          * @param v1 initial sparse vector
          * @param v2 initial sparse vector
          * @param res_v scalar that represents dot product result
          */
          template<typename SparseVectorType, typename ScalarType>
          void sparse_inner_prod(const SparseVectorType& v1, const SparseVectorType& v2, ScalarType& res_v){
              typename SparseVectorType::const_iterator v_it1 = v1.begin();
              typename SparseVectorType::const_iterator v_it2 = v2.begin();
              while((v_it1 != v1.end())&&(v_it2 != v2.end())){
                  if(v_it1->first == v_it2->first){
                      res_v += (v_it1->second)*(v_it2->second);
                      ++v_it1;
                      ++v_it2;
                  }
                  else if(v_it1->first < v_it2->first){
                      ++v_it1;
                  }
                  else 
                      ++v_it2;
                      
                  
              }
          }
          
          /** @brief Building a new set of column indices J_u, cf. Kallischko dissertation p.31
          * @param A_v_c vectorized column-wise initial matrix 
          * @param res residual vector 
          * @param J set of column indices
          * @param J_u set of new column indices 
          * @param tag SPAI tag with parameters
          */
          template <typename SparseVectorType, typename ScalarType>
          bool buildAugmentedIndexSet(const std::vector<SparseVectorType>& A_v_c, 
                                      const SparseVectorType& res,
                                      std::vector<unsigned int>& J,
                                      std::vector<unsigned int>& J_u,
                                      const spai_tag& tag){
              std::vector<std::pair<unsigned int, ScalarType> > p;
              size_t cur_size = 0;
              ScalarType inprod, norm2;
              //print_sparse_vector(res);
              for(typename SparseVectorType::const_iterator res_it = res.begin(); res_it != res.end(); ++res_it){
                  if(!isInIndexSet(J, res_it->first) && (std::abs(res_it->second) > tag.getResidualThreshold())){
                      inprod = norm2 = 0;
                      sparse_inner_prod(res, A_v_c[res_it->first], inprod);
                      sparse_norm_2(A_v_c[res_it->first], norm2);
                      p.push_back(std::pair<size_t, ScalarType>(res_it->first, (inprod*inprod)/(norm2*norm2)));
                  }
              }
              
              std::sort(p.begin(), p.end(), CompareSecond());
              while ((cur_size < J.size())&&(p.size() > 0)) {
                  J_u.push_back(p[0].first);
                  p.erase(p.begin());
                  cur_size++;
              }
              p.clear();
              return (cur_size > 0);
          }
          
          /** @brief Building a new indices to current set of row indices I_n, cf. Kallischko dissertation p.32
          * @param A_v_c vectorized column-wise initial matrix
          * @param I set of previous determined row indices
          * @param J_n set of new column indices
          * @param I_n set of new indices 
          */
          template<typename SparseVectorType>
          void buildNewRowSet(const std::vector<SparseVectorType>& A_v_c, const std::vector<unsigned int>& I, 
                              const std::vector<unsigned int>& J_n, std::vector<unsigned int>& I_n){
              for(size_t i = 0; i < J_n.size(); ++i){
                  for(typename SparseVectorType::const_iterator col_it = A_v_c[J_n[i]].begin(); col_it!=A_v_c[J_n[i]].end(); ++col_it){
                      if(!isInIndexSet(I, col_it->first)&&!isInIndexSet(I_n, col_it->first)){
                          I_n.push_back(col_it->first);
                      }
                  }
              }
          }
          
          /** @brief Composition of new block for QR factorization cf. Kallischko dissertation p.82, figure 4.7
          * @param A_I_J previously composed block
          * @param A_I_J_u matrix Q'*A(I, \\tilde J), where \\tilde J - set of new column indices
          * @param A_I_u_J_u is composition of lower part A(I, \\tilde J) and  A(\\tilde I, \\tilde J) - new block for QR decomposition
          */
          template<typename MatrixType>
          void QRBlockComposition(const MatrixType& A_I_J, const MatrixType& A_I_J_u, MatrixType& A_I_u_J_u){
              typedef typename MatrixType::value_type ScalarType;
              size_t row_n1 = A_I_J_u.size1() - A_I_J.size2();
              size_t row_n2 = A_I_u_J_u.size1();
              size_t row_n = row_n1 + row_n2;
              size_t col_n = A_I_J_u.size2();
              MatrixType C = boost::numeric::ublas::zero_matrix<ScalarType>(row_n, col_n);
              boost::numeric::ublas::project(C, boost::numeric::ublas::range(0, row_n1), boost::numeric::ublas::range(0, col_n)) += 
              boost::numeric::ublas::project(A_I_J_u, boost::numeric::ublas::range(A_I_J.size2(), A_I_J_u.size1()),
                                            boost::numeric::ublas::range(0, col_n));
                                            
              boost::numeric::ublas::project(C, boost::numeric::ublas::range(row_n1, row_n1 + row_n2),
                                            boost::numeric::ublas::range(0, col_n)) += A_I_u_J_u;
              A_I_u_J_u = C;
          }
          
          /** @brief CPU-based dynamic update for SPAI preconditioner
          * @param A initial sparse matrix
          * @param A_v_c vectorized column-wise initial matrix
          * @param g_res container of residuals for all columns
          * @param g_is_update container with identificators that shows which block should be modified 
          * @param g_I container of row index sets for all columns
          * @param g_J container of column index sets for all columns
          * @param g_b_v container of vectors of beta for Q recovery(cf. Golub Van Loan "Matrix Computations", 3rd edition p.211)
          * @param g_A_I_J container of block matrices from previous update  
          * @param tag  SPAI configuration tag
          */
          template<typename SparseMatrixType, typename SparseVectorType, typename DenseMatrixType, typename VectorType>
          void block_update(const SparseMatrixType& A, const std::vector<SparseVectorType>& A_v_c,
                            std::vector<SparseVectorType>& g_res,
                            std::vector<bool>& g_is_update,
                            std::vector<std::vector<unsigned int> >& g_I,
                            std::vector<std::vector<unsigned int> >& g_J, 
                            std::vector<VectorType>& g_b_v,
                            std::vector<DenseMatrixType>& g_A_I_J,
                            spai_tag const & tag){
              typedef typename DenseMatrixType::value_type ScalarType;
              //set of new column indices
              std::vector<std::vector<unsigned int> > g_J_u(g_J.size());
              //set of new row indices 
              std::vector<std::vector<unsigned int> > g_I_u(g_J.size());
              //matrix A(I, \tilde J), cf. Kallischko p.31-32
              std::vector<DenseMatrixType> g_A_I_J_u(g_J.size());
              //matrix A(\tilde I, \tilde J), cf. Kallischko
              std::vector<DenseMatrixType> g_A_I_u_J_u(g_J.size());
              //new vector of beta coefficients from QR factorization
              std::vector<VectorType> g_b_v_u(g_J.size());
#ifdef _OPENMP
              #pragma omp parallel for
#endif              
              for(std::size_t i = 0; i < g_J.size(); ++i){
                  if(g_is_update[i]){
                      if(buildAugmentedIndexSet<SparseVectorType, ScalarType>(A_v_c, g_res[i], g_J[i], g_J_u[i], tag)){
                          //initialize matrix A_I_\hatJ
                          initProjectSubMatrix(A, g_J_u[i], g_I[i], g_A_I_J_u[i]);
                          //multiplication of Q'*A_I_\hatJ
                          apply_q_trans_mat(g_A_I_J[i], g_b_v[i], g_A_I_J_u[i]);
                          //building new rows index set \hatI
                          buildNewRowSet(A_v_c, g_I[i], g_J_u[i], g_I_u[i]);
                          initProjectSubMatrix(A, g_J_u[i], g_I_u[i], g_A_I_u_J_u[i]);
                          //composition of block for new QR factorization
                          QRBlockComposition(g_A_I_J[i], g_A_I_J_u[i], g_A_I_u_J_u[i]);
                          //QR factorization
                          single_qr(g_A_I_u_J_u[i], g_b_v_u[i]);
                          //composition of new R and new vector b_v
                          composeNewR(g_A_I_J_u[i], g_A_I_u_J_u[i], g_A_I_J[i]);
                          composeNewVector(g_b_v_u[i], g_b_v[i]);
                          //composition of new sets: I and J
                          g_J[i].insert(g_J[i].end(), g_J_u[i].begin(), g_J_u[i].end());
                          g_I[i].insert(g_I[i].end(), g_I_u[i].begin(), g_I_u[i].end());
                      }else{
                          g_is_update[i] = false;
                      }
                  }
              }
          }
          /**************************************************** GPU SPAI Update ****************************************************************/
          
          
          //performs Q'*A(I, \tilde J) on GPU
          /** @brief Performs multiplication Q'*A(I, \\tilde J) on GPU
          * @param g_J_u container of sets of new column indices
          * @param g_I container of row indices
          * @param g_A_I_J_vcl block matrix composed from previous blocks, they are blocks of R
          * @param g_bv_vcl block of beta vectors 
          * @param g_A_I_J_u_vcl block of matrices A(I, \\tilde J) 
          * @param g_is_update indicators, that show if a certain block should be processed
          * @param cur_iter current iteration, used to make sure that kernel compiles just once
          */
          template<typename ScalarType>
          void block_q_multiplication(const std::vector<std::vector<unsigned int> >& g_J_u, 
                                      const std::vector<std::vector<unsigned int> >& g_I, 
                                      block_matrix& g_A_I_J_vcl, 
                                      block_vector& g_bv_vcl, 
                                      block_matrix& g_A_I_J_u_vcl, 
                                      std::vector<cl_uint>& g_is_update,
                                      const unsigned int cur_iter){
              unsigned int local_r_n, local_c_n, sz_blocks;
              get_max_block_size(g_I, local_r_n);
              get_max_block_size(g_J_u, local_c_n);
              //for debug 
              std::vector<cl_uint> matrix_dims(g_I.size()*2, static_cast<cl_uint>(0));
              std::vector<cl_uint> blocks_ind(g_I.size() + 1, static_cast<cl_uint>(0));
              compute_blocks_size(g_I, g_J_u, sz_blocks, blocks_ind, matrix_dims);
              std::vector<ScalarType> con_A_I_J(sz_blocks, static_cast<ScalarType>(0));
              /*if(cur_iter == 1){
                  //if first run - compile the program
                  std::string block_q_kernel_file_name = "kernels/spai/block_q.cl";
                  std::string block_q_kernel_source;
                  read_kernel_from_file(block_q_kernel_file_name, block_q_kernel_source);
                  viennacl::ocl::program & block_q_prog = viennacl::ocl::current_context().add_program(block_q_kernel_source.c_str(), "block_q_kernel_source");
                  block_q_prog.add_kernel("block_q_mult");
                  //
              }*/
              
              viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                  static_cast<unsigned int>(sizeof(cl_uint)*(g_is_update.size())),
                                                                                                            &(g_is_update[0]));
              viennacl::ocl::kernel& block_q_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_q_mult");
              block_q_kernel.local_work_size(0, local_c_n);
              block_q_kernel.global_work_size(0, 256);
              viennacl::ocl::enqueue(block_q_kernel(g_A_I_J_vcl.handle(), g_A_I_J_vcl.handle2(), g_A_I_J_u_vcl.handle(), g_A_I_J_u_vcl.handle2(), 
                                                    g_bv_vcl.handle(),   
                                                    g_bv_vcl.handle1(), g_A_I_J_vcl.handle1(), g_A_I_J_u_vcl.handle1(), g_is_update_vcl,
                                                    viennacl::ocl::local_mem(static_cast<unsigned int>(sizeof(ScalarType)*(local_r_n*local_c_n))),
                                                    static_cast<cl_uint>(g_I.size())));
          }
          
          /** @brief Assembly of container of index row sets: I_q, row indices for new "QR block"
          * @param g_I container of row indices
          * @param g_J container of column indices
          * @param g_I_u container of new row indices 
          * @param g_I_q container of row indices for new QR blocks
          */
          void assemble_qr_row_inds(const std::vector<std::vector<unsigned int> >& g_I, const std::vector<std::vector<unsigned int> > g_J, 
                                    const std::vector<std::vector<unsigned int> >& g_I_u, 
                                    std::vector<std::vector<unsigned int> >& g_I_q){
#ifdef _OPENMP
              #pragma omp parallel for
#endif              
              for(std::size_t i = 0; i < g_I.size(); ++i){
                  for(std::size_t j = g_J[i].size(); j < g_I[i].size(); ++j){
                      g_I_q[i].push_back(g_I[i][j]);
                  }
                  
                  for(std::size_t j = 0; j < g_I_u[i].size(); ++j){
                      g_I_q[i].push_back(g_I_u[i][j]);
                  }
              }
          }

          /** @brief Performs assembly for new QR block
          * @param g_J container of column indices 
          * @param g_I container of row indices
          * @param g_J_u container of new column indices
          * @param g_I_u container of new row indices 
          * @param g_I_q container of row indices for new QR blocks 
          * @param g_A_I_J_u_vcl blocks of Q'*A(I, \\tilde J)
          * @param matrix_dimensions array with matrix dimensions for all blocks
          * @param g_A_I_u_J_u_vcl blocks A(\\tilde I, \\tilde J)
          * @param g_is_update container with update indicators
          * @param is_empty_block indicator if all previous blocks A(\\tilde I, \\tilde J) - are empty, in case if they are empty kernel with smaller number of
                                  arguments is used
          * @param cur_iter current iteration, used to make sure that kernel compiles just once 
          */
          template<typename ScalarType>
          void assemble_qr_block(
                                const std::vector<std::vector<unsigned int> >& g_J, 
                                const std::vector<std::vector<unsigned int> >& g_I, 
                                const std::vector<std::vector<unsigned int> >& g_J_u,
                                const std::vector<std::vector<unsigned int> >& g_I_u, 
                                std::vector<std::vector<unsigned int> >& g_I_q,
                                block_matrix& g_A_I_J_u_vcl, 
                                viennacl::ocl::handle<cl_mem>& matrix_dimensions, 
                                block_matrix& g_A_I_u_J_u_vcl, 
                                std::vector<cl_uint>& g_is_update,
                                const bool is_empty_block,
                                const unsigned int cur_iter){
              //std::vector<std::vector<unsigned int> > g_I_q(g_I.size());
              assemble_qr_row_inds(g_I, g_J, g_I_u, g_I_q);
              unsigned int sz_blocks;
              std::vector<cl_uint> matrix_dims(g_I.size()*2, static_cast<cl_uint>(0));
              std::vector<cl_uint> blocks_ind(g_I.size() + 1, static_cast<cl_uint>(0));
              compute_blocks_size(g_I_q, g_J_u, sz_blocks, blocks_ind, matrix_dims);
              std::vector<ScalarType> con_A_I_J_q(sz_blocks, static_cast<ScalarType>(0));
              
              /*if(cur_iter == 1){
                  std::string qr_block_asm_file_name = "kernels/spai/qr_block_assembly_g.cl";
                  std::string qr_block_asm_source;
                  read_kernel_from_file(qr_block_asm_file_name, qr_block_asm_source);
                  viennacl::ocl::program & qr_block_asm_prog = viennacl::ocl::current_context().add_program(qr_block_asm_source.c_str(), 
                                                                                                        "qr_block_assembly_kernel_source");
                  
                  qr_block_asm_prog.add_kernel("block_qr_assembly");
                  
                  
                  //extra kernel in case of empty block A_I_u_J_u
                  std::string qr_block_asm_file_name_1 = "kernels/spai/qr_block_assembly_1_g.cl";
                  std::string qr_block_asm_source_1;
                  read_kernel_from_file(qr_block_asm_file_name_1, qr_block_asm_source_1);
                  viennacl::ocl::program & qr_block_asm_prog_1 = viennacl::ocl::current_context().add_program(qr_block_asm_source_1.c_str(), 
                                                                                                            "qr_block_assembly_kernel_source_1");
                  
                  qr_block_asm_prog_1.add_kernel("block_qr_assembly_1");
                  
              }*/
              block_matrix g_A_I_J_q_vcl;
              //need to allocate memory for QR block
              g_A_I_J_q_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                      static_cast<unsigned int>(sizeof(ScalarType)*sz_blocks),
                                                                                      &(con_A_I_J_q[0]));
              g_A_I_J_q_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                static_cast<unsigned int>(sizeof(cl_uint)*2*static_cast<unsigned int>(g_I.size())), 
                                                                                &(matrix_dims[0]));
              g_A_I_J_q_vcl.handle2() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                  static_cast<unsigned int>(sizeof(cl_uint)*2*static_cast<unsigned int>(g_I.size() + 1)),
                                                                                      &(blocks_ind[0]));
              viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                                            static_cast<unsigned int>(sizeof(cl_uint)*(g_is_update.size())),
                                                                                                            &(g_is_update[0]));
              
              if(!is_empty_block){
                  viennacl::ocl::kernel& qr_assembly_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_qr_assembly");
                  qr_assembly_kernel.local_work_size(0, 1);
                  qr_assembly_kernel.global_work_size(0, 256);
                  viennacl::ocl::enqueue(qr_assembly_kernel(matrix_dimensions, 
                                                            g_A_I_J_u_vcl.handle(), 
                                                            g_A_I_J_u_vcl.handle2(), 
                                                            g_A_I_J_u_vcl.handle1(), 
                                                            g_A_I_u_J_u_vcl.handle(),
                                                            g_A_I_u_J_u_vcl.handle2(), 
                                                            g_A_I_u_J_u_vcl.handle1(), 
                                                            g_A_I_J_q_vcl.handle(), 
                                                            g_A_I_J_q_vcl.handle2(), 
                                                            g_A_I_J_q_vcl.handle1(),
                                                            g_is_update_vcl,
                                                            static_cast<unsigned int>(g_I.size())));
              }else{
                  viennacl::ocl::kernel& qr_assembly_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_qr_assembly_1");
                  qr_assembly_kernel.local_work_size(0, 1);
                  qr_assembly_kernel.global_work_size(0, 256);
                  viennacl::ocl::enqueue(qr_assembly_kernel(matrix_dimensions, g_A_I_J_u_vcl.handle(), g_A_I_J_u_vcl.handle2(), 
                                                            g_A_I_J_u_vcl.handle1(),
                                                            g_A_I_J_q_vcl.handle(), 
                                                            g_A_I_J_q_vcl.handle2(), g_A_I_J_q_vcl.handle1(),
                                                            g_is_update_vcl,
                                                            static_cast<unsigned int>(g_I.size())));
              }
              g_A_I_u_J_u_vcl.handle() = g_A_I_J_q_vcl.handle();
              g_A_I_u_J_u_vcl.handle1() = g_A_I_J_q_vcl.handle1();
              g_A_I_u_J_u_vcl.handle2() = g_A_I_J_q_vcl.handle2();
          }

          /** @brief Performs assembly for new R matrix on GPU
          * @param g_I container of row indices
          * @param g_J container of column indices
          * @param g_A_I_J_vcl container of block matrices from previous update  
          * @param g_A_I_J_u_vcl container of block matrices Q'*A(I, \\tilde J)
          * @param g_A_I_u_J_u_vcl container of block matrices QR factored on current iteration
          * @param g_bv_vcl block of beta vectors from previous iteration
          * @param g_bv_vcl_u block of updated beta vectors got after recent QR factorization
          * @param g_is_update container with identificators that shows which block should be modified 
          * @param cur_iter current iteration, used to make sure that kernel compiles just once 
          */ 
          template<typename ScalarType>
          void assemble_r(std::vector<std::vector<unsigned int> >& g_I, std::vector<std::vector<unsigned int> >& g_J,
                          block_matrix& g_A_I_J_vcl, 
                          block_matrix& g_A_I_J_u_vcl,
                          block_matrix& g_A_I_u_J_u_vcl, 
                          block_vector& g_bv_vcl, 
                          block_vector& g_bv_vcl_u,
                          std::vector<cl_uint>& g_is_update,
                          const unsigned int cur_iter){
              std::vector<cl_uint> matrix_dims(g_I.size()*2, static_cast<cl_uint>(0));
              std::vector<cl_uint> blocks_ind(g_I.size() + 1, static_cast<cl_uint>(0));
              std::vector<cl_uint> start_bv_r_inds(g_I.size() + 1, 0);
              unsigned int sz_blocks, bv_size;
              compute_blocks_size(g_I, g_J, sz_blocks, blocks_ind, matrix_dims);
              get_size(g_J, bv_size);
              init_start_inds(g_J, start_bv_r_inds);
              std::vector<ScalarType> con_A_I_J_r(sz_blocks, static_cast<ScalarType>(0));
              std::vector<ScalarType> b_v_r(bv_size, static_cast<ScalarType>(0));
              /*if(cur_iter == 1){
                  std::string r_block_asm_file_name = "kernels/spai/r_block_assembly_g.cl";
                  std::string r_block_asm_source;
                  read_kernel_from_file(r_block_asm_file_name, r_block_asm_source);
                  viennacl::ocl::program & r_block_asm_prog = viennacl::ocl::current_context().add_program(r_block_asm_source.c_str(), 
                                                                                                            "r_block_assembly_kernel_source");
                  r_block_asm_prog.add_kernel("block_r_assembly");
                  
                  std::string bv_block_asm_file_name = "kernels/spai/bv_block_assembly_g.cl";
                  std::string bv_block_asm_source;
                  read_kernel_from_file(bv_block_asm_file_name, bv_block_asm_source);
                  viennacl::ocl::program & bv_block_asm_prog = viennacl::ocl::current_context().add_program(bv_block_asm_source.c_str(), 
                                                                                                          "bv_block_assembly_kernel_source");
                  bv_block_asm_prog.add_kernel("block_bv_assembly");
              }*/
              block_matrix g_A_I_J_r_vcl;
              block_vector g_bv_r_vcl;
              g_A_I_J_r_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                      static_cast<unsigned int>(sizeof(ScalarType)*sz_blocks),
                                                                                      &(con_A_I_J_r[0]));
              g_A_I_J_r_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                          static_cast<unsigned int>(sizeof(cl_uint)*2*static_cast<unsigned int>(g_I.size())), 
                                                                                      &(matrix_dims[0]));
              g_A_I_J_r_vcl.handle2() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                          static_cast<unsigned int>(sizeof(cl_uint)*2*static_cast<unsigned int>(g_I.size() + 1)),
                                                                                      &(blocks_ind[0]));
              g_bv_r_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,  
                                                                        static_cast<unsigned int>(sizeof(ScalarType)*bv_size), 
                                                                        &(b_v_r[0]));
              g_bv_r_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                static_cast<unsigned int>(sizeof(cl_uint)*(g_I.size() + 1)), 
                                                                                &(start_bv_r_inds[0]));
              viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                          static_cast<unsigned int>(sizeof(cl_uint)*(g_is_update.size())),
                                                                                                            &(g_is_update[0]));
              viennacl::ocl::kernel& r_assembly_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_r_assembly");
              r_assembly_kernel.local_work_size(0, 1);
              r_assembly_kernel.global_work_size(0, 256);
              
              viennacl::ocl::enqueue(r_assembly_kernel(g_A_I_J_vcl.handle(), g_A_I_J_vcl.handle2(), g_A_I_J_vcl.handle1(), 
                                                      g_A_I_J_u_vcl.handle(), g_A_I_J_u_vcl.handle2(), g_A_I_J_u_vcl.handle1(), 
                                                      g_A_I_u_J_u_vcl.handle(), g_A_I_u_J_u_vcl.handle2(), g_A_I_u_J_u_vcl.handle1(), 
                                                      g_A_I_J_r_vcl.handle(), g_A_I_J_r_vcl.handle2(), g_A_I_J_r_vcl.handle1(),
                                                      g_is_update_vcl, static_cast<cl_uint>(g_I.size())));
              
              viennacl::ocl::kernel & bv_assembly_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_bv_assembly");
              bv_assembly_kernel.local_work_size(0, 1);
              bv_assembly_kernel.global_work_size(0, 256);
              viennacl::ocl::enqueue(bv_assembly_kernel(g_bv_vcl.handle(), g_bv_vcl.handle1(), g_A_I_J_vcl.handle1(), g_bv_vcl_u.handle(),
                                                        g_bv_vcl_u.handle1(), g_A_I_J_u_vcl.handle1(),
                                                        g_bv_r_vcl.handle(), g_bv_r_vcl.handle1(), g_A_I_J_r_vcl.handle1(), g_is_update_vcl,
                                                        static_cast<cl_uint>(g_I.size())));
              g_bv_vcl.handle() = g_bv_r_vcl.handle();
              g_bv_vcl.handle1() = g_bv_r_vcl.handle1();
              
              g_A_I_J_vcl.handle() = g_A_I_J_r_vcl.handle();
              g_A_I_J_vcl.handle2() = g_A_I_J_r_vcl.handle2();
              g_A_I_J_vcl.handle1() = g_A_I_J_r_vcl.handle1();
          }
          
          /** @brief GPU-based block update 
          * @param A sparse matrix
          * @param A_v_c vectorized column-wise initial matrix 
          * @param g_is_update container with identificators that shows which block should be modified 
          * @param g_res container of residuals for all columns
          * @param g_J container of column index sets for all columns
          * @param g_I container of row index sets for all columns
          * @param g_A_I_J_vcl container of block matrices from previous update  
          * @param g_bv_vcl block of beta vectors from previous iteration
          * @param tag SPAI configuration tag
          * @param cur_iter current iteration, used to make sure that kernel compiles just once 
          */ 
          template<typename ScalarType, unsigned int MAT_ALIGNMENT, typename SparseVectorType>
          void block_update(const viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT>& A, const std::vector<SparseVectorType>& A_v_c,
                            std::vector<cl_uint>& g_is_update,
                            std::vector<SparseVectorType>& g_res,
                            std::vector<std::vector<unsigned int> >& g_J, 
                            std::vector<std::vector<unsigned int> >& g_I,
                            block_matrix& g_A_I_J_vcl, 
                            block_vector& g_bv_vcl,
                            spai_tag const & tag,
                            const unsigned int cur_iter){
              //updated index set for columns
              std::vector<std::vector<unsigned int> > g_J_u(g_J.size());
              //updated index set for rows
              std::vector<std::vector<unsigned int> > g_I_u(g_J.size());
              //mixed index set of old and updated indices for rows 
              std::vector<std::vector<unsigned int> > g_I_q(g_J.size());
              //GPU memory for A_I_\hatJ
              block_matrix g_A_I_J_u_vcl;
              //GPU memory for A_\hatI_\hatJ
              block_matrix g_A_I_u_J_u_vcl;
              bool is_empty_block;
              //GPU memory for new b_v
              block_vector g_bv_u_vcl;
#ifdef _OPENMP
              #pragma omp parallel for
#endif              
              for(std::size_t i = 0; i < g_J.size(); ++i){
                  if(g_is_update[i]){
                      if(buildAugmentedIndexSet<SparseVectorType, ScalarType>(A_v_c, g_res[i], g_J[i], g_J_u[i], tag)){
                          buildNewRowSet(A_v_c, g_I[i], g_J_u[i], g_I_u[i]);
                      }
                  }
              }
              //assemble new A_I_J_u blocks on GPU and multiply them with Q'
              block_assembly(A, g_J_u, g_I, g_A_I_J_u_vcl, g_is_update, is_empty_block, cur_iter);
              //I have matrix A_I_J_u ready..
              block_q_multiplication<ScalarType>(g_J_u, g_I, g_A_I_J_vcl, g_bv_vcl, g_A_I_J_u_vcl, g_is_update, cur_iter);
              //assemble A_\hatI_\hatJ
              block_assembly(A, g_J_u, g_I_u, g_A_I_u_J_u_vcl, g_is_update, is_empty_block, cur_iter);
              assemble_qr_block<ScalarType>(g_J, g_I, g_J_u, g_I_u, g_I_q, g_A_I_J_u_vcl, g_A_I_J_vcl.handle1(),
                                            g_A_I_u_J_u_vcl, g_is_update, is_empty_block, cur_iter);
              
              block_qr<ScalarType>(g_I_q, g_J_u, g_A_I_u_J_u_vcl, g_bv_u_vcl, g_is_update, cur_iter);
              //concatanation of new and old indices
#ifdef _OPENMP
              #pragma omp parallel for
#endif              
              for(std::size_t i = 0; i < g_J.size(); ++i){
                  g_J[i].insert(g_J[i].end(), g_J_u[i].begin(), g_J_u[i].end());
                  g_I[i].insert(g_I[i].end(), g_I_u[i].begin(), g_I_u[i].end());
              }
              assemble_r<ScalarType>(g_I, g_J, g_A_I_J_vcl, g_A_I_J_u_vcl, g_A_I_u_J_u_vcl,  g_bv_vcl,  g_bv_u_vcl, g_is_update, cur_iter);
          }
        
        }        
      }        
    }
}
#endif