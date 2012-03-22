#ifndef VIENNACL_LINALG_DETAIL_SPAI_QR_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_QR_HPP

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

/** @file viennacl/linalg/detail/spai/qr.hpp
    @brief Implementation of a simultaneous QR factorization of multiple matrices. Experimental.
    
    SPAI code contributed by Nikolay Lukash
*/

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

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"

#include "viennacl/linalg/detail/spai/block_matrix.hpp"
#include "viennacl/linalg/detail/spai/block_vector.hpp"
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
        
          
          
          //********** DEBUG FUNCTIONS *****************//
          template< typename T, typename InputIterator>
          void Print(std::ostream& ostr, InputIterator it_begin, InputIterator it_end){
              //std::ostream_iterator<int> it_os(ostr, delimiter);
              std::string delimiters = " ";
              std::copy(it_begin, it_end, std::ostream_iterator<T>(ostr, delimiters.c_str()));
              ostr<<std::endl;
          }
          
          template<typename VectorType, typename MatrixType>
          void write_to_block(VectorType& con_A_I_J, unsigned int start_ind,  const std::vector<unsigned int>& I, const std::vector<unsigned int>& J, MatrixType& m){
              m.resize(I.size(), J.size(), false);
              for(size_t i = 0; i < J.size(); ++i){
                  for(size_t j = 0; j < I.size(); ++j){
                      m(j,i) = con_A_I_J[start_ind + i*I.size() + j];
                  }
              }
          }
          
          template<typename VectorType>
          void print_continious_matrix(VectorType& con_A_I_J, std::vector<cl_uint>& blocks_ind,
                                      const std::vector<std::vector<unsigned int> >& g_I, const std::vector<std::vector<unsigned int> >& g_J){
              typedef typename VectorType::value_type ScalarType;
              std::vector<boost::numeric::ublas::matrix<ScalarType> > com_A_I_J(g_I.size());
              for(size_t i = 0; i < g_I.size(); ++i){
                  write_to_block( con_A_I_J, blocks_ind[i], g_I[i], g_J[i], com_A_I_J[i]);
                  std::cout<<com_A_I_J[i]<<std::endl;
              }
          }
          template<typename VectorType>
          void print_continious_vector(VectorType& con_v, std::vector<cl_uint>& block_ind, const std::vector<std::vector<unsigned int> >& g_J){
              typedef typename VectorType::value_type ScalarType;
              std::vector<boost::numeric::ublas::vector<ScalarType> > com_v(g_J.size());
              //Print<ScalarType>(std::cout, con_v.begin(), con_v.end());
              for(size_t i = 0; i < g_J.size(); ++i){
                  com_v[i].resize(g_J[i].size());
                  for(size_t j = 0; j < g_J[i].size(); ++j){
                      com_v[i](j) = con_v[block_ind[i] + j];
                  }
                  std::cout<<com_v[i]<<std::endl;
              }
          }
          
          ///**************************************** BLOCK FUNCTIONS ************************************//
          /** @brief Computes size of elements, start indices and matrix dimensions for a certain block
          * @param g_I container of row indices 
          * @param g_J container of column indices 
          * @param sz general size for all elements in a certain block
          * @param blocks_ind start indices in a certain
          * @param matrix_dims matrix dimensions for each block
          */ 
          void compute_blocks_size(const std::vector<std::vector<unsigned int> >& g_I, const std::vector<std::vector<unsigned int> >& g_J, 
                                  unsigned int& sz, std::vector<cl_uint>& blocks_ind, std::vector<cl_uint>& matrix_dims){
              sz = 0;
              for(size_t i = 0; i < g_I.size(); ++i){
                  sz += static_cast<unsigned int>(g_I[i].size()*g_J[i].size());
                  matrix_dims[2*i] = static_cast<cl_uint>(g_I[i].size());
                  matrix_dims[2*i + 1] = static_cast<cl_uint>(g_J[i].size());
                  blocks_ind[i+1] = blocks_ind[i] + static_cast<cl_uint>(g_I[i].size()*g_J[i].size());
                  
              }
          }
          /** @brief Computes size of particular container of index set 
          * @param inds container of index sets 
          * @param size output size 
          */ 
          void get_size(const std::vector<std::vector<unsigned int> >& inds, unsigned int& size){
              size = 0;
              for (size_t i = 0; i < inds.size(); ++i) {
                  size += static_cast<unsigned int>(inds[i].size());
              }
          }
          
          /** @brief Initializes start indices of particular index set 
          * @param inds container of index sets 
          * @param start_inds output index set 
          */
          void init_start_inds(const std::vector<std::vector<unsigned int> >& inds, std::vector<cl_uint>& start_inds){
              for(size_t i = 0; i < inds.size(); ++i){
                  start_inds[i+1] = start_inds[i] + static_cast<cl_uint>(inds[i].size());
              }
          }

          //*************************************  QR FUNCTIONS  ***************************************//
          /** @brief Dot prod of particular column of martix A with it's self starting at a certain index beg_ind 
          * @param A init matrix 
          * @param beg_ind starting index
          * @param res result of dot product
          */
          template<typename MatrixType, typename ScalarType>
          void dot_prod(const MatrixType& A,  unsigned int beg_ind, ScalarType& res){
              res = static_cast<ScalarType>(0);
              for(size_t i = beg_ind; i < A.size1(); ++i){
                  res += A(i, beg_ind-1)*A(i, beg_ind-1);
              }
          }
          /** @brief Dot prod of particular matrix column with arbitrary vector: A(:, col_ind)
          * @param A init matrix
          * @param v input vector
          * @param col_ind starting column index
          * @param start_ind starting index inside column
          * @param res result of dot product
          */
          template<typename MatrixType, typename VectorType, typename ScalarType>
          void custom_inner_prod(const MatrixType& A, const VectorType& v, unsigned int col_ind, unsigned int start_ind, ScalarType& res){
              res = static_cast<ScalarType>(0);
              for(unsigned int i = start_ind; i < static_cast<unsigned int>(A.size1()); ++i){
                  res += A(i, col_ind)*v(i);  
              }
          }
          
          /** @brief Copying part of matrix column 
          * @param A init matrix
          * @param v output vector
          * @param beg_ind start index for copying
          */
          template<typename MatrixType, typename VectorType>
          void copy_vector(const MatrixType & A, VectorType & v, const unsigned int beg_ind){
              for(unsigned int i = beg_ind; i < static_cast<unsigned int>(A.size1()); ++i){
                  v(i) = A( i, beg_ind-1);
              }
          }
          
          //householder reflection c.f. Gene H. Golub, Charles F. Van Loan "Matrix Computations" 3rd edition p.210
          /** @brief Coputation of Householder vector, householder reflection c.f. Gene H. Golub, Charles F. Van Loan "Matrix Computations" 3rd edition p.210
          * @param A init matrix
          * @param j start index for computations 
          * @param v output Householder vector 
          * @param b beta
          */
          template<typename MatrixType, typename VectorType, typename ScalarType>
          void householder_vector(const MatrixType& A, unsigned int j, VectorType& v, ScalarType& b){
              ScalarType sg;
              //
              dot_prod(A, j+1, sg); 
              copy_vector(A, v, j+1);
              ScalarType mu;
              v(j) = static_cast<ScalarType>(1.0);
              if(sg == 0){
                  b = 0;
              }
              else{
                  mu = std::sqrt(A(j,j)*A(j, j) + sg);
                  if(A(j, j) <= 0){
                      v(j) = A(j, j) - mu;
                  }else{
                      v(j) = -sg/(A(j, j) + mu);
                  }
                  b = 2*(v(j)*v(j))/(sg + v(j)*v(j));
                  v = v/v(j);
              }
          }
          /** @brief Inplace application of Householder vector to a matrix A
          * @param A init matrix
          * @param iter_cnt current iteration
          * @param v Householder vector 
          * @param b beta
          */
          template<typename MatrixType, typename VectorType, typename ScalarType>
          void apply_householder_reflection(MatrixType& A, unsigned int iter_cnt, VectorType& v, ScalarType b){
              //update every column of matrix A
              ScalarType in_prod_res;
              for(unsigned int i = iter_cnt; i < static_cast<unsigned int>(A.size2()); ++i){
                  //update each column in a fashion: ai = ai - b*v*(v'*ai)
                  custom_inner_prod(A, v, i, iter_cnt, in_prod_res);
                  for(unsigned int j = iter_cnt; j < static_cast<unsigned int>(A.size1()); ++j){
                      A(j, i) -= b*in_prod_res*v(j);
                  }
              }
          }
          
          /** @brief Storage of vector v in column(A, ind), starting from ind-1 index of a column
          * @param A init matrix
          * @param ind index of a column
          * @param v vector that should be stored
          */
          template<typename MatrixType, typename VectorType>
          void store_householder_vector(MatrixType& A, unsigned int ind, VectorType& v){
              for(unsigned int i = ind; i < static_cast<unsigned int>(A.size1()); ++i){
                  A(i, ind-1) = v(i);
              }
          }
          
          
          //QR algorithm 
          /** @brief Inplace QR factorization via Householder reflections c.f. Gene H. Golub, Charles F. Van Loan "Matrix Computations" 3rd edition p.224
          * @param R input matrix 
          * @param b_v vector of betas
          */
          template<typename MatrixType, typename VectorType>
          void single_qr(MatrixType& R, VectorType& b_v){
              typedef typename MatrixType::value_type ScalarType;
              if((R.size1() > 0) && (R.size2() > 0)){
                  VectorType v = (VectorType)boost::numeric::ublas::zero_vector<ScalarType>(R.size1());
                  b_v = (VectorType)boost::numeric::ublas::zero_vector<ScalarType>(R.size2());
                  for(unsigned int i = 0; i < static_cast<unsigned int>(R.size2()); ++i){
                      householder_vector(R, i, v, b_v[i]);
                      apply_householder_reflection(R, i, v, b_v[i]);
                      if(i < R.size1()) store_householder_vector(R, i+1, v);
                  }
              }
          }
          
          //********************** HELP FUNCTIONS FOR GPU-based QR factorization *************************//
          /** @brief Reading from text file into string
          * @param file_name file name
          * @param kernel_source string that contains file
          */
          void read_kernel_from_file(std::string& file_name, std::string& kernel_source){
              std::ifstream ifs(file_name.c_str(), std::ifstream::in);
              
              if (!ifs)
                std::cerr << "WARNING: Cannot open file " << file_name << std::endl;
              
              std::string line;
              std::ostringstream ost;
              while (std::getline(ifs, line)) {
                  ost<<line<<std::endl;
              }
              kernel_source = ost.str();
          }
          
          /** @brief Getting max size of rows/columns from container of index set
          * @param inds container of index set
          * @param max_size max size that corresponds to that container
          */
          void get_max_block_size(const std::vector<std::vector<unsigned int> >& inds, unsigned int& max_size){
              max_size = 0;
              for(unsigned int i = 0; i < inds.size(); ++i){
                  if(inds[i].size() > max_size){
                      max_size = static_cast<unsigned int>(inds[i].size());
                  }
              }
          }
          
          /** @brief Dot_prod(column(A, ind), v) starting from index ind+1 
          * @param A input matrix
          * @param v input vector
          * @param ind index
          * @param res result value
          */
          template<typename MatrixType, typename VectorType, typename ScalarType>
          void custom_dot_prod(const MatrixType& A, const VectorType& v, unsigned int ind, ScalarType& res){
              res = static_cast<ScalarType>(0);
              for(unsigned int j = ind; j < A.size1(); ++j){
                  if(j == ind){
                      res += v(j);
                  }else{
                      res += A(j, ind)*v(j);
                  }
              }
          }
          
          /** @brief Recovery Q from matrix R and vector of betas b_v
          * @param R input matrix
          * @param b_v vector of betas
          * @param y output vector
          */
          template<typename MatrixType, typename VectorType>
          void apply_q_trans_vec(const MatrixType& R, const VectorType& b_v, VectorType& y){
              typedef typename MatrixType::value_type ScalarType;
              ScalarType inn_prod = static_cast<ScalarType>(0);
              for(size_t i = 0; i < R.size2(); ++i){
                  custom_dot_prod(R, y, static_cast<unsigned int>(i), inn_prod);
                  for(size_t j = i; j < R.size1(); ++j){
                      if(i == j){
                          y(j) -= b_v(i)*inn_prod;
                      }
                      else{
                          y(j) -= b_v(i)*inn_prod*R(j,i);
                      }
                  }
              }
          }
          
          /** @brief Multiplication of Q'*A, where Q is in implicit for lower part of R and vector of betas - b_v
          * @param R input matrix
          * @param b_v vector of betas
          * @param A output matrix
          */
          template<typename MatrixType, typename VectorType>
          void apply_q_trans_mat(const MatrixType& R, const VectorType& b_v, MatrixType& A){
              VectorType tmp_v;
              for(size_t i = 0; i < A.size2(); ++i){
                  tmp_v = (VectorType)column(A,i);
                  apply_q_trans_vec(R, b_v, tmp_v);
                  column(A,i) = tmp_v;
              }
          }
          
          //parallel QR for GPU
          /** @brief Inplace QR factorization via Householder reflections c.f. Gene H. Golub, Charles F. Van Loan "Matrix Computations" 3rd edition p.224 performed
                      on GPU
          * @param g_I container of row indices
          * @param g_J container of column indices
          * @param g_A_I_J_vcl contigious matrices, GPU memory is used
          * @param g_bv_vcl contigios vectors beta, GPU memory is used 
          * @param g_is_update container of indicators that show active blocks  
          * @param cur_iter current iteration
          */
          template<typename ScalarType>
          void block_qr(std::vector<std::vector<unsigned int> >& g_I, 
                        std::vector<std::vector<unsigned int> >& g_J, 
                        block_matrix& g_A_I_J_vcl,
                        block_vector& g_bv_vcl,
                        std::vector<cl_uint>& g_is_update,
                        const unsigned int cur_iter){
              //typedef typename MatrixType::value_type ScalarType;
              unsigned int bv_size;
              unsigned int v_size;
              //set up arguments for GPU
              //find maximum size of rows/columns
              unsigned int local_r_n, local_c_n;
              //find max size for blocks
              get_max_block_size(g_I, local_r_n);
              get_max_block_size(g_J, local_c_n);
              //get size
              get_size(g_J, bv_size);
              get_size(g_I, v_size);
              //get start indices
              std::vector<cl_uint> start_bv_inds(g_I.size() + 1, 0);
              std::vector<cl_uint> start_v_inds(g_I.size() + 1, 0);
              init_start_inds(g_J, start_bv_inds);
              init_start_inds(g_I, start_v_inds);
              //init arrays
              std::vector<ScalarType> b_v(bv_size, static_cast<ScalarType>(0));
              std::vector<ScalarType> v(v_size, static_cast<ScalarType>(0));
              //call qr program
              block_vector v_vcl;
              /*if(cur_iter == 0)
              {
                  //if first run - compile the program
                  std::string qr_kernel_file_name = "kernels/spai/qr3_a_n.cl";
                  std::string qr_kernel_source;
                  read_kernel_from_file(qr_kernel_file_name, qr_kernel_source);
                  viennacl::ocl::program & qr_prog = viennacl::ocl::current_context().add_program(qr_kernel_source.c_str(), "qr_kernel_source");
                  qr_prog.add_kernel("block_qr");
                  //
              }*/
              
              g_bv_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,  
                                                                      static_cast<unsigned int>(sizeof(ScalarType)*bv_size), 
                                                                      &(b_v[0]));
              
              v_vcl.handle() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,  
                                                                                                  static_cast<unsigned int>(sizeof(ScalarType)*v_size), 
                                                                                                  &(v[0]));
              //the same as j_start_inds
              g_bv_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                static_cast<unsigned int>(sizeof(cl_uint)*g_I.size()), 
                                                                                &(start_bv_inds[0]));
              
              v_vcl.handle1() = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                                static_cast<unsigned int>(sizeof(cl_uint)*g_I.size()), 
                                                                                &(start_v_inds[0]));
              viennacl::ocl::handle<cl_mem> g_is_update_vcl = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                                static_cast<unsigned int>(sizeof(cl_uint)*g_is_update.size()),
                                                                                &(g_is_update[0]));
              //local memory
              //viennacl::ocl::enqueue(k(vcl_vec, size, viennacl::ocl::local_mem(sizeof(SCALARTYPE) * k.local_work_size()), temp));
              viennacl::ocl::kernel& qr_kernel = viennacl::ocl::get_kernel(viennacl::linalg::kernels::spai<ScalarType, 1>::program_name(), "block_qr");
              qr_kernel.local_work_size(0, local_c_n);
              qr_kernel.global_work_size(0, 256);
              viennacl::ocl::enqueue(qr_kernel(g_A_I_J_vcl.handle(), g_A_I_J_vcl.handle1(), g_bv_vcl.handle(), 
                                              v_vcl.handle(), g_A_I_J_vcl.handle2(), 
                                              g_bv_vcl.handle1(), v_vcl.handle1(), g_is_update_vcl,
                                              viennacl::ocl::local_mem(static_cast<unsigned int>(sizeof(ScalarType)*(local_r_n*local_c_n))),
                                              static_cast<cl_uint>(g_I.size())));
              
          }
        }
      }
    }
}
#endif