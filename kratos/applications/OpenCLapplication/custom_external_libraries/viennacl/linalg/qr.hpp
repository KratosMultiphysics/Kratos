#ifndef VIENNACL_LINALG_QR_HPP
#define VIENNACL_LINALG_QR_HPP

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

/** @file viennacl/linalg/qr.hpp
    @brief Proivdes a QR factorization using a block-based approach.  Experimental in 1.2.x.
*/

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "boost/numeric/ublas/vector_proxy.hpp"
#include "boost/numeric/ublas/io.hpp"
#include "boost/numeric/ublas/matrix_expression.hpp"

#include "viennacl/matrix.hpp"
#include "viennacl/matrix_proxy.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/range.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      
      // orthogonalises j-th column of A
      template <typename MatrixType, typename VectorType>
      typename MatrixType::value_type setup_householder_vector(MatrixType const & A, VectorType & v, std::size_t j)
      {
        typedef typename MatrixType::value_type   ScalarType;
        
        //compute norm of column below diagonal:
        ScalarType sigma = 0;
        ScalarType beta = 0;
        for (std::size_t k = j+1; k<A.size1(); ++k)
          sigma += A(k, j) * A(k, j);

        //get v from A:
        v[j] = 1;
        //ScalarType scaling = sqrt(sigma + A(j,j)*A(j,j));
        //ScalarType scaling = sqrt(sigma);
        ScalarType scaling = 1.0;
        for (std::size_t k = j+1; k<A.size1(); ++k)
          v[k] = A(k, j) / scaling;
        sigma = sigma / (scaling * scaling);
        ScalarType A_jj = A(j,j) / scaling;
        
        std::cout << "sigma: " << sigma << std::endl;
        assert( sigma >= 0.0  && "sigma must be non-negative!");

        
        if (sigma == 0)
          return 0;
        else
        {
          ScalarType mu = sqrt(sigma + A_jj*A_jj);
          std::cout << "mu: " << mu << std::endl;
          std::cout << "sigma: " << sigma << std::endl;
          
          ScalarType v1;
          if (A_jj <= 0)
            v1 = A_jj - mu;
          else
            v1 = -sigma / (A_jj + mu);
          
          beta = 2.0 * v1 * v1 / (sigma + v1 * v1);
          
          //divide v by its diagonal element v[j]
          v[j] = 1;
          std::cout << "v1: " << v1 << std::endl;
          for (std::size_t k = j+1; k<A.size1(); ++k)
            v[k] /= v1;
        }
          
        return beta;
      }
      

      template <typename MatrixType, typename VectorType>
      typename MatrixType::value_type setup_householder_vector_ublas(MatrixType const & A, VectorType & v, MatrixType & matrix_1x1, std::size_t j)
      {
        using boost::numeric::ublas::range;
        using boost::numeric::ublas::project;
        
        typedef typename MatrixType::value_type   ScalarType;
        
        //compute norm of column below diagonal:
        //ScalarType sigma = 0;
        //for (std::size_t k = j+1; k<A.size1(); ++k)
        //  sigma += A(k, j) * A(k, j);
        matrix_1x1 = prod( trans(project(A, range(j+1, A.size1()), range(j, j+1))),
                                 project(A, range(j+1, A.size1()), range(j, j+1))
                         );
        ScalarType sigma = matrix_1x1(0,0);
        ScalarType beta = 0;
        ScalarType A_jj = A(j,j);
        
        assert( sigma >= 0.0  && "sigma must be non-negative!");

        //get v from A:
        //for (std::size_t k = j+1; k<A.size1(); ++k)
        //  v[k] = A(k, j);
        v(j,0) = 1.0;
        project(v, range(j+1, A.size1()), range(0,1)) = project(A, range(j+1, A.size1()), range(j,j+1));
        
        if (sigma == 0)
          return 0;
        else
        {
          ScalarType mu = sqrt(sigma + A_jj*A_jj);
          //std::cout << "mu: " << mu << std::endl;
          //std::cout << "sigma: " << sigma << std::endl;
          
          ScalarType v1;
          if (A_jj <= 0)
            v1 = A_jj - mu;
          else
            v1 = -sigma / (A_jj + mu);
          
          beta = 2.0 * v1 * v1 / (sigma + v1 * v1);
          
          //divide v by its diagonal element v[j]
          //v[j] = 1;
          //for (std::size_t k = j+1; k<A.size1(); ++k)
          //  v[k] /= v1;
          project(v, range(j+1, A.size1()), range(0,1)) /= v1;
        }
          
        return beta;
      }


      template <typename MatrixType, typename VectorType>
      typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type 
      setup_householder_vector_viennacl(MatrixType const & A, VectorType & v, MatrixType & matrix_1x1, std::size_t j)
      {
        //using boost::numeric::ublas::range;
        //using boost::numeric::ublas::project;
        using viennacl::range;
        using viennacl::project;
        
        typedef typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type   ScalarType;
        
        //compute norm of column below diagonal:
        //ScalarType sigma = 0;
        //for (std::size_t k = j+1; k<A.size1(); ++k)
        //  sigma += A(k, j) * A(k, j);
        matrix_1x1 = viennacl::linalg::prod( trans(project(A, range(j+1, A.size1()), range(j, j+1))),
                                                   project(A, range(j+1, A.size1()), range(j, j+1))
                                           );
        ScalarType sigma = matrix_1x1(0,0);
        ScalarType beta = 0;
        ScalarType A_jj = A(j,j);

        //std::cout << "sigma: " << sigma << std::endl;
        assert( sigma >= 0.0  && "sigma must be non-negative!");


        //get v from A:
        //for (std::size_t k = j+1; k<A.size1(); ++k)
        //  v[k] = A(k, j);
        v(j,0) = 1.0;
        project(v, range(j+1, A.size1()), range(0,1)) = project(A, range(j+1, A.size1()), range(j,j+1));
        
        if (sigma == 0)
          return 0;
        else
        {
          ScalarType mu = sqrt(sigma + A_jj*A_jj);
          //std::cout << "mu: " << mu << std::endl;
          //std::cout << "sigma: " << sigma << std::endl;
          
          ScalarType v1;
          if (A_jj <= 0)
            v1 = A_jj - mu;
          else
            v1 = -sigma / (A_jj + mu);
          
          beta = 2.0 * v1 * v1 / (sigma + v1 * v1);
          
          //divide v by its diagonal element v[j]
          //v[j] = 1;
          //for (std::size_t k = j+1; k<A.size1(); ++k)
          //  v[k] /= v1;
          //v(j,0) = 1.0;
          project(v, range(j+1, A.size1()), range(0,1)) /= v1;
        }
          
        return beta;
      }


      // Apply (I - beta v v^T) to the k-th column of A, where v is the reflector starting at j-th row/column
      template <typename MatrixType, typename VectorType, typename ScalarType>
      void householder_reflect(MatrixType & A, VectorType & v, ScalarType beta, std::size_t j, std::size_t k)
      {
        ScalarType v_in_col = A(j,k);
        for (std::size_t i=j+1; i<A.size1(); ++i)
          v_in_col += v[i] * A(i,k);

        assert(v[j] == 1.0);
        //std::cout << "v[]: " << v[0] << ", " << v[1] << ", " << v[2] << std::endl;
        //std::cout << "v_in_col: " << v_in_col << std::endl;
        
        for (std::size_t i=j; i<A.size1(); ++i)
          A(i,k) -= beta * v_in_col * v[i];
      }

      template <typename MatrixType, typename VectorType, typename ScalarType>
      void householder_reflect_ublas(MatrixType & A, VectorType & v, MatrixType & matrix_1x1, ScalarType beta, std::size_t j, std::size_t k)
      {
        using boost::numeric::ublas::range;
        using boost::numeric::ublas::project;
        
        ScalarType v_in_col = A(j,k);
        //for (std::size_t i=j+1; i<A.size1(); ++i)
        //  v_in_col += v[i] * A(i,k);

        matrix_1x1 = prod(trans(project(v, range(j+1, A.size1()), range(0, 1))),
                         project(A, range(j+1, A.size1()), range(k,k+1)));
        v_in_col += matrix_1x1(0,0);
                         
        //for (std::size_t i=j; i<A.size1(); ++i)
        //  A(i,k) -= beta * v_in_col * v[i];
        
        project(A, range(j, A.size1()), range(k, k+1)) -= (beta * v_in_col) * project(v, range(j, A.size1()), range(0, 1));
      }

      template <typename MatrixType, typename VectorType, typename ScalarType>
      void householder_reflect_viennacl(MatrixType & A, VectorType & v, MatrixType & matrix_1x1, ScalarType beta, std::size_t j, std::size_t k)
      {
        //using boost::numeric::ublas::range;
        //using boost::numeric::ublas::project;
        using viennacl::range;
        using viennacl::project;
        
        ScalarType v_in_col = A(j,k);
        //for (std::size_t i=j+1; i<A.size1(); ++i)
        //  v_in_col += v[i] * A(i,k);

        matrix_1x1 = viennacl::linalg::prod(trans(project(v, range(j+1, A.size1()), range(0, 1))),
                                                  project(A, range(j+1, A.size1()), range(k,k+1)));
        v_in_col += matrix_1x1(0,0);
                         
        //for (std::size_t i=j; i<A.size1(); ++i)
        //  A(i,k) -= beta * v_in_col * v[i];
        
        if ( beta * v_in_col != 0.0)
        {
          VectorType temp = project(v, range(j, A.size1()), range(0, 1));
          project(v, range(j, A.size1()), range(0, 1)) *= (beta * v_in_col);
          project(A, range(j, A.size1()), range(k, k+1)) -= project(v, range(j, A.size1()), range(0, 1));
          project(v, range(j, A.size1()), range(0, 1)) = temp;
        }
      }


      // Apply (I - beta v v^T) to A, where v is the reflector starting at j-th row/column
      template <typename MatrixType, typename VectorType, typename ScalarType>
      void householder_reflect(MatrixType & A, VectorType & v, ScalarType beta, std::size_t j)
      {
        std::size_t column_end = A.size2();
        
        for (std::size_t k=j; k<column_end; ++k) //over columns
          householder_reflect(A, v, beta, j, k);
      }
      
      
      template <typename MatrixType, typename VectorType>
      void write_householder_to_A(MatrixType & A, VectorType const & v, std::size_t j)
      {
        for (std::size_t i=j+1; i<A.size1(); ++i)
          A(i,j) = v[i];
      }
      
      template <typename MatrixType, typename VectorType>
      void write_householder_to_A_ublas(MatrixType & A, VectorType const & v, std::size_t j)
      {
        //for (std::size_t i=j+1; i<A.size1(); ++i)
        //  A(i,j) = v[i];
        using boost::numeric::ublas::range;
        using boost::numeric::ublas::project;
        
        //VectorType temp = project(v, range(j+1, A.size1()));
        project( A, range(j+1, A.size1()), range(j, j+1) ) = project(v, range(j+1, A.size1()), range(0, 1) );;
      }

      template <typename MatrixType, typename VectorType>
      void write_householder_to_A_viennacl(MatrixType & A, VectorType const & v, std::size_t j)
      {
        //for (std::size_t i=j+1; i<A.size1(); ++i)
        //  A(i,j) = v[i];
        //using boost::numeric::ublas::range;
        //using boost::numeric::ublas::project;
        using viennacl::range;
        using viennacl::project;
        
        //VectorType temp = project(v, range(j+1, A.size1()));
        project( A, range(j+1, A.size1()), range(j, j+1) ) = project(v, range(j+1, A.size1()), range(0, 1) );;
      }

      
      /*template<typename MatrixType>
      std::vector<typename MatrixType::value_type> qr(MatrixType & A)
      {
        typedef typename MatrixType::value_type   ScalarType;
        
        std::vector<ScalarType> betas(A.size2());
        std::vector<ScalarType> v(A.size1());

        //copy A to Q:
        for (size_t j=0; j<A.size2(); ++j)
        {
            betas[j] = setup_householder_vector(A, v, j);
            householder_reflect(A, v, betas[j], j);
            write_householder_to_A(A, v, j);
        }
        
        return betas;
      }*/
      
      
      
      
      class range
      {
        public:
          range(std::size_t start, std::size_t end) : start_(start), end_(end) {}
          
          std::size_t lower() const { return start_; }
          std::size_t upper() const { return end_; }
          
        private:
          std::size_t start_;
          std::size_t end_;
      };

      template <typename MatrixType>
      class sub_matrix
      {
        public:
          typedef typename MatrixType::value_type value_type;
          
          sub_matrix(MatrixType & mat,
                      range row_range,
                      range col_range) : mat_(mat), row_range_(row_range), col_range_(col_range) {}
                      
          value_type operator()(size_t row, size_t col) const
          {
            assert(row < size1());
            assert(col < size2());
            return mat_(row + row_range_.lower(), col + col_range_.lower()); 
          }
                      
          std::size_t size1() const { return row_range_.upper() - row_range_.lower(); }
          std::size_t size2() const { return col_range_.upper() - col_range_.lower(); }
          
        private:
          MatrixType & mat_;
          range row_range_;
          range col_range_;
      };


      //computes C = prod(A, B)
      template <typename MatrixTypeA, typename MatrixTypeB, typename MatrixTypeC>
      void prod_AA(MatrixTypeA const & A, MatrixTypeB const & B, MatrixTypeC & C)
      {
        assert(C.size1() == A.size1());
        assert(A.size2() == B.size1());
        assert(B.size2() == C.size2());
        
        typedef typename MatrixTypeC::value_type   ScalarType;
        
        for (std::size_t i=0; i<C.size1(); ++i)
        {
          for (std::size_t j=0; j<C.size2(); ++j)
          {
            ScalarType val = 0;
            for (std::size_t k=0; k<A.size2(); ++k)
              val += A(i, k) * B(k, j);
            C(i, j) = val;
          }
        }
      }
      
      //computes C = prod(A^T, B)
      template <typename MatrixTypeA, typename MatrixTypeB, typename MatrixTypeC>
      void prod_TA(MatrixTypeA const & A, MatrixTypeB const & B, MatrixTypeC & C)
      {
        assert(C.size1() == A.size2());
        assert(A.size1() == B.size1());
        assert(B.size2() == C.size2());
        
        typedef typename MatrixTypeC::value_type   ScalarType;
        
        for (std::size_t i=0; i<C.size1(); ++i)
        {
          for (std::size_t j=0; j<C.size2(); ++j)
          {
            ScalarType val = 0;
            for (std::size_t k=0; k<A.size1(); ++k)
              val += A(k, i) * B(k, j);
            C(i, j) = val;
          }
        }
      }
      

    } //namespace detail
        


    //takes an inplace QR matrix A and generates Q and R explicitly
    template <typename MatrixType, typename VectorType>
    void recoverQ(MatrixType const & A, VectorType const & betas, MatrixType & Q, MatrixType & R)
    {
      typedef typename MatrixType::value_type   ScalarType;
      
      std::vector<ScalarType> v(A.size1());

      Q.clear();
      R.clear();

      //
      // Recover R from upper-triangular part of A:
      //
      std::size_t i_max = std::min(R.size1(), R.size2());
      for (std::size_t i=0; i<i_max; ++i)
        for (std::size_t j=i; j<R.size2(); ++j)
          R(i,j) = A(i,j);
      
      //
      // Recover Q by applying all the Householder reflectors to the identity matrix:
      //
      for (std::size_t i=0; i<Q.size1(); ++i)
        Q(i,i) = 1.0;

      std::size_t j_max = std::min(A.size1(), A.size2());
      for (std::size_t j=0; j<j_max; ++j)
      {
        std::size_t col_index = j_max - j - 1;
        v[col_index] = 1.0;
        for (std::size_t i=col_index+1; i<A.size1(); ++i)
          v[i] = A(i, col_index);
        
        /*std::cout << "Recovery with beta = " << betas[col_index] << ", j=" << col_index << std::endl;
        std::cout << "v: ";
        for (size_t i=0; i<v.size(); ++i)
          std::cout << v[i] << ", ";
        std::cout << std::endl;*/

        if (betas[col_index] != 0)
          detail::householder_reflect(Q, v, betas[col_index], col_index);
      }
    }


    /** @brief Implementation of inplace-QR factorization for a general Boost.uBLAS compatible matrix A 
     * 
     * @param A            A dense compatible to Boost.uBLAS
     * @param block_size   The block size to be used. The number of columns of A must be a multiple of block_size
     */
    template<typename MatrixType>
    std::vector<typename MatrixType::value_type> inplace_qr_ublas(MatrixType & A, std::size_t block_size = 32)
    {
      typedef typename MatrixType::value_type   ScalarType;
      typedef boost::numeric::ublas::matrix_range<MatrixType>  MatrixRange;
      
      using boost::numeric::ublas::range;
      using boost::numeric::ublas::project;
      
      std::vector<ScalarType> betas(A.size2());
      //boost::numeric::ublas::vector<ScalarType> v(A.size1());
      MatrixType v(A.size1(), 1);
      MatrixType matrix_1x1(1,1);

      MatrixType Y(A.size1(), block_size); Y.clear(); Y.resize(A.size1(), block_size);
      MatrixType W(A.size1(), block_size); W.clear(); W.resize(A.size1(), block_size);
        
      //run over A in a block-wise manner:
      for (std::size_t j = 0; j < std::min(A.size1(), A.size2()); j += block_size)
      {
        //determine Householder vectors:
        for (std::size_t k = 0; k < block_size; ++k)
        {
          betas[j+k] = detail::setup_householder_vector_ublas(A, v, matrix_1x1, j+k);
          
          for (std::size_t l = k; l < block_size; ++l)
            detail::householder_reflect_ublas(A, v, matrix_1x1, betas[j+k], j+k, j+l);

          detail::write_householder_to_A_ublas(A, v, j+k);
        }

        //
        // Setup Y:
        //
        Y.clear();  Y.resize(A.size1(), block_size);
        for (std::size_t k = 0; k < block_size; ++k)
        {
          //write Householder to Y:
          Y(j+k,k) = 1.0;
          project(Y, range(j+k+1, A.size1()), range(k, k+1)) = project(A, range(j+k+1, A.size1()), range(j+k, j+k+1));
        }
        
        //
        // Setup W:
        //
        
        //first vector:
        W.clear();  W.resize(A.size1(), block_size);
        W(j, 0) = -betas[j];
        project(W, range(j+1, A.size1()), range(0, 1)) = -betas[j] * project(A, range(j+1, A.size1()), range(j, j+1));
        
        
        //k-th column of W is given by -beta * (Id + W*Y^T) v_k, where W and Y have k-1 columns
        for (std::size_t k = 1; k < block_size; ++k)
        {
          MatrixRange Y_old = project(Y, range(j, A.size1()), range(0, k));
          MatrixRange v_k   = project(Y, range(j, A.size1()), range(k, k+1));
          MatrixRange W_old = project(W, range(j, A.size1()), range(0, k));
          MatrixRange z     = project(W, range(j, A.size1()), range(k, k+1));
          
          MatrixType YT_prod_v = boost::numeric::ublas::prod(boost::numeric::ublas::trans(Y_old), v_k);
          z = - betas[j+k] * (v_k + prod(W_old, YT_prod_v));
        }

        //
        //apply (I+WY^T)^T = I + Y W^T to the remaining columns of A:
        //
        
        if (A.size2() - j - block_size > 0)
        {
          
          MatrixRange A_part(A, range(j, A.size1()), range(j+block_size, A.size2()));
          MatrixRange W_part(W, range(j, A.size1()), range(0, block_size));
          MatrixType temp = boost::numeric::ublas::prod(trans(W_part), A_part);
          
          A_part += prod(project(Y, range(j, A.size1()), range(0, Y.size2())),
                         temp);
        }
      }

      return betas;
    }


    /** @brief Implementation of a OpenCL-only QR factorization for GPUs (or multi-core CPU) 
     * 
     * Performance is rather poor at small matrix sizes.
     * Prefer the use of the hybrid version, which is automatically chosen using the interface function inplace_qr()
     * 
     * @param A            A dense ViennaCL matrix to be factored
     * @param block_size   The block size to be used. The number of columns of A must be a multiple of block_size
     */
    template<typename MatrixType>
    std::vector< typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type > 
    inplace_qr_viennacl(MatrixType & A, std::size_t block_size = 16)
    {
      typedef typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type   ScalarType;
      typedef viennacl::matrix_range<MatrixType>  MatrixRange;
      
      //using boost::numeric::ublas::range;
      //using boost::numeric::ublas::project;
      using viennacl::range;
      using viennacl::project;
      
      std::vector<ScalarType> betas(A.size2());
      //boost::numeric::ublas::vector<ScalarType> v(A.size1());
      MatrixType v(A.size1(), 1);
      MatrixType matrix_1x1(1,1);

      MatrixType Y(A.size1(), block_size); Y.clear();
      MatrixType W(A.size1(), block_size); W.clear();

      MatrixType YT_prod_v(block_size, 1);
      MatrixType z(A.size1(), 1);      
      
      //run over A in a block-wise manner:
      for (std::size_t j = 0; j < std::min(A.size1(), A.size2()); j += block_size)
      {
        
        //determine Householder vectors:
        for (std::size_t k = 0; k < block_size; ++k)
        {
          betas[j+k] = detail::setup_householder_vector_viennacl(A, v, matrix_1x1, j+k);
          for (std::size_t l = k; l < block_size; ++l)
            detail::householder_reflect_viennacl(A, v, matrix_1x1, betas[j+k], j+k, j+l);

          detail::write_householder_to_A_viennacl(A, v, j+k);
        }

        //
        // Setup Y:
        //
        Y.clear();
        for (std::size_t k = 0; k < block_size; ++k)
        {
          //write Householder to Y:
          Y(j+k,k) = 1.0;
          project(Y, range(j+k+1, A.size1()), range(k, k+1)) = project(A, range(j+k+1, A.size1()), range(j+k, j+k+1));
        }
        
        //
        // Setup W:
        //
        
        //first vector:
        W.clear();
        W(j, 0) = -betas[j];
        //project(W, range(j+1, A.size1()), range(0, 1)) = -betas[j] * project(A, range(j+1, A.size1()), range(j, j+1));
        project(W, range(j+1, A.size1()), range(0, 1)) = project(A, range(j+1, A.size1()), range(j, j+1));
        project(W, range(j+1, A.size1()), range(0, 1)) *= -betas[j];
        
        
        //k-th column of W is given by -beta * (Id + W*Y^T) v_k, where W and Y have k-1 columns
        for (std::size_t k = 1; k < block_size; ++k)
        {
          MatrixRange Y_old = project(Y, range(j, A.size1()), range(0, k));
          MatrixRange v_k   = project(Y, range(j, A.size1()), range(k, k+1));
          MatrixRange W_old = project(W, range(j, A.size1()), range(0, k));
          //MatrixRange z     = project(W, range(0, A.size1()), range(k, k+1));
         
          //std::cout << "should: " << k << std::endl;
          project(YT_prod_v, range(0, k), range(0,1)) = prod(trans(Y_old), v_k);
          project(z, range(j, A.size1()), range(0,1)) = prod(W_old, project(YT_prod_v, range(0, k), range(0,1)));
          //project(W, range(0, A.size1()), range(k, k+1)) = - betas[j+k] * (v_k + prod(W_old, YT_prod_v));
          project(W, range(j, A.size1()), range(k, k+1)) = project(z, range(j, A.size1()), range(0,1));
          project(W, range(j, A.size1()), range(k, k+1)) += v_k;
          project(W, range(j, A.size1()), range(k, k+1)) *= - betas[j+k];
        }

        //
        //apply (I+WY^T)^T = I + Y W^T to the remaining columns of A:
        //
        
        if (A.size2() - j - block_size > 0)
        {
          
          MatrixRange A_part(A, range(j, A.size1()), range(j+block_size, A.size2()));
          MatrixRange W_part(W, range(j, A.size1()), range(0, block_size));
          MatrixType temp = prod(trans(W_part), A_part);
          
          A_part += prod(project(Y, range(j, A.size1()), range(0, Y.size2())),
                         temp);
        }
      }

      return betas;
    }






    //MatrixType is ViennaCL-matrix
    /** @brief Implementation of a hybrid QR factorization using uBLAS on the CPU and ViennaCL for GPUs (or multi-core CPU) 
     * 
     * Prefer the use of the convenience interface inplace_qr()
     * 
     * @param A            A dense ViennaCL matrix to be factored
     * @param block_size   The block size to be used. The number of columns of A must be a multiple of block_size
     */
    template<typename MatrixType>
    std::vector< typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type > 
    inplace_qr_hybrid(MatrixType & A, std::size_t block_size = 16)
    {
      typedef typename viennacl::result_of::cpu_value_type< typename MatrixType::value_type >::type   ScalarType;

      typedef viennacl::matrix_range<MatrixType>                    VCLMatrixRange;
      typedef boost::numeric::ublas::matrix<ScalarType>             UblasMatrixType;
      typedef boost::numeric::ublas::matrix_range<UblasMatrixType>  UblasMatrixRange;
      
      //using boost::numeric::ublas::range;
      //using boost::numeric::ublas::project;
      
      std::vector<ScalarType> betas(A.size2());
      UblasMatrixType v(A.size1(), 1);
      UblasMatrixType matrix_1x1(1,1);

      UblasMatrixType ublasW(A.size1(), block_size); ublasW.clear(); ublasW.resize(A.size1(), block_size);
      UblasMatrixType ublasY(A.size1(), block_size); ublasY.clear(); ublasY.resize(A.size1(), block_size);
      
      UblasMatrixType ublasA(A.size1(), A.size1());
      
      MatrixType vclW(ublasW.size1(), ublasW.size2());
      MatrixType vclY(ublasY.size1(), ublasY.size2());
      
        
      //run over A in a block-wise manner:
      for (std::size_t j = 0; j < std::min(A.size1(), A.size2()); j += block_size)
      {
        UblasMatrixRange ublasA_part = boost::numeric::ublas::project(ublasA,
                                                                      boost::numeric::ublas::range(0, A.size1()),
                                                                      boost::numeric::ublas::range(j, j+block_size));
        viennacl::copy(viennacl::project(A,
                                         viennacl::range(0, A.size1()),
                                         viennacl::range(j, j+block_size)),
                       ublasA_part
                      );
        
        //determine Householder vectors:
        for (std::size_t k = 0; k < block_size; ++k)
        {
          betas[j+k] = detail::setup_householder_vector_ublas(ublasA, v, matrix_1x1, j+k);
          
          for (std::size_t l = k; l < block_size; ++l)
            detail::householder_reflect_ublas(ublasA, v, matrix_1x1, betas[j+k], j+k, j+l);

          detail::write_householder_to_A_ublas(ublasA, v, j+k);
        }

        //
        // Setup Y:
        //
        ublasY.clear();  ublasY.resize(A.size1(), block_size);
        for (std::size_t k = 0; k < block_size; ++k)
        {
          //write Householder to Y:
          ublasY(j+k,k) = 1.0;
          boost::numeric::ublas::project(ublasY, 
                                         boost::numeric::ublas::range(j+k+1, A.size1()), 
                                         boost::numeric::ublas::range(k, k+1)) 
            = boost::numeric::ublas::project(ublasA, 
                                             boost::numeric::ublas::range(j+k+1, A.size1()),
                                             boost::numeric::ublas::range(j+k, j+k+1));
        }
        
        //
        // Setup W:
        //
        
        //first vector:
        ublasW.clear();  ublasW.resize(A.size1(), block_size);
        ublasW(j, 0) = -betas[j];
        boost::numeric::ublas::project(ublasW, 
                                       boost::numeric::ublas::range(j+1, A.size1()), 
                                       boost::numeric::ublas::range(0, 1)) 
           = -betas[j] * boost::numeric::ublas::project(ublasA, 
                                                        boost::numeric::ublas::range(j+1, A.size1()), 
                                                        boost::numeric::ublas::range(j, j+1));
        
        
        //k-th column of W is given by -beta * (Id + W*Y^T) v_k, where W and Y have k-1 columns
        for (std::size_t k = 1; k < block_size; ++k)
        {
          UblasMatrixRange Y_old = boost::numeric::ublas::project(ublasY,
                                                                  boost::numeric::ublas::range(j, A.size1()),
                                                                  boost::numeric::ublas::range(0, k));
          UblasMatrixRange v_k   = boost::numeric::ublas::project(ublasY,
                                                                  boost::numeric::ublas::range(j, A.size1()),
                                                                  boost::numeric::ublas::range(k, k+1));
          UblasMatrixRange W_old = boost::numeric::ublas::project(ublasW, 
                                                                  boost::numeric::ublas::range(j, A.size1()), 
                                                                  boost::numeric::ublas::range(0, k));
          UblasMatrixRange z     = boost::numeric::ublas::project(ublasW, 
                                                                  boost::numeric::ublas::range(j, A.size1()), 
                                                                  boost::numeric::ublas::range(k, k+1));
          
          UblasMatrixType YT_prod_v = boost::numeric::ublas::prod(boost::numeric::ublas::trans(Y_old), v_k);
          z = - betas[j+k] * (v_k + prod(W_old, YT_prod_v));
        }
        
        

        //
        //apply (I+WY^T)^T = I + Y W^T to the remaining columns of A:
        //
        
        VCLMatrixRange A_part = viennacl::project(A,
                                                  viennacl::range(0, A.size1()),
                                                  viennacl::range(j, j+block_size));
        
        viennacl::copy(boost::numeric::ublas::project(ublasA,
                                                      boost::numeric::ublas::range(0, A.size1()),
                                                      boost::numeric::ublas::range(j, j+block_size)),
                       A_part);
        
        viennacl::copy(ublasW, vclW);
        viennacl::copy(ublasY, vclY);
        
        if (A.size2() - j - block_size > 0)
        {
          
          VCLMatrixRange A_part(A, range(j, A.size1()), range(j+block_size, A.size2()));
          VCLMatrixRange W_part(vclW, range(j, A.size1()), range(0, block_size));
          MatrixType temp = viennacl::linalg::prod(trans(W_part), A_part);
          
          A_part += prod(viennacl::project(vclY, 
                                           viennacl::range(j, A.size1()), 
                                           viennacl::range(0, vclY.size2())),
                         temp);
        }
      }

      return betas;
    }



    /** @brief Overload of inplace-QR factorization of a ViennaCL matrix A 
     * 
     * @param A            A dense ViennaCL matrix to be factored
     * @param block_size   The block size to be used. The number of columns of A must be a multiple of block_size
     */
    template<typename T, typename F, unsigned int ALIGNMENT>
    std::vector<T> inplace_qr(viennacl::matrix<T, F, ALIGNMENT> & A, std::size_t block_size = 16)
    {
      if (A.size2() % block_size != 0)
        std::cerr << "ViennaCL: Warning in inplace_qr(): Number of columns is not a multiple of the block size" << std::endl;
      
      return inplace_qr_hybrid(A, block_size);
    }

    /** @brief Overload of inplace-QR factorization for a general Boost.uBLAS compatible matrix A 
     * 
     * @param A            A dense compatible to Boost.uBLAS
     * @param block_size   The block size to be used. The number of columns of A must be a multiple of block_size
     */
    template<typename MatrixType>
    std::vector<typename MatrixType::value_type> inplace_qr(MatrixType & A, std::size_t block_size = 16)
    {
      if (A.size2() % block_size != 0)
        std::cerr << "ViennaCL: Warning in inplace_qr(): Number of columns is not a multiple of the block size" << std::endl;
      
      return inplace_qr_ublas(A, block_size);
    }


        
  } //linalg
} //viennacl


#endif
