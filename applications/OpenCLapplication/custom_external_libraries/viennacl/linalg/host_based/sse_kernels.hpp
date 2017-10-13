#ifndef VIENNACL_LINALG_HOST_BASED_SSE_KERNELS_HPP_
#define VIENNACL_LINALG_HOST_BASED_SSE_KERNELS_HPP_

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

/** @file viennacl/linalg/host_based/sse_kernels.hpp
*   @brief optimized linear algebra operations for the CPU
*
*   Contributed by Alex Christensen.
*/

#ifdef VIENNACL_WITH_OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <vector>

//for std::min
#include <algorithm>

#include "viennacl/linalg/host_based/sse_blas.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace host_based
    {
      namespace detail
      {

        // returns true if the matrix is hermitian (or real symmetric), false otherwise
        template <typename ScalarType>
        bool isHermitian(ScalarType ** const A, vcl_size_t n)
        {
          for(vcl_size_t i=0;i<n;i++)
            for(vcl_size_t j=i;j<n;j++)
              if(A[i][j] != conjIfComplex(A[j][i]))
                return false;
          return true;
        }

        // returns the bandwidth of a hermitian (or real symmetric) matrix
        template <typename ScalarType>
        vcl_size_t getHermitianBandwidth(ScalarType ** const A, vcl_size_t n)
        {
          for(vcl_size_t i=n-1;i>=0;i--)
            for(vcl_size_t j=0;j<n-i;j++)
              if(A[i+j][j]!=ScalarType(0))
                return 2*i+1;
          return 0;
        }

        // helper for tridiagonalizeBandedMatrix
        // does a householder similarity transform to eliminate a range of nonzeros in a row of a hermitian matrix
        template <typename ScalarType>
        void eliminateHermitian(ScalarType ** A, vcl_size_t row, vcl_size_t from, vcl_size_t to, vcl_size_t width, ScalarType * ss)
        {
          if(from>=to)
            return;

          ScalarType norm=_nrm2(&A[row][row+from],to-from);

          if(norm != ScalarType(0))
          {

            //pick the better of two reflectors, to 1 or -1
            //this is wierd syntax that also works with std::complex
            if(std::abs(A[row][row+from]-ScalarType(1))>std::abs(A[row][row+from]+ScalarType(1)))
              norm=-norm;
            for(vcl_size_t i=row+from;i<row+to;i++)
              A[row][i]/=norm;
            A[row][row+from]+=ScalarType(1);

            //apply the similarity transformation

            //left transformation
            for(vcl_size_t j=row+1;j<row+width;j++)
            {
              ScalarType s=_dotc(to-from,&A[row][row+from],&A[j][row+from]);
              s=-s/A[row][row+from];
              _axpy(&A[row][row+from],&A[j][row+from],to-from,s);
            }

            //conjugate householder reflector for right transformation
            for(vcl_size_t i=row+from;i<row+to;i++)
              A[row][i]=conjIfComplex(A[row][i]);

            //right transformation (cache aligned)
            for(vcl_size_t i=0;i<width;i++)
              ss[i]=ScalarType(0);
            for(vcl_size_t i=from;i<to;i++)
              _axpy(&A[row+i][row],ss,width,conjIfComplex(A[row][row+i]));
            for(vcl_size_t i=0;i<width;i++)
              ss[i]=-ss[i]/A[row][row+from];
            for(vcl_size_t i=from;i<to;i++)
              _axpy(ss,&A[row+i][row],width,A[row][row+i]);

            //clean up the householder reflector
            for(vcl_size_t col=row+from;col<row+to;col++)
              A[row][col]=conjIfComplex(A[col][row]);

          }
        }

        // reduces a hermitian (or symmetric real) banded matrix to a hermitian (or symmetric real) tridiagonal matrix,
        // using householder similarity transforms, so eigenvalues are preserved.
        // bandwidth should be an odd integer, such as 3 for an already tridiagonal matrix
        // based on http://www.netlib.org/lapack/lawnspdf/lawn208.pdf
        template<typename ScalarType>
        void tridiagonalizeHermitianBandedMatrix(ScalarType ** A, vcl_size_t n, vcl_size_t bandwidth)
        {
          if(bandwidth<=3)
            return;

          vcl_size_t belowDiagonal=(bandwidth-1)/2;
          ScalarType *ss=new ScalarType[bandwidth+belowDiagonal];

          //eliminate and chase bulges where the elimination makes a bulge
          vcl_size_t k=0;
          for(;k<n-belowDiagonal;k++)
          {

              //eliminate below the diagonal
              eliminateHermitian(A,k,1,1+belowDiagonal,std::min(n-k,2*belowDiagonal+1),ss);

              //chase the bulge
              for(vcl_size_t bulgeStart=k+1;bulgeStart<n-belowDiagonal;bulgeStart+=belowDiagonal)
                  for(vcl_size_t i=0;i<belowDiagonal-1;i++)
                      eliminateHermitian(A,bulgeStart+i,belowDiagonal,std::min(n-bulgeStart-i,belowDiagonal*2-i),std::min(n-bulgeStart-i,bandwidth+belowDiagonal),ss);
          }

          //eliminate beyond where elimination makes bulges
          for(;k<n-2;k++)
              eliminateHermitian(A,k,1,n-k,n-k,ss);

          delete [] ss;
        }

        // reduces a hermitian (or symmetric real) matrix to a hermitian (or symmetric real) banded matrix with bandwidth 2*block_size+1
        // using householder similarity transformations, so eigenvalues are preserved. reduceToBandedMatrix(A,1) reduces the matrix to tridiagonal
        template<typename ScalarType>
        void reduceHermitianToBandedMatrix(ScalarType ** A, vcl_size_t n, vcl_size_t block_size, vcl_size_t num_threads)
        {
          ScalarType* norms=new ScalarType[block_size];
          ScalarType* ss=new ScalarType[n];

          for (vcl_size_t k=0;k<n-block_size;k+=block_size)
          {
            for(vcl_size_t bi=0;bi<std::min(block_size,n-k-block_size);bi++)
            {

              //this is the same as the norm of the column, since it's hermetian
              norms[bi]=_nrm2(&A[k+bi][k+bi+block_size],n-k-bi-block_size);

              if(norms[bi]!=ScalarType(0))
              {

                //pick the better of two reflectors, to 1 or -1
                //this is wierd syntax that also works with std::complex
                if(std::abs(A[k+bi][k+bi+block_size]-ScalarType(1))>std::abs(A[k+bi][k+bi+block_size]+ScalarType(1)))
                    norms[bi]=-norms[bi];
                for(vcl_size_t i=k+bi+block_size;i<n;i++)
                    A[k+bi][i]/=norms[bi];
                A[k+bi][k+bi+block_size]+=ScalarType(1);

                // Apply transformation to remaining rows within the block
                for(vcl_size_t j=k+bi+1;j<k+block_size;j++)
                {
                    ScalarType s=_dotc(n-k-bi-block_size,&A[k+bi][k+bi+block_size],&A[j][k+bi+block_size]);
                    s=-s/A[k+bi][k+bi+block_size];
                    _axpy(&A[k+bi][k+bi+block_size],&A[j][k+bi+block_size],n-k-bi-block_size,s);
                }
              }
            }

            //apply transformations from block to remaining rows and columns the block in parallel

            //left transformations
  #ifdef VIENNACL_WITH_OPENMP
  #pragma omp parallel for
            for(int j=k+block_size;j<(int)n;j++)
  #else
            for(vcl_size_t j=k+block_size;j<n;j++)
  #endif
            {
              for(vcl_size_t bi=0;bi<std::min(block_size,n-k-block_size);bi++)
              {
                if(norms[bi]!=ScalarType(0))
                {
                  ScalarType s=_dotc(n-k-bi-block_size,&A[k+bi][k+bi+block_size],&A[j][k+bi+block_size]);
                  s=-s/A[k+bi][k+bi+block_size];
                  _axpy(&A[k+bi][k+bi+block_size],&A[j][k+bi+block_size],n-k-bi-block_size,s);
                }
              }
            }

            //conjugate householder reflectors for right transformations
            for(vcl_size_t bi=0;bi<block_size;bi++)
              for(vcl_size_t i=k+bi+block_size;i<n;i++)
                A[k+bi][i]=conjIfComplex(A[k+bi][i]);

            //right transformations (cache aligned)
  #ifdef VIENNACL_WITH_OPENMP
  #pragma omp parallel for
            for(int section=0;section<(int)num_threads;section++)
  #else
            for(vcl_size_t section=0;section<num_threads;section++)
  #endif
            {
              vcl_size_t start=((n-k)*(section+0))/num_threads+k;
              vcl_size_t end  =((n-k)*(section+1))/num_threads+k;
              vcl_size_t length=end-start;
              for(vcl_size_t bi=0;bi<std::min(block_size,n-k-block_size);bi++)
              {
                if(norms[bi]!=ScalarType(0))
                {
                  for(vcl_size_t i=start;i<end;i++)
                    ss[i]=ScalarType(0);
                  for(vcl_size_t i=k+bi+block_size;i<n;i++)
                    _axpy(&A[i][start],ss+start,length,conjIfComplex(A[k+bi][i]));
                  for(vcl_size_t i=start;i<end;i++)
                    ss[i]=-ss[i]/A[k+bi][k+bi+block_size];
                  for(vcl_size_t i=k+bi+block_size;i<n;i++)
                    _axpy(ss+start,&A[i][start],length,A[k+bi][i]);
                }
              }
            }

            //clean up householder reflectors
            for(vcl_size_t row=k;row<k+block_size;row++)
              for(vcl_size_t col=row+block_size;col<n;col++)
                A[row][col]=conjIfComplex(A[col][row]);
          }
          delete [] norms;
          delete [] ss;
        }

      } //namespace detail

      /** @brief Inplace reduction of a dense n x n row-major or column-major hermitian (or real symmetric) matrix
      *         to tridiagonal form using householder similarity transforms (preserving eigenvalues)
      *
      * @param A            A dense hermitian matrix to be tridiagonalized
      * @param n            The height and width of the hermitian matrix
      * @param block_size   The block size to be used
      * @param num_threads  The number of threads to be used with OpenMP
      */
      template<typename ScalarType>
      void inplace_tred2(ScalarType ** A, vcl_size_t n, vcl_size_t block_size = 1, vcl_size_t num_threads = 1)
      {
        if(!detail::isHermitian(A,n))
          std::cerr << "ViennaCL: Warning in inplace_tred2(): Matrix is not hermitian (or real symmetric)" << std::endl;

        // Don't touch the whole matrix if the bandwidth is already small.
        // There's nothing numerically significant about n*4,
        // it's just a point I chose to switch to assuming the matrix is full.
        vcl_size_t bandwidth=detail::getHermitianBandwidth(A,n);
        if(bandwidth*bandwidth*num_threads<n*4 || 2*block_size+1>bandwidth)
          detail::tridiagonalizeHermitianBandedMatrix(A,n,bandwidth);
        else
        {
          detail::reduceHermitianToBandedMatrix(A,n,block_size,num_threads);
          detail::tridiagonalizeHermitianBandedMatrix(A,n,2*block_size+1);
        }
      }

      /** @brief Inplace lu factorization of an m x n dense row-major matrix with optional partial pivoting,
      *         returning true for an even number of pivots, false for an odd number of pivots.  Factorization
      *         is successful if there are no nonzero values on the diagonal.
      *
      * @param A            A dense row-major matrix to be factorized
      * @param m            The height of the matrix
      * @param n            The width of the matrix
      * @param piv          The optional pivot vector to store the pivot indices.  If piv is NULL, no partial pivoting will be performed.
      * @param block_size   The block size to be used
      */
      template <typename ScalarType>
      bool lu_factorize_row_major(ScalarType ** A, vcl_size_t m, vcl_size_t n, vcl_size_t * piv = NULL, vcl_size_t block_size = 8)
      {
        // Use a parallel "left-looking", row-operation-based, block Crout/Doolittle algorithm.
        if(piv)
          for(vcl_size_t i=0; i<m; i++)
            piv[i]=i;
        bool pivsign=true;

        // Outer loop.
            for(vcl_size_t j=0; j<std::min(m,n); j+=block_size)
        {
                  block_size=std::min(std::min(m-j,n-j),block_size);

          //do Gaussian elimination with partial pivoting in the block
          //(in the first few columns of the matrix)
          for(vcl_size_t bi=0;bi<block_size;bi++)
          {
            // Find pivot and exchange if necessary.
            vcl_size_t p=j+bi;
            if(piv)
            {
              for(vcl_size_t i=j+bi+1; i<m; i++)
                if(std::abs(A[i][j+bi])>std::abs(A[p][j+bi]))
                  p=i;

              if (p!=j+bi)
              {
                for(vcl_size_t k=0; k<n; k++)
                {
                  ScalarType t=A[p][k];
                  A[p][k]=A[j+bi][k];
                  A[j+bi][k]=t;
                }

                //swap pivot vector
                vcl_size_t k = piv[p];
                piv[p] = piv[j+bi];
                piv[j+bi] = k;
                pivsign = !pivsign;
              }
            }

            //eliminate below the diagonal in the block
            ScalarType elimVal=A[j+bi][j+bi];
            if(elimVal==ScalarType(0))
            {
              //apply previous transformations from the block to the top of the submatrix
              for(vcl_size_t row=j+1;row<j+bi;row++)
                for(vcl_size_t bi_=0;bi_<row-j;bi_++)
                  if(A[row][j+bi_]!=ScalarType(0))
                    _axpy(&(A[j+bi_][j+block_size]),&(A[row][j+block_size]),n-j-block_size,-A[row][j+bi_]);
              return pivsign;
            }
            for(vcl_size_t row=j+bi+1;row<m;row++)
            {
              ScalarType multiplier=A[row][j+bi]/elimVal;
                for(vcl_size_t col=j+bi;col<j+block_size;col++)
                  A[row][col]-=multiplier*A[j+bi][col];
                    A[row][j+bi]=multiplier;
            }
          }

          //at this point, the matrix looks something like this (if block size were 4)
          //
          //U U U U * * * *
          //L U U U * * * *
          //L L U U * * * *
          //L L L U * * * *
          //L L L L * * * *
          //L L L L * * * *
          //L L L L * * * *
          //L L L L * * * *

          //apply previous transformations from the block to the top of the submatrix
          for(vcl_size_t row=j+1;row<j+block_size;row++)
            for(vcl_size_t bi=0;bi<row-j;bi++)
              if(A[row][j+bi]!=ScalarType(0))
                _axpy(&(A[j+bi][j+block_size]),&(A[row][j+block_size]),n-j-block_size,-A[row][j+bi]);

          //at this point, the matrix looks something like this (if block size were 4)
          //
          //U U U U U U U U
          //L U U U U U U U
          //L L U U U U U U
          //L L L U U U U U
          //L L L L * * * *
          //L L L L * * * *
          //L L L L * * * *
          //L L L L * * * *

          //apply previous transformations from the block in parallel to the rest of the submatrix
  #ifdef VIENNACL_OPENMP
  #pragma omp parallel for
          for(int row=j+block_size;row<(int)m;row++)
  #else
          for(vcl_size_t row=j+block_size;row<m;row++)
  #endif
          for(vcl_size_t bi=0;bi<block_size;bi++)
            if(A[row][j+bi]!=ScalarType(0))
                _axpy(&(A[j+bi][j+block_size]),&(A[row][j+block_size]),n-j-block_size,-A[row][j+bi]);
        }
        return pivsign;
      }

      /** @brief Inplace qr factorization of an m x n dense column-major matrix, returning the householder normalization coefficients
      *
      * @param A            A dense column-major matrix to be factorized
      * @param m            The height of the matrix
      * @param n            The width of the matrix
      * @param block_size   The block size to be used
      */
      template <typename ScalarType>
          std::vector<ScalarType> inplace_qr_col_major(ScalarType ** A, vcl_size_t m, vcl_size_t n, vcl_size_t block_size = 8)
      {
        std::vector<ScalarType> betas(std::min(m,n));
        ScalarType* norms=new ScalarType[block_size];

        for(vcl_size_t k=0; k<std::min(m,n); k+=block_size)
        {
          block_size=std::min(std::min(m-k,n-k),block_size);

          for(vcl_size_t bi=0;bi<block_size;bi++)
          {

            // Compute 2-norm of k+bi-th column below the diagonal
            norms[bi]=_nrm2(&A[k+bi][k+bi],m-k-bi);

            if(norms[bi]!=ScalarType(0))
            {
              //pick the better of two reflectors, to 1 or -1,
              //this is wierd syntax that also works with std::complex
              if(std::abs(A[k+bi][k+bi]-ScalarType(1))>std::abs(A[k+bi][k+bi]+ScalarType(1)))
                norms[bi]*=-1;
              for(vcl_size_t i=k+bi;i<m;i++)
                A[k+bi][i]/=norms[bi];
              A[k+bi][k+bi]+=ScalarType(1);

              // Apply transformation to columns within the block
              for(vcl_size_t j=k+bi+1; j<k+block_size; j++)
              {
                ScalarType s=_dotc(m-k-bi,&A[k+bi][k+bi],&A[j][k+bi]);
                s = -s/A[k+bi][k+bi];
                _axpy(&A[k+bi][k+bi],&A[j][k+bi],m-k-bi,s);
              }
            }
            //temporarily store the diagonal value of R in betas
            betas[k+bi]=-norms[bi];
          }

          //apply transformations from block to remaining columns to the right of the block in parallel
  #ifdef VIENNACL_OPENMP
  #pragma omp parallel for
          for(int j=k+block_size; j<(int)n; j++)
  #else
          for(vcl_size_t j=k+block_size; j<n; j++)
  #endif
          {
            for(vcl_size_t bi=0;bi<block_size;bi++)
            {
              if(norms[bi]!=ScalarType(0))
                          {
                ScalarType s=_dotc(m-k-bi,&A[k+bi][k+bi],&A[j][k+bi]);
                s = -s/A[k+bi][k+bi];
                _axpy(&A[k+bi][k+bi],A[j]+k+bi,m-k-bi,s);
              }
            }
          }
        }

        //normalize the householder reflectors and store the betas
        for(vcl_size_t j=0;j<std::min(m,n);j++)
        {
          ScalarType beta=A[j][j];
          for(vcl_size_t i=j+1;i<m;i++)
            A[j][i]/=beta;
          A[j][j]=betas[j];//R diagonal values were stored temporarily in betas
          betas[j]=beta;
        }

        delete [] norms;
        return betas;
      }

      /** @brief Inplace qr factorization of an m x n dense row-major matrix, returning the householder normalization coefficients
      *
      * @param A            A dense row-major matrix to be factorized
      * @param m            The height of the matrix
      * @param n            The width of the matrix
      * @param block_size   The block size to be used
      * @param num_threads  Number of threads to be used
      */
      template <typename ScalarType>
      std::vector<ScalarType> inplace_qr_row_major(ScalarType ** A, vcl_size_t m, vcl_size_t n, vcl_size_t block_size = 8, vcl_size_t num_threads = 1)
      {
        std::vector<ScalarType> betas(std::min(m,n));
        ScalarType* norms=new ScalarType[block_size];
        ScalarType* ss=new ScalarType[n];

        //allocate O(m) memory for temporary column-major storage of the block for blas functions
        ScalarType** block_cols=new ScalarType*[block_size];
        for(vcl_size_t i=0;i<block_size;i++)
          block_cols[i]=new ScalarType[m];

        for(vcl_size_t k=0; k<std::min(m,n); k+=block_size)
        {
          block_size=std::min(std::min(m-k,n-k),block_size);

          //copy the block to column-major storage for cache alignment (necessary for _nrm2)
          for(vcl_size_t i=0;i<m-k;i++)
            for(vcl_size_t bi=0;bi<block_size;bi++)
              block_cols[bi][i]=A[k+i][k+bi];

          for(vcl_size_t bi=0;bi<block_size;bi++)
          {

            // Compute 2-norm of k+bi-th column below the diagonal
            norms[bi]=_nrm2(&block_cols[bi][bi],m-k-bi);

            if(norms[bi]!=ScalarType(0))
            {
              //pick the better of two reflectors, to 1 or -1,
              //this is wierd syntax that also works with std::complex
              if(std::abs(block_cols[bi][bi]-ScalarType(1))>std::abs(block_cols[bi][bi]+ScalarType(1)))
                norms[bi]*=-1;
              for(vcl_size_t i=bi;i<m-k;i++)
                block_cols[bi][i]/=norms[bi];
              block_cols[bi][bi]+=ScalarType(1);

              // Apply transformation to columns within the block
              for(vcl_size_t j=bi+1; j<block_size; j++)
              {
                ScalarType s=_dotc(m-k-bi,&block_cols[bi][bi],&block_cols[j][bi]);
                s = -s/block_cols[bi][bi];
                _axpy(&block_cols[bi][bi],&block_cols[j][bi],m-k-bi,s);
              }
            }
            //temporarily store the diagonal value of R in betas
            betas[k+bi]=-norms[bi];
          }

          //copy the block back to row-major storage
          for(vcl_size_t i=0;i<m-k;i++)
            for(vcl_size_t bi=0;bi<block_size;bi++)
              A[k+i][k+bi]=block_cols[bi][i];

          //apply transformations from block to remaining rows to the right of the block in parallel
  #ifdef VIENNACL_OPENMP
  #pragma omp parallel for
          for(int section=0;section<(int)num_threads;section++)
  #else
          for(vcl_size_t section=0;section<num_threads;section++)
  #endif
          {
            vcl_size_t start=((n-k-block_size)*(section+0))/num_threads+k+block_size;
            vcl_size_t end  =((n-k-block_size)*(section+1))/num_threads+k+block_size;
            vcl_size_t length=end-start;
            for(vcl_size_t bi=0;bi<block_size;bi++)
            {
              if(norms[bi]!=ScalarType(0))
              {
                for(vcl_size_t i=start;i<end;i++)
                  ss[i]=ScalarType(0);
                for(vcl_size_t i=k+bi;i<m;i++)
                  _axpy(&A[i][start],ss+start,length,A[i][k+bi]);
                for(vcl_size_t i=start;i<end;i++)
                  ss[i]=-ss[i]/A[k+bi][k+bi];
                for(vcl_size_t i=k+bi;i<m;i++)
                  _axpy(ss+start,&A[i][start],length,A[i][k+bi]);
              }
            }
          }
        }

        //normalize the householder reflectors and store the betas
        for(vcl_size_t j=0;j<std::min(m,n);j++)
        {
          ScalarType beta=A[j][j];
          for(vcl_size_t i=j+1;i<m;i++)
            A[i][j]/=beta;
          A[j][j]=betas[j];//R diagonal values were stored temporarily in betas
          betas[j]=beta;
        }

        delete [] norms;
        for(vcl_size_t i=0;i<block_size;i++)
          delete [] block_cols[i];
        delete [] block_cols;
        delete [] ss;

        return betas;
      }

    } //namespace host_based
  } //namespace linalg
} //namespace viennacl
#endif
