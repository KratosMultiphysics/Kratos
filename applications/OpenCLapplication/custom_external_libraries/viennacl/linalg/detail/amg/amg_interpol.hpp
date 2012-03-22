#ifndef VIENNACL_LINALG_DETAIL_AMG_AMG_INTERPOL_HPP
#define VIENNACL_LINALG_DETAIL_AMG_AMG_INTERPOL_HPP

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

/** @file amg_interpol.hpp
    @brief Implementations of several variants of the AMG interpolation operators (setup phase). Experimental.
*/

#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include "viennacl/linalg/amg.hpp"

#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "amg_debug.hpp"

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      namespace amg
      {
    
    /** @brief Calls the right function to build interpolation matrix
     * @param level    Coarse level identifier
     * @param A      Operator matrix on all levels
     * @param P      Prolongation matrices. P[level] is constructed
     * @param Pointvector  Vector of points on all levels
     * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_interpol(unsigned int level, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
  switch (tag.get_interpol())
  {
    case VIENNACL_AMG_INTERPOL_DIRECT: amg_interpol_direct (level, A, P, Pointvector, tag); break;
    case VIENNACL_AMG_INTERPOL_CLASSIC: amg_interpol_classic (level, A, P, Pointvector, tag); break;
    case VIENNACL_AMG_INTERPOL_AG: amg_interpol_ag (level, A, P, Pointvector, tag); break;
    case VIENNACL_AMG_INTERPOL_SA: amg_interpol_sa (level, A, P, Pointvector, tag); break;
  }
    } 
    /** @brief Direct interpolation. Multi-threaded! (VIENNACL_AMG_INTERPOL_DIRECT)
     * @param level    Coarse level identifier
     * @param A      Operator matrix on all levels
     * @param P      Prolongation matrices. P[level] is constructed
     * @param Pointvector  Vector of points on all levels
     * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_interpol_direct(unsigned int level, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      ScalarType temp_res;
      ScalarType row_sum, c_sum, diag;
      //int diag_sign;
      unsigned int x, y;
      amg_point *pointx, *pointy;
      unsigned int c_points = Pointvector[level].get_cpoints();

      // Setup Prolongation/Interpolation matrix
      P[level] = SparseMatrixType(A[level].size1(),c_points);
      P[level].clear();
      
      // Assign indices to C points
      Pointvector[level].build_index();
      
      // Direct Interpolation (Yang, p.14)
#ifdef _OPENMP
      #pragma omp parallel for private (pointx,pointy,row_sum,c_sum,temp_res,y,x,diag) shared (P,A,Pointvector,tag)
#endif      
      for (x=0; x < Pointvector[level].size(); ++x)
      {
  pointx = Pointvector[level][x];
  /*if (A[level](x,x) > 0) 
    diag_sign = 1;
  else
    diag_sign = -1;*/
  
  // When the current line corresponds to a C point then the diagonal coefficient is 1 and the rest 0
  if (pointx->is_cpoint())
    P[level](x,pointx->get_coarse_index()) = 1;
  
  // When the current line corresponds to a F point then the diagonal is 0 and the rest has to be computed (Yang, p.14)
  if (pointx->is_fpoint())
  {
    // Jump to row x
    InternalRowIterator row_iter = A[level].begin1();
    row_iter += x;
    
    // Row sum of coefficients (without diagonal) and sum of influencing C point coefficients has to be computed
    row_sum = c_sum = diag = 0;
    for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
    {
      y = col_iter.index2();
      if (x == y)// || *col_iter * diag_sign > 0)
      {
        diag += *col_iter;
        continue;
      }
      
      // Sum all other coefficients in line x
      row_sum += *col_iter;

      pointy = Pointvector[level][y];
      // Sum all coefficients that correspond to a strongly influencing C point
      if (pointy->is_cpoint())
        if (pointx->is_influencing(pointy))
    c_sum += *col_iter;        
    }
    temp_res = -row_sum/(c_sum*diag);

    // Iterate over all strongly influencing points of point x
    for (amg_point::iterator iter = pointx->begin_influencing(); iter != pointx->end_influencing(); ++iter)
    {    
      pointy = *iter;
      // The value is only non-zero for columns that correspond to a C point
      if (pointy->is_cpoint())
      {
        if (temp_res != 0)
    P[level](x, pointy->get_coarse_index()) = temp_res * A[level](x,pointy->get_index());
      }
    }
    
    //Truncate interpolation if chosen
    if (tag.get_interpolweight() != 0)
      amg_truncate_row(P[level], x, tag);
  }
      }
      
      // P test
      //test_interpolation(A[level], P[level], Pointvector[level]);
      
      #ifdef DEBUG
      std::cout << "Prolongation Matrix:" << std::endl;
      printmatrix (P[level]);
      #endif  
    }
    /** @brief Classical interpolation. Don't use with onepass classical coarsening or RS0 (Yang, p.14)! Multi-threaded! (VIENNACL_AMG_INTERPOL_CLASSIC)
     * @param level    Coarse level identifier
     * @param A      Operator matrix on all levels
     * @param P      Prolongation matrices. P[level] is constructed
     * @param Pointvector  Vector of points on all levels
     * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_interpol_classic(unsigned int level, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      ScalarType temp_res;
      ScalarType weak_sum, strong_sum;
      int diag_sign;
      amg_sparsevector<ScalarType> c_sum_row;
      amg_point *pointx, *pointy, *pointk, *pointm;
      unsigned int x, y, k, m;
      
      unsigned int c_points = Pointvector[level].get_cpoints();
      
      // Setup Prolongation/Interpolation matrix
      P[level] = SparseMatrixType(A[level].size1(), c_points);
      P[level].clear();
      
      // Assign indices to C points
      Pointvector[level].build_index();
      
      // Classical Interpolation (Yang, p.13-14)
#ifdef _OPENMP
      #pragma omp parallel for private (pointx,pointy,pointk,pointm,weak_sum,strong_sum,c_sum_row,temp_res,x,y,k,m,diag_sign) shared (A,P,Pointvector)
#endif      
      for (x=0; x < Pointvector[level].size(); ++x)
      {
  pointx = Pointvector[level][x];
  if (A[level](x,x) > 0) 
    diag_sign = 1;
  else
    diag_sign = -1;
  
  // When the current line corresponds to a C point then the diagonal coefficient is 1 and the rest 0
  if (pointx->is_cpoint())
    P[level](x,pointx->get_coarse_index()) = 1;

  // When the current line corresponds to a F point then the diagonal is 0 and the rest has to be computed (Yang, p.14)
  if (pointx->is_fpoint())
  {  
    // Jump to row x
    InternalRowIterator row_iter = A[level].begin1();
    row_iter += x;
    
    weak_sum = 0;
    c_sum_row = amg_sparsevector<ScalarType>(A[level].size1());
    c_sum_row.clear();
    for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
    {
      k = col_iter.index2();
      pointk = Pointvector[level][k];
      
      // Sum of weakly influencing neighbors + diagonal coefficient
      if (x == k || !pointx->is_influencing(pointk))// || *col_iter * diag_sign > 0)
      {
        weak_sum += *col_iter;
        continue;
      }
        
      // Sums of coefficients in row k (strongly influening F neighbors) of C point neighbors of x are calculated
      if (pointk->is_fpoint() && pointx->is_influencing(pointk))
      {
        for (amg_point::iterator iter = pointx->begin_influencing(); iter != pointx->end_influencing(); ++iter)
        {
    pointm = *iter;
    m = pointm->get_index();
    
    if (pointm->is_cpoint())
      // Only use coefficients that have opposite sign of diagonal.
      if (A[level](k,m) * diag_sign < 0)
        c_sum_row[k] += A[level](k,m);
        }
        continue;
      }
    }
    
    // Iterate over all strongly influencing points of point x
    for (amg_point::iterator iter = pointx->begin_influencing(); iter != pointx->end_influencing(); ++iter)
    {    
      pointy = *iter;
      y = pointy->get_index();
      
      // The value is only non-zero for columns that correspond to a C point
      if (pointy->is_cpoint())
      {
        strong_sum = 0;
        // Calculate term for strongly influencing F neighbors
        for (typename amg_sparsevector<ScalarType>::iterator iter2 = c_sum_row.begin(); iter2 != c_sum_row.end(); ++iter2)
        {
    k = iter2.index();
    // Only use coefficients that have opposite sign of diagonal.
    if (A[level](k,y) * diag_sign < 0)
      strong_sum += (A[level](x,k) * A[level](k,y)) / (*iter2);
        }
        
        // Calculate coefficient
        temp_res = - (A[level](x,y) + strong_sum) / (weak_sum);
        if (temp_res != 0)
    P[level](x,pointy->get_coarse_index()) = temp_res;   
      }
    }
    
    //Truncate iteration if chosen
    if (tag.get_interpolweight() != 0)
      amg_truncate_row(P[level], x, tag);
  }
      }
      
      #ifdef DEBUG
      std::cout << "Prolongation Matrix:" << std::endl;
      printmatrix (P[level]);
      #endif  
    }
    
    /** @brief Interpolation truncation (for VIENNACL_AMG_INTERPOL_DIRECT and VIENNACL_AMG_INTERPOL_CLASSIC)
    *
    * @param P    Interpolation matrix
    * @param row  Row which has to be truncated
    * @param tag  AMG preconditioner tag
    */
    template <typename SparseMatrixType>
    void amg_truncate_row(SparseMatrixType & P, unsigned int row, amg_tag & tag)
    {
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      ScalarType row_max, row_min, row_sum_pos, row_sum_neg, row_sum_pos_scale, row_sum_neg_scale;
      
      InternalRowIterator row_iter = P.begin1();
      row_iter += row;
      
      row_max = 0;
      row_min = 0;
      row_sum_pos = 0;
      row_sum_neg = 0;
      
      // Truncate interpolation by making values to zero that are a lot smaller than the biggest value in a row
      // Determine max entry and sum of row (seperately for negative and positive entries)
      for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
      {
  if (*col_iter > row_max)
    row_max = *col_iter;
  if (*col_iter < row_min)
    row_min = *col_iter;
  if (*col_iter > 0)
    row_sum_pos += *col_iter;
  if (*col_iter < 0)
    row_sum_neg += *col_iter;
      }
      
      row_sum_pos_scale = row_sum_pos;
      row_sum_neg_scale = row_sum_neg;
      
      // Make certain values to zero (seperately for negative and positive entries)
      for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
      {
  if (*col_iter > 0 && *col_iter < tag.get_interpolweight() * row_max)
  {
    row_sum_pos_scale -= *col_iter;
    *col_iter = 0;
  }
  if (*col_iter < 0 && *col_iter > tag.get_interpolweight() * row_min)
  {
    row_sum_pos_scale -= *col_iter;
    *col_iter = 0;
  }
      }
      
      // Scale remaining values such that row sum is unchanged
      for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
      {
  if (*col_iter > 0)
    *col_iter = *col_iter *(row_sum_pos/row_sum_pos_scale);
  if (*col_iter < 0)
    *col_iter = *col_iter *(row_sum_neg/row_sum_neg_scale);
      }
    }
    
    /** @brief AG (aggregation based) interpolation. Multi-Threaded! (VIENNACL_INTERPOL_SA)
     * @param level    Coarse level identifier
     * @param A      Operator matrix on all levels
     * @param P      Prolongation matrices. P[level] is constructed
     * @param Pointvector  Vector of points on all levels
     * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_interpol_ag(unsigned int level, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      unsigned int x;
      amg_point *pointx, *pointy;
      unsigned int c_points = Pointvector[level].get_cpoints();
      
      P[level] = SparseMatrixType(A[level].size1(), c_points);
      P[level].clear();
      
      // Assign indices to C points
      Pointvector[level].build_index();
      
      // Set prolongation such that F point is interpolated (weight=1) by the aggregate it belongs to (Vanek et al p.6)
#ifdef _OPENMP
      #pragma omp parallel for private (x,pointx) shared (P)
#endif      
      for (x=0; x<Pointvector[level].size(); ++x)
      {
  pointx = Pointvector[level][x];
  pointy = Pointvector[level][pointx->get_aggregate()];
  // Point x belongs to aggregate y.
  P[level](x,pointy->get_coarse_index()) = 1;
      }
      
      #ifdef DEBUG
      std::cout << "Aggregation based Prolongation:" << std::endl;
      printmatrix(P[level]);
      #endif
    }
      
    /** @brief SA (smoothed aggregate) interpolation. Multi-Threaded! (VIENNACL_INTERPOL_SA)
     * @param level    Coarse level identifier
     * @param A      Operator matrix on all levels
     * @param P      Prolongation matrices. P[level] is constructed
     * @param Pointvector  Vector of points on all levels
     * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_interpol_sa(unsigned int level, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      unsigned int x,y;
      ScalarType diag;
      unsigned int c_points = Pointvector[level].get_cpoints();
           
      InternalType1 P_tentative = InternalType1(P.size());
      SparseMatrixType Jacobi = SparseMatrixType(A[level].size1(), A[level].size2());
      Jacobi.clear();
      P[level] = SparseMatrixType(A[level].size1(), c_points);
      P[level].clear();      
           
      // Build Jacobi Matrix via filtered A matrix (Vanek et al. p.6)
#ifdef _OPENMP
      #pragma omp parallel for private (x,y,diag) shared (A,Pointvector)
#endif      
      for (x=0; x<A[level].size1(); ++x)
      {
  diag = 0;
  InternalRowIterator row_iter = A[level].begin1();
  row_iter += x;
  for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
  {
    y = col_iter.index2();
    // Determine the structure of the Jacobi matrix by using a filtered matrix of A:
    // The diagonal consists of the diagonal coefficient minus all coefficients of points not in the neighborhood of x.
    // All other coefficients are the same as in A.
    // Already use Jacobi matrix to save filtered A matrix to speed up computation.
    if (x == y)
      diag += *col_iter;
    else if (!Pointvector[level][x]->is_influencing(Pointvector[level][y]))
      diag += -*col_iter;
    else
      Jacobi (x,y) = *col_iter;      
  }
  InternalRowIterator row_iter2 = Jacobi.begin1();
  row_iter2 += x;
  // Traverse through filtered A matrix and compute the Jacobi filtering
  for (InternalColIterator col_iter2 = row_iter2.begin(); col_iter2 != row_iter2.end(); ++col_iter2)
  {
      *col_iter2 = - tag.get_interpolweight()/diag * *col_iter2;
  }
  // Diagonal can be computed seperately.
  Jacobi (x,x) = 1 - tag.get_interpolweight();
      }
          
      #ifdef DEBUG
      std::cout << "Jacobi Matrix:" << std::endl;
      printmatrix(Jacobi);
      #endif
      
      // Use AG interpolation as tentative prolongation
      amg_interpol_ag(level, A, P_tentative, Pointvector, tag);
      
      #ifdef DEBUG
      std::cout << "Tentative Prolongation:" << std::endl;
      printmatrix(P_tentative[level]);
      #endif
      
      // Multiply Jacobi matrix with tentative prolongation to get actual prolongation
      amg_mat_prod(Jacobi,P_tentative[level],P[level]);
      
      #ifdef DEBUG
      std::cout << "Prolongation Matrix:" << std::endl;
      printmatrix (P[level]);
      #endif    
    }
      } //namespace amg
    }
  }
}

#endif
