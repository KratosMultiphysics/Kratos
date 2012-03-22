#ifndef VIENNACL_LINALG_DETAIL_AMG_AMG_COARSE_HPP
#define VIENNACL_LINALG_DETAIL_AMG_AMG_COARSE_HPP

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

/** @file amg_coarse.hpp
    @brief Implementations of several variants of the AMG coarsening procedure (setup phase). Experimental.
*/

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
    
    /** @brief Calls the right coarsening procedure
      * @param level    Coarse level identifier
      * @param A    Operator matrix on all levels
      * @param Pointvector   Vector of points on all levels
      * @param Slicing    Partitioning of the system matrix to different processors (only used in RS0 and RS3)
      * @param tag    AMG preconditioner tag
      */
    template <typename InternalType1, typename InternalType2, typename InternalType3>
    void amg_coarse(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, InternalType3 & Slicing, amg_tag & tag)
    {
  switch (tag.get_coarse())
  {
    case VIENNACL_AMG_COARSE_RS: amg_coarse_classic (level, A, Pointvector, tag); break;
    case VIENNACL_AMG_COARSE_ONEPASS: amg_coarse_classic_onepass (level, A, Pointvector, tag); break;
    case VIENNACL_AMG_COARSE_RS0: amg_coarse_rs0 (level, A, Pointvector, Slicing, tag); break;
    case VIENNACL_AMG_COARSE_RS3: amg_coarse_rs3 (level, A, Pointvector, Slicing, tag); break;
    case VIENNACL_AMG_COARSE_AG:   amg_coarse_ag (level, A, Pointvector, tag); break;
  }
    } 
    
    /** @brief Determines strong influences in system matrix, classical approach (RS). Multithreaded!
    * @param level    Coarse level identifier
    * @param A      Operator matrix on all levels
    * @param Pointvector   Vector of points on all levels
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_influence(unsigned int level, InternalType1 const & A, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::value_type ScalarType;
      typedef typename SparseMatrixType::const_iterator1 ConstRowIterator;
      typedef typename SparseMatrixType::const_iterator2 ConstColIterator;
      
      ScalarType max;
      int diag_sign;
      //unsigned int i;
        
#ifdef _OPENMP
      #pragma omp parallel for private (max,diag_sign) shared (A,Pointvector)
#endif      
      for (unsigned int i=0; i<A[level].size1(); ++i)
      {  
  diag_sign = 1;
  if (A[level](i,i) < 0)
    diag_sign = -1;
  
  ConstRowIterator row_iter = A[level].begin1();
  row_iter += i;
  // Find greatest non-diagonal negative value (positive if diagonal is negative) in row
  max = 0;
  for (ConstColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
  {
      if (i == (unsigned int) col_iter.index2()) continue;
      if (diag_sign == 1)
        if (max > *col_iter)  max = *col_iter;
      if (diag_sign == -1)
        if (max < *col_iter)  max = *col_iter;
  }
  
  // If maximum is 0 then the row is independent of the others
  if (max == 0)
    continue;
  
  // Find all points that strongly influence current point (Yang, p.5)
  for (ConstColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
  {
    unsigned int j = col_iter.index2();  
    if (i == j) continue;
    if (diag_sign * (-*col_iter) >= tag.get_threshold() * (diag_sign * (-max)))
    {
      // Strong influence from j to i found, save information
      Pointvector[level][i]->add_influencing_point(Pointvector[level][j]);
    }
  }
      }
      
      #ifdef DEBUG
      std::cout << "Influence Matrix: " << std::endl;
      boost::numeric::ublas::matrix<bool> mat;
      Pointvector[level].get_influence_matrix(mat);
      printmatrix (mat);
      #endif
      
      // Save influenced points
      for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
      {
  for (typename amg_point::iterator iter2 = (*iter)->begin_influencing(); iter2 != (*iter)->end_influencing(); ++iter2)
  {
    (*iter2)->add_influenced_point(*iter);
  }
      }
        
      #ifdef DEBUG
      std::cout << "Influence Measures: " << std::endl;
      boost::numeric::ublas::vector<unsigned int> temp;
      Pointvector[level].get_influence(temp);
      printvector (temp);
      std::cout << "Point Sorting: " << std::endl;
      Pointvector[level].get_sorting(temp);
      printvector (temp);
      #endif 
    }
        
    /** @brief Classical (RS) one-pass coarsening. Single-Threaded! (VIENNACL_AMG_COARSE_CLASSIC_ONEPASS)
    * @param level     Course level identifier
    * @param A      Operator matrix on all levels
    * @param Pointvector   Vector of points on all levels
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_coarse_classic_onepass(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      amg_point* c_point, *point1, *point2;
      unsigned int i;
        
      // Check and save all strong influences
      amg_influence (level, A, Pointvector, tag);    
      
      // Traverse through points and calculate initial influence measure
#ifdef _OPENMP
      #pragma omp parallel for private (i) shared (Pointvector)
#endif      
      for (i=0; i<Pointvector[level].size(); ++i)
  Pointvector[level][i]->calc_influence();
      
       // Do initial sorting
      Pointvector[level].sort();
      
      // Get undecided point with highest influence measure
      while ((c_point = Pointvector[level].get_nextpoint()) != NULL)
      {    
  // Make this point C point
  Pointvector[level].make_cpoint(c_point);
  
  // All strongly influenced points become F points
  for (typename amg_point::iterator iter = c_point->begin_influenced(); iter != c_point->end_influenced(); ++iter)
  {
    point1 = *iter;
    // Found strong influence from C point (c_point influences point1), check whether point is still undecided, otherwise skip
    if (!point1->is_undecided()) continue;
    // Make this point F point if it is still undecided point
    Pointvector[level].make_fpoint(point1);
    
    // Add +1 to influence measure for all undecided points that strongly influence new F point
    for (typename amg_point::iterator iter2 = point1->begin_influencing(); iter2 != point1->end_influencing(); ++iter2)
    {
      point2 = *iter2;
      // Found strong influence to F point (point2 influences point1)
      if (point2->is_undecided())
        Pointvector[level].add_influence(point2,1);
    }
  }
      }
      
      // If a point is neither C nor F point but is nevertheless influenced by other points make it F point
      // (this situation can happen when this point does not influence other points and the points that influence this point became F points already)
      /*#pragma omp parallel for private (i,point1)
      for (i=0; i<Pointvector[level].size(); ++i)
      {
  point1 = Pointvector[level][i];
  if (point1->is_undecided())
  {
    // Undecided point found. Check whether it is influenced by other point and if so: Make it F point.
    if (point1->number_influencing() > 0)
    {
      #pragma omp critical
      Pointvector[level].make_fpoint(point1);
    }
  }
      }*/

      #if defined (DEBUG)//  or defined (DEBUGBENCH)
      unsigned int c_points = Pointvector[level].get_cpoints();
      unsigned int f_points = Pointvector[level].get_fpoints();
      std::cout << "1st pass: Level " << level << ": ";
      std::cout << "No of C points = " << c_points << ", ";
      std::cout << "No of F points = " << f_points << std::endl;
      #endif

      #ifdef DEBUG
      std::cout << "Coarse Points:" << std::endl;
      boost::numeric::ublas::vector<bool> C;
      Pointvector[level].get_C(C);
      printvector (C);
      std::cout << "Fine Points:" << std::endl;
      boost::numeric::ublas::vector<bool> F;
      Pointvector[level].get_F(F);
      printvector (F);
      #endif
    }    
        
    /** @brief Classical (RS) two-pass coarsening. Single-Threaded! (VIENNACL_AMG_COARSE_CLASSIC)
    * @param level    Coarse level identifier
    * @param A      Operator matrix on all levels
    * @param Pointvector   Vector of points on all levels
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_coarse_classic(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      bool add_C;
      amg_point *c_point, *point1, *point2;     
      
      // Use one-pass-coarsening as first pass.
      amg_coarse_classic_onepass(level, A, Pointvector, tag);
    
      // 2nd pass: Add more C points if F-F connection does not have a common C point.
      for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
      {
  point1 = *iter;
  // If point is F point, check for strong connections.
  if (point1->is_fpoint())
  {
    // Check for strong connections from influencing and influenced points.
    amg_point::iterator iter2 = point1->begin_influencing();
    amg_point::iterator iter3 = point1->begin_influenced();
    
    // Iterate over both lists at once. This makes sure that points are no checked twice when influence relation is symmetric (which is often the case).
    // Note: Only works because influencing and influenced lists are sorted by point-index.
    while(iter2 != point1->end_influencing() || iter3 != point1->end_influenced())
    {     
      if (iter2 == point1->end_influencing())
      {
        point2 = *iter3;
        ++iter3;
      }
      else if (iter3 == point1->end_influenced())
      {
        point2 = *iter2;
        ++iter2;
      }
      else
      {      
        if ((*iter2)->get_index() == (*iter3)->get_index())   
        {
    point2 = *iter2;
    ++iter2;
    ++iter3;
        }
        else if ((*iter2)->get_index() < (*iter3)->get_index())
        {
    point2 = *iter2;
    ++iter2;
        }
        else
        {
    point2 = *iter3;
    ++iter3;
        }
      }
      // Only check points with higher index as points with lower index have been checked already.
      if (point2->get_index() < point1->get_index())
        continue;
      
      // If there is a strong connection then it has to either be a C point or a F point with common C point.
      // C point? Then skip as everything is ok.
      if (point2->is_cpoint())
        continue;
      // F point? Then check whether F points point1 and point2 have a common C point.
      if (point2->is_fpoint())
      {
        add_C = true;
        // C point is common for two F points if they are both strongly influenced by that C point.
        // Compare strong influences for point1 and point2.
        for (amg_point::iterator iter3 = point1->begin_influencing(); iter3 != point1 -> end_influencing(); ++iter3)
        {
    c_point = *iter3;
    // Stop search when strong common influence is found via c_point.
    if (c_point->is_cpoint())
    {
      if (point2->is_influencing(c_point))
      {
        add_C = false;
        break;            
      }
    }
        }
        // No common C point found? Then make second F point to C point.
        if (add_C == true)
    Pointvector[level].switch_ftoc(point2);
      }
    }
  }
      }
      
      #ifdef DEBUG
      std::cout << "After 2nd pass:" << std::endl;
      std::cout << "Coarse Points:" << std::endl;
      boost::numeric::ublas::vector<bool> C;
      Pointvector[level].get_C(C);
      printvector (C);
      std::cout << "Fine Points:" << std::endl;
      boost::numeric::ublas::vector<bool> F;
      Pointvector[level].get_F(F);
      printvector (F);
      #endif

      #ifdef DEBUG
#ifdef _OPENMP
      #pragma omp critical
#endif      
      {
      std::cout << "No C and no F point: ";
      for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
  if ((*iter)->is_undecided())
    std::cout << (*iter)->get_index() << " ";
      std::cout << std::endl;
      }
      #endif
    }

    /** @brief Parallel classical RS0 coarsening. Multi-Threaded! (VIENNACL_AMG_COARSE_RS0 || VIENNACL_AMG_COARSE_RS3)
    * @param level    Coarse level identifier
    * @param A      Operator matrix on all level
    * @param Pointvector   Vector of points on all levels
    * @param Slicing    Partitioning of the system matrix and the other data structures to different processors
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2, typename InternalType3>
    void amg_coarse_rs0(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, InternalType3 & Slicing, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      unsigned int total_points;
      
      // Slice matrix into parts such that points are distributed among threads
      Slicing.slice(level, A, Pointvector);     
      
      // Run classical coarsening in parallel
      total_points = 0;
#ifdef _OPENMP
      #pragma omp parallel for shared (total_points,Slicing,level)
#endif      
      for (unsigned int i=0; i<Slicing._threads; ++i)
      {
  amg_coarse_classic(level,Slicing.A_slice[i],Slicing.Pointvector_slice[i],tag);
  
  // Save C points (using Slicing.Offset on the next level as temporary memory)
  // Note: Number of C points for point i is saved in i+1!! (makes it easier later to compute offset)
  Slicing.Offset[i+1][level+1] = Slicing.Pointvector_slice[i][level].get_cpoints();
#ifdef _OPENMP
  #pragma omp critical
#endif  
  total_points += Slicing.Pointvector_slice[i][level].get_cpoints();
      }      
      
      // If no coarser level can be found on any level then resume and coarsening will stop in amg_coarse()
      if (total_points != 0)
      {    
#ifdef _OPENMP
  #pragma omp parallel for shared (Slicing)
#endif  
  for (unsigned int i=0; i<Slicing._threads; ++i)
  {
    // If no higher coarse level can be found on slice i (saved in Slicing.Offset[i+1][level+1]) then pull C point(s) to the next level
    if (Slicing.Offset[i+1][level+1] == 0)
    {
      // All points become C points
      for (unsigned int j=0; j<Slicing.A_slice[i][level].size1(); ++j)
        Slicing.Pointvector_slice[i][level].make_cpoint(Slicing.Pointvector_slice[i][level][j]);
      Slicing.Offset[i+1][level+1] = Slicing.A_slice[i][level].size1();
    }
  }
    
  // Build slicing offset from number of C points (offset = total sum of C points on threads with lower number)
  for (unsigned int i=2; i<=Slicing._threads; ++i)
    Slicing.Offset[i][level+1] += Slicing.Offset[i-1][level+1];
      
  // Join C and F points
  Slicing.join(level, Pointvector);
      }
      
      // Calculate global influence measures for interpolation and/or RS3.
      amg_influence(level, A, Pointvector, tag); 
      
      #if defined(DEBUG)// or defined (DEBUGBENCH)
      for (unsigned int i=0; i<Slicing._threads; ++i)
      {
  unsigned int c_points = Slicing.Pointvector_slice[i][level].get_cpoints();
  unsigned int f_points = Slicing.Pointvector_slice[i][level].get_fpoints();
  std::cout << "Thread " << i << ": ";
  std::cout << "No of C points = " << c_points << ", ";
  std::cout << "No of F points = " << f_points << std::endl;
      }
      #endif
    }
    
    /** @brief RS3 coarsening. Single-Threaded! (VIENNACL_AMG_COARSE_RS3)
    * @param level    Coarse level identifier
    * @param A      Operator matrix on all levels
    * @param Pointvector   Vector of points on all levels
    * @param Slicing    Partitioning of the system matrix and the other data structures to different processors
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2, typename InternalType3>
    void amg_coarse_rs3(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, InternalType3 & Slicing, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      amg_point *c_point, *point1, *point2;
      bool add_C;
      unsigned int i, j;
            
      // Run RS0 first (parallel).
      amg_coarse_rs0(level, A, Pointvector, Slicing, tag);
      
      // Save slicing offset
      boost::numeric::ublas::vector<unsigned int> Offset = boost::numeric::ublas::vector<unsigned int> (Slicing.Offset.size());
      for (i=0; i<Slicing.Offset.size(); ++i)
  Offset[i] = Slicing.Offset[i][level];
      
      // Correct the coarsening with a third pass: Don't allow strong F-F connections without common C point
      for (i=0; i<Slicing._threads; ++i)
      {
  //for (j=Slicing.Offset[i][level]; j<Slicing.Offset[i+1][level]; ++j)
  for (j=Offset[i]; j<Offset[i+1]; ++j)
  {
    point1 = Pointvector[level][j];
    // If point is F point, check for strong connections.
    if (point1->is_fpoint())
    {
      // Check for strong connections from influencing and influenced points.
      amg_point::iterator iter2 = point1->begin_influencing();
      amg_point::iterator iter3 = point1->begin_influenced();
      
      // Iterate over both lists at once. This makes sure that points are no checked twice when influence relation is symmetric (which is often the case).
      // Note: Only works because influencing and influenced lists are sorted by point-index.
      while(iter2 != point1->end_influencing() || iter3 != point1->end_influenced())
      {     
        if (iter2 == point1->end_influencing())
        {
    point2 = *iter3;
    ++iter3;
        }
        else if (iter3 == point1->end_influenced())
        {
    point2 = *iter2;
    ++iter2;
        }
        else
        {      
    if ((*iter2)->get_index() == (*iter3)->get_index())   
    {
      point2 = *iter2;
      ++iter2;
      ++iter3;
    }
    else if ((*iter2)->get_index() < (*iter3)->get_index())
    {
      point2 = *iter2;
      ++iter2;
    }
    else
    {
      point2 = *iter3;
      ++iter3;
    }
        }
              
        // Only check points with higher index as points with lower index have been checked already.
        if (point2->get_index() < point1->get_index())
    continue;
                
        // Only check points that are outside the slicing boundaries (interior F-F connections have already been checked in second pass)
        //if (point2->get_index() >= Slicing.Offset[i][level] || point2->get_index() < Slicing.Offset[i+1][level])
        if (point2->get_index() >= Offset[i] && point2->get_index() < Offset[i+1])
    continue;
        
        // If there is a strong connection then it has to either be a C point or a F point with common C point.
        // C point? Then skip as everything is ok.
        if (point2->is_cpoint())
    continue;
        // F point? Then check whether F points point1 and point2 have a common C point.
        if (point2->is_fpoint())
        {
    add_C = true;
    // C point is common for two F points if they are both strongly influenced by that C point.
    // Compare strong influences for point1 and point2.
    for (amg_point::iterator iter3 = point1->begin_influencing(); iter3 != point1 -> end_influencing(); ++iter3)
    {
      c_point = *iter3;
      // Stop search when strong common influence is found via c_point.
      if (c_point->is_cpoint())
      {
        if (point2->is_influencing(c_point))
        {
          add_C = false;
          break;            
        }
      }
    }
    // No common C point found? Then make second F point to C point.
    if (add_C == true)
    {
      Pointvector[level].switch_ftoc(point2);
      // Add +1 to offsets as one C point has been added.
      for (unsigned int j=i+1; j<=Slicing._threads; ++j)
        Slicing.Offset[j][level+1]++;
    }
        }
      }
    }
  }
      }
      
      #ifdef DEBUG
      std::cout << "After 3rd pass:" << std::endl;
      std::cout << "Coarse Points:" << std::endl;
      boost::numeric::ublas::vector<bool> C;
      Pointvector[level].get_C(C);
      printvector (C);
      std::cout << "Fine Points:" << std::endl;
      boost::numeric::ublas::vector<bool> F;
      Pointvector[level].get_F(F);
      printvector (F);
      #endif

      #ifdef DEBUG
      unsigned int i;
#ifdef _OPENMP
      #pragma omp critical
#endif      
      {
      std::cout << "No C and no F point: ";
      for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
  if ((*iter)->is_undecided())
    std::cout << i << " ";
      std::cout << std::endl;
      }
      #endif
    }
    
    /** @brief AG (aggregation based) coarsening. Single-Threaded! (VIENNACL_AMG_COARSE_SA)
    *
    * @param level    Coarse level identifier
    * @param A      Operator matrix on all levels
    * @param Pointvector   Vector of points on all levels
    * @param tag    AMG preconditioner tag
    */
    template <typename InternalType1, typename InternalType2>
    void amg_coarse_ag(unsigned int level, InternalType1 & A, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;
      typedef typename SparseMatrixType::value_type ScalarType;
      
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      unsigned int x,y;
      ScalarType diag;
      amg_point *pointx, *pointy;
    
      // Cannot determine aggregates if size == 1 as then a new aggregate would always consist of this point (infinite loop)
      if (A[level].size1() == 1) return;
      
      // SA algorithm (Vanek et al. p.6)     
      // Build neighborhoods
#ifdef _OPENMP
      #pragma omp parallel for private (x,y,diag) shared (A)
#endif      
      for (x=0; x<A[level].size1(); ++x)
      {
  InternalRowIterator row_iter = A[level].begin1();
  row_iter += x;
  diag = A[level](x,x);
  for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
  {
    y = col_iter.index2();
    if (y == x || (std::abs(*col_iter) >= tag.get_threshold()*pow(0.5,level-1) * sqrt(std::abs(diag*A[level](y,y)))))
    {
      // Neighborhood x includes point y
      Pointvector[level][x]->add_influencing_point(Pointvector[level][y]);
    }
  }
      }
      
      #ifdef DEBUG
      std::cout << "Neighborhoods:" << std::endl;
      boost::numeric::ublas::matrix<bool> mat;
      Pointvector[level].get_influence_matrix(mat);
      printmatrix (mat);
      #endif

      // Build aggregates from neighborhoods  
      for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
      {
  pointx = (*iter);
  
  if (pointx->is_undecided())
  {
    // Make center of aggregate to C point and include it to aggregate x.
    Pointvector[level].make_cpoint(pointx);
    pointx->set_aggregate (pointx->get_index());
    for (amg_point::iterator iter2 = pointx->begin_influencing(); iter2 != pointx->end_influencing(); ++iter2)
    {
     pointy = (*iter2);
      
      if (pointy->is_undecided())
      {
        // Make neighbor y to F point and include it to aggregate x.
        Pointvector[level].make_fpoint(pointy);
        pointy->set_aggregate (pointx->get_index());
      }
    }
  }
      }
      
      #ifdef DEBUG
      std::cout << "After aggregation:" << std::endl;
      std::cout << "Coarse Points:" << std::endl;
      boost::numeric::ublas::vector<bool> C;
      Pointvector[level].get_C(C);
      printvector (C);
      std::cout << "Fine Points:" << std::endl;
      boost::numeric::ublas::vector<bool> F;
      Pointvector[level].get_F(F);
      printvector (F);
      std::cout << "Aggregates:" << std::endl;
      printvector (Aggregates[level]);          
      #endif
    }
      } //namespace amg
    }
  }
}

#endif
