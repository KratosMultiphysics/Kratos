#ifndef VIENNACL_LINALG_AMG_HPP_
#define VIENNACL_LINALG_AMG_HPP_

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

/** @file viennacl/linalg/amg.hpp
    @brief Main include file for algebraic multigrid (AMG) preconditioners.  Experimental in 1.2.x.
    
    Implementation contributed by Markus Wagner
*/

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/triangular.hpp> 
#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/tools/tools.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/direct_solve.hpp"

#include "viennacl/linalg/detail/amg/amg_base.hpp"
#include "viennacl/linalg/detail/amg/amg_coarse.hpp"
#include "viennacl/linalg/detail/amg/amg_interpol.hpp"

#include <map>

#ifdef _OPENMP
 #include <omp.h>
#endif

#include "viennacl/linalg/detail/amg/amg_debug.hpp"

#define VIENNACL_AMG_COARSE_LIMIT 50
#define VIENNACL_AMG_MAX_LEVELS 100

namespace viennacl
{
  namespace linalg
  {    
    typedef detail::amg::amg_tag          amg_tag;
    
    
    
    /** @brief Setup AMG preconditioner
    *
    * @param A      Operator matrices on all levels
    * @param P      Prolongation/Interpolation operators on all levels
    * @param Pointvector  Vector of points on all levels
    * @param tag    AMG preconditioner tag 
    */
    template <typename InternalType1, typename InternalType2>
    void amg_setup(InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename InternalType1::value_type SparseMatrixType;
      typedef typename InternalType2::value_type PointVectorType;     
      typedef typename SparseMatrixType::value_type ScalarType;   
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      unsigned int i, iterations, c_points, f_points;
      detail::amg::amg_slicing<InternalType1,InternalType2> Slicing;
      
      // Set number of iterations. If automatic coarse grid construction is chosen (0), then set a maximum size and stop during the process.
      iterations = tag.get_coarselevels();
      if (iterations == 0)
        iterations = VIENNACL_AMG_MAX_LEVELS;
             
      // For parallel coarsenings build data structures (number of threads set automatically).
      if (tag.get_coarse() == VIENNACL_AMG_COARSE_RS0 || tag.get_coarse() == VIENNACL_AMG_COARSE_RS3)
        Slicing.init(iterations);
      
      for (i=0; i<iterations; ++i)
      {  
        // Initialize Pointvector on level i and construct points.
        Pointvector[i] = PointVectorType(A[i].size1());
        Pointvector[i].init_points();
        
        // Construct C and F points on coarse level (i is fine level, i+1 coarse level).
        detail::amg::amg_coarse (i, A, Pointvector, Slicing, tag);

        // Calculate number of C and F points on level i.
        c_points = Pointvector[i].get_cpoints();
        f_points = Pointvector[i].get_fpoints();

        #if defined (DEBUG) //or defined(DEBUGBENCH) 
        std::cout << "Level " << i << ": ";
        std::cout << "No of C points = " << c_points << ", ";
        std::cout << "No of F points = " << f_points << std::endl;
        #endif
        
        // Stop routine when the maximal coarse level is found (no C or F point). Coarsest level is level i.
        if (c_points == 0 || f_points == 0)
          break;
          
        // Construct interpolation matrix for level i.
        detail::amg::amg_interpol (i, A, P, Pointvector, tag);
        
        // Compute coarse grid operator (A[i+1] = R * A[i] * P) with R = trans(P).
        detail::amg::amg_galerkin_prod(A[i], P[i], A[i+1]);
        
        // Test triple matrix product. Very slow for large matrix sizes (ublas).
        // test_triplematprod(A[i],P[i],A[i+1]);
        
        Pointvector[i].delete_points();
        
        #ifdef DEBUG
        std::cout << "Coarse Grid Operator Matrix:" << std::endl;
        printmatrix (A[i+1]);
        #endif  
        
        // If Limit of coarse points is reached then stop. Coarsest level is level i+1.
        if (tag.get_coarselevels() == 0 && c_points <= VIENNACL_AMG_COARSE_LIMIT)
        {
          tag.set_coarselevels(i+1);
          return;
        }
      }
      tag.set_coarselevels(i);  
    }
    
    /** @brief Initialize AMG preconditioner
    *
    * @param mat    System matrix
    * @param A      Operator matrices on all levels
    * @param P      Prolongation/Interpolation operators on all levels
    * @param Pointvector  Vector of points on all levels
    * @param tag    AMG preconditioner tag 
    */
    template <typename MatrixType, typename InternalType1, typename InternalType2>
    void amg_init(MatrixType const & mat, InternalType1 & A, InternalType1 & P, InternalType2 & Pointvector, amg_tag & tag)
    {
      typedef typename MatrixType::value_type ScalarType;
      typedef typename InternalType1::value_type SparseMatrixType;
      
      if (tag.get_coarselevels() > 0)
      {
        A.resize(tag.get_coarselevels()+1);
        P.resize(tag.get_coarselevels());
        Pointvector.resize(tag.get_coarselevels());
      }
      else
      {
        A.resize(VIENNACL_AMG_MAX_LEVELS+1);
        P.resize(VIENNACL_AMG_MAX_LEVELS);
        Pointvector.resize(VIENNACL_AMG_MAX_LEVELS);
      }
      
      // Insert operator matrix as operator for finest level.
      SparseMatrixType A0 (mat);
      A.insert_element (0, A0);  
    }
    
    /** @brief Save operators after setup phase for CPU computation.
    *
    * @param A      Operator matrices on all levels on the CPU
    * @param P      Prolongation/Interpolation operators on all levels on the CPU
    * @param R      Restriction operators on all levels on the CPU
    * @param A_setup    Operators matrices on all levels from setup phase
    * @param P_setup    Prolongation/Interpolation operators on all levels from setup phase
    * @param tag    AMG preconditioner tag 
    */
    template <typename InternalType1, typename InternalType2>
    void amg_transform_cpu (InternalType1 & A, InternalType1 & P, InternalType1 & R, InternalType2 & A_setup, InternalType2 & P_setup, amg_tag & tag)
    {  
      typedef typename InternalType1::value_type MatrixType;
      
      // Resize internal data structures to actual size.
      A.resize(tag.get_coarselevels()+1);
      P.resize(tag.get_coarselevels());
      R.resize(tag.get_coarselevels());
      
      // Transform into matrix type.
      for (unsigned int i=0; i<tag.get_coarselevels()+1; ++i)
      {
        A[i].resize(A_setup[i].size1(),A_setup[i].size2(),false);
        A[i] = A_setup[i];
      }
      for (unsigned int i=0; i<tag.get_coarselevels(); ++i)
      {
        P[i].resize(P_setup[i].size1(),P_setup[i].size2(),false);
        P[i] = P_setup[i];
      }
      for (unsigned int i=0; i<tag.get_coarselevels(); ++i)
      {
        R[i].resize(P_setup[i].size2(),P_setup[i].size1(),false);
        P_setup[i].set_trans(true);
        R[i] = P_setup[i];
        P_setup[i].set_trans(false);
      }
    }
    
    /** @brief Save operators after setup phase for GPU computation.
    *
    * @param A      Operator matrices on all levels on the GPU
    * @param P      Prolongation/Interpolation operators on all levels on the GPU
    * @param R      Restriction operators on all levels on the GPU
    * @param A_setup    Operators matrices on all levels from setup phase
    * @param P_setup    Prolongation/Interpolation operators on all levels from setup phase
    * @param tag    AMG preconditioner tag 
    */
    template <typename InternalType1, typename InternalType2>
    void amg_transform_gpu (InternalType1 & A, InternalType1 & P, InternalType1 & R, InternalType2 & A_setup, InternalType2 & P_setup, amg_tag & tag)
    {
      typedef typename InternalType1::value_type MatrixType;
      typedef typename InternalType2::value_type::value_type ScalarType;
      
      // Resize internal data structures to actual size.
      A.resize(tag.get_coarselevels()+1);
      P.resize(tag.get_coarselevels());
      R.resize(tag.get_coarselevels());
      
      // Copy to GPU using the internal sparse matrix structure: std::vector<std::map>.
      for (unsigned int i=0; i<tag.get_coarselevels()+1; ++i)
      {
        A[i].resize(A_setup[i].size1(),A_setup[i].size2(),false);
        viennacl::copy(*(A_setup[i].get_internal_pointer()),A[i]);
      }
      for (unsigned int i=0; i<tag.get_coarselevels(); ++i)
      {
        P[i].resize(P_setup[i].size1(),P_setup[i].size2(),false);
        viennacl::copy(*(P_setup[i].get_internal_pointer()),P[i]);
        //viennacl::copy((boost::numeric::ublas::compressed_matrix<ScalarType>)P_setup[i],P[i]);
      }
      for (unsigned int i=0; i<tag.get_coarselevels(); ++i)
      {
        R[i].resize(P_setup[i].size2(),P_setup[i].size1(),false);
        P_setup[i].set_trans(true);
        viennacl::copy(*(P_setup[i].get_internal_pointer()),R[i]);
        P_setup[i].set_trans(false);
      }
    }
    
    /** @brief Setup data structures for precondition phase.
    *
    * @param result    Result vector on all levels
    * @param rhs    RHS vector on all levels
    * @param residual    Residual vector on all levels
    * @param A      Operators matrices on all levels from setup phase
    * @param tag    AMG preconditioner tag 
    */
    template <typename InternalVectorType, typename SparseMatrixType>
    void amg_setup_apply (InternalVectorType & result, InternalVectorType & rhs, InternalVectorType & residual, SparseMatrixType const & A, amg_tag const & tag)
    {        
      typedef typename InternalVectorType::value_type VectorType;
      
      result.resize(tag.get_coarselevels()+1);
      rhs.resize(tag.get_coarselevels()+1);
      residual.resize(tag.get_coarselevels());
            
      for (unsigned int level=0; level < tag.get_coarselevels()+1; ++level)
      {
        result[level] = VectorType(A[level].size1());
        result[level].clear();
        rhs[level] = VectorType(A[level].size1());
        rhs[level].clear();
      }
      for (unsigned int level=0; level < tag.get_coarselevels(); ++level)
      {
        residual[level] = VectorType(A[level].size1());
        residual[level].clear();
      }
    }
    /** @brief Pre-compute LU factorization for direct solve (ublas library).
     *  @brief Speeds up precondition phase as this is computed only once overall instead of once per iteration.
    *
    * @param op      Operator matrix for direct solve
    * @param Permutation  Permutation matrix which saves the factorization result
    * @param A      Operator matrix on coarsest level
    */
    template <typename ScalarType, typename SparseMatrixType>
    void amg_lu(boost::numeric::ublas::compressed_matrix<ScalarType> & op, boost::numeric::ublas::permutation_matrix<ScalarType> & Permutation, SparseMatrixType const & A)
    {
      typedef typename SparseMatrixType::const_iterator1 ConstRowIterator;
      typedef typename SparseMatrixType::const_iterator2 ConstColIterator;
    
      // Copy to operator matrix. Needed 
      op.resize(A.size1(),A.size2(),false);
      for (ConstRowIterator row_iter = A.begin1(); row_iter != A.end1(); ++row_iter)
        for (ConstColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
          op (col_iter.index1(), col_iter.index2()) = *col_iter;
      
      // Permutation matrix has to be reinitialized with actual size. Do not clear() or resize()!
      Permutation = boost::numeric::ublas::permutation_matrix<ScalarType> (op.size1());
      boost::numeric::ublas::lu_factorize(op,Permutation);
    }
    
    /** @brief AMG preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class amg_precond
    {
      typedef typename MatrixType::value_type ScalarType;
      typedef boost::numeric::ublas::vector<ScalarType> VectorType;
      typedef detail::amg::amg_sparsematrix<ScalarType> SparseMatrixType;
      typedef detail::amg::amg_pointvector PointVectorType;
      
      typedef typename SparseMatrixType::const_iterator1 InternalConstRowIterator;
      typedef typename SparseMatrixType::const_iterator2 InternalConstColIterator;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      boost::numeric::ublas::vector <SparseMatrixType> A_setup;
      boost::numeric::ublas::vector <SparseMatrixType> P_setup;
      boost::numeric::ublas::vector <MatrixType> A;
      boost::numeric::ublas::vector <MatrixType> P;
      boost::numeric::ublas::vector <MatrixType> R;
      boost::numeric::ublas::vector <PointVectorType> Pointvector;
       
      mutable boost::numeric::ublas::compressed_matrix<ScalarType> op;
      mutable boost::numeric::ublas::permutation_matrix<ScalarType> Permutation;  
      
      mutable boost::numeric::ublas::vector <VectorType> result;
      mutable boost::numeric::ublas::vector <VectorType> rhs;
      mutable boost::numeric::ublas::vector <VectorType> residual;
      
      mutable bool done_init_apply;
          
      amg_tag _tag;
    public:
    
      amg_precond(): Permutation(0) {}
      /** @brief The constructor. Saves system matrix, tag and builds data structures for setup.
      *
      * @param mat  System matrix
      * @param tag  The AMG tag
      */
      amg_precond(MatrixType const & mat, amg_tag const & tag): Permutation(0)
      {
        _tag = tag;
        // Initialize data structures.
        amg_init (mat,A_setup,P_setup,Pointvector,_tag);
        
        done_init_apply = false;
      }
      
      /** @brief Start setup phase for this class and copy data structures.
      */
      void setup()
      {
        // Start setup phase.
        amg_setup(A_setup,P_setup,Pointvector,_tag);
        // Transform to CPU-Matrixtype for precondition phase.
        amg_transform_cpu(A,P,R,A_setup,P_setup,_tag);
        
        done_init_apply = false;
      }
      
      /** @brief Prepare data structures for preconditioning:
       *  Build data structures for precondition phase.
       *  Do LU factorization on coarsest level.
      */
      void init_apply() const
      {
        // Setup precondition phase (Data structures).
        amg_setup_apply(result,rhs,residual,A_setup,_tag);
        // Do LU factorization for direct solve.
        amg_lu(op,Permutation,A_setup[_tag.get_coarselevels()]);
        
        done_init_apply = true;
      }
      
      /** @brief Returns complexity measures.
      *
      * @param avgstencil  Average stencil sizes on all levels
      * @return     Operator complexity of AMG method
      */
      template <typename VectorType>
      ScalarType calc_complexity(VectorType & avgstencil)
      {
        avgstencil = VectorType (_tag.get_coarselevels()+1);
        unsigned int nonzero=0, systemmat_nonzero=0, level_coefficients=0;
          
        for (unsigned int level=0; level < _tag.get_coarselevels()+1; ++level)
        {
          level_coefficients = 0;
          for (InternalRowIterator row_iter = A_setup[level].begin1(); row_iter != A_setup[level].end1(); ++row_iter)
          { 
            for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {
              if (level == 0)
                systemmat_nonzero++;
              nonzero++;
              level_coefficients++;
            }
          }
          avgstencil[level] = level_coefficients/(double)A_setup[level].size1();
        }
        return nonzero/static_cast<double>(systemmat_nonzero);
      }

      /** @brief Precondition Operation
      *
      * @param vec The vector to which preconditioning is applied to (ublas version)
      */
      template <typename VectorType>
      void apply(VectorType & vec) const
      {   
        // Build data structures and do lu factorization before first iteration step.
        if (!done_init_apply)
          init_apply();
        
        int level;
        
        // Precondition operation (Yang, p.3)
        rhs[0] = vec;
        for (level=0; level <(signed)_tag.get_coarselevels(); level++)
        {    
          result[level].clear();
          
          // Apply Smoother _presmooth times.
          smooth_jacobi (level, _tag.get_presmooth(), result[level], rhs[level]);    
          
          #ifdef DEBUG
          std::cout << "After presmooth:" << std::endl;
          printvector(result[level]);
          #endif

          // Compute residual.
          residual[level] = rhs[level] - boost::numeric::ublas::prod (A[level],result[level]);
          
          #ifdef DEBUG          
          std::cout << "Residual:" << std::endl;
          printvector(residual[level]);
          #endif
          
          // Restrict to coarse level. Restricted residual is RHS of coarse level.
          rhs[level+1] = boost::numeric::ublas::prod (R[level],residual[level]);
          
          #ifdef DEBUG
          std::cout << "Restricted Residual: " << std::endl;
          printvector(rhs[level+1]);
          #endif
        }
          
        // On highest level use direct solve to solve equation.
        result[level] = rhs[level];
        boost::numeric::ublas::lu_substitute(op,Permutation,result[level]);

        #ifdef DEBUG
        std::cout << "After direct solve: " << std::endl;
        printvector (result[level]);
        #endif
          
        for (level=_tag.get_coarselevels()-1; level >= 0; level--)
        {       
          #ifdef DEBUG
          std::cout << "Coarse Error: " << std::endl;
          printvector(result[level+1]);
          #endif
          
          // Interpolate error to fine level. Correct solution by adding error.
          result[level] += boost::numeric::ublas::prod (P[level], result[level+1]);
              
          #ifdef DEBUG
          std::cout << "Corrected Result: " << std::endl;
          printvector (result[level]);
          #endif
          
          // Apply Smoother _postsmooth times.
          smooth_jacobi (level, _tag.get_postsmooth(), result[level], rhs[level]);
          
          #ifdef DEBUG
          std::cout << "After postsmooth: " << std::endl;
          printvector (result[level]);
          #endif
        }    
        vec = result[0];
      }
      
      /** @brief (Weighted) Jacobi Smoother (CPU version)
      * @param level    Coarse level to which smoother is applied to
      * @param iterations  Number of smoother iterations    
      * @param x     The vector smoothing is applied to
      * @param rhs    The right hand side of the equation for the smoother
      */
      template <typename VectorType>
      void smooth_jacobi(int level, int const iterations, VectorType & x, VectorType const & rhs) const
      {
        VectorType old_result (x.size());
        unsigned int index;
        ScalarType sum = 0, diag = 1;
        
        for (int i=0; i<iterations; ++i)
        {
          old_result = x;
          x.clear();
#ifdef _OPENMP
          #pragma omp parallel for private (sum,diag) shared (rhs,x)
#endif          
          for (index=0; index<A_setup[level].size1(); ++index)  
          {
            InternalConstRowIterator row_iter = A_setup[level].begin1();
            row_iter += index; 
            sum = 0;
            diag = 1;
            for (InternalConstColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {
              if (col_iter.index1() == col_iter.index2())
                diag = *col_iter;
              else
                sum += *col_iter * old_result[col_iter.index2()];
            }
            x[index]= _tag.get_jacobiweight() * (rhs[index] - sum) / diag + (1-_tag.get_jacobiweight()) * old_result[index];
          }
        }
      }
      
      amg_tag & tag() { return _tag; }
    };
    
    /** @brief AMG preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT>
    class amg_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT> >
    {
      typedef viennacl::compressed_matrix<ScalarType, MAT_ALIGNMENT> MatrixType;
      typedef viennacl::vector<ScalarType> VectorType;
      typedef detail::amg::amg_sparsematrix<ScalarType> SparseMatrixType;
      typedef detail::amg::amg_pointvector PointVectorType;
      
      typedef typename SparseMatrixType::const_iterator1 InternalConstRowIterator;
      typedef typename SparseMatrixType::const_iterator2 InternalConstColIterator;
      typedef typename SparseMatrixType::iterator1 InternalRowIterator;
      typedef typename SparseMatrixType::iterator2 InternalColIterator;
      
      boost::numeric::ublas::vector <SparseMatrixType> A_setup;
      boost::numeric::ublas::vector <SparseMatrixType> P_setup;
      boost::numeric::ublas::vector <MatrixType> A;
      boost::numeric::ublas::vector <MatrixType> P;
      boost::numeric::ublas::vector <MatrixType> R;
      boost::numeric::ublas::vector <PointVectorType> Pointvector;
      
      mutable boost::numeric::ublas::compressed_matrix<ScalarType> op;
      mutable boost::numeric::ublas::permutation_matrix<ScalarType> Permutation;  
      
      mutable boost::numeric::ublas::vector <VectorType> result;
      mutable boost::numeric::ublas::vector <VectorType> rhs;
      mutable boost::numeric::ublas::vector <VectorType> residual;
          
      mutable bool done_init_apply;
    
      amg_tag _tag;
      
    public:
      
      amg_precond(): Permutation(0) {}
      
      /** @brief The constructor. Builds data structures.
      *
      * @param mat  System matrix
      * @param tag  The AMG tag
      */
      amg_precond(compressed_matrix<ScalarType, MAT_ALIGNMENT> const & mat, amg_tag const & tag): Permutation(0)
      {
        _tag = tag;
        
        // Copy to CPU. Internal structure of sparse matrix is used for copy operation.
        std::vector<std::map<unsigned int, ScalarType> > mat2 = std::vector<std::map<unsigned int, ScalarType> >(mat.size1());
        viennacl::copy(mat, mat2);
        
        // Initialize data structures.
        amg_init (mat2,A_setup,P_setup,Pointvector,_tag);
          
        done_init_apply = false;
      }
      
      /** @brief Start setup phase for this class and copy data structures.
      */
      void setup()
      {
        // Start setup phase.
        amg_setup(A_setup,P_setup,Pointvector,_tag);  
        // Transform to GPU-Matrixtype for precondition phase.
        amg_transform_gpu(A,P,R,A_setup,P_setup,_tag);  
        
        done_init_apply = false;
      }
      
      /** @brief Prepare data structures for preconditioning:
       *  Build data structures for precondition phase.
       *  Do LU factorization on coarsest level.
      */
      void init_apply() const
      {
        // Setup precondition phase (Data structures).
        amg_setup_apply(result,rhs,residual,A_setup,_tag);
        // Do LU factorization for direct solve.
        amg_lu(op,Permutation,A_setup[_tag.get_coarselevels()]);
        
        done_init_apply = true;
      }

      /** @brief Returns complexity measures
      *
      * @param avgstencil  Average stencil sizes on all levels
      * @return     Operator complexity of AMG method
      */
      template <typename VectorType>
      ScalarType calc_complexity(VectorType & avgstencil)
      {
        avgstencil = VectorType (_tag.get_coarselevels()+1);
        unsigned int nonzero=0, systemmat_nonzero=0, level_coefficients=0;
        
        for (unsigned int level=0; level < _tag.get_coarselevels()+1; ++level)
        {
          level_coefficients = 0;
          for (InternalRowIterator row_iter = A_setup[level].begin1(); row_iter != A_setup[level].end1(); ++row_iter)
          { 
            for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {
              if (level == 0)
                systemmat_nonzero++;
              nonzero++;
              level_coefficients++;
            }
          }
          avgstencil[level] = level_coefficients/(double)A[level].size1();
        }
        return nonzero/static_cast<double>(systemmat_nonzero);
      }

      /** @brief Precondition Operation
      *
      * @param vec The vector to which preconditioning is applied to 
      */
      template <typename VectorType>
      void apply(VectorType & vec) const
      {
        if (!done_init_apply)
          init_apply();  
        
        int level;
        
        // Precondition operation (Yang, p.3).
        rhs[0] = vec;
        for (level=0; level <(signed)_tag.get_coarselevels(); level++)
        {    
          result[level].clear();
          
          // Apply Smoother _presmooth times.
          smooth_jacobi (level, _tag.get_presmooth(), result[level], rhs[level]);

          #ifdef DEBUG
          std::cout << "After presmooth: " << std::endl;
          printvector(result[level]);
          #endif
          
          // Compute residual.
          residual[level] = rhs[level] - viennacl::linalg::prod (A[level],result[level]);
          
          #ifdef DEBUG             
          std::cout << "Residual: " << std::endl;
          printvector(residual[level]);
          #endif
          
          // Restrict to coarse level. Result is RHS of coarse level equation.
          //residual_coarse[level] = viennacl::linalg::prod(R[level],residual[level]);
          rhs[level+1] = viennacl::linalg::prod(R[level],residual[level]);
          
          #ifdef DEBUG
          std::cout << "Restricted Residual: " << std::endl;
          printvector(rhs[level+1]);
          #endif
        }
          
        // On highest level use direct solve to solve equation (on the CPU)
        //TODO: Use GPU direct solve!
        result[level] = rhs[level];
        boost::numeric::ublas::vector <ScalarType> result_cpu (result[level].size());

        copy (result[level],result_cpu);
        boost::numeric::ublas::lu_substitute(op,Permutation,result_cpu);
        copy (result_cpu, result[level]);
        
        #ifdef DEBUG
        std::cout << "After direct solve: " << std::endl;
        printvector (result[level]);
        #endif
          
        for (level=_tag.get_coarselevels()-1; level >= 0; level--)
        {   
          #ifdef DEBUG
          std::cout << "Coarse Error: " << std::endl;
          printvector(result[level+1]);
          #endif
          
          // Interpolate error to fine level and correct solution.
          result[level] += viennacl::linalg::prod(P[level],result[level+1]);
              
          #ifdef DEBUG
          std::cout << "Corrected Result: " << std::endl;
          printvector (result[level]);
          #endif
          
          // Apply Smoother _postsmooth times.
          smooth_jacobi (level, _tag.get_postsmooth(), result[level], rhs[level]);
          
          #ifdef DEBUG
          std::cout << "After postsmooth: " << std::endl;
          printvector (result[level]);
          #endif
        }    
        vec = result[0];
      }
      
      /** @brief Jacobi Smoother (GPU version)
      * @param level       Coarse level to which smoother is applied to
      * @param iterations  Number of smoother iterations    
      * @param x           The vector smoothing is applied to
      * @param rhs         The right hand side of the equation for the smoother
      */
      template <typename VectorType>
      void smooth_jacobi(int level, unsigned int iterations, VectorType & x, VectorType const & rhs) const
      {     
        VectorType old_result (x.size());
  
  //viennacl::ocl::program & p = viennacl::ocl::current_context().add_program
  //          (viennacl::tools::make_double_kernel(jacobi_kernel,viennacl::ocl::current_device().info()), "jacobi_kernel");
  //viennacl::ocl::kernel & k = p.add_kernel("jacobi");
  
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::compressed_matrix<ScalarType, MAT_ALIGNMENT>::program_name(),
                    "jacobi");
  
        for (unsigned int i=0; i<iterations; ++i)
        {
          old_result = x;    
          x.clear();
          viennacl::ocl::enqueue(k(A[level].handle1(), A[level].handle2(), A[level].handle(),
                                  static_cast<ScalarType>(_tag.get_jacobiweight()), 
                                  old_result,
                                  x,
                                  rhs,
                                  static_cast<cl_uint>(rhs.size()))); 
          
        }
      }
      
      amg_tag & tag() { return _tag; }
    };

  }
}



#endif

