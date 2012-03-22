#ifndef VIENNACL_JACOBI_PRECOND_HPP_
#define VIENNACL_JACOBI_PRECOND_HPP_

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

/** @file jacobi_precond.hpp
    @brief Implementation of a simple Jacobi preconditioner
*/

#include <vector>
#include <cmath>
#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/tools/tools.hpp"

#include <map>

namespace viennacl
{
  namespace linalg
  {
    
    /** @brief A tag for a jacobi preconditioner
    */
    class jacobi_tag {};
    

    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class jacobi_precond
    {
      typedef typename MatrixType::value_type      ScalarType;
      
      public:
        jacobi_precond(MatrixType const & mat, jacobi_tag const & tag) : system_matrix(mat)
        {
          assert(mat.size1() == mat.size2());
          diag_A_inv.resize(mat.size1());  //resize without preserving values
          
          for (typename MatrixType::const_iterator1 row_it = system_matrix.begin1();
                row_it != system_matrix.end1();
                ++row_it)
          {
            bool diag_found = false;
            for (typename MatrixType::const_iterator2 col_it = row_it.begin();
                  col_it != row_it.end();
                  ++col_it)
            {
              if (col_it.index1() == col_it.index2())
              {
                diag_A_inv[col_it.index1()] = static_cast<ScalarType>(1.0) / *col_it;
                diag_found = true;
              }
            }
            if (!diag_found)
              throw "ViennaCL: Zero in diagonal encountered while setting up Jacobi preconditioner!";
          }
        }
        
        
        /** @brief Apply to res = b - Ax, i.e. jacobi applied vec (right hand side),  */
        template <typename VectorType>
        void apply(VectorType & vec) const
        {
          assert(vec.size() == diag_A_inv.size());
          for (size_t i=0; i<vec.size(); ++i)
          {
            vec[i] *= diag_A_inv[i];
          }
        }
        
      private:
        MatrixType const & system_matrix;
        std::vector<ScalarType> diag_A_inv;
    };

    
    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT>
    class jacobi_precond< compressed_matrix<ScalarType, MAT_ALIGNMENT> >
    {
      typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>   MatrixType;
      
      public:
        jacobi_precond(MatrixType const & mat, jacobi_tag const & tag) : system_matrix(mat), diag_A_inv(mat.size1())
        {
          assert(system_matrix.size1() == system_matrix.size2());

          init_gpu();
        }
          
        /*void init_cpu()
        {
          
          std::vector< std::map<unsigned int, ScalarType> > cpu_check;
          std::vector<ScalarType> diag_A_inv_cpu(system_matrix.size1());
          
          copy(system_matrix, cpu_check);
          viennacl::tools::const_sparse_matrix_adapter<ScalarType> cpu_check_adapter(cpu_check);
          
          for (typename viennacl::tools::const_sparse_matrix_adapter<ScalarType>::const_iterator1 row_it = cpu_check_adapter.begin1();
                row_it != cpu_check_adapter.end1();
                ++row_it)
          {
            bool diag_found = false;
            for (typename viennacl::tools::const_sparse_matrix_adapter<ScalarType>::const_iterator2 col_it = row_it.begin();
                  col_it != row_it.end();
                  ++col_it)
            {
              if (col_it.index1() == col_it.index2())
              {
                diag_found = true;
                diag_A_inv_cpu[col_it.index1()] = static_cast<ScalarType>(1.0) / *col_it;
              }
            }
            if (!diag_found)
              throw "ViennaCL: Zero in diagonal encountered while setting up Jacobi preconditioner!";
          }
          
          diag_A_inv.resize(system_matrix.size1(), false);
          viennacl::fast_copy(diag_A_inv_cpu, diag_A_inv);
        }*/
        
        void init_gpu()
        {
          viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(
                                              viennacl::linalg::kernels::compressed_matrix<ScalarType, MAT_ALIGNMENT>::program_name(),
                                              "jacobi_precond");

          viennacl::ocl::enqueue( k(system_matrix.handle1(), system_matrix.handle2(), system_matrix.handle(), 
                                    diag_A_inv, static_cast<cl_uint>(diag_A_inv.size())) );        
        }
        
        
        template <unsigned int ALIGNMENT>
        void apply(viennacl::vector<ScalarType, ALIGNMENT> & vec) const
        {
          assert(viennacl::traits::size1(system_matrix) == viennacl::traits::size(vec));
          
          //run kernel:
          viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<ScalarType, ALIGNMENT>::program_name(),
                                                                "diag_precond");

          viennacl::ocl::enqueue(
             k(viennacl::traits::handle(diag_A_inv), cl_uint(viennacl::traits::start(diag_A_inv)), cl_uint(viennacl::traits::size(diag_A_inv)),
               viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)) )
                                );        
        }
        
      private:
        MatrixType const & system_matrix;
        viennacl::vector<ScalarType> diag_A_inv;
    };

  }
}




#endif



