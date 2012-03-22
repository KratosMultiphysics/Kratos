#ifndef VIENNACL_LINALG_ROW_SCALING_HPP_
#define VIENNACL_LINALG_ROW_SCALING_HPP_

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

/** @file row_scaling.hpp
    @brief A row normalization preconditioner is implemented here
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
    
    /** @brief A tag for a row preconditioner
    */
    class row_scaling_tag
    {
      public:
        /** @brief Constructor
        *
        * @param p   Integer selecting the desired row norm.
        */
        row_scaling_tag(unsigned int p = 2) : norm_(p) {}
        
        /** @brief Returns the index p of the l^p-norm (1... sum(abs(x)), 2... sqrt(sum(x_i^2))). Currently only p=1 and p=2 supported*/
        unsigned int norm() const { return norm_; }
        
      private:
        unsigned int norm_;
    };
    

    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines
    */
    template <typename MatrixType>
    class row_scaling
    {
      typedef typename MatrixType::value_type      ScalarType;
      
      public:
        /** @brief Constructor for the preconditioner
        *
        * @param mat   The system matrix
        * @param tag   A row scaling tag holding the desired norm.
        */
        row_scaling(MatrixType const & mat, row_scaling_tag const & tag) : system_matrix(mat), tag_(tag)
        {
          assert(mat.size1() == mat.size2());
          diag_M_inv.resize(mat.size1());  //resize without preserving values
          
          for (typename MatrixType::const_iterator1 row_it = system_matrix.begin1();
                row_it != system_matrix.end1();
                ++row_it)
          {
            diag_M_inv[row_it.index1()];
            for (typename MatrixType::const_iterator2 col_it = row_it.begin();
                  col_it != row_it.end();
                  ++col_it)
            {
              if (tag_.norm() == 1)
                diag_M_inv[col_it.index1()] += std::fabs(*col_it);
              else
                diag_M_inv[col_it.index1()] += (*col_it) * (*col_it);
            }
            if (diag_M_inv[row_it.index1()] == 0)
              throw "ViennaCL: Zero row encountered while setting up row scaling preconditioner!";
            
            if (tag_.norm() == 1)
              diag_M_inv[row_it.index1()] = static_cast<ScalarType>(1.0) / diag_M_inv[row_it.index1()];
            else
              diag_M_inv[row_it.index1()] = static_cast<ScalarType>(1.0) / std::sqrt(diag_M_inv[row_it.index1()]);
          }
        }
        
        
        /** @brief Apply to res = b - Ax, i.e. row applied vec (right hand side),  */
        template <typename VectorType>
        void apply(VectorType & vec) const
        {
          assert(vec.size() == diag_M_inv.size());
          for (size_t i=0; i<vec.size(); ++i)
          {
            vec[i] *= diag_M_inv[i];
          }
        }
        
      private:
        MatrixType const & system_matrix;
        row_scaling_tag const & tag_;
        std::vector<ScalarType> diag_M_inv;
    };

    
    /** @brief Jacobi preconditioner class, can be supplied to solve()-routines.
    *
    *  Specialization for compressed_matrix
    */
    template <typename ScalarType, unsigned int MAT_ALIGNMENT>
    class row_scaling< compressed_matrix<ScalarType, MAT_ALIGNMENT> >
    {
      typedef compressed_matrix<ScalarType, MAT_ALIGNMENT>   MatrixType;
      
      public:
        /** @brief Constructor for the preconditioner
        *
        * @param mat   The system matrix
        * @param tag   A row scaling tag holding the desired norm.
        */
        row_scaling(MatrixType const & mat, row_scaling_tag const & tag) : system_matrix(mat), tag_(tag), diag_M_inv(mat.size1())
        {
          assert(system_matrix.size1() == system_matrix.size2());
          
          init_gpu();
        }
        
        /*
        void init_cpu()
        {
          std::vector< std::map<unsigned int, ScalarType> > cpu_check;
          std::vector<ScalarType> diag_M_inv_cpu(system_matrix.size1());
          
          copy(system_matrix, cpu_check);
          viennacl::tools::const_sparse_matrix_adapter<ScalarType> cpu_check_adapter(cpu_check);
          
          for (typename viennacl::tools::const_sparse_matrix_adapter<ScalarType>::const_iterator1 row_it = cpu_check_adapter.begin1();
                row_it != cpu_check_adapter.end1();
                ++row_it)
          {
            diag_M_inv_cpu[row_it.index1()] = 0;
            for (typename viennacl::tools::const_sparse_matrix_adapter<ScalarType>::const_iterator2 col_it = row_it.begin();
                  col_it != row_it.end();
                  ++col_it)
            {
              if (tag_.norm() == 1)
                diag_M_inv_cpu[col_it.index1()] += std::fabs(*col_it);
              else
                diag_M_inv_cpu[col_it.index1()] += (*col_it) * (*col_it);
            }
            if (diag_M_inv_cpu[row_it.index1()] == 0)
              throw "ViennaCL: Zero row encountered while setting up row scaling preconditioner!";
            
            if (tag_.norm() == 1)
              diag_M_inv_cpu[row_it.index1()] = static_cast<ScalarType>(1.0) / diag_M_inv_cpu[row_it.index1()];
            else
              diag_M_inv_cpu[row_it.index1()] = static_cast<ScalarType>(1.0) / std::sqrt(diag_M_inv_cpu[row_it.index1()]);
          }
          
          diag_M_inv.resize(system_matrix.size1(), false);
          viennacl::fast_copy(diag_M_inv_cpu, diag_M_inv);
        } */
        
        void init_gpu()
        {
          if (tag_.norm() == 1)
          {
            viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(
                                                viennacl::linalg::kernels::compressed_matrix<ScalarType, MAT_ALIGNMENT>::program_name(),
                                                "row_scaling_1");

            viennacl::ocl::enqueue( k(system_matrix.handle1(), system_matrix.handle2(), system_matrix.handle(), 
                                      diag_M_inv, static_cast<cl_uint>(diag_M_inv.size())) );        
          }
          else
          {
            viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(
                                                viennacl::linalg::kernels::compressed_matrix<ScalarType, MAT_ALIGNMENT>::program_name(),
                                                "row_scaling_2");

            viennacl::ocl::enqueue( k(system_matrix.handle1(), system_matrix.handle2(), system_matrix.handle(), 
                                      diag_M_inv, static_cast<cl_uint>(diag_M_inv.size())) );        
          }
        }
        
        template <unsigned int ALIGNMENT>
        void apply(viennacl::vector<ScalarType, ALIGNMENT> & vec) const
        {
          assert(viennacl::traits::size1(system_matrix) == viennacl::traits::size(vec));
          
          //run kernel:
          viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<ScalarType, ALIGNMENT>::program_name(),
                                                                "diag_precond");

          viennacl::ocl::enqueue(
             k(viennacl::traits::handle(diag_M_inv), cl_uint(viennacl::traits::start(diag_M_inv)), cl_uint(viennacl::traits::size(diag_M_inv)),
               viennacl::traits::handle(vec), cl_uint(viennacl::traits::start(vec)), cl_uint(viennacl::traits::size(vec)) )
                                );        
          
        }
        
      private:
        MatrixType const & system_matrix;
        row_scaling_tag const & tag_;
        viennacl::vector<ScalarType> diag_M_inv;
    };

  }
}




#endif



