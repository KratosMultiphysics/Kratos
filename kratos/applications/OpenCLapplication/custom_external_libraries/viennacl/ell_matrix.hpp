#ifndef VIENNACL_ELL_MATRIX_HPP_
#define VIENNACL_ELL_MATRIX_HPP_

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

/** @file ell_matrix.hpp
    @brief Implementation of the ell_matrix class
    
    Contributed by Volodymyr Kysenko.
*/


#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/backend.hpp"

#include "viennacl/linalg/kernels/ell_matrix_kernels.h"

namespace viennacl
{
    template<typename SCALARTYPE, unsigned int ALIGNMENT /* see forwards.h for default argument */>
    class ell_matrix
    {

      public:
        ell_matrix() 
        {
          viennacl::linalg::kernels::ell_matrix<SCALARTYPE, ALIGNMENT>::init();
        }
        
        ell_matrix(std::size_t row_num, std::size_t col_num) 
        {
          viennacl::linalg::kernels::ell_matrix<SCALARTYPE, ALIGNMENT>::init();
        }
    
      public:
        std::size_t internal_size1() const { return viennacl::tools::roundUpToNextMultiple<std::size_t>(rows_, ALIGNMENT); }
        std::size_t internal_size2() const { return viennacl::tools::roundUpToNextMultiple<std::size_t>(cols_, ALIGNMENT); }

        std::size_t size1() const { return rows_; }
        std::size_t size2() const { return cols_; }
        
        std::size_t internal_maxnnz() const {return viennacl::tools::roundUpToNextMultiple<std::size_t>(maxnnz_, ALIGNMENT); }
        std::size_t maxnnz() const { return maxnnz_; }

        std::size_t nnz() const { return rows_ * maxnnz_; }
        std::size_t internal_nnz() const { return internal_size1() * internal_maxnnz(); }

        const viennacl::ocl::handle<cl_mem>& handle1( ) const { return elements_; } 
        const viennacl::ocl::handle<cl_mem>& handle2() const { return coords_; }

      #if defined(_MSC_VER) && _MSC_VER < 1500          //Visual Studio 2005 needs special treatment
        template <typename CPU_MATRIX>
        friend void copy(const CPU_MATRIX & cpu_matrix, ell_matrix & gpu_matrix );
      #else
        template <typename CPU_MATRIX, typename T, unsigned int ALIGN>
        friend void copy(const CPU_MATRIX & cpu_matrix, ell_matrix<T, ALIGN> & gpu_matrix );
      #endif        
        
      private:
        std::size_t rows_;
        std::size_t cols_;
        std::size_t maxnnz_;

        viennacl::ocl::handle<cl_mem> coords_;
        viennacl::ocl::handle<cl_mem> elements_;        
    };

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX& cpu_matrix, ell_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix )
    {
      if(cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0)
      {
        //determine max capacity for row
        std::size_t max_entries_per_row = 0;
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
          std::size_t num_entries = 0;
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
          {
              ++num_entries;
          }

          max_entries_per_row = std::max(max_entries_per_row, num_entries);
        }

        //setup GPU matrix
        gpu_matrix.maxnnz_ = max_entries_per_row;
        gpu_matrix.rows_ = cpu_matrix.size1();
        gpu_matrix.cols_ = cpu_matrix.size2();

        std::size_t nnz = gpu_matrix.internal_nnz();

        std::vector<cl_uint> coords(nnz, 0);
        std::vector<SCALARTYPE> elements(nnz, 0);

        // std::cout << "ELL_MATRIX copy " << gpu_matrix.maxnnz_ << " " << gpu_matrix.rows_ << " " << gpu_matrix.cols_ << " " 
        //             << gpu_matrix.internal_maxnnz() << "\n";

        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
          std::size_t data_index = 0;
          
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
          {
            coords[gpu_matrix.internal_size1() * data_index + col_it.index1()]   = col_it.index2();
            elements[gpu_matrix.internal_size1() * data_index + col_it.index1()] = *col_it;
            //std::cout << *col_it << "\n";
              data_index++;
          }
        }


        gpu_matrix.coords_   = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY, coords);
        gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_ONLY, elements);
      }
    }

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const ell_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix, CPU_MATRIX& cpu_matrix)
    {
      if(gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0)
      {
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2());

        std::vector<SCALARTYPE> elements(gpu_matrix.internal_nnz());
        std::vector<cl_uint> coords(gpu_matrix.internal_nnz());

        cl_int err;

        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle1(), CL_TRUE, 0, sizeof(SCALARTYPE) * elements.size(), &(elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle2(), CL_TRUE, 0, sizeof(cl_uint) * coords.size(), &(coords[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);

        viennacl::ocl::get_queue().finish();

        for(std::size_t row = 0; row < gpu_matrix.size1(); row++)
        {
          for(std::size_t ind = 0; ind < gpu_matrix.internal_maxnnz(); ind++)
          {
            std::size_t offset = gpu_matrix.internal_size1() * ind + row;
            
            if(elements[offset] == static_cast<SCALARTYPE>(0.0))
            {
                continue;
            }

            if(coords[offset] >= gpu_matrix.size2())
            {
                std::cerr << "ViennaCL encountered invalid data " << offset << " " << ind << " " << row << " " << coords[offset] << " " << gpu_matrix.size2() << std::endl;
                return;
            }

            cpu_matrix(row, coords[offset]) = elements[offset];
          }
        }
      }
    }

    namespace linalg
    {
      /** @brief Returns a proxy class that represents matrix-vector multiplication with a hyb_matrix
      *
      * This is used for the convenience expression result = prod(mat, vec);
      *
      * @param mat    The matrix
      * @param vec    The vector
      */
      template<class SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
      vector_expression<const ell_matrix<SCALARTYPE, ALIGNMENT>,
                        const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                        op_prod > prod_impl(const ell_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                                            const vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec)
      {
        return vector_expression<const ell_matrix<SCALARTYPE, ALIGNMENT>,
                                 const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                 op_prod >(mat, vec);
      }
      
      template<class TYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
      void prod_impl(
                      const viennacl::ell_matrix<TYPE, ALIGNMENT>& mat, 
                      const viennacl::vector<TYPE, VECTOR_ALIGNMENT>& vec,
                      viennacl::vector<TYPE, VECTOR_ALIGNMENT>& result)
      {
        assert(mat.size1() == result.size());
        assert(mat.size2() == vec.size());

        result.clear();

        std::stringstream ss;
        ss << "vec_mul_" << 1;//(ALIGNMENT != 1?4:1);
        viennacl::ocl::kernel& k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::ell_matrix<TYPE, ALIGNMENT>::program_name(), "vec_mul");

        unsigned int thread_num = 128;
        unsigned int group_num = 256;

        k.local_work_size(0, thread_num);
        k.global_work_size(0, thread_num * group_num);

        viennacl::ocl::enqueue(k(mat.handle2(), 
                                 mat.handle1(),
                                 vec,
                                 result,
                                 cl_uint(mat.size1()),
                                 cl_uint(mat.size2()),
                                 cl_uint(mat.internal_size1()),
                                 cl_uint(mat.maxnnz()),
                                 cl_uint(mat.internal_maxnnz())
                                ) 
        );


      }
    }
    
    
    /** @brief Implementation of the operation v1 = A * v2, where A is a matrix
    *
    * @param proxy  An expression template proxy class.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    template <unsigned int MAT_ALIGNMENT>
    viennacl::vector<SCALARTYPE, ALIGNMENT> & 
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const ell_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                                          const viennacl::vector<SCALARTYPE, ALIGNMENT>,
                                                                                          viennacl::op_prod> & proxy) 
    {
      // check for the special case x = A * x
      if (proxy.rhs().handle().get() == this->handle().get())
      {
        viennacl::vector<SCALARTYPE, ALIGNMENT> result(proxy.rhs().size());
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
        *this = result;
        return *this;
      }
      else
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      return *this;
    }
    
}

#endif

















