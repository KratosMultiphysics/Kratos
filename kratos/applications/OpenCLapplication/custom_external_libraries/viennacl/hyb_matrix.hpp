#ifndef VIENNACL_HYB_MATRIX_HPP_
#define VIENNACL_HYB_MATRIX_HPP_

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

/** @file hyb_matrix.hpp
    @brief Implementation of the hyb_matrix class
    
    Contributed by Volodymyr Kysenko.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"

#include "viennacl/tools/tools.hpp"
#include "viennacl/ocl/backend.hpp"

#include "viennacl/linalg/kernels/hyb_matrix_kernels.h"

namespace viennacl
{
    template<typename SCALARTYPE, unsigned int ALIGNMENT  /* see forwards.h for default argument */>
    class hyb_matrix
    {

      public:
        hyb_matrix() : csr_threshold_(0.8), rows_(0), cols_(0) 
        {
          viennacl::linalg::kernels::hyb_matrix<SCALARTYPE, ALIGNMENT>::init();
        }
        
        hyb_matrix(std::size_t row_num, std::size_t col_num) : csr_threshold_(0.8), rows_(row_num), cols_(col_num)
        {
          viennacl::linalg::kernels::hyb_matrix<SCALARTYPE, ALIGNMENT>::init();
        }

        SCALARTYPE  csr_threshold()  const { return csr_threshold_; }
        void csr_threshold(SCALARTYPE thr) { csr_threshold_ = thr; }

        std::size_t internal_size1() const { return viennacl::tools::roundUpToNextMultiple<std::size_t>(rows_, ALIGNMENT); }
        std::size_t internal_size2() const { return viennacl::tools::roundUpToNextMultiple<std::size_t>(cols_, ALIGNMENT); }

        std::size_t size1() const { return rows_; }
        std::size_t size2() const { return cols_; }

        std::size_t internal_ellnnz() const {return viennacl::tools::roundUpToNextMultiple<std::size_t>(ellnnz_, ALIGNMENT); }
        std::size_t ell_nnz() const { return ellnnz_; }
        std::size_t csr_nnz() const { return csrnnz_; }

        const viennacl::ocl::handle<cl_mem>& handle1() const { return ell_elements_; } 
        const viennacl::ocl::handle<cl_mem>& handle2() const { return ell_coords_; }
        const viennacl::ocl::handle<cl_mem>& handle3() const { return csr_rows_; } 
        const viennacl::ocl::handle<cl_mem>& handle4() const { return csr_cols_; } 
        const viennacl::ocl::handle<cl_mem>& handle5() const { return csr_elements_; }  
    
      public:    
      #if defined(_MSC_VER) && _MSC_VER < 1500          //Visual Studio 2005 needs special treatment
        template <typename CPU_MATRIX>
        friend void copy(const CPU_MATRIX & cpu_matrix, hyb_matrix & gpu_matrix );
      #else
        template <typename CPU_MATRIX, typename T, unsigned int ALIGN>
        friend void copy(const CPU_MATRIX & cpu_matrix, hyb_matrix<T, ALIGN> & gpu_matrix );
      #endif
        
      private:
        SCALARTYPE  csr_threshold_;
        std::size_t rows_;
        std::size_t cols_;
        std::size_t ellnnz_;
        std::size_t csrnnz_;

        viennacl::ocl::handle<cl_mem> ell_coords_; // ell coords
        viennacl::ocl::handle<cl_mem> ell_elements_; // ell elements
        
        viennacl::ocl::handle<cl_mem> csr_rows_;
        viennacl::ocl::handle<cl_mem> csr_cols_;
        viennacl::ocl::handle<cl_mem> csr_elements_;
    };

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX& cpu_matrix, hyb_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix )
    {
      if(cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0)
      {
        //determine max capacity for row
        std::size_t max_entries_per_row = 0;
        std::vector<std::size_t> hist_entries(cpu_matrix.size1(), 0);

        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
            std::size_t num_entries = 0;
            for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
            {
                ++num_entries;
            }

            hist_entries[num_entries] += 1;
            max_entries_per_row = std::max(max_entries_per_row, num_entries);
        }
        
        std::size_t sum = 0;
        for(std::size_t ind = 0; ind <= max_entries_per_row; ind++)
        {
            sum += hist_entries[ind];
            
            if(sum >= gpu_matrix.csr_threshold() * cpu_matrix.size1())
            {
                max_entries_per_row = ind;
                break;
            }
            }

        //setup GPU matrix
        gpu_matrix.ellnnz_ = max_entries_per_row;
        gpu_matrix.rows_ = cpu_matrix.size1();
        gpu_matrix.cols_ = cpu_matrix.size2();

        std::size_t nnz = gpu_matrix.internal_size1() * gpu_matrix.internal_ellnnz();

        std::vector<cl_uint> ell_coords(nnz, 0);
        std::vector<cl_uint> csr_rows(cpu_matrix.size1() + 1, 0);
        std::vector<cl_uint> csr_cols;

        std::vector<SCALARTYPE> ell_elements(nnz, 0.0f);
        std::vector<SCALARTYPE> csr_elements;

        std::size_t csr_index = 0;

        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1(); row_it != cpu_matrix.end1(); ++row_it)
        {
          std::size_t data_index = 0;
  
          csr_rows[row_it.index1()] = csr_index;
          
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin(); col_it != row_it.end(); ++col_it)
          {
            if(data_index < max_entries_per_row)
            {
                ell_coords[gpu_matrix.internal_size1() * data_index + col_it.index1()]   = col_it.index2();
                ell_elements[gpu_matrix.internal_size1() * data_index + col_it.index1()] = *col_it;                        
            }
            else
            {
                csr_cols.push_back(col_it.index2());
                csr_elements.push_back(*col_it);

                csr_index++;
            }

            data_index++;
          }

        }

        if(csr_cols.empty())
        {
          csr_cols.push_back(0);
          csr_elements.push_back(0);
        }

        csr_rows[csr_rows.size() - 1] = csr_index;

        gpu_matrix.csrnnz_ = csr_cols.size();

        gpu_matrix.ell_coords_   = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, ell_coords);
        gpu_matrix.ell_elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, ell_elements);

        gpu_matrix.csr_rows_   = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, csr_rows);
        gpu_matrix.csr_cols_   = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, csr_cols);
        gpu_matrix.csr_elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, csr_elements);

      }
    }

    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const hyb_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix, CPU_MATRIX& cpu_matrix)
    {
      if(gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0)
      {
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2());

        std::vector<SCALARTYPE> ell_elements(gpu_matrix.internal_size1() * gpu_matrix.internal_ellnnz());
        std::vector<cl_uint> ell_coords(gpu_matrix.internal_size1() * gpu_matrix.internal_ellnnz());

        std::vector<SCALARTYPE> csr_elements(gpu_matrix.csr_nnz());
        std::vector<cl_uint> csr_rows(gpu_matrix.size1() + 1);
        std::vector<cl_uint> csr_cols(gpu_matrix.csr_nnz());

        cl_int err;

        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle1(), CL_TRUE, 0, sizeof(SCALARTYPE) * ell_elements.size(), &(ell_elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle2(), CL_TRUE, 0, sizeof(cl_uint) * ell_coords.size(), &(ell_coords[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle3(), CL_TRUE, 0, sizeof(cl_uint) * csr_rows.size(), &(csr_rows[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle4(), CL_TRUE, 0, sizeof(cl_uint) * csr_cols.size(), &(csr_cols[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle5(), CL_TRUE, 0, sizeof(SCALARTYPE) * csr_elements.size(), &(csr_elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);

        viennacl::ocl::get_queue().finish();

        for(std::size_t row = 0; row < gpu_matrix.size1(); row++)
        {
          for(std::size_t ind = 0; ind < gpu_matrix.internal_ellnnz(); ind++)
          {
            std::size_t offset = gpu_matrix.internal_size1() * ind + row;
            
            if(ell_elements[offset] == static_cast<SCALARTYPE>(0.0))
            {
              continue;
            }

            if(ell_coords[offset] >= gpu_matrix.size2())
            {
              std::cerr << "ViennaCL encountered invalid data " << offset << " " << ind << " " << row << " " << ell_coords[offset] << " " << gpu_matrix.size2() << std::endl;
              return;
            }

            cpu_matrix(row, ell_coords[offset]) = ell_elements[offset];
          }

          for(std::size_t ind = csr_rows[row]; ind < csr_rows[row+1]; ind++)
          {
            if(csr_elements[ind] == static_cast<SCALARTYPE>(0.0))
            {
              continue;
            }

            if(csr_cols[ind] >= gpu_matrix.size2())
            {
              std::cerr << "ViennaCL encountered invalid data " << std::endl;
              return;
            }

            cpu_matrix(row, csr_cols[ind]) = csr_elements[ind];
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
      vector_expression<const hyb_matrix<SCALARTYPE, ALIGNMENT>,
                        const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                        op_prod > prod_impl(const hyb_matrix<SCALARTYPE, ALIGNMENT> & mat, 
                                      const vector<SCALARTYPE, VECTOR_ALIGNMENT> & vec)
      {
        return vector_expression<const hyb_matrix<SCALARTYPE, ALIGNMENT>,
                                const vector<SCALARTYPE, VECTOR_ALIGNMENT>, 
                                op_prod >(mat, vec);
      }
      
      template<class TYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
      void prod_impl( const viennacl::hyb_matrix<TYPE, ALIGNMENT>& mat, 
                      const viennacl::vector<TYPE, VECTOR_ALIGNMENT>& vec,
                      viennacl::vector<TYPE, VECTOR_ALIGNMENT>& result)
      {
        assert(mat.size1() == result.size());
        assert(mat.size2() == vec.size());

        result.clear();

        viennacl::ocl::kernel& k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::hyb_matrix<TYPE, ALIGNMENT>::program_name(), "vec_mul");

        unsigned int thread_num = 256;
        unsigned int group_num = 32;

        k.local_work_size(0, thread_num);
        k.global_work_size(0, thread_num * group_num);

        viennacl::ocl::enqueue(k(mat.handle2(), 
                                mat.handle1(),
                                mat.handle3(),
                                mat.handle4(),
                                mat.handle5(),
                                vec,
                                result,
                                cl_uint(mat.size1()),
                                cl_uint(mat.internal_size1()),
                                cl_uint(mat.ell_nnz()),
                                cl_uint(mat.internal_ellnnz())
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
    viennacl::vector<SCALARTYPE, ALIGNMENT>::operator=(const viennacl::vector_expression< const hyb_matrix<SCALARTYPE, MAT_ALIGNMENT>,
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