#ifndef VIENNACL_COORDINATE_MATRIX_HPP_
#define VIENNACL_COORDINATE_MATRIX_HPP_

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

/** @file coordinate_matrix.hpp
    @brief Implementation of the coordinate_matrix class
*/

#include <map>
#include <vector>
#include <list>

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/vector.hpp"

#include "viennacl/linalg/coordinate_matrix_operations.hpp"

namespace viennacl
{
  
    
    //provide copy-operation:
    /** @brief Copies a sparse matrix from the host to the OpenCL device (either GPU or multi-core CPU)
    *
    * For the requirements on the CPU_MATRIX type, see the documentation of the function copy(CPU_MATRIX, compressed_matrix<>)
    *
    * @param cpu_matrix   A sparse matrix on the host.
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX & cpu_matrix,
                     coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
    {
      size_t group_num = 64;
      
      // Step 1: Determine nonzeros:
      if ( cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0 )
      {
        std::size_t num_entries = 0;
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            ++num_entries;
          }
        }
        
        // Step 2: Set up matrix data:
        std::cout << "Number of entries: " << num_entries << std::endl;
        gpu_matrix.nonzeros_ = num_entries;
        gpu_matrix.rows_ = cpu_matrix.size1();
        gpu_matrix.cols_ = cpu_matrix.size2();

        std::vector<cl_uint> coord_buffer(2*gpu_matrix.internal_nnz());
        std::vector<cl_uint> group_boundaries(group_num + 1);
        std::vector<SCALARTYPE> elements(gpu_matrix.internal_nnz());
        
        std::size_t data_index = 0;
        std::size_t current_fraction = 0;
        
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            coord_buffer[2*data_index] = static_cast<cl_uint>(col_it.index1());
            coord_buffer[2*data_index + 1] = static_cast<cl_uint>(col_it.index2());
            elements[data_index] = *col_it;
            ++data_index;
          }
          
          if (data_index > (current_fraction + 1) / static_cast<double>(group_num) * num_entries)    //split data equally over 64 groups
            group_boundaries[++current_fraction] = data_index;
        }
        
        //write end of last group:
        group_boundaries[group_num] = data_index;
        //group_boundaries[1] = data_index; //for one compute unit
        
        /*std::cout << "Group boundaries: " << std::endl;
        for (size_t i=0; i<group_boundaries.size(); ++i)
          std::cout << group_boundaries[i] << std::endl;*/
        
        gpu_matrix.coord_buffer_     = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, coord_buffer);
        gpu_matrix.elements_         = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, elements);
        gpu_matrix.group_boundaries_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, group_boundaries);
      }
    }

    /** @brief Copies a sparse matrix in the std::vector< std::map < > > format to an OpenCL device.
    *
    * @param cpu_matrix   A sparse square matrix on the host.
    * @param gpu_matrix   A coordinate_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const std::vector< std::map<unsigned int, SCALARTYPE> > & cpu_matrix,
                     coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
    {
      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(cpu_matrix, cpu_matrix.size(), cpu_matrix.size()), gpu_matrix);
    }
    
    //gpu to cpu:
    /** @brief Copies a sparse matrix from the OpenCL device (either GPU or multi-core CPU) to the host.
    *
    * There are two type requirements on the CPU_MATRIX type (fulfilled by e.g. boost::numeric::ublas):
    * - resize(rows, cols)  A resize function to bring the matrix into the correct size
    * - operator(i,j)       Write new entries via the parenthesis operator
    *
    * @param gpu_matrix   A coordinate_matrix from ViennaCL
    * @param cpu_matrix   A sparse matrix on the host.
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
                     CPU_MATRIX & cpu_matrix )
    {
      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2(), false);
        
        //get raw data from memory:
        std::vector<cl_uint> coord_buffer(2*gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        //std::cout << "GPU nonzeros: " << gpu_matrix.nnz() << std::endl;
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), gpu_matrix.handle12().get(), CL_TRUE, 0, sizeof(cl_uint)* 2 *gpu_matrix.nnz(), &(coord_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
        
        //fill the cpu_matrix:
        for (std::size_t index = 0; index < gpu_matrix.nnz(); ++index)
        {
          cpu_matrix(coord_buffer[2*index], coord_buffer[2*index+1]) = elements[index];
        }
      }
    }

    /** @brief Copies a sparse matrix from an OpenCL device to the host. The host type is the std::vector< std::map < > > format .
    *
    * @param gpu_matrix   A coordinate_matrix from ViennaCL
    * @param cpu_matrix   A sparse matrix on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              std::vector< std::map<unsigned int, SCALARTYPE> > & cpu_matrix)
    {
      tools::sparse_matrix_adapter<SCALARTYPE> temp(cpu_matrix, gpu_matrix.size1(), gpu_matrix.size2());
      copy(gpu_matrix, temp);
    }


    //////////////////////// coordinate_matrix //////////////////////////
    /** @brief A sparse square matrix, where entries are stored as triplets (i,j, val), where i and j are the row and column indices and val denotes the entry.
    *
    * The present implementation of coordinate_matrix suffers from poor runtime efficiency. Users are adviced to use compressed_matrix in the meanwhile.
    *
    * @tparam SCALARTYPE    The floating point type (either float or double, checked at compile time)
    * @tparam ALIGNMENT     The internal memory size for the arrays, given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT /* see VCLForwards.h */ >
    class coordinate_matrix
    {
    public:
      typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
      
      /** @brief Default construction of a coordinate matrix. No memory is allocated */
      coordinate_matrix() : rows_(0), cols_(0), nonzeros_(0) { viennacl::linalg::kernels::coordinate_matrix<SCALARTYPE, ALIGNMENT>::init(); }
      
      /** @brief Construction of a coordinate matrix with the supplied number of rows and columns. If the number of nonzeros is positive, memory is allocated
      *
      * @param rows     Number of rows
      * @param cols     Number of columns
      * @param nonzeros Optional number of nonzeros for memory preallocation
      */
      coordinate_matrix(std::size_t rows, std::size_t cols, std::size_t nonzeros = 0) : 
        rows_(rows), cols_(cols), nonzeros_(nonzeros)
      {
        viennacl::linalg::kernels::coordinate_matrix<SCALARTYPE, ALIGNMENT>::init();
        if (nonzeros > 0)
        {
          coord_buffer_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(cl_uint) * 2 * internal_nnz());
          elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * internal_nnz());
          group_boundaries_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(cl_uint) * (rows + 1));
        }
      }
        
      /** @brief Allocate memory for the supplied number of nonzeros in the matrix. Old values are preserved. */
      void reserve(std::size_t new_nonzeros)
      {
        if (new_nonzeros > nonzeros_)
        {
          viennacl::ocl::handle<cl_mem> coord_buffer_old = coord_buffer_;
          viennacl::ocl::handle<cl_mem> elements_old = elements_;
          coord_buffer_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(cl_uint) * 2 * internal_nnz());
          elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * internal_nnz());
          
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), coord_buffer_old.get(), coord_buffer_.get(), 0, 0, sizeof(cl_uint) * 2 * nonzeros_, 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), elements_old.get(), elements_.get(), 0, 0, sizeof(SCALARTYPE)*nonzeros_, 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);

          //new memory must be padded with zeros:
          std::vector<long> temp(internal_nnz() - nonzeros_);
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), coord_buffer_old.get(), coord_buffer_.get(), 0, nonzeros_, sizeof(cl_uint) * 2 * temp.size(), 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), elements_old.get(), elements_.get(), 0, nonzeros_, sizeof(SCALARTYPE)*temp.size(), 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
        }
      }

      /** @brief Resize the matrix.
      *
      * @param new_size1    New number of rows
      * @param new_size2    New number of columns
      * @param preserve     If true, the old values are preserved. At present, old values are always discarded.
      */
      void resize(std::size_t new_size1, std::size_t new_size2, bool preserve = true)
      {
        assert (new_size1 > 0 && new_size2 > 0);
                
        if (new_size1 < rows_ || new_size2 < cols_) //enlarge buffer
        {
          std::vector<std::map<unsigned int, SCALARTYPE> > stl_sparse_matrix;
          if (rows_ > 0)
            stl_sparse_matrix.resize(rows_);
          
          if (preserve && rows_ > 0)
            viennacl::copy(*this, stl_sparse_matrix);
            
          stl_sparse_matrix.resize(new_size1);
          
          std::cout << "Cropping STL matrix of size " << stl_sparse_matrix.size() << std::endl;
          if (new_size2 < cols_ && rows_ > 0)
          {
            for (std::size_t i=0; i<stl_sparse_matrix.size(); ++i)
            {
              std::list<unsigned int> to_delete;
              for (typename std::map<unsigned int, SCALARTYPE>::iterator it = stl_sparse_matrix[i].begin();
                   it != stl_sparse_matrix[i].end();
                  ++it)
              {
                if (it->first >= new_size2)
                  to_delete.push_back(it->first);
              }
              
              for (std::list<unsigned int>::iterator it = to_delete.begin(); it != to_delete.end(); ++it)
                stl_sparse_matrix[i].erase(*it);
            }
            std::cout << "Cropping done..." << std::endl;
          }
          
          rows_ = new_size1;
          cols_ = new_size2;
          viennacl::copy(stl_sparse_matrix, *this);
        }
          
        rows_ = new_size1;
        cols_ = new_size2;
      }


      /** @brief  Returns the number of rows */
      std::size_t size1() const { return rows_; }
      /** @brief  Returns the number of columns */
      std::size_t size2() const { return cols_; }
      /** @brief  Returns the number of nonzero entries */
      std::size_t nnz() const { return nonzeros_; }
      /** @brief  Returns the number of internal nonzero entries */
      std::size_t internal_nnz() const { return viennacl::tools::roundUpToNextMultiple<std::size_t>(nonzeros_, ALIGNMENT);; }
      
      /** @brief  Returns the OpenCL handle to the (row, column) index array */
      const viennacl::ocl::handle<cl_mem> & handle12() const { return coord_buffer_; }
      /** @brief  Returns the OpenCL handle to the matrix entry array */
      const viennacl::ocl::handle<cl_mem> & handle() const { return elements_; }
      /** @brief  Returns the OpenCL handle to the group start index array */
      const viennacl::ocl::handle<cl_mem> & handle3() const { return group_boundaries_; }
      
      #if defined(_MSC_VER) && _MSC_VER < 1500      //Visual Studio 2005 needs special treatment
      template <typename CPU_MATRIX>
      friend void copy(const CPU_MATRIX & cpu_matrix, coordinate_matrix & gpu_matrix );
      #else
      template <typename CPU_MATRIX, typename SCALARTYPE2, unsigned int ALIGNMENT2>
      friend void copy(const CPU_MATRIX & cpu_matrix, coordinate_matrix<SCALARTYPE2, ALIGNMENT2> & gpu_matrix );
      #endif

    private:
      /** @brief Copy constructor is by now not available. */
      coordinate_matrix(coordinate_matrix const &);
      
      /** @brief Assignment is by now not available. */
      coordinate_matrix & operator=(coordinate_matrix const &);
      
      
      std::size_t rows_;
      std::size_t cols_;
      std::size_t nonzeros_;
      viennacl::ocl::handle<cl_mem> coord_buffer_;
      viennacl::ocl::handle<cl_mem> elements_;
      viennacl::ocl::handle<cl_mem> group_boundaries_;
    };


}

#endif
