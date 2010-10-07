/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_COMPRESSED_MATRIX_HPP_
#define _VIENNACL_COMPRESSED_MATRIX_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/kernel.hpp"

#include "viennacl/linalg/compressed_matrix_operations.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl
{
    //////////////////////// compressed_matrix //////////////////////////
    /** @brief A sparse square matrix in compressed sparse rows format.
    *
    * @tparam SCALARTYPE    The floating point type (either float or double, checked at compile time)
    * @tparam ALIGNMENT     The internal memory size for the entries in each row is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT /* see VCLForwards.h */>
    class compressed_matrix
    {
    public:
      typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
      
      /** @brief Default construction of a compressed matrix. No memory is allocated */
      compressed_matrix() : _rows(0), _cols(0), _nonzeros(0) { viennacl::linalg::kernels::compressed_matrix<SCALARTYPE, ALIGNMENT>::init(); }
      
      /** @brief Construction of a compressed matrix with the supplied number of rows and columns. If the number of nonzeros is positive, memory is allocated
      *
      * @param rows     Number of rows
      * @param cols     Number of columns
      * @param nonzeros Optional number of nonzeros for memory preallocation
      */
      explicit compressed_matrix(unsigned int rows, unsigned int cols, unsigned int nonzeros = 0) : 
        _rows(rows), _cols(cols), _nonzeros(nonzeros)
      {
        viennacl::linalg::kernels::compressed_matrix<SCALARTYPE, ALIGNMENT>::init();
        
        if (rows > 0)
          _row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * rows);
        if (nonzeros > 0)
        {
          _col_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * nonzeros);
          _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * nonzeros);
        }
      }
      
      /** @brief Sets the row, column and value arrays of the compressed matrix
      *
      * @param row_jumper     Pointer to an array holding the indices of the first element of each row (starting with zero). E.g. row_jumper[10] returns the index of the first entry of the 11th row. The array length is 'cols + 1'
      * @param col_buffer     Pointer to an array holding the column index of each entry. The array length is 'nonzeros'
      * @param elements       Pointer to an array holding the entries of the sparse matrix. The array length is 'elements'
      * @param cols           Number of columns (and rows) of the sparse matrix
      * @param nonzeros       Number of nonzeros
      */
      void set(unsigned int * row_jumper, unsigned int * col_buffer, SCALARTYPE * elements, unsigned int cols, unsigned int nonzeros)
      {
        assert(cols > 0);
        assert(nonzeros > 0);
        _row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * (cols + 1), row_jumper);
        _col_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * nonzeros, col_buffer);
        _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * nonzeros, elements);
        _nonzeros = nonzeros;
        _rows = cols;
        _cols = cols;
      }
        
      /** @brief Allocate memory for the supplied number of nonzeros in the matrix. Old values are preserved. */
      void reserve(unsigned int new_nonzeros)
      {
        if (new_nonzeros > _nonzeros)
        {
          viennacl::ocl::handle<cl_mem> _col_buffer_old = _col_buffer;
          viennacl::ocl::handle<cl_mem> _elements_old = _elements;
          _col_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * new_nonzeros);
          _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * new_nonzeros);
          
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _col_buffer_old.get(), _col_buffer.get(), 0, 0, sizeof(unsigned int)*_nonzeros, 0, NULL, NULL);
          CL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _elements_old.get(), _elements.get(), 0, 0, sizeof(SCALARTYPE)*_nonzeros, 0, NULL, NULL);
          CL_ERR_CHECK(err);

          _nonzeros = new_nonzeros;
        }
      }

      /** @brief Resize the matrix.
      *
      * @param new_size1    New number of rows
      * @param new_size2    New number of columns
      * @param preserve     If true, the old values are preserved. At present, old values are always discarded.
      */
      void resize(unsigned int new_size1, unsigned int new_size2, bool preserve = true)
      {
        if (new_size1 > 0 && new_size2 > 0)
        {
          if (new_size1 > _rows) //enlarge buffer
          {
            if (_rows == 0)
            {
              _row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int)*(new_size1 + 1));
              _col_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int)*(new_size1 + 1));
              _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*(new_size1 + 1));
              //set new memory to zero:
              std::vector<unsigned int> coord_temp(new_size1 + 1);
              std::vector<SCALARTYPE> temp(new_size1 + 1);
              cl_int err;
              err = clEnqueueWriteBuffer(viennacl::ocl::device().queue().get(), _row_buffer.get(), CL_TRUE, 0, sizeof(unsigned int)*coord_temp.size(), &(coord_temp[0]), 0, NULL, NULL);
              CL_ERR_CHECK(err);
              err = clEnqueueWriteBuffer(viennacl::ocl::device().queue().get(), _col_buffer.get(), CL_TRUE, 0, sizeof(unsigned int)*coord_temp.size(), &(coord_temp[0]), 0, NULL, NULL);
              CL_ERR_CHECK(err);
              err = clEnqueueWriteBuffer(viennacl::ocl::device().queue().get(), _elements.get(), CL_TRUE, 0, sizeof(SCALARTYPE)*temp.size(), &(temp[0]), 0, NULL, NULL);
              CL_ERR_CHECK(err);
            }
            else //enlarge only row array, because no entries are added to the matrix
            {
              viennacl::ocl::handle<cl_mem> _row_buffer_old = _row_buffer;
              _row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int)*(new_size1 + 1));
              cl_int err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _row_buffer_old.get(), _row_buffer.get(), 0, 0, sizeof(unsigned int)* (_rows + 1), 0, NULL, NULL);
              CL_ERR_CHECK(err);
              
              //set new memory to zero:
              std::vector<SCALARTYPE> temp(new_size1 - _rows + 1);
              err = clEnqueueWriteBuffer(viennacl::ocl::device().queue().get(), _elements.get(), CL_TRUE, sizeof(SCALARTYPE)*(_rows + 1), sizeof(SCALARTYPE)*temp.size(), &(temp[0]), 0, NULL, NULL);
              CL_ERR_CHECK(err);
            }
          }
          else if (new_size1 < _rows) //reduce buffer
          {
              viennacl::ocl::handle<cl_mem> _row_buffer_old = _row_buffer;
              _row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int)*(new_size1 + 1));
              cl_int err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _row_buffer_old.get(), _row_buffer.get(), 0, 0, sizeof(unsigned int) * (new_size1 + 1), 0, NULL, NULL);
              CL_ERR_CHECK(err);
            
            //TODO: discard entries in the matrix that are beyond the allowed sizes
          }
          
          _rows = new_size1;
          _cols = new_size2;
        }
      }

      /** @brief  Returns the number of rows */
      const unsigned int & size1() const { return _rows; }
      /** @brief  Returns the number of columns */
      const unsigned int & size2() const { return _cols; }
      /** @brief  Returns the number of nonzero entries */
      const unsigned int & nnz() const { return _nonzeros; }
      
      /** @brief  Returns the OpenCL handle to the row index array */
      const viennacl::ocl::handle<cl_mem> & handle1() const { return _row_buffer; }
      /** @brief  Returns the OpenCL handle to the column index array */
      const viennacl::ocl::handle<cl_mem> & handle2() const { return _col_buffer; }
      /** @brief  Returns the OpenCL handle to the matrix entry array */
      const viennacl::ocl::handle<cl_mem> & handle() const { return _elements; }
      
    private:
      unsigned int _rows;
      unsigned int _cols;
      unsigned int _nonzeros;
      viennacl::ocl::handle<cl_mem> _row_buffer;
      viennacl::ocl::handle<cl_mem> _col_buffer;
      viennacl::ocl::handle<cl_mem> _elements;
    };

    
    
    

    //provide copy-operation:
    /** @brief Copies a sparse matrix from the host to the OpenCL device (either GPU or multi-core CPU)
    *
    * There are some type requirements on the CPU_MATRIX type (fulfilled by e.g. boost::numeric::ublas):
    * - .size1() returns the number of rows
    * - .size2() returns the number of columns
    * - const_iterator1    is a type definition for an iterator along increasing row indices
    * - const_iterator2    is a type definition for an iterator along increasing columns indices
    * - The const_iterator1 type provides an iterator of type const_iterator2 via members .begin() and .end() that iterates along column indices in the current row.
    * - The types const_iterator1 and const_iterator2 provide members functions .index1() and .index2() that return the current row and column indices respectively.
    * - Dereferenciation of an object of type const_iterator2 returns the entry.
    *
    * @param cpu_matrix   A sparse matrix on the host.
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX & cpu_matrix,
                     compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
    {
      if ( cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0 )
      {
        gpu_matrix.resize(static_cast<unsigned int>(cpu_matrix.size1()), static_cast<unsigned int>(cpu_matrix.size2()), false);
        
        //determine nonzeros:
        long num_entries = 0;
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          unsigned int entries_per_row = 0;
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            ++entries_per_row;
          }
          num_entries += viennacl::tools::roundUpToNextMultiple<unsigned int>(entries_per_row, ALIGNMENT);
        }
        
        //std::cout << "CPU->GPU, Number of entries: " << num_entries << std::endl;
        
        //set up matrix entries:
        std::vector<unsigned int> row_buffer(cpu_matrix.size1() + 1);
        std::vector<unsigned int> col_buffer(num_entries);
        std::vector<SCALARTYPE> elements(num_entries);
        
        unsigned int row_index = 0;
        unsigned int data_index = 0;
        
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          row_buffer[row_index] = data_index;
          ++row_index;
          
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            col_buffer[data_index] = static_cast<unsigned int>(col_it.index2());
            elements[data_index] = *col_it;
            ++data_index;
          }
          data_index = viennacl::tools::roundUpToNextMultiple<unsigned int>(data_index, ALIGNMENT); //take care of alignment
        }
        row_buffer[row_index] = data_index;
        
        /*gpu_matrix._row_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, row_buffer);
        gpu_matrix._col_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, col_buffer);
        gpu_matrix._elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, elements);
        
        gpu_matrix._nonzeros = num_entries;*/
        gpu_matrix.set(&row_buffer[0], &col_buffer[0], &elements[0], static_cast<unsigned int>(cpu_matrix.size1()), num_entries);
      }
    }
    
    
    //adapted for std::vector< std::map < > > argument:
    /** @brief Copies a sparse matrix in the std::vector< std::map < > > format to an OpenCL device.
    *
    * @param cpu_matrix   A sparse matrix on the host.
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const std::vector< std::map<unsigned int, SCALARTYPE> > & cpu_matrix,
                     compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
    {
      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(cpu_matrix), gpu_matrix);
    }
    
    //gpu to cpu:
    /** @brief Copies a sparse matrix from the OpenCL device (either GPU or multi-core CPU) to the host.
    *
    * There are two type requirements on the CPU_MATRIX type (fulfilled by e.g. boost::numeric::ublas):
    * - resize(rows, cols)  A resize function to bring the matrix into the correct size
    * - operator(i,j)       Write new entries via the parenthesis operator
    *
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    * @param cpu_matrix   A sparse matrix on the host.
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
                     CPU_MATRIX & cpu_matrix )
    {
      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2());
        
        //get raw data from memory:
        std::vector<unsigned int> row_buffer(gpu_matrix.size1() + 1);
        std::vector<unsigned int> col_buffer(gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        //std::cout << "GPU->CPU, nonzeros: " << gpu_matrix.nnz() << std::endl;
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle1().get(), CL_TRUE, 0, sizeof(unsigned int)*(gpu_matrix.size1() + 1), &(row_buffer[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle2().get(), CL_TRUE, 0, sizeof(unsigned int)*gpu_matrix.nnz(), &(col_buffer[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        viennacl::ocl::finish();
        
        //fill the cpu_matrix:
        unsigned int data_index = 0;
        for (unsigned int row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            if (col_buffer[data_index] >= gpu_matrix.size1())
            {
              std::cerr << "ViennaCL encountered invalid data at colbuffer[" << data_index << "]: " << col_buffer[data_index] << std::endl;
              return;
            }
            
            if (elements[data_index] != static_cast<SCALARTYPE>(0.0))
              cpu_matrix(row-1, col_buffer[data_index]) = elements[data_index];
            ++data_index;
          }
        }
      }
    }
    
    /** @brief Copies a sparse matrix from an OpenCL device to the host. The host type is the std::vector< std::map < > > format .
    *
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    * @param cpu_matrix   A sparse matrix on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              std::vector< std::map<unsigned int, SCALARTYPE> > & cpu_matrix)
    {
      tools::sparse_matrix_adapter<SCALARTYPE> temp(cpu_matrix);
      copy(gpu_matrix, temp);
    }
    

}

#endif
