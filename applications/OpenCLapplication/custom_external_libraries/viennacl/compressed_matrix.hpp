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
======================================================================= */

#ifndef _VIENNACL_COMPRESSED_MATRIX_HPP_
#define _VIENNACL_COMPRESSED_MATRIX_HPP_

/** @file compressed_matrix.hpp
    @brief Implementation of the compressed_matrix class
*/

#include <vector>
#include <list>
#include <map>
#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/vector.hpp"

#include "viennacl/linalg/compressed_matrix_operations.hpp"

#include "viennacl/tools/tools.hpp"

namespace viennacl
{
    

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
      //std::cout << "copy for (" << cpu_matrix.size1() << ", " << cpu_matrix.size2() << ", " << cpu_matrix.nnz() << ")" << std::endl;
      
      if ( cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0 )
      {
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
        
        if (num_entries == 0) //we copy an empty matrix
        {
          num_entries = 1;
        }
        
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
    
    #ifdef VIENNACL_HAVE_EIGEN
    template <typename SCALARTYPE, int flags, unsigned int ALIGNMENT>
    void copy(const Eigen::SparseMatrix<SCALARTYPE, flags> & eigen_matrix,
              compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix)
    {
      std::vector< std::map<unsigned int, SCALARTYPE> >  stl_matrix(eigen_matrix.rows());
      
      for (int k=0; k < eigen_matrix.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<SCALARTYPE, flags>::InnerIterator it(eigen_matrix, k); it; ++it)
          stl_matrix[it.row()][it.col()] = it.value();
        
      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(stl_matrix), gpu_matrix);
    }
    #endif
    
    
    #ifdef VIENNACL_HAVE_MTL4
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const mtl::compressed2D<SCALARTYPE> & cpu_matrix,
              compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix)
    {
      typedef mtl::compressed2D<SCALARTYPE>  MatrixType;
      
      std::vector< std::map<unsigned int, SCALARTYPE> >  stl_matrix(cpu_matrix.num_rows());
      
      using mtl::traits::range_generator;
      using mtl::traits::range::min;

      // Choose between row and column traversal
      typedef typename min<range_generator<mtl::tag::row, MatrixType>,
                           range_generator<mtl::tag::col, MatrixType> >::type   range_type;
      range_type                                                      my_range;

      // Type of outer cursor
      typedef typename range_type::type                               c_type;
      // Type of inner cursor
      typedef typename mtl::traits::range_generator<mtl::tag::nz, c_type>::type ic_type;

      // Define the property maps
      typename mtl::traits::row<MatrixType>::type                              row(cpu_matrix); 
      typename mtl::traits::col<MatrixType>::type                              col(cpu_matrix);
      typename mtl::traits::const_value<MatrixType>::type                      value(cpu_matrix); 

      // Now iterate over the matrix    
      for (c_type cursor(my_range.begin(cpu_matrix)), cend(my_range.end(cpu_matrix)); cursor != cend; ++cursor)
        for (ic_type icursor(mtl::begin<mtl::tag::nz>(cursor)), icend(mtl::end<mtl::tag::nz>(cursor)); icursor != icend; ++icursor)
          stl_matrix[row(*icursor)][col(*icursor)] = value(*icursor);
      
      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(stl_matrix), gpu_matrix);
    }
    #endif
    
    
    
    
    
    
    
    //
    // gpu to cpu:
    //
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
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2(), false);
        
        //get raw data from memory:
        std::vector<unsigned int> row_buffer(gpu_matrix.size1() + 1);
        std::vector<unsigned int> col_buffer(gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        //std::cout << "GPU->CPU, nonzeros: " << gpu_matrix.nnz() << std::endl;
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle1(), CL_TRUE, 0, sizeof(unsigned int)*(gpu_matrix.size1() + 1), &(row_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle2(), CL_TRUE, 0, sizeof(unsigned int)*gpu_matrix.nnz(), &(col_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
        
        //fill the cpu_matrix:
        unsigned int data_index = 0;
        for (unsigned int row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            if (col_buffer[data_index] >= gpu_matrix.size2())
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
    
    
    #ifdef VIENNACL_HAVE_EIGEN
    template <typename SCALARTYPE, int flags, unsigned int ALIGNMENT>
    void copy(compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              Eigen::SparseMatrix<SCALARTYPE, flags> & eigen_matrix)
    {
      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        assert(static_cast<unsigned int>(eigen_matrix.rows()) >= gpu_matrix.size1()
               && static_cast<unsigned int>(eigen_matrix.cols()) >= gpu_matrix.size2()
               && "Provided Eigen compressed matrix is too small!");
        
        //get raw data from memory:
        std::vector<unsigned int> row_buffer(gpu_matrix.size1() + 1);
        std::vector<unsigned int> col_buffer(gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle1(),
                                  CL_TRUE, 0, sizeof(unsigned int)*(gpu_matrix.size1() + 1), &(row_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle2(),
                                  CL_TRUE, 0, sizeof(unsigned int)*gpu_matrix.nnz(), &(col_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle(),
                                  CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
        
        eigen_matrix.setZero();
        eigen_matrix.startFill();
        unsigned int data_index = 0;
        for (unsigned int row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            assert(col_buffer[data_index] < gpu_matrix.size2() && "ViennaCL encountered invalid data at col_buffer");
            if (elements[data_index] != static_cast<SCALARTYPE>(0.0))
              eigen_matrix.fill(row-1, col_buffer[data_index]) = elements[data_index];
            ++data_index;
          }
        }
        eigen_matrix.endFill();
      }
    }
    #endif
    
    
    
    #ifdef VIENNACL_HAVE_MTL4
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              mtl::compressed2D<SCALARTYPE> & mtl4_matrix)
    {
      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        assert(mtl4_matrix.num_rows() >= gpu_matrix.size1()
               && mtl4_matrix.num_cols() >= gpu_matrix.size2()
               && "Provided MTL4 compressed matrix is too small!");
        
        //get raw data from memory:
        std::vector<unsigned int> row_buffer(gpu_matrix.size1() + 1);
        std::vector<unsigned int> col_buffer(gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle1(),
                                  CL_TRUE, 0, sizeof(unsigned int)*(gpu_matrix.size1() + 1), &(row_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle2(),
                                  CL_TRUE, 0, sizeof(unsigned int)*gpu_matrix.nnz(), &(col_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), gpu_matrix.handle(),
                                  CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
        
        //set_to_zero(mtl4_matrix);  
        //mtl4_matrix.change_dim(gpu_matrix.size1(), gpu_matrix.size2());
        
        mtl::matrix::inserter< mtl::compressed2D<SCALARTYPE> >  ins(mtl4_matrix);
        unsigned int data_index = 0;
        for (unsigned int row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            assert(col_buffer[data_index] < gpu_matrix.size2() && "ViennaCL encountered invalid data at col_buffer");
            if (elements[data_index] != static_cast<SCALARTYPE>(0.0))
              ins(row-1, col_buffer[data_index]) << typename mtl::Collection< mtl::compressed2D<SCALARTYPE> >::value_type(elements[data_index]);
            ++data_index;
          }
        }
      }
    }
    #endif
    
    
    
    
    
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
          _row_buffer = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(unsigned int) * rows);
        if (nonzeros > 0)
        {
          _col_buffer = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(unsigned int) * nonzeros);
          _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * nonzeros);
        }
      }
      
      explicit compressed_matrix(cl_mem mem_row_buffer, cl_mem mem_col_buffer, cl_mem mem_elements, unsigned int rows, unsigned int cols, unsigned int nonzeros) : 
        _rows(rows), _cols(cols), _nonzeros(nonzeros)
      {
          _row_buffer = mem_row_buffer;
          _row_buffer.inc();             //prevents that the user-provided memory is deleted once the matrix object is destroyed.
          _col_buffer = mem_col_buffer;
          _col_buffer.inc();             //prevents that the user-provided memory is deleted once the matrix object is destroyed.
          _elements = mem_elements;
          _elements.inc();               //prevents that the user-provided memory is deleted once the matrix object is destroyed.
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
        //std::cout << "Setting memory: " << cols + 1 << ", " << nonzeros << std::endl;
        _row_buffer = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(unsigned int) * (cols + 1), row_jumper);
        _col_buffer = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(unsigned int) * nonzeros, col_buffer);
        _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * nonzeros, elements);
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
          _col_buffer = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(unsigned int) * new_nonzeros);
          _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * new_nonzeros);
          
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), _col_buffer_old, _col_buffer, 0, 0, sizeof(unsigned int)*_nonzeros, 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), _elements_old, _elements, 0, 0, sizeof(SCALARTYPE)*_nonzeros, 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);

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
        assert(new_size1 > 0 && new_size2 > 0);
        //std::cout << "Resizing from (" << _rows << ", " << _cols << ") to (" << new_size1 << ", " << new_size2 << ")" << std::endl;
        
        if (new_size1 != _rows || new_size2 != _cols)
        {
          std::vector<std::map<unsigned int, SCALARTYPE> > stl_sparse_matrix;
          if (_rows > 0)
            stl_sparse_matrix.resize(_rows);
          
          if (preserve && _rows > 0)
            viennacl::copy(*this, stl_sparse_matrix);
            
          stl_sparse_matrix.resize(new_size1);
          
          //discard entries with column index larger than new_size2
          if (new_size2 < _cols && _rows > 0)
          {
            for (size_t i=0; i<stl_sparse_matrix.size(); ++i)
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
          }
          
          copy(stl_sparse_matrix, *this);
          
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

    
    

}

#endif
