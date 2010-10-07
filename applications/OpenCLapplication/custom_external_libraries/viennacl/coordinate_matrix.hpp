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

#ifndef _VIENNACL_COORDINATE_MATRIX_HPP_
#define _VIENNACL_COORDINATE_MATRIX_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/kernel.hpp"

#include "viennacl/linalg/coordinate_matrix_operations.hpp"

namespace viennacl
{
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
      coordinate_matrix() : _rows(0), _cols(0), _nonzeros(0) { viennacl::linalg::kernels::coordinate_matrix<SCALARTYPE, ALIGNMENT>::init(); }
      
      /** @brief Construction of a coordinate matrix with the supplied number of rows and columns. If the number of nonzeros is positive, memory is allocated
      *
      * @param rows     Number of rows
      * @param cols     Number of columns
      * @param nonzeros Optional number of nonzeros for memory preallocation
      */
      explicit coordinate_matrix(unsigned int rows, unsigned int cols, unsigned int nonzeros = 0) : 
        _rows(rows), _cols(cols), _nonzeros(nonzeros), _internal_nonzeros(viennacl::tools::roundUpToNextMultiple<unsigned int>(nonzeros, ALIGNMENT))
      {
        viennacl::linalg::kernels::coordinate_matrix<SCALARTYPE, ALIGNMENT>::init();
        if (nonzeros > 0)
        {
          _coord_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * 2 * _internal_nonzeros);
          _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * _internal_nonzeros);
        }
      }
        
      /** @brief Allocate memory for the supplied number of nonzeros in the matrix. Old values are preserved. */
      void reserve(unsigned int new_nonzeros)
      {
        if (new_nonzeros > _nonzeros)
        {
          _internal_nonzeros = viennacl::tools::roundUpToNextMultiple<unsigned int>(new_nonzeros, ALIGNMENT);
          viennacl::ocl::handle<cl_mem> _coord_buffer_old = _coord_buffer;
          viennacl::ocl::handle<cl_mem> _elements_old = _elements;
          _coord_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(unsigned int) * 2 * _internal_nonzeros);
          _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE) * _internal_nonzeros);
          
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _coord_buffer_old.get(), _coord_buffer.get(), 0, 0, sizeof(unsigned int) * 2 * _nonzeros, 0, NULL, NULL);
          CL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _elements_old.get(), _elements.get(), 0, 0, sizeof(SCALARTYPE)*_nonzeros, 0, NULL, NULL);
          CL_ERR_CHECK(err);

          //new memory must be padded with zeros:
          std::vector<long> temp(_internal_nonzeros - _nonzeros);
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _coord_buffer_old.get(), _coord_buffer.get(), 0, _nonzeros, sizeof(unsigned int) * 2 * temp.size(), 0, NULL, NULL);
          CL_ERR_CHECK(err);
          err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), _elements_old.get(), _elements.get(), 0, _nonzeros, sizeof(SCALARTYPE)*temp.size(), 0, NULL, NULL);
          CL_ERR_CHECK(err);
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
            //no need to do anything...
          }
          else if (new_size1 < _rows) //reduce buffer
          {
            //no need to reduce buffer
            
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
      /** @brief  Returns the number of internal nonzero entries */
      const unsigned int & internal_nnz() const { return _internal_nonzeros; }
      
      /** @brief  Returns the OpenCL handle to the (row, column) index array */
      const viennacl::ocl::handle<cl_mem> & handle12() const { return _coord_buffer; }
      /** @brief  Returns the OpenCL handle to the matrix entry array */
      const viennacl::ocl::handle<cl_mem> & handle() const { return _elements; }
      
      template <typename CPU_MATRIX, typename SCALARTYPE2, unsigned int ALIGNMENT2>
      friend void copy(const CPU_MATRIX & cpu_matrix, coordinate_matrix<SCALARTYPE2, ALIGNMENT2> & gpu_matrix );

    private:
      unsigned int _rows;
      unsigned int _cols;
      unsigned int _nonzeros;
      unsigned int _internal_nonzeros;
      viennacl::ocl::handle<cl_mem> _coord_buffer;
      viennacl::ocl::handle<cl_mem> _elements;
    };

    
    
    
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
      if ( cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0 )
      {
        gpu_matrix.resize(static_cast<unsigned int>(cpu_matrix.size1()), static_cast<unsigned int>(cpu_matrix.size2()), false);
        
        //determine nonzeros:
        unsigned int num_entries = 0;
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
        
        //set up matrix entries:
        gpu_matrix._internal_nonzeros = num_entries;
        if (num_entries % ALIGNMENT != 0)
          gpu_matrix._internal_nonzeros = ((num_entries / ALIGNMENT) + 1) * ALIGNMENT;

        std::vector<unsigned int> coord_buffer(2*gpu_matrix._internal_nonzeros);
        std::vector<SCALARTYPE> elements(gpu_matrix._internal_nonzeros);
        
        unsigned int data_index = 0;
        
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            coord_buffer[2*data_index] = static_cast<unsigned int>(col_it.index1());
            coord_buffer[2*data_index + 1] = static_cast<unsigned int>(col_it.index2());
            elements[data_index] = *col_it;
            ++data_index;
          }
        }
        
        gpu_matrix._coord_buffer = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, coord_buffer);
        gpu_matrix._elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, elements);
        
        gpu_matrix._nonzeros = num_entries;
      }
    }

    /** @brief Copies a sparse matrix in the std::vector< std::map < > > format to an OpenCL device.
    *
    * @param cpu_matrix   A sparse matrix on the host.
    * @param gpu_matrix   A coordinate_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const std::vector< std::map<unsigned int, SCALARTYPE> > & cpu_matrix,
                     coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
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
    * @param gpu_matrix   A coordinate_matrix from ViennaCL
    * @param cpu_matrix   A sparse matrix on the host.
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const coordinate_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
                     CPU_MATRIX & cpu_matrix )
    {
      if ( cpu_matrix.size1() > 0 && cpu_matrix.size2() > 0 )
      {
        cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2());
        
        //get raw data from memory:
        std::vector<unsigned int> coord_buffer(2*gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());
        
        //std::cout << "GPU nonzeros: " << gpu_matrix.nnz() << std::endl;
        
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle12().get(), CL_TRUE, 0, sizeof(unsigned int)* 2 *gpu_matrix.nnz(), &(coord_buffer[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.nnz(), &(elements[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        viennacl::ocl::finish();
        
        //fill the cpu_matrix:
        for (unsigned int index = 0; index < gpu_matrix.nnz(); ++index)
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
      tools::sparse_matrix_adapter<SCALARTYPE> temp(cpu_matrix);
      copy(gpu_matrix, temp);
    }

}

#endif
