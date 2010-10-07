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

#ifndef _VIENNACL_MATRIX_HPP_
#define _VIENNACL_MATRIX_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/ocl/kernel.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/matrix_operations.hpp"

namespace viennacl
{
    /** @brief A tag for row-major storage of a dense matrix. */
    struct row_major
    {
      unsigned int operator()(unsigned int i, unsigned int j, unsigned int num_rows, unsigned int num_cols)
      {
        return i * num_rows + j;
      }
    };

    
    /* not provided yet 
    struct column_major
    {
      unsigned int operator()(unsigned int i, unsigned int j, unsigned int num_rows, unsigned int num_cols)
      {
        return i + j * num_cols;
      }
    }; */
    
    /** @brief A proxy class for outer products of two vectors */
    template <typename SCALARTYPE, unsigned int A>
    class outer_prod_proxy
    {
      public:
        outer_prod_proxy( vector<SCALARTYPE, A> const & vec1, 
                          vector<SCALARTYPE, A> const & vec2) : _vector1(vec1), _vector2(vec2) {}
        
        vector<SCALARTYPE, A> const & get_vector1() const { return _vector1; }
        vector<SCALARTYPE, A> const & get_vector2() const { return _vector2; }
      
      private:
        vector<SCALARTYPE, A> const & _vector1;
        vector<SCALARTYPE, A> const & _vector2;
    };

    /** @brief A proxy class for scaled outer products of two vectors */
    template <typename SCALARTYPE, unsigned int A>
    class scalar_times_outer_prod_proxy
    {
      public:
        scalar_times_outer_prod_proxy( SCALARTYPE val,
                                       vector<SCALARTYPE, A> const & vec1, 
                                       vector<SCALARTYPE, A> const & vec2) : _val(val), _vector1(vec1), _vector2(vec2) {}
        
        SCALARTYPE get_scalar() const { return _val; }
        vector<SCALARTYPE, A> const & get_vector1() const { return _vector1; }
        vector<SCALARTYPE, A> const & get_vector2() const { return _vector2; }
      
      private:
        SCALARTYPE _val;
        vector<SCALARTYPE, A> const & _vector1;
        vector<SCALARTYPE, A> const & _vector2;
    };
    

    /** @brief A proxy class for a transposed matrix. Allows several operations without rearranging the underlying values in the GPU memory. */
    template <typename SCALARTYPE, typename F, unsigned int A>
    class transposed_matrix_proxy
    {
      public:
        transposed_matrix_proxy(matrix<SCALARTYPE, F, A> const & mat) : _matrix(mat) {}
        
        matrix<SCALARTYPE, F, A> const & get_matrix() const { return _matrix; }
        
        /** @brief The number of rows of the transposed matrix. Equals the number of columns in the original matrix */
        unsigned int size1() const { return _matrix.size2(); }
        /** @brief The number of columns of the transposed matrix. Equals the number of rows in the original matrix */
        unsigned int size2() const { return _matrix.size1(); }
      
      private:
        matrix<SCALARTYPE, F, A> const & _matrix;
    };

    /** @brief A tag indicating iteration along increasing row index of a matrix */
    struct row_iteration {};
    
    /** @brief A tag indicating iteration along increasing columns index of a matrix */
    struct col_iteration {};

    //STL-like iterator. TODO: STL-compliance...
    template <typename ROWCOL, typename MATRIXTYPE>
    class matrix_iterator
    {
        typedef matrix_iterator<ROWCOL, MATRIXTYPE>    self_type;
      public:
        typedef typename MATRIXTYPE::value_type       value_type;
        
        matrix_iterator(MATRIXTYPE & mat, unsigned int start_row, unsigned int start_col) : _mat(mat), _row(start_row), _col(start_col) {};
        
        value_type operator*(void) { return _mat(_row, _col); }
        self_type & operator++(void) { viennacl::tools::MATRIX_ITERATOR_INCREMENTER<ROWCOL, MATRIXTYPE>::apply(_mat, _row, _col); return *this; }
        self_type & operator++(int) { self_type tmp = *this; ++(*this); return tmp; }
        
        bool operator==(self_type const & other) { return (_row == other._row) && (_col == other._col); }
        bool operator!=(self_type const & other) { return !(*this == other); }
        
        unsigned int index1() { return _row; }
        unsigned int index2() { return _col; }
        
        MATRIXTYPE & operator()(void) const { return _mat; }
      
      private:
        MATRIXTYPE & _mat;
        unsigned int _row;
        unsigned int _col;
    };

    /** @brief A dense matrix class
    *
    * @tparam SCALARTYPE   The underlying scalar type (either float or double)
    * @tparam F            Storage layout: Either row_major of column_major (at present only row_major is possible)
    * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    class matrix
    {
      
    public:
      
      typedef matrix_iterator<row_iteration, matrix<SCALARTYPE, F, ALIGNMENT> >   iterator1;
      typedef matrix_iterator<col_iteration, matrix<SCALARTYPE, F, ALIGNMENT> >   iterator2;
      
      /** @brief The default constructor. Does not allocate any memory. */
      matrix() : _rows(0), _columns(0) { viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::init();  };
      
      /** @brief Creates the matrix with the given dimensions
      *
      * @param rows     Number of rows
      * @param columns  Number of columns
      */
      explicit matrix(unsigned int rows, unsigned int columns) :
        _rows(rows), _columns(columns), _internal_size(internal_size1() * internal_size2())
      {
        viennacl::linalg::kernels::matrix<SCALARTYPE, ALIGNMENT>::init(); 
        _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
      }

      //copy constructor:
      matrix(const matrix<SCALARTYPE, F, ALIGNMENT> & mat) :
        _rows(mat.rows()), _columns(mat.columns()), _internal_size(mat.internal_size()),
        _elements(viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size()))
      {
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), mat.handle().get(), handle().get(), 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
        CL_ERR_CHECK(err);
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator=(const matrix<SCALARTYPE, F, ALIGNMENT> & mat)
      {
        resize(mat.size());
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), mat.handle().get(), handle().get(), 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        return *this;
      }
      

      //TODO: preserve entries if requested!
      /** @brief Resizes the matrix.
      *
      * @param rows       New number of rows
      * @param columns    New number of columns
      * @param preserve   If true, existing values are preserved. However, in the current version of ViennaCL entries are always discarded.
      */
      void resize(unsigned int rows, unsigned int columns, bool preserve = true)
      {
        _rows = rows;
        _columns = columns;
        if (internal_size1() * internal_size2() > internal_size())
        {
          _internal_size = internal_size1() * internal_size2();
          _elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
        }
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const outer_prod_proxy<SCALARTYPE, ALIGNMENT> & proxy) 
      {
        viennacl::linalg::rank_1_update(*this, proxy.get_vector1(), proxy.get_vector2());
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator -= (const outer_prod_proxy<SCALARTYPE, ALIGNMENT> & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, static_cast<SCALARTYPE>(-1.0), proxy.get_vector1(), proxy.get_vector2());
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const scalar_times_outer_prod_proxy<SCALARTYPE, ALIGNMENT> & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, proxy.get_scalar(), proxy.get_vector1(), proxy.get_vector2());
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator -= (const scalar_times_outer_prod_proxy<SCALARTYPE, ALIGNMENT> & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, static_cast<SCALARTYPE>(-1.0) * proxy.get_scalar(), proxy.get_vector1(), proxy.get_vector2());
        return *this;
      }


      /** @brief Returns the number of rows */
      const unsigned int & size1() const { return _rows;}
      /** @brief Returns the number of columns */
      const unsigned int & size2() const { return _columns; }
      
      /** @brief Resets all entries to zero */
      void clear()
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::init();
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::clear.setArgument(pos++, _elements);
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::clear.setArgument(pos++, _internal_size);
        
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::clear.start1D();
      }
      
      
      //const unsigned int row_stride() const { return roundUpToNextMultiple<unsigned int>(columns(), ALIGNMENT); }
      /** @brief Returns the internal number of rows. Usually required for launching OpenCL kernels only */
      const unsigned int internal_size1() const { return _rows; }
      /** @brief Returns the internal number of columns. Usually required for launching OpenCL kernels only */
      const unsigned int internal_size2() const { return viennacl::tools::roundUpToNextMultiple<unsigned int>(size2(), ALIGNMENT); }
      /** @brief Returns the total amount of allocated memory in multiples of sizeof(SCALARTYPE) */
      const unsigned int internal_size() const { return internal_size1() * internal_size2(); }
      
      /** @brief Returns the OpenCL handle */
      const viennacl::ocl::handle<cl_mem> & handle() const { return _elements; }
      
      template <typename CPU_MATRIX, typename SCALARTYPE2, typename F2, unsigned int ALIGNMENT2>
      friend void copy(const CPU_MATRIX & cpu_matrix, matrix<SCALARTYPE2, F2, ALIGNMENT2> & gpu_matrix );

    private:
      unsigned int _rows;
      unsigned int _columns;
      unsigned int _internal_size;
      viennacl::ocl::handle<cl_mem> _elements;
    };

    /*
    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT, typename CPUMatrixType>
    void copy(const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix, CPUMatrixType & cpu_matrix)
    {
      std::vector<SCALARTYPE> tmp(gpu_matrix.internal_size());
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::device().queue(), gpu_matrix.elements(), CL_TRUE, 0, sizeof(SCALARTYPE) * gpu_matrix.internal_size(), &tmp[0], 0, NULL, NULL);
      CL_ERR_CHECK(err);
      clFinish(viennacl::ocl::device().queue());
      
      cpu_matrix.resize(gpu_matrix.size1(), gpu_matrix.size2());

      for (unsigned int i = 0; i < gpu_matrix.size1(); ++i)
        for (unsigned int j = 0; j < gpu_matrix.size2(); ++j)
          cpu_matrix(i,j) = tmp[ i * gpu_matrix.internal_size2() + j ];
    }*/
    
    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream & s, const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix)
    {
      std::vector<SCALARTYPE> tmp(gpu_matrix.internal_size());
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE) * gpu_matrix.internal_size(), &tmp[0], 0, NULL, NULL);
      CL_ERR_CHECK(err);
      viennacl::ocl::finish();
      
      s << "[" << gpu_matrix.size1() << "," << gpu_matrix.size2() << "]";
      
      s << "(";
      for (unsigned int i = 0; i < gpu_matrix.size1(); ++i)
      {
        s << "(";
        for (unsigned int j = 0; j < gpu_matrix.size2(); ++j)
        {
          s << tmp[ i * gpu_matrix.internal_size2() + j ];
          if (j < gpu_matrix.size2() - 1)
            s << ",";
        }
        s << ")";
        if (i < gpu_matrix.size1() - 1)
          s << ",";
      }
      s << ")";
      return s;
    }
    
    /** @brief Returns an expression template class representing a transposed matrix */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT> trans(matrix<SCALARTYPE, F, ALIGNMENT> & mat)
    {
      return transposed_matrix_proxy<SCALARTYPE, F, ALIGNMENT>(mat);
    }
    
    
    //transfer operations:
    //gpu to cpu:
    /** @brief Copies a dense matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense memory on the host. Type requirements: .size1() returns number of rows, .size2() returns number of columns. Access to entries via operator()
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX & cpu_matrix,
              matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix )
    {
      gpu_matrix.resize(static_cast<unsigned int>(cpu_matrix.size1()), static_cast<unsigned int>(cpu_matrix.size2()), false);

      std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
      typename std::vector<SCALARTYPE>::iterator dst = data.begin();
      for (unsigned int i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (unsigned int j = 0; j < gpu_matrix.size2(); ++j) 
          (*dst++) = cpu_matrix(i,j);
        for (unsigned int j = gpu_matrix.size2(); j < gpu_matrix.internal_size2(); ++j) 
          (*dst++) = static_cast<SCALARTYPE>(0);
      }
      
      //in case that rows and columns are padded:
      for (unsigned int i = gpu_matrix.size1(); i < gpu_matrix.internal_size1(); ++i)
        for (unsigned int j = 0; j < gpu_matrix.internal_size2(); ++j)
          (*dst++) = static_cast<SCALARTYPE>(0);
      
      gpu_matrix._elements = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, data);
    }
    
    //gpu to cpu:
    /** @brief Copies a dense matrix from the OpenCL device (GPU or multi-core CPU) to the host (CPU). 
    *
    * @param gpu_matrix   A dense ViennaCL matrix
    * @param cpu_matrix   A dense memory on the host. Must have at least as many rows and columns as the gpu_matrix! Type requirement: Access to entries via operator()
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void copy(const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix,
                     CPU_MATRIX & cpu_matrix )
    {
      if ( (gpu_matrix.size1() > 0) && (gpu_matrix.size2() > 0) )
      {
        std::vector<SCALARTYPE> temp_buffer(gpu_matrix.internal_size());
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), &(temp_buffer[0]), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        viennacl::ocl::finish();
        
        //now copy entries to cpu_matrix:
        for (unsigned int i = 0; i < gpu_matrix.size1(); ++i)
          for (unsigned int j = 0; j < gpu_matrix.size2(); ++j) 
            cpu_matrix(i,j) = temp_buffer[i * gpu_matrix.internal_size2() + j];
      }
    }
   
   
    template<class SCALARTYPE,unsigned int VECTOR_ALIGNMENT>
    viennacl::scalar_times_outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT> operator*(viennacl::outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT> const & proxy, SCALARTYPE val)
    {
      return scalar_times_outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT>(val, proxy.get_vector1(), proxy.get_vector2());
    }

    template<class SCALARTYPE, unsigned int VECTOR_ALIGNMENT>
    viennacl::scalar_times_outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT> operator*(SCALARTYPE val, viennacl::outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT> const & proxy)
    {
      return scalar_times_outer_prod_proxy<SCALARTYPE, VECTOR_ALIGNMENT>(val, proxy.get_vector1(), proxy.get_vector2());
    }
    
   

} //namespace viennacl

#endif
