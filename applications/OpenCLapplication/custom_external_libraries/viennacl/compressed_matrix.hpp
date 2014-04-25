#ifndef VIENNACL_COMPRESSED_MATRIX_HPP_
#define VIENNACL_COMPRESSED_MATRIX_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file viennacl/compressed_matrix.hpp
    @brief Implementation of the compressed_matrix class
*/

#include <vector>
#include <list>
#include <map>
#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"

#include "viennacl/linalg/sparse_matrix_operations.hpp"

#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/entry_proxy.hpp"

namespace viennacl
{
    namespace detail
    {
      template <typename CPU_MATRIX, typename SCALARTYPE, unsigned int ALIGNMENT>
      void copy_impl(const CPU_MATRIX & cpu_matrix,
                     compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
                     vcl_size_t nonzeros)
      {
        assert( (gpu_matrix.size1() == 0 || viennacl::traits::size1(cpu_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
        assert( (gpu_matrix.size2() == 0 || viennacl::traits::size2(cpu_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), cpu_matrix.size1() + 1);
        viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), nonzeros);
        std::vector<SCALARTYPE> elements(nonzeros);

        vcl_size_t row_index  = 0;
        vcl_size_t data_index = 0;

        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          row_buffer.set(row_index, data_index);
          ++row_index;

          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            col_buffer.set(data_index, col_it.index2());
            elements[data_index] = *col_it;
            ++data_index;
          }
          data_index = viennacl::tools::align_to_multiple<vcl_size_t>(data_index, ALIGNMENT); //take care of alignment
        }
        row_buffer.set(row_index, data_index);

        gpu_matrix.set(row_buffer.get(),
                       col_buffer.get(),
                       &elements[0],
                       cpu_matrix.size1(),
                       cpu_matrix.size2(),
                       nonzeros);
      }
    }

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
        //determine nonzeros:
        vcl_size_t num_entries = 0;
        for (typename CPU_MATRIX::const_iterator1 row_it = cpu_matrix.begin1();
              row_it != cpu_matrix.end1();
              ++row_it)
        {
          vcl_size_t entries_per_row = 0;
          for (typename CPU_MATRIX::const_iterator2 col_it = row_it.begin();
                col_it != row_it.end();
                ++col_it)
          {
            ++entries_per_row;
          }
          num_entries += viennacl::tools::align_to_multiple<vcl_size_t>(entries_per_row, ALIGNMENT);
        }

        if (num_entries == 0) //we copy an empty matrix
          num_entries = 1;

        //set up matrix entries:
        viennacl::detail::copy_impl(cpu_matrix, gpu_matrix, num_entries);
      }
    }


    //adapted for std::vector< std::map < > > argument:
    /** @brief Copies a sparse square matrix in the std::vector< std::map < > > format to an OpenCL device. Use viennacl::tools::sparse_matrix_adapter for non-square matrices.
    *
    * @param cpu_matrix   A sparse square matrix on the host using STL types
    * @param gpu_matrix   A compressed_matrix from ViennaCL
    */
    template <typename SizeType, typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const std::vector< std::map<SizeType, SCALARTYPE> > & cpu_matrix,
              compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix )
    {
      vcl_size_t nonzeros = 0;
      vcl_size_t max_col = 0;
      for (vcl_size_t i=0; i<cpu_matrix.size(); ++i)
      {
        if (cpu_matrix[i].size() > 0)
        nonzeros += ((cpu_matrix[i].size() - 1) / ALIGNMENT + 1) * ALIGNMENT;
        if (cpu_matrix[i].size() > 0)
          max_col = std::max<vcl_size_t>(max_col, (cpu_matrix[i].rbegin())->first);
      }

      viennacl::detail::copy_impl(tools::const_sparse_matrix_adapter<SCALARTYPE, SizeType>(cpu_matrix, cpu_matrix.size(), max_col + 1),
                                  gpu_matrix,
                                  nonzeros);
    }

#ifdef VIENNACL_WITH_UBLAS
    template <typename ScalarType, typename F, vcl_size_t IB, typename IA, typename TA>
    void copy(const boost::numeric::ublas::compressed_matrix<ScalarType, F, IB, IA, TA> & ublas_matrix,
              viennacl::compressed_matrix<ScalarType, 1> & gpu_matrix)
    {
      assert( (gpu_matrix.size1() == 0 || viennacl::traits::size1(ublas_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (gpu_matrix.size2() == 0 || viennacl::traits::size2(ublas_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

      //we just need to copy the CSR arrays:
      viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), ublas_matrix.size1() + 1);
      for (vcl_size_t i=0; i<=ublas_matrix.size1(); ++i)
        row_buffer.set(i, ublas_matrix.index1_data()[i]);

      viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), ublas_matrix.nnz());
      for (vcl_size_t i=0; i<ublas_matrix.nnz(); ++i)
        col_buffer.set(i, ublas_matrix.index2_data()[i]);

      gpu_matrix.set(row_buffer.get(),
                     col_buffer.get(),
                     &(ublas_matrix.value_data()[0]),
                     ublas_matrix.size1(),
                     ublas_matrix.size2(),
                     ublas_matrix.nnz());

    }
#endif

    #ifdef VIENNACL_WITH_EIGEN
    template <typename SCALARTYPE, int flags, unsigned int ALIGNMENT>
    void copy(const Eigen::SparseMatrix<SCALARTYPE, flags> & eigen_matrix,
              compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix)
    {
      assert( (gpu_matrix.size1() == 0 || static_cast<vcl_size_t>(eigen_matrix.rows()) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (gpu_matrix.size2() == 0 || static_cast<vcl_size_t>(eigen_matrix.cols()) == gpu_matrix.size2()) && bool("Size mismatch") );

      std::vector< std::map<unsigned int, SCALARTYPE> >  stl_matrix(eigen_matrix.rows());

      for (int k=0; k < eigen_matrix.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<SCALARTYPE, flags>::InnerIterator it(eigen_matrix, k); it; ++it)
          stl_matrix[it.row()][it.col()] = it.value();

      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(stl_matrix, eigen_matrix.rows(), eigen_matrix.cols()), gpu_matrix);
    }
#endif


#ifdef VIENNACL_WITH_MTL4
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(const mtl::compressed2D<SCALARTYPE> & cpu_matrix,
              compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix)
    {
      assert( (gpu_matrix.size1() == 0 || static_cast<vcl_size_t>(cpu_matrix.num_rows()) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (gpu_matrix.size2() == 0 || static_cast<vcl_size_t>(cpu_matrix.num_cols()) == gpu_matrix.size2()) && bool("Size mismatch") );

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

      copy(tools::const_sparse_matrix_adapter<SCALARTYPE>(stl_matrix, cpu_matrix.num_rows(), cpu_matrix.num_cols()), gpu_matrix);
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
      assert( (viennacl::traits::size1(cpu_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (viennacl::traits::size2(cpu_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        //get raw data from memory:
        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), cpu_matrix.size1() + 1);
        viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());

        //std::cout << "GPU->CPU, nonzeros: " << gpu_matrix.nnz() << std::endl;

        viennacl::backend::memory_read(gpu_matrix.handle1(), 0, row_buffer.raw_size(), row_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle2(), 0, col_buffer.raw_size(), col_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle(),  0, sizeof(SCALARTYPE)* gpu_matrix.nnz(), &(elements[0]));

        //fill the cpu_matrix:
        vcl_size_t data_index = 0;
        for (vcl_size_t row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            if (col_buffer[data_index] >= gpu_matrix.size2())
            {
              std::cerr << "ViennaCL encountered invalid data at colbuffer[" << data_index << "]: " << col_buffer[data_index] << std::endl;
              return;
            }

            if (elements[data_index] != static_cast<SCALARTYPE>(0.0))
              cpu_matrix(row-1, static_cast<vcl_size_t>(col_buffer[data_index])) = elements[data_index];
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
      tools::sparse_matrix_adapter<SCALARTYPE> temp(cpu_matrix, cpu_matrix.size(), cpu_matrix.size());
      copy(gpu_matrix, temp);
    }

#ifdef VIENNACL_WITH_UBLAS
    template <typename ScalarType, unsigned int ALIGNMENT, typename F, vcl_size_t IB, typename IA, typename TA>
    void copy(viennacl::compressed_matrix<ScalarType, ALIGNMENT> const & gpu_matrix,
              boost::numeric::ublas::compressed_matrix<ScalarType> & ublas_matrix)
    {
      assert( (viennacl::traits::size1(ublas_matrix) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (viennacl::traits::size2(ublas_matrix) == gpu_matrix.size2()) && bool("Size mismatch") );

      viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), gpu_matrix.size1() + 1);
      viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), gpu_matrix.nnz());

      viennacl::backend::memory_read(gpu_matrix.handle1(), 0, row_buffer.raw_size(), row_buffer.get());
      viennacl::backend::memory_read(gpu_matrix.handle2(), 0, col_buffer.raw_size(), col_buffer.get());

      ublas_matrix.clear();
      ublas_matrix.reserve(gpu_matrix.nnz());

      ublas_matrix.set_filled(gpu_matrix.size1() + 1, gpu_matrix.nnz());

      for (vcl_size_t i=0; i<ublas_matrix.size1() + 1; ++i)
        ublas_matrix.index1_data()[i] = row_buffer[i];

      for (vcl_size_t i=0; i<ublas_matrix.nnz(); ++i)
        ublas_matrix.index2_data()[i] = col_buffer[i];

      viennacl::backend::memory_read(gpu_matrix.handle(),  0, sizeof(ScalarType) * gpu_matrix.nnz(), &(ublas_matrix.value_data()[0]));

    }
#endif

#ifdef VIENNACL_WITH_EIGEN
    template <typename SCALARTYPE, int flags, unsigned int ALIGNMENT>
    void copy(compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              Eigen::SparseMatrix<SCALARTYPE, flags> & eigen_matrix)
    {
      assert( (static_cast<vcl_size_t>(eigen_matrix.rows()) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (static_cast<vcl_size_t>(eigen_matrix.cols()) == gpu_matrix.size2()) && bool("Size mismatch") );

      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {
        //get raw data from memory:
        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), gpu_matrix.size1() + 1);
        viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());

        viennacl::backend::memory_read(gpu_matrix.handle1(), 0, row_buffer.raw_size(), row_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle2(), 0, col_buffer.raw_size(), col_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle(),  0, sizeof(SCALARTYPE)* gpu_matrix.nnz(),        &(elements[0]));

        eigen_matrix.setZero();
        vcl_size_t data_index = 0;
        for (vcl_size_t row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            assert(col_buffer[data_index] < gpu_matrix.size2() && bool("ViennaCL encountered invalid data at col_buffer"));
            if (elements[data_index] != static_cast<SCALARTYPE>(0.0))
              eigen_matrix.insert(row-1, col_buffer[data_index]) = elements[data_index];
            ++data_index;
          }
        }
      }
    }
#endif



#ifdef VIENNACL_WITH_MTL4
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(compressed_matrix<SCALARTYPE, ALIGNMENT> & gpu_matrix,
              mtl::compressed2D<SCALARTYPE> & mtl4_matrix)
    {
      assert( (static_cast<vcl_size_t>(mtl4_matrix.num_rows()) == gpu_matrix.size1()) && bool("Size mismatch") );
      assert( (static_cast<vcl_size_t>(mtl4_matrix.num_cols()) == gpu_matrix.size2()) && bool("Size mismatch") );

      if ( gpu_matrix.size1() > 0 && gpu_matrix.size2() > 0 )
      {

        //get raw data from memory:
        viennacl::backend::typesafe_host_array<unsigned int> row_buffer(gpu_matrix.handle1(), gpu_matrix.size1() + 1);
        viennacl::backend::typesafe_host_array<unsigned int> col_buffer(gpu_matrix.handle2(), gpu_matrix.nnz());
        std::vector<SCALARTYPE> elements(gpu_matrix.nnz());

        viennacl::backend::memory_read(gpu_matrix.handle1(), 0, row_buffer.raw_size(), row_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle2(), 0, col_buffer.raw_size(), col_buffer.get());
        viennacl::backend::memory_read(gpu_matrix.handle(),  0, sizeof(SCALARTYPE)* gpu_matrix.nnz(), &(elements[0]));

        //set_to_zero(mtl4_matrix);
        //mtl4_matrix.change_dim(gpu_matrix.size1(), gpu_matrix.size2());

        mtl::matrix::inserter< mtl::compressed2D<SCALARTYPE> >  ins(mtl4_matrix);
        vcl_size_t data_index = 0;
        for (vcl_size_t row = 1; row <= gpu_matrix.size1(); ++row)
        {
          while (data_index < row_buffer[row])
          {
            assert(col_buffer[data_index] < gpu_matrix.size2() && bool("ViennaCL encountered invalid data at col_buffer"));
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
        typedef viennacl::backend::mem_handle                                                              handle_type;
        typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
        typedef vcl_size_t                                                                                 size_type;

        /** @brief Default construction of a compressed matrix. No memory is allocated */
        compressed_matrix() : rows_(0), cols_(0), nonzeros_(0) {}

        /** @brief Construction of a compressed matrix with the supplied number of rows and columns. If the number of nonzeros is positive, memory is allocated
        *
        * @param rows     Number of rows
        * @param cols     Number of columns
        * @param nonzeros Optional number of nonzeros for memory preallocation
        * @param ctx      Optional context in which the matrix is created (one out of multiple OpenCL contexts, CUDA, host)
        */
        explicit compressed_matrix(vcl_size_t rows, vcl_size_t cols, vcl_size_t nonzeros = 0, viennacl::context ctx = viennacl::context())
          : rows_(rows), cols_(cols), nonzeros_(nonzeros)
        {
          row_buffer_.switch_active_handle_id(ctx.memory_type());
          col_buffer_.switch_active_handle_id(ctx.memory_type());
            elements_.switch_active_handle_id(ctx.memory_type());

#ifdef VIENNACL_WITH_OPENCL
          if (ctx.memory_type() == OPENCL_MEMORY)
          {
            row_buffer_.opencl_handle().context(ctx.opencl_context());
            col_buffer_.opencl_handle().context(ctx.opencl_context());
              elements_.opencl_handle().context(ctx.opencl_context());
          }
#endif
          if (rows > 0)
          {
            viennacl::backend::memory_create(row_buffer_, viennacl::backend::typesafe_host_array<unsigned int>().element_size() * (rows + 1), ctx);
          }
          if (nonzeros > 0)
          {
            viennacl::backend::memory_create(col_buffer_, viennacl::backend::typesafe_host_array<unsigned int>().element_size() * nonzeros, ctx);
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE) * nonzeros, ctx);
          }
        }

        /** @brief Construction of a compressed matrix with the supplied number of rows and columns. If the number of nonzeros is positive, memory is allocated
        *
        * @param rows     Number of rows
        * @param cols     Number of columns
        * @param ctx      Context in which to create the matrix
        */
        explicit compressed_matrix(vcl_size_t rows, vcl_size_t cols, viennacl::context ctx)
          : rows_(rows), cols_(cols), nonzeros_(0)
        {
          row_buffer_.switch_active_handle_id(ctx.memory_type());
          col_buffer_.switch_active_handle_id(ctx.memory_type());
            elements_.switch_active_handle_id(ctx.memory_type());

#ifdef VIENNACL_WITH_OPENCL
          if (ctx.memory_type() == OPENCL_MEMORY)
          {
            row_buffer_.opencl_handle().context(ctx.opencl_context());
            col_buffer_.opencl_handle().context(ctx.opencl_context());
              elements_.opencl_handle().context(ctx.opencl_context());
          }
#endif
          if (rows > 0)
          {
            viennacl::backend::memory_create(row_buffer_, viennacl::backend::typesafe_host_array<unsigned int>().element_size() * (rows + 1), ctx);
          }
        }

        explicit compressed_matrix(viennacl::context ctx) : rows_(0), cols_(0), nonzeros_(0)
        {
          row_buffer_.switch_active_handle_id(ctx.memory_type());
          col_buffer_.switch_active_handle_id(ctx.memory_type());
            elements_.switch_active_handle_id(ctx.memory_type());

#ifdef VIENNACL_WITH_OPENCL
          if (ctx.memory_type() == OPENCL_MEMORY)
          {
            row_buffer_.opencl_handle().context(ctx.opencl_context());
            col_buffer_.opencl_handle().context(ctx.opencl_context());
              elements_.opencl_handle().context(ctx.opencl_context());
          }
#endif
        }


#ifdef VIENNACL_WITH_OPENCL
        explicit compressed_matrix(cl_mem mem_row_buffer, cl_mem mem_col_buffer, cl_mem mem_elements,
                                  vcl_size_t rows, vcl_size_t cols, vcl_size_t nonzeros) :
          rows_(rows), cols_(cols), nonzeros_(nonzeros)
        {
            row_buffer_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
            row_buffer_.opencl_handle() = mem_row_buffer;
            row_buffer_.opencl_handle().inc();             //prevents that the user-provided memory is deleted once the matrix object is destroyed.
            row_buffer_.raw_size(sizeof(cl_uint) * (rows + 1));

            col_buffer_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
            col_buffer_.opencl_handle() = mem_col_buffer;
            col_buffer_.opencl_handle().inc();             //prevents that the user-provided memory is deleted once the matrix object is destroyed.
            col_buffer_.raw_size(sizeof(cl_uint) * nonzeros);

            elements_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
            elements_.opencl_handle() = mem_elements;
            elements_.opencl_handle().inc();               //prevents that the user-provided memory is deleted once the matrix object is destroyed.
            elements_.raw_size(sizeof(SCALARTYPE) * nonzeros);
        }
#endif


        /** @brief Assignment a compressed matrix from possibly another memory domain. */
        compressed_matrix & operator=(compressed_matrix const & other)
        {
          assert( (rows_ == 0 || rows_ == other.size1()) && bool("Size mismatch") );
          assert( (cols_ == 0 || cols_ == other.size2()) && bool("Size mismatch") );

          rows_ = other.size1();
          cols_ = other.size2();
          nonzeros_ = other.nnz();

          viennacl::backend::typesafe_memory_copy<unsigned int>(other.row_buffer_, row_buffer_);
          viennacl::backend::typesafe_memory_copy<unsigned int>(other.col_buffer_, col_buffer_);
          viennacl::backend::typesafe_memory_copy<SCALARTYPE>(other.elements_, elements_);

          return *this;
        }


        /** @brief Sets the row, column and value arrays of the compressed matrix
        *
        * @param row_jumper     Pointer to an array holding the indices of the first element of each row (starting with zero). E.g. row_jumper[10] returns the index of the first entry of the 11th row. The array length is 'cols + 1'
        * @param col_buffer     Pointer to an array holding the column index of each entry. The array length is 'nonzeros'
        * @param elements       Pointer to an array holding the entries of the sparse matrix. The array length is 'elements'
        * @param rows           Number of rows of the sparse matrix
        * @param cols           Number of columns of the sparse matrix
        * @param nonzeros       Number of nonzeros
        */
        void set(const void * row_jumper,
                 const void * col_buffer,
                 const SCALARTYPE * elements,
                 vcl_size_t rows,
                 vcl_size_t cols,
                 vcl_size_t nonzeros)
        {
          assert( (rows > 0)     && bool("Error in compressed_matrix::set(): Number of rows must be larger than zero!"));
          assert( (cols > 0)     && bool("Error in compressed_matrix::set(): Number of columns must be larger than zero!"));
          assert( (nonzeros > 0) && bool("Error in compressed_matrix::set(): Number of nonzeros must be larger than zero!"));
          //std::cout << "Setting memory: " << cols + 1 << ", " << nonzeros << std::endl;

          //row_buffer_.switch_active_handle_id(viennacl::backend::OPENCL_MEMORY);
          viennacl::backend::memory_create(row_buffer_, viennacl::backend::typesafe_host_array<unsigned int>(row_buffer_).element_size() * (rows + 1), viennacl::traits::context(row_buffer_), row_jumper);

          //col_buffer_.switch_active_handle_id(viennacl::backend::OPENCL_MEMORY);
          viennacl::backend::memory_create(col_buffer_, viennacl::backend::typesafe_host_array<unsigned int>(col_buffer_).element_size() * nonzeros, viennacl::traits::context(col_buffer_), col_buffer);

          //elements_.switch_active_handle_id(viennacl::backend::OPENCL_MEMORY);
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE) * nonzeros, viennacl::traits::context(elements_), elements);

          nonzeros_ = nonzeros;
          rows_ = rows;
          cols_ = cols;
        }

        /** @brief Allocate memory for the supplied number of nonzeros in the matrix. Old values are preserved. */
        void reserve(vcl_size_t new_nonzeros)
        {
          if (new_nonzeros > nonzeros_)
          {
            handle_type col_buffer_old;
            handle_type elements_old;
            viennacl::backend::memory_shallow_copy(col_buffer_, col_buffer_old);
            viennacl::backend::memory_shallow_copy(elements_,   elements_old);

            viennacl::backend::typesafe_host_array<unsigned int> size_deducer(col_buffer_);
            viennacl::backend::memory_create(col_buffer_, size_deducer.element_size() * new_nonzeros, viennacl::traits::context(col_buffer_));
            viennacl::backend::memory_create(elements_,   sizeof(SCALARTYPE) * new_nonzeros,          viennacl::traits::context(elements_));

            viennacl::backend::memory_copy(col_buffer_old, col_buffer_, 0, 0, size_deducer.element_size() * nonzeros_);
            viennacl::backend::memory_copy(elements_old,   elements_,   0, 0, sizeof(SCALARTYPE)* nonzeros_);

            nonzeros_ = new_nonzeros;
          }
        }

        /** @brief Resize the matrix.
        *
        * @param new_size1    New number of rows
        * @param new_size2    New number of columns
        * @param preserve     If true, the old values are preserved. At present, old values are always discarded.
        */
        void resize(vcl_size_t new_size1, vcl_size_t new_size2, bool preserve = true)
        {
          assert(new_size1 > 0 && new_size2 > 0 && bool("Cannot resize to zero size!"));

          if (new_size1 != rows_ || new_size2 != cols_)
          {
            std::vector<std::map<unsigned int, SCALARTYPE> > stl_sparse_matrix;
            if (rows_ > 0)
            {
              if (preserve)
              {
                stl_sparse_matrix.resize(rows_);
                viennacl::copy(*this, stl_sparse_matrix);
              } else
                stl_sparse_matrix[0][0] = 0;
            } else {
              stl_sparse_matrix.resize(new_size1);
              stl_sparse_matrix[0][0] = 0;      //enforces nonzero array sizes if matrix was initially empty
            }

            stl_sparse_matrix.resize(new_size1);

            //discard entries with column index larger than new_size2
            if (new_size2 < cols_ && rows_ > 0)
            {
              for (vcl_size_t i=0; i<stl_sparse_matrix.size(); ++i)
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

            viennacl::copy(stl_sparse_matrix, *this);

            rows_ = new_size1;
            cols_ = new_size2;
          }
        }

        /** @brief Returns a reference to the (i,j)-th entry of the sparse matrix. If (i,j) does not exist (zero), it is inserted (slow!) */
        entry_proxy<SCALARTYPE> operator()(vcl_size_t i, vcl_size_t j)
        {
          assert( (i < rows_) && (j < cols_) && bool("compressed_matrix access out of bounds!"));

          vcl_size_t index = element_index(i, j);

          // check for element in sparsity pattern
          if (index < nonzeros_)
            return entry_proxy<SCALARTYPE>(index, elements_);

          // Element not found. Copying required. Very slow, but direct entry manipulation is painful anyway...
          std::vector< std::map<unsigned int, SCALARTYPE> > cpu_backup(rows_);
          tools::sparse_matrix_adapter<SCALARTYPE> adapted_cpu_backup(cpu_backup, rows_, cols_);
          viennacl::copy(*this, adapted_cpu_backup);
          cpu_backup[i][static_cast<unsigned int>(j)] = 0.0;
          viennacl::copy(adapted_cpu_backup, *this);

          index = element_index(i, j);

          assert(index < nonzeros_);

          return entry_proxy<SCALARTYPE>(index, elements_);
        }

        /** @brief  Returns the number of rows */
        const vcl_size_t & size1() const { return rows_; }
        /** @brief  Returns the number of columns */
        const vcl_size_t & size2() const { return cols_; }
        /** @brief  Returns the number of nonzero entries */
        const vcl_size_t & nnz() const { return nonzeros_; }

        /** @brief  Returns the OpenCL handle to the row index array */
        const handle_type & handle1() const { return row_buffer_; }
        /** @brief  Returns the OpenCL handle to the column index array */
        const handle_type & handle2() const { return col_buffer_; }
        /** @brief  Returns the OpenCL handle to the matrix entry array */
        const handle_type & handle() const { return elements_; }

        /** @brief  Returns the OpenCL handle to the row index array */
        handle_type & handle1() { return row_buffer_; }
        /** @brief  Returns the OpenCL handle to the column index array */
        handle_type & handle2() { return col_buffer_; }
        /** @brief  Returns the OpenCL handle to the matrix entry array */
        handle_type & handle() { return elements_; }

        void switch_memory_context(viennacl::context new_ctx)
        {
          viennacl::backend::switch_memory_context<unsigned int>(row_buffer_, new_ctx);
          viennacl::backend::switch_memory_context<unsigned int>(col_buffer_, new_ctx);
          viennacl::backend::switch_memory_context<SCALARTYPE>(elements_, new_ctx);
        }

        viennacl::memory_types memory_context() const
        {
          return row_buffer_.get_active_handle_id();
        }

      private:

        vcl_size_t element_index(vcl_size_t i, vcl_size_t j)
        {
          //read row indices
          viennacl::backend::typesafe_host_array<unsigned int> row_indices(row_buffer_, 2);
          viennacl::backend::memory_read(row_buffer_, row_indices.element_size()*i, row_indices.element_size()*2, row_indices.get());

          //get column indices for row i:
          viennacl::backend::typesafe_host_array<unsigned int> col_indices(col_buffer_, row_indices[1] - row_indices[0]);
          viennacl::backend::memory_read(col_buffer_, col_indices.element_size()*row_indices[0], row_indices.element_size()*col_indices.size(), col_indices.get());

          //get entries for row i:
          viennacl::backend::typesafe_host_array<SCALARTYPE> row_entries(elements_, row_indices[1] - row_indices[0]);
          viennacl::backend::memory_read(elements_, sizeof(SCALARTYPE)*row_indices[0], sizeof(SCALARTYPE)*row_entries.size(), row_entries.get());

          for (vcl_size_t k=0; k<col_indices.size(); ++k)
          {
            if (col_indices[k] == j)
              return row_indices[0] + k;
          }

          // if not found, return index past the end of the matrix (cf. matrix.end() in the spirit of the STL)
          return nonzeros_;
        }

        // /** @brief Copy constructor is by now not available. */
        //compressed_matrix(compressed_matrix const &);


        vcl_size_t rows_;
        vcl_size_t cols_;
        vcl_size_t nonzeros_;
        handle_type row_buffer_;
        handle_type col_buffer_;
        handle_type elements_;
    };



    //
    // Specify available operations:
    //

    /** \cond */

    namespace linalg
    {
      namespace detail
      {
        // x = A * y
        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              // check for the special case x = A * x
              if (viennacl::traits::handle(lhs) == viennacl::traits::handle(rhs.rhs()))
              {
                viennacl::vector<T> temp(lhs);
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
                lhs = temp;
              }
              else
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs);
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs += temp;
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs -= temp;
            }
        };


        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const compressed_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::linalg::prod_impl(rhs.lhs(), temp, lhs);
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const compressed_matrix<T, A>, vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs += temp_result;
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const compressed_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const compressed_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs(), viennacl::traits::context(rhs));
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs -= temp_result;
            }
        };

     } // namespace detail
   } // namespace linalg

   /** \endcond */
}

#endif
