#ifndef VIENNACL_MATRIX_HPP_
#define VIENNACL_MATRIX_HPP_

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

/** @file matrix.hpp
    @brief Implementation of the dense matrix class
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/matrix_size_deducer.hpp"
#include "viennacl/tools/matrix_kernel_class_deducer.hpp"
#include "viennacl/meta/result_of.hpp"
#include "viennacl/meta/enable_if.hpp"

namespace viennacl
{
    /** @brief A tag for row-major storage of a dense matrix. */
    struct row_major
    {
      /** @brief Returns the memory offset for entry (i,j) of a dense matrix.
      *
      * @param i   row index
      * @param j   column index
      * @param num_rows  number of entries per row (including alignment)
      * @param num_cols  number of entries per column (including alignment)
      */
      static vcl_size_t mem_index(vcl_size_t i, vcl_size_t j, vcl_size_t num_rows, vcl_size_t num_cols)
      {
        return i * num_cols + j;
      }
      
      static vcl_size_t internal_size1(vcl_size_t rows, vcl_size_t alignment)
      {
        return viennacl::tools::roundUpToNextMultiple<vcl_size_t>(rows, alignment);;
      }
      
      static vcl_size_t internal_size2(vcl_size_t cols, vcl_size_t alignment)
      {
        return viennacl::tools::roundUpToNextMultiple<vcl_size_t>(cols, alignment);
      }
    };

    struct column_major
    {
      /** @brief Returns the memory offset for entry (i,j) of a dense matrix.
      *
      * @param i   row index
      * @param j   column index
      * @param num_rows  number of entries per row (including alignment)
      * @param num_cols  number of entries per column (including alignment)
      */
      static vcl_size_t mem_index(vcl_size_t i, vcl_size_t j, vcl_size_t num_rows, vcl_size_t num_cols)
      {
        return i + j * num_rows;
      }
      
      static vcl_size_t internal_size1(vcl_size_t rows, vcl_size_t alignment)
      {
        return viennacl::tools::roundUpToNextMultiple<vcl_size_t>(rows, alignment);
      }
      
      static vcl_size_t internal_size2(vcl_size_t cols, vcl_size_t alignment)
      {
        return viennacl::tools::roundUpToNextMultiple<vcl_size_t>(cols, alignment);
      }
    };
    
    template <typename LHS, typename RHS, typename OP>
    class matrix_expression
    {
      public:
        ///** @brief Extracts the vector type from the two operands.
        //*/
        //typedef typename viennacl::tools::VECTOR_EXTRACTOR<LHS, RHS>::ResultType    VectorType;
      
        matrix_expression(LHS & lhs, RHS & rhs) : _lhs(lhs), _rhs(rhs) {}
        
        /** @brief Get left hand side operand
        */
        LHS & lhs() const { return _lhs; }
        /** @brief Get right hand side operand
        */
        RHS & rhs() const { return _rhs; }
        
        /** @brief Returns the size of the result vector */
        std::size_t size1() const { return viennacl::tools::MATRIX_SIZE_DEDUCER<LHS, RHS, OP>::size1(_lhs, _rhs); }
        std::size_t size2() const { return viennacl::tools::MATRIX_SIZE_DEDUCER<LHS, RHS, OP>::size2(_lhs, _rhs); }
        
      private:
        /** @brief The left hand side operand */
        typename result_of::matrix_expression_internal_storage<LHS>::type _lhs;
        /** @brief The right hand side operand */
        typename result_of::matrix_expression_internal_storage<RHS>::type _rhs;
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
        
        matrix_iterator(MATRIXTYPE & mat, 
                        std::size_t start_row,
                        std::size_t start_col) : mat_(mat), row_(start_row), col_(start_col) {};
        
        value_type operator*(void) { return mat_(row_, col_); }
        self_type & operator++(void) { viennacl::tools::MATRIX_ITERATOR_INCREMENTER<ROWCOL, MATRIXTYPE>::apply(mat_, row_, col_); return *this; }
        self_type & operator++(int) { self_type tmp = *this; ++(*this); return tmp; }
        
        bool operator==(self_type const & other) { return (row_ == other.row_) && (col_ == other.col_); }
        bool operator!=(self_type const & other) { return !(*this == other); }
        
        vcl_size_t index1() { return row_; }
        vcl_size_t index2() { return col_; }
        
        MATRIXTYPE & operator()(void) const { return mat_; }
      
      private:
        MATRIXTYPE & mat_;
        vcl_size_t row_;
        vcl_size_t col_;
    };

    /** @brief A dense matrix class
    *
    * @tparam SCALARTYPE   The underlying scalar type (either float or double)
    * @tparam F            Storage layout: Either row_major or column_major (at present only row_major is supported)
    * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    class matrix
    {
      typedef matrix<SCALARTYPE, F, ALIGNMENT>          self_type;
    public:
      
      typedef matrix_iterator<row_iteration, matrix<SCALARTYPE, F, ALIGNMENT> >   iterator1;
      typedef matrix_iterator<col_iteration, matrix<SCALARTYPE, F, ALIGNMENT> >   iterator2;
      typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
      typedef vcl_size_t                                                          size_type;
      
      /** @brief The default constructor. Does not allocate any memory. */
      matrix() : rows_(0), columns_(0)
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
      };
      
      /** @brief Creates the matrix with the given dimensions
      *
      * @param rows     Number of rows
      * @param columns  Number of columns
      */
      explicit matrix(size_type rows, size_type columns) :
        rows_(rows), columns_(columns)
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
        elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
      }

      explicit matrix(cl_mem mem, size_type rows, size_type columns) :
        rows_(rows), columns_(columns)
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
        elements_ = mem;
        elements_.inc(); //prevents that the user-provided memory is deleted once the matrix object is destroyed.
      }

      template <typename LHS, typename RHS, typename OP>
      matrix(matrix_expression< LHS, RHS, OP> const & proxy) : rows_(proxy.size1()), columns_(proxy.size2())
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
        elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
        
        *this = proxy;
      }

      matrix(matrix_range<self_type> const & proxy) : rows_(proxy.size1()), columns_(proxy.size2())
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
        elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
        
        *this = proxy;
      }

      matrix(matrix_range<const self_type> const & proxy) : rows_(proxy.size1()), columns_(proxy.size2())
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        KernelClass::init();
        elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
        
        *this = proxy;
      }



      //copy constructor:
      matrix(const matrix<SCALARTYPE, F, ALIGNMENT> & mat) :
        rows_(mat.size1()), columns_(mat.size2()),
        elements_(viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size()))
      {
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), mat.handle().get(), elements_.get(), 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator=(const matrix<SCALARTYPE, F, ALIGNMENT> & mat)
      {
        resize(mat.size1(), mat.size2(), false);
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), mat.handle().get(), elements_.get(), 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        return *this;
      }
      
      matrix<SCALARTYPE, F, ALIGNMENT> & operator=(const matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                            const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                            op_trans> & proxy)
      {
        assert(elements_.get() != proxy.lhs().handle().get() && "Self-assignment of matrix transpose not implemented");
        assert(proxy.lhs().size1() == size2() && "Matrix dimensions do not match!");
        assert(proxy.lhs().size2() == size1() && "Matrix dimensions do not match!");

        resize(proxy.lhs().size2(), proxy.lhs().size1(), false);
        
        std::vector<SCALARTYPE> temp(proxy.lhs().internal_size());
        
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                         proxy.lhs().handle().get(), CL_TRUE, 0,
                                         sizeof(SCALARTYPE)*proxy.lhs().internal_size(),
                                         &(temp[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();

        /*
        for (size_t i=0; i<proxy.lhs().size1(); ++i)
        {
          for (size_t j=0; j<proxy.lhs().size2(); ++j)
            std::cout << temp[F::mem_index(i,j, proxy.lhs().internal_size1(), proxy.lhs().internal_size2())] << ", ";
        }*/
        
        std::vector<SCALARTYPE> temp_trans(internal_size());

        for (vcl_size_t i=0; i<proxy.lhs().size1(); ++i)
          for (vcl_size_t j=0; j<proxy.lhs().size2(); ++j)
            temp_trans[F::mem_index(j,i, internal_size1(), internal_size2())] 
             = temp[F::mem_index(i,j, proxy.lhs().internal_size1(), proxy.lhs().internal_size2())];

        /*     
        for (size_t i=0; i<proxy.lhs().size1(); ++i)
        {
          for (size_t j=0; j<proxy.lhs().size2(); ++j)
            std::cout << temp_trans[F::mem_index(i,j, proxy.lhs().internal_size1(), proxy.lhs().internal_size2())] << ", ";
        }*/
        
        elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, 
                                                                   sizeof(SCALARTYPE)*internal_size(),
                                                                   &(temp_trans[0]));
          
        return *this;
      }


      matrix<SCALARTYPE, F, ALIGNMENT> & operator=(const matrix_range<self_type> & mat)
      {
        resize(mat.size1(), mat.size2(), false);
        
        // clear matrix:
        clear();
        
        // use inplace_add:
        viennacl::linalg::inplace_add(*this, mat);
        
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator=(const matrix_range<const self_type> & mat)
      {
        resize(mat.size1(), mat.size2(), false);
        
        // clear matrix:
        clear();
        
        // use inplace_add:
        viennacl::linalg::inplace_add(*this, mat);
        
        return *this;
      }




      /** @brief Resizes the matrix.
      *   Existing entries can be preserved, but 
      *
      * @param rows       New number of rows
      * @param columns    New number of columns
      * @param preserve   If true, existing values are preserved. 
      */
      void resize(size_type rows, size_type columns, bool preserve = true)
      {
        assert(rows > 0 && columns > 0);
        if (preserve)
        {
          //get old entries:
          std::vector< SCALARTYPE > old_entries(internal_size());
          cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), //src
                                           elements_.get(), //dest
                                           CL_TRUE, //blocking
                                           0, //offset
                                           sizeof(SCALARTYPE)*internal_size(), //size
                                           &(old_entries[0]), //destination
                                           0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
          
          //set up entries of new matrix:
          std::vector< SCALARTYPE > new_entries(F::internal_size1(rows, ALIGNMENT) * F::internal_size2(columns, ALIGNMENT));
          for (size_type i=0; i<rows; ++i)
          {
            if (i >= rows_)
              continue;
              
            for (size_type j=0; j<columns; ++j)
            {
              if (j >= columns_)
                continue;
              new_entries[F::mem_index(i, j, F::internal_size1(rows, ALIGNMENT), F::internal_size2(columns, ALIGNMENT))] 
                 = old_entries[F::mem_index(i, j, internal_size1(), internal_size2())];
            }
          }
          
          //copy new entries to GPU:
          elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, new_entries);
          rows_ = rows;
          columns_ = columns;
        }
        else //discard old entries:
        {
          rows_ = rows;
          columns_ = columns;
          
          std::vector< SCALARTYPE > new_entries(F::internal_size1(rows, ALIGNMENT) * F::internal_size2(columns, ALIGNMENT));
          elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, new_entries);
        }
      }
      
      
      //read-write access to an element of the vector
      /** @brief Read-write access to a single element of the vector
      */
      entry_proxy<SCALARTYPE> operator()(size_type row_index, size_type col_index)
      {
        return entry_proxy<SCALARTYPE>(F::mem_index(row_index, col_index, internal_size1(), internal_size2()), elements_);
      }
      
      /** @brief Read access to a single element of the vector
      */
      scalar<SCALARTYPE> operator()(size_type row_index, size_type col_index) const
      {
        scalar<SCALARTYPE> tmp;
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(),
                                  elements_.get(),
                                  tmp.handle().get(),
                                  sizeof(SCALARTYPE) * F::mem_index(row_index, col_index, internal_size1(), internal_size2()),
                                  0,
                                  sizeof(SCALARTYPE),
                                  0,
                                  NULL,
                                  NULL);
        //assert(err == CL_SUCCESS);
        VIENNACL_ERR_CHECK(err);
        return tmp;
      }
      

      matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                         const matrix<SCALARTYPE, F, ALIGNMENT>,
                         op_add >
      operator + (const matrix< SCALARTYPE, F, ALIGNMENT> & other) 
      {
        return matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                  const matrix<SCALARTYPE, F, ALIGNMENT>,
                                  op_add > (*this, other);
      }


      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const matrix< SCALARTYPE, F, ALIGNMENT> & other) 
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const matrix_range< matrix<SCALARTYPE, F, ALIGNMENT> > & other) 
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }

      matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                         const matrix<SCALARTYPE, F, ALIGNMENT>,
                         op_sub >
      operator - (const matrix< SCALARTYPE, F, ALIGNMENT> & other) 
      {
        return matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                  const matrix<SCALARTYPE, F, ALIGNMENT>,
                                  op_sub > (*this, other);
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator -= (const matrix< SCALARTYPE, F, ALIGNMENT> & other) 
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }

      template <unsigned int A1, unsigned int A2>
      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const matrix_expression< const vector<SCALARTYPE, A1>,
                                                                               const vector<SCALARTYPE, A2>,
                                                                               op_prod > & proxy) 
      {
        viennacl::linalg::rank_1_update(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      template <unsigned int A1, unsigned int A2>
      matrix<SCALARTYPE, F, ALIGNMENT> & operator -= (const matrix_expression< const vector<SCALARTYPE, A1>,
                                                                               const vector<SCALARTYPE, A2>,
                                                                               op_prod > & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, static_cast<SCALARTYPE>(-1.0), proxy.lhs(), proxy.rhs());
        return *this;
      }

      template <unsigned int A1, unsigned int A2>
      matrix<SCALARTYPE, F, ALIGNMENT> & operator += (const matrix_expression< const matrix_expression< const vector<SCALARTYPE, A1>,
                                                                                                        const vector<SCALARTYPE, A2>,
                                                                                                        op_prod >,
                                                                               const SCALARTYPE,
                                                                               op_prod > & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, proxy.rhs(), proxy.lhs().lhs(), proxy.lhs().rhs());
        return *this;
      }

      template <unsigned int A1, unsigned int A2>
      matrix<SCALARTYPE, F, ALIGNMENT> & operator -= (const matrix_expression< const matrix_expression< const vector<SCALARTYPE, A1>,
                                                                                                        const vector<SCALARTYPE, A2>,
                                                                                                        op_prod >,
                                                                               const SCALARTYPE,
                                                                               op_prod > & proxy) 
      {
        viennacl::linalg::scaled_rank_1_update(*this, static_cast<SCALARTYPE>(-1.0) * proxy.rhs(), proxy.lhs().lhs(), proxy.lhs().rhs());
        return *this;
      }
      
      
      matrix<SCALARTYPE, F, ALIGNMENT> & operator *= (SCALARTYPE val) 
      {
        viennacl::linalg::inplace_mult(*this, val);
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator *= (scalar<SCALARTYPE> const & val) 
      {
        viennacl::linalg::inplace_mult(*this, val);
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator /= (SCALARTYPE val) 
      {
        viennacl::linalg::inplace_mult(*this, SCALARTYPE(1.0) / val);
        return *this;
      }

      matrix<SCALARTYPE, F, ALIGNMENT> & operator /= (scalar<SCALARTYPE> const & val) 
      {
        viennacl::linalg::inplace_divide(*this, val);
        return *this;
      }


      //this = A * B and related (with trans())
      template <typename MatrixType1, typename MatrixType2>
      matrix<SCALARTYPE, F, ALIGNMENT> & operator = (const matrix_expression< MatrixType1,
                                                                              MatrixType2,
                                                                              op_prod > & proxy) 
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      //this = A + B
      template <typename T1, typename T2>
      matrix<SCALARTYPE, F, ALIGNMENT> &
      operator = (const matrix_expression< const T1,
                                           const T2,
                                           op_add > & proxy) 
      {
        viennacl::linalg::add(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      
      //this = A - B
      template <typename T1, typename T2>
      matrix<SCALARTYPE, F, ALIGNMENT> &
      operator = (const matrix_expression< const T1,
                                           const T2,
                                           op_sub > & proxy) 
      {
        viennacl::linalg::add(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      
      
      

      //this = A - B
      matrix<SCALARTYPE, F, ALIGNMENT> & operator = (const matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                               const matrix<SCALARTYPE, F, ALIGNMENT>,
                                                                               op_sub > & proxy) 
      {
        viennacl::linalg::sub(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }


      /** @brief Returns the number of rows */
      const size_type & size1() const { return rows_;}
      /** @brief Returns the number of columns */
      const size_type & size2() const { return columns_; }
      
      /** @brief Resets all entries to zero */
      void clear()
      {
        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< matrix<SCALARTYPE, F, ALIGNMENT> >::ResultType    KernelClass;
        
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "clear");
        viennacl::ocl::enqueue(k(elements_,
                                 cl_uint(0), cl_uint(0),
                                 cl_uint(size1()), cl_uint(size2()),
                                 cl_uint(internal_size1()), cl_uint(internal_size2())
                                )
                              );
      }
      
      
      //const unsigned int row_stride() const { return roundUpToNextMultiple<unsigned int>(columns(), ALIGNMENT); }
      /** @brief Returns the internal number of rows. Usually required for launching OpenCL kernels only */
      const size_type internal_size1() const { return F::internal_size1(size1(), ALIGNMENT); }
      /** @brief Returns the internal number of columns. Usually required for launching OpenCL kernels only */
      const size_type internal_size2() const { return F::internal_size2(size2(), ALIGNMENT); }
      /** @brief Returns the total amount of allocated memory in multiples of sizeof(SCALARTYPE) */
      const size_type internal_size() const { return internal_size1() * internal_size2(); }
      
      /** @brief Returns the OpenCL handle */
      const viennacl::ocl::handle<cl_mem> & handle() const { return elements_; }
      
      #if defined(_MSC_VER) && _MSC_VER < 1500          //Visual Studio 2005 needs special treatment
      template <typename CPU_MATRIX>
      friend void copy(const CPU_MATRIX & cpu_matrix,
                      matrix & gpu_matrix );
      
      template <typename SCALARTYPE2, typename A1, typename A2>
      friend void copy(const std::vector< std::vector<SCALARTYPE2, A1>, A2> & cpu_matrix,
                      matrix & gpu_matrix );
      
      template <typename SCALARTYPE2>
      friend void fast_copy(SCALARTYPE2 * cpu_matrix_begin,
                            SCALARTYPE2 * cpu_matrix_end,
                            matrix & gpu_matrix);
      
      #ifdef VIENNACL_HAVE_EIGEN
      friend void copy(const Eigen::MatrixXf & cpu_matrix,
                       matrix & gpu_matrix);
      
      friend void copy(const Eigen::MatrixXd & cpu_matrix,
                       matrix & gpu_matrix);
      #endif
      
      #ifdef VIENNACL_HAVE_MTL4
      template <typename SCALARTYPE2, typename T>
      friend void copy(const mtl::dense2D<SCALARTYPE2, T>& cpu_matrix,
                       matrix & gpu_matrix);
      #endif
      #else
      template <typename CPU_MATRIX, typename SCALARTYPE2, typename F2, unsigned int ALIGNMENT2>
      friend void copy(const CPU_MATRIX & cpu_matrix,
                      matrix<SCALARTYPE2, F2, ALIGNMENT2> & gpu_matrix );
                      
      template <typename SCALARTYPE2, typename A1, typename A2, typename F2, unsigned int ALIGNMENT2>
      friend void copy(const std::vector< std::vector<SCALARTYPE2, A1>, A2> & cpu_matrix,
                       matrix<SCALARTYPE2, F2, ALIGNMENT2> & gpu_matrix );
      
      template <typename SCALARTYPE2, typename F2, unsigned int ALIGNMENT2>
      friend void fast_copy(SCALARTYPE2 * cpu_matrix_begin,
                            SCALARTYPE2 * cpu_matrix_end,
                            matrix<SCALARTYPE2, F2, ALIGNMENT2> & gpu_matrix);
      
      #ifdef VIENNACL_HAVE_EIGEN
      template <typename F2, unsigned int ALIGNMENT2>
      friend void copy(const Eigen::MatrixXf & cpu_matrix,
                matrix<float, F2, ALIGNMENT2> & gpu_matrix);
      
      template <typename F2, unsigned int ALIGNMENT2>
      friend void copy(const Eigen::MatrixXd & cpu_matrix,
                matrix<double, F2, ALIGNMENT2> & gpu_matrix);
      #endif
      
      #ifdef VIENNACL_HAVE_MTL4
      template <typename SCALARTYPE2, typename T, typename F2, unsigned int ALIGNMENT2>
      friend void copy(const mtl::dense2D<SCALARTYPE2, T>& cpu_matrix,
                       matrix<SCALARTYPE2, F2, ALIGNMENT2> & gpu_matrix);
      #endif
      #endif                 
      
    private:
      size_type rows_;
      size_type columns_;
      viennacl::ocl::handle<cl_mem> elements_;
    }; //matrix

    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream & s, const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix)
    {
      typedef typename matrix<SCALARTYPE, F, ALIGNMENT>::size_type      size_type;
      
      std::vector<SCALARTYPE> tmp(gpu_matrix.internal_size());
      cl_int err;
      err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE) * gpu_matrix.internal_size(), &tmp[0], 0, NULL, NULL);
      VIENNACL_ERR_CHECK(err);
      viennacl::ocl::get_queue().finish();
      
      s << "[" << gpu_matrix.size1() << "," << gpu_matrix.size2() << "]";
      
      s << "(";
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        s << "(";
        for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        {
          s << tmp[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
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

    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param expr         A matrix expression
    */
    template<typename LHS, typename RHS, typename OP>
    std::ostream & operator<<(std::ostream & s, const matrix_expression<LHS, RHS, OP> & expr)
    {
      typedef typename viennacl::tools::CPU_SCALAR_TYPE_DEDUCER< typename tools::CONST_REMOVER<LHS>::ResultType >::ResultType     ScalarType;

      matrix<ScalarType> temp = expr;
      s << temp;
      return s;
    }
    
    /** @brief Returns an expression template class representing a transposed matrix */
    template<class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                       const matrix<SCALARTYPE, F, ALIGNMENT>,
                       op_trans> trans(const matrix<SCALARTYPE, F, ALIGNMENT> & mat)
    {
      return matrix_expression< const matrix<SCALARTYPE, F, ALIGNMENT>,
                                const matrix<SCALARTYPE, F, ALIGNMENT>,
                                op_trans>(mat, mat);
    }
    
    
    /////////////////////// transfer operations: //////////////////////////////////////

    //
    //cpu to gpu, generic type:
    //
    /** @brief Copies a dense matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense matrix on the host. Type requirements: .size1() returns number of rows, .size2() returns number of columns. Access to entries via operator()
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void copy(const CPU_MATRIX & cpu_matrix,
              matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix )
    {
      typedef typename matrix<SCALARTYPE, F, ALIGNMENT>::size_type      size_type;
      
      //std::cout << "Copying CPU_MATRIX!" << std::endl;
      //std::cout << "Size at begin: " << gpu_matrix.size1() << ", " << gpu_matrix.size2() << std::endl;
      if (gpu_matrix.size1() == 0 || gpu_matrix.size2() == 0)
      {
        gpu_matrix.resize(cpu_matrix.size1(),
                          cpu_matrix.size2(), false);
      }
      else
      {
        assert( (gpu_matrix.size1() == cpu_matrix.size1()) 
               && (gpu_matrix.size2() == cpu_matrix.size2())
              );
      }

      std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
          data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
      }
      
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
      //std::cout << "Size at end: " << gpu_matrix.size1() << ", " << gpu_matrix.size2() << std::endl;
    }
    
    //
    //cpu to gpu, STL type:
    //
    /** @brief Copies a dense STL-type matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense matrix on the host of type std::vector< std::vector<> >. cpu_matrix[i][j] returns the element in the i-th row and j-th columns (both starting with zero)
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename SCALARTYPE, typename A1, typename A2, typename F, unsigned int ALIGNMENT>
    void copy(const std::vector< std::vector<SCALARTYPE, A1>, A2> & cpu_matrix,
              matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix )
    {
      typedef typename matrix<SCALARTYPE, F, ALIGNMENT>::size_type      size_type;
      
      if (gpu_matrix.size1() == 0 || gpu_matrix.size2() == 0)
      {
        gpu_matrix.resize(cpu_matrix.size(),
                          cpu_matrix[0].size(),
                          false);
      }
      else
      {
        assert( (gpu_matrix.size1() == cpu_matrix.size()) 
               && (gpu_matrix.size2() == cpu_matrix[0].size())
              );
      }

      std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
          data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix[i][j];
      }
      
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
    }
    
    
    //
    //cpu to gpu, another STL type:
    //
    /** @brief Copies a dense matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU) without temporary. Matrix-Layout on CPU must be equal to the matrix-layout on the GPU.
    *
    * @param cpu_matrix_begin   Pointer to the first matrix entry. Cf. iterator concept in STL
    * @param cpu_matrix_end     Pointer past the last matrix entry. Cf. iterator concept in STL
    * @param gpu_matrix         A dense ViennaCL matrix
    */
    template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void fast_copy(SCALARTYPE * cpu_matrix_begin,
                   SCALARTYPE * cpu_matrix_end,
                   matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix)
    {
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                            sizeof(SCALARTYPE) * (cpu_matrix_end - cpu_matrix_begin),
                                                                            cpu_matrix_begin);
    }
    
   
    #ifdef VIENNACL_HAVE_EIGEN
    /** @brief Copies a dense Eigen matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense MTL matrix. cpu_matrix(i, j) returns the element in the i-th row and j-th columns (both starting with zero)
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename F, unsigned int ALIGNMENT>
    void copy(const Eigen::MatrixXf & cpu_matrix,
              matrix<float, F, ALIGNMENT> & gpu_matrix)
    {
      typedef typename matrix<float, F, ALIGNMENT>::size_type      size_type;
      
      if (gpu_matrix.size1() == 0 || gpu_matrix.size2() == 0)
      {
        gpu_matrix.resize(cpu_matrix.rows(),
                          cpu_matrix.cols(),
                          false);
      }
      else
      {
        assert( (gpu_matrix.size1() == static_cast<std::size_t>(cpu_matrix.rows())) 
               && (gpu_matrix.size2() == static_cast<std::size_t>(cpu_matrix.cols()))
              );
      }

      std::vector<float> data(gpu_matrix.internal_size());
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
          data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
      }
      
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
    }
    
    /** @brief Copies a dense Eigen matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense MTL matrix. cpu_matrix(i, j) returns the element in the i-th row and j-th columns (both starting with zero)
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename F, unsigned int ALIGNMENT>
    void copy(const Eigen::MatrixXd & cpu_matrix,
              matrix<double, F, ALIGNMENT> & gpu_matrix)
    {
      typedef typename matrix<double, F, ALIGNMENT>::size_type      size_type;
      
      if (gpu_matrix.size1() == 0 || gpu_matrix.size2() == 0)
      {
        gpu_matrix.resize(cpu_matrix.rows(),
                          cpu_matrix.cols(),
                          false);
      }
      else
      {
        assert( (gpu_matrix.size1() == static_cast<std::size_t>(cpu_matrix.rows())) 
               && (gpu_matrix.size2() == static_cast<std::size_t>(cpu_matrix.cols()))
              );
      }

      std::vector<double> data(gpu_matrix.internal_size());
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
          data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
      }
      
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
    }
    #endif
    
    #ifdef VIENNACL_HAVE_MTL4
    /** @brief Copies a dense MTL matrix from the host (CPU) to the OpenCL device (GPU or multi-core CPU)
    *
    * @param cpu_matrix   A dense MTL matrix. cpu_matrix(i, j) returns the element in the i-th row and j-th columns (both starting with zero)
    * @param gpu_matrix   A dense ViennaCL matrix
    */
    template <typename SCALARTYPE, typename T, typename F, unsigned int ALIGNMENT>
    void copy(const mtl::dense2D<SCALARTYPE, T>& cpu_matrix,
              matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix)
    {
      typedef typename matrix<SCALARTYPE, F, ALIGNMENT>::size_type      size_type;
      
      if (gpu_matrix.size1() == 0 || gpu_matrix.size2() == 0)
      {
        gpu_matrix.resize(cpu_matrix.num_rows(),
                          cpu_matrix.num_cols(),
                          false);
      }
      else
      {
        assert( (gpu_matrix.size1() == cpu_matrix.num_rows()) 
               && (gpu_matrix.size2() == cpu_matrix.num_cols())
              );
      }

      std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
          data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix[i][j];
      }
      
      gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
    }
    #endif
    
    
    
    
    //
    //gpu to cpu, generic type
    //
    /** @brief Copies a dense matrix from the OpenCL device (GPU or multi-core CPU) to the host (CPU). 
    *
    * @param gpu_matrix   A dense ViennaCL matrix
    * @param cpu_matrix   A dense memory on the host. Must have at least as many rows and columns as the gpu_matrix! Type requirement: Access to entries via operator()
    */
    template <typename CPU_MATRIX, typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void copy(const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix,
              CPU_MATRIX & cpu_matrix )
    {
      typedef typename matrix<float, F, ALIGNMENT>::size_type      size_type;
      
      if ( (gpu_matrix.size1() > 0) && (gpu_matrix.size2() > 0) )
      {
        std::vector<SCALARTYPE> temp_buffer(gpu_matrix.internal_size());
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), &(temp_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        //now copy entries to cpu_matrix:
        for (size_type i = 0; i < gpu_matrix.size1(); ++i)
          for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
            cpu_matrix(i,j) = temp_buffer[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
      }
    }

    //gpu to cpu, STL type
    /** @brief Copies a dense matrix from the OpenCL device (GPU or multi-core CPU) to the host (CPU). 
    *
    * @param gpu_matrix   A dense ViennaCL matrix
    * @param cpu_matrix   A dense memory on the host using STL types, typically std::vector< std::vector<> > Must have at least as many rows and columns as the gpu_matrix! Type requirement: Access to entries via operator()
    */
    template <typename SCALARTYPE, typename A1, typename A2, typename F, unsigned int ALIGNMENT>
    void copy(const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix,
              std::vector< std::vector<SCALARTYPE, A1>, A2> & cpu_matrix)
    {
      typedef typename matrix<float, F, ALIGNMENT>::size_type      size_type;
      
      if ( (gpu_matrix.size1() > 0) && (gpu_matrix.size2() > 0) 
         && (cpu_matrix.size() >= gpu_matrix.size1()) && (cpu_matrix[0].size() >= gpu_matrix.size2()))
      {
        std::vector<SCALARTYPE> temp_buffer(gpu_matrix.internal_size());
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), gpu_matrix.handle().get(), CL_TRUE, 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), &(temp_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        //now copy entries to cpu_matrix:
        for (size_type i = 0; i < gpu_matrix.size1(); ++i)
          for (size_type j = 0; j < gpu_matrix.size2(); ++j) 
            cpu_matrix[i][j] = temp_buffer[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
      }
    }

    //gpu to cpu, STL type
    /** @brief Copies a dense matrix from the OpenCL device (GPU or multi-core CPU) to the host (CPU). 
    *
    * @param gpu_matrix         A dense ViennaCL matrix
    * @param cpu_matrix_begin   Pointer to the output memory on the CPU. User must ensure that provided memory is large enough.
    */
    template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void fast_copy(const matrix<SCALARTYPE, F, ALIGNMENT> & gpu_matrix,
                   SCALARTYPE * cpu_matrix_begin)
    {
      cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                       gpu_matrix.handle().get(), 
                                       CL_TRUE, 0,
                                       sizeof(SCALARTYPE)*gpu_matrix.internal_size(),
                                       cpu_matrix_begin, 0, NULL, NULL);
      VIENNACL_ERR_CHECK(err);
    }









    // outer_prod(v1, v2) * val;
    template<typename CPU_SCALAR, typename SCALARTYPE,unsigned int VECTOR_ALIGNMENT>
    viennacl::matrix_expression< const viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                    const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                    op_prod>,
                                 const SCALARTYPE,
                                 op_prod>  operator*(const viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                                        const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                                        op_prod> & proxy,
                                                     CPU_SCALAR const & val)
    {
      return viennacl::matrix_expression< const viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                             const viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT>,
                                                                             op_prod>,
                                          const SCALARTYPE,
                                          op_prod>(proxy, static_cast<SCALARTYPE>(val));
    }

    // val * outer_prod(v1, v2);
    template <typename CPU_SCALAR, typename SCALARTYPE, unsigned int VA1, unsigned int VA2>
    viennacl::matrix_expression< const viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VA1>,
                                                                    const viennacl::vector<SCALARTYPE, VA2>,
                                                                    op_prod>,
                                 const SCALARTYPE,
                                 op_prod>  operator*(CPU_SCALAR const & val,
                                                     viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VA1>,
                                                                                  const viennacl::vector<SCALARTYPE, VA2>,
                                                                                  op_prod> const & proxy)
    {
      return viennacl::matrix_expression< const viennacl::matrix_expression< const viennacl::vector<SCALARTYPE, VA1>,
                                                                             const viennacl::vector<SCALARTYPE, VA2>,
                                                                             op_prod>,
                                          const SCALARTYPE,
                                          op_prod>(proxy, static_cast<SCALARTYPE>(val));
    }
    
   

} //namespace viennacl

#endif
