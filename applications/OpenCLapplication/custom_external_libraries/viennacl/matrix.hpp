#ifndef VIENNACL_MATRIX_HPP_
#define VIENNACL_MATRIX_HPP_

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

/** @file viennacl/matrix.hpp
    @brief Implementation of the dense matrix class
*/

#include "viennacl/forwards.h"
#include "viennacl/scalar.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/linalg/matrix_operations.hpp"
#include "viennacl/linalg/sparse_matrix_operations.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/matrix_size_deducer.hpp"
#include "viennacl/meta/result_of.hpp"
#include "viennacl/meta/enable_if.hpp"
//#include "viennacl/rand/utils.hpp"
#include "viennacl/traits/handle.hpp"

namespace viennacl
{
  /** @brief Base class for representing matrices where the individual entries are not all stored explicitly, e.g. identity_matrix<>
    *
    * Examples are identity_matrix, scalar_matrix, and zero_matrix.
    */
  template<typename SCALARTYPE>
  class implicit_matrix_base
  {
    protected:
      typedef vcl_size_t        size_type;
      implicit_matrix_base(size_type size1, size_type size2, std::pair<SCALARTYPE, bool> value, bool diag) : size1_(size1), size2_(size2), value_(value), diag_(diag){ }
    public:
      typedef SCALARTYPE const & const_reference;
      typedef SCALARTYPE cpu_value_type;

      size_type size1() const { return size1_; }
      size_type size2() const { return size2_; }

      SCALARTYPE  value() const { return value_.first; }
      bool is_value_static( ) const { return value_.second; }
      bool diag() const { return diag_; }

      const_reference operator()(size_type i, size_type j) const {
        if(diag_) return (i == j) ? value_.first : 0;
        return value_.first;
      }

    protected:
      size_type size1_;
      size_type size2_;
      std::pair<SCALARTYPE, bool> value_;
      bool diag_;
  };

  //
  // Initializer types
  //
  /** @brief Represents a vector consisting of 1 at a given index and zeros otherwise. To be used as an initializer for viennacl::vector, vector_range, or vector_slize only. */
  template <typename SCALARTYPE>
  class identity_matrix
  {
    public:
      typedef vcl_size_t         size_type;
      typedef SCALARTYPE const & const_reference;

      identity_matrix(size_type s, viennacl::context ctx = viennacl::context()) : size_(s), diag_(1), off_diag_(0), ctx_(ctx) {}

      size_type size1() const { return size_; }
      size_type size2() const { return size_; }
      const_reference operator()(size_type i, size_type j) const { return (i == j) ? diag_ : off_diag_; }

      viennacl::context context() const { return ctx_; }

    private:
      size_type size_;
      SCALARTYPE diag_;
      SCALARTYPE off_diag_;
      viennacl::context ctx_;
  };


  /** @brief Represents a vector consisting of zeros only. To be used as an initializer for viennacl::vector, vector_range, or vector_slize only. */
  template <typename SCALARTYPE>
  class zero_matrix
  {
    public:
      typedef vcl_size_t         size_type;
      typedef SCALARTYPE const & const_reference;

      zero_matrix(size_type s1, size_type s2, viennacl::context ctx = viennacl::context()) : size1_(s1), size2_(s2), val_(0), ctx_(ctx) {}

      size_type size1() const { return size1_; }
      size_type size2() const { return size2_; }
      const_reference operator()(size_type /*i*/, size_type /*j*/) const { return val_; }

      viennacl::context context() const { return ctx_; }

    private:
      size_type size1_;
      size_type size2_;
      SCALARTYPE val_;
      viennacl::context ctx_;
  };


  /** @brief Represents a vector consisting of scalars 's' only, i.e. v[i] = s for all i. To be used as an initializer for viennacl::vector, vector_range, or vector_slize only. */
  template <typename SCALARTYPE>
  class scalar_matrix
  {
    public:
      typedef vcl_size_t         size_type;
      typedef SCALARTYPE const & const_reference;

      scalar_matrix(size_type s1, size_type s2, const_reference val, viennacl::context ctx = viennacl::context()) : size1_(s1), size2_(s2), value_(val), ctx_(ctx) {}

      size_type size1() const { return size1_; }
      size_type size2() const { return size2_; }
      const_reference operator()(size_type /*i*/, size_type /*j*/) const { return value_; }

      viennacl::context context() const { return ctx_; }

    private:
      size_type size1_;
      size_type size2_;
      SCALARTYPE value_;
      viennacl::context ctx_;
  };



//#ifdef VIENNACL_WITH_OPENCL
//  template<class SCALARTYPE, class DISTRIBUTION>
//  rand::random_matrix_t<SCALARTYPE, DISTRIBUTION> random_matrix(unsigned int size1, unsigned int size2, DISTRIBUTION const & distribution){
//      return rand::random_matrix_t<SCALARTYPE,DISTRIBUTION>(size1,size2,distribution);
//  }
//#endif

  /** @brief Expression template class for representing a tree of expressions which ultimately result in a matrix.
    *
    * @tparam LHS   The left hand side of the expression tree
    * @tparam RHS   The right hand side of the expression tree
    * @tparam OP    The operator to apply to LHS and RHS to obtain the result.
    */
  template <typename LHS, typename RHS, typename OP>
  class matrix_expression
  {
      typedef typename viennacl::result_of::reference_if_nonscalar<LHS>::type     lhs_reference_type;
      typedef typename viennacl::result_of::reference_if_nonscalar<RHS>::type     rhs_reference_type;

    public:
      typedef vcl_size_t       size_type;

      matrix_expression(LHS & lhs, RHS & rhs) : lhs_(lhs), rhs_(rhs) {}

      /** @brief Get left hand side operand
      */
      LHS & lhs() const { return lhs_; }
      /** @brief Get right hand side operand
      */
      RHS & rhs() const { return rhs_; }

      /** @brief Returns the size of the result vector */
      vcl_size_t size1() const { return viennacl::tools::MATRIX_SIZE_DEDUCER<LHS, RHS, OP>::size1(lhs_, rhs_); }
      vcl_size_t size2() const { return viennacl::tools::MATRIX_SIZE_DEDUCER<LHS, RHS, OP>::size2(lhs_, rhs_); }

    private:
      /** @brief The left hand side operand */
      lhs_reference_type lhs_;
      /** @brief The right hand side operand */
      rhs_reference_type rhs_;
  };


  /** @brief A tag indicating iteration along increasing row index of a matrix */
  struct row_iteration {};

  /** @brief A tag indicating iteration along increasing columns index of a matrix */
  struct col_iteration {};

  //STL-like iterator. TODO: STL-compliance...
  /** @brief uBLAS-like iterator class for iterating over the entries of a dense matrix. */
  template <typename ROWCOL, typename MATRIXTYPE>
  class matrix_iterator
  {
      typedef matrix_iterator<ROWCOL, MATRIXTYPE>    self_type;
    public:
      typedef typename MATRIXTYPE::value_type       value_type;

      matrix_iterator(MATRIXTYPE & mat,
                      vcl_size_t start_row,
                      vcl_size_t start_col) : mat_(mat), row_(start_row), col_(start_col) {}

      value_type operator*(void) { return mat_(row_, col_); }
      self_type & operator++(void) { viennacl::tools::MATRIX_ITERATOR_INCREMENTER<ROWCOL, MATRIXTYPE>::apply(mat_, row_, col_); return *this; }
      self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }

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
  template <class SCALARTYPE, typename F, typename SizeType /* see forwards.h for default type */, typename DistanceType /* see forwards.h for default type */>
  class matrix_base
  {
      typedef matrix_base<SCALARTYPE, F, SizeType, DistanceType>          self_type;
    public:

      typedef matrix_iterator<row_iteration, self_type >   iterator1;
      typedef matrix_iterator<col_iteration, self_type >   iterator2;
      typedef scalar<SCALARTYPE>                                                  value_type;
      typedef SCALARTYPE                                                          cpu_value_type;
      typedef SizeType                                                            size_type;
      typedef DistanceType                                                        difference_type;
      typedef viennacl::backend::mem_handle                                       handle_type;
      typedef F                                                                   orientation_functor;
      typedef typename F::orientation_category                                    orientation_category;

      static const size_type alignment = 128;


      /** @brief The default constructor. Does not allocate any memory. */
      explicit matrix_base() : size1_(0), size2_(0), start1_(0), start2_(0), stride1_(1), stride2_(1), internal_size1_(0), internal_size2_(0) {}

      /** @brief Creates the matrix with the given dimensions
      *
      * @param rows     Number of rows
      * @param columns  Number of columns
      * @param ctx      Optional context in which the matrix is created (one out of multiple OpenCL contexts, CUDA, host)
      */
      explicit matrix_base(size_type rows, size_type columns, viennacl::context ctx = viennacl::context())
          : size1_(rows), size2_(columns), start1_(0), start2_(0), stride1_(1), stride2_(1),
            internal_size1_(viennacl::tools::align_to_multiple<size_type>(rows, alignment)),
            internal_size2_(viennacl::tools::align_to_multiple<size_type>(columns, alignment))
      {
        if (rows > 0 && columns > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), ctx);
          clear();
        }
      }


      /** @brief Constructor for creating a matrix_range or matrix_stride from some other matrix/matrix_range/matrix_stride */
      explicit matrix_base(viennacl::backend::mem_handle & h,
                           size_type mat_size1, size_type mat_start1, difference_type mat_stride1, size_type mat_internal_size1,
                           size_type mat_size2, size_type mat_start2, difference_type mat_stride2, size_type mat_internal_size2)
        : size1_(mat_size1), size2_(mat_size2),
          start1_(mat_start1), start2_(mat_start2),
          stride1_(mat_stride1), stride2_(mat_stride2),
          internal_size1_(mat_internal_size1), internal_size2_(mat_internal_size2),
          elements_(h) {}

      template <typename LHS, typename RHS, typename OP>
      explicit matrix_base(matrix_expression<const LHS, const RHS, OP> const & proxy) :
        size1_(viennacl::traits::size1(proxy)), size2_(viennacl::traits::size2(proxy)), start1_(0), start2_(0), stride1_(1), stride2_(1),
        internal_size1_(viennacl::tools::align_to_multiple<size_type>(size1_, alignment)),
        internal_size2_(viennacl::tools::align_to_multiple<size_type>(size2_, alignment))
      {
        elements_.switch_active_handle_id(viennacl::traits::active_handle_id(proxy));
        if (internal_size() > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(proxy));
          clear();
          self_type::operator=(proxy);
        }
      }

      // CUDA or host memory:
      explicit matrix_base(SCALARTYPE * ptr_to_mem, viennacl::memory_types mem_type,
                           size_type mat_size1, size_type mat_start1, difference_type mat_stride1, size_type mat_internal_size1,
                           size_type mat_size2, size_type mat_start2, difference_type mat_stride2, size_type mat_internal_size2)
        : size1_(mat_size1), size2_(mat_size2),
          start1_(mat_start1), start2_(mat_start2),
          stride1_(mat_stride1), stride2_(mat_stride2),
          internal_size1_(mat_internal_size1), internal_size2_(mat_internal_size2)
      {
        if (mem_type == viennacl::CUDA_MEMORY)
        {
#ifdef VIENNACL_WITH_CUDA
          elements_.switch_active_handle_id(viennacl::CUDA_MEMORY);
          elements_.cuda_handle().reset(reinterpret_cast<char*>(ptr_to_mem));
          elements_.cuda_handle().inc(); //prevents that the user-provided memory is deleted once the vector object is destroyed.
#else
          throw cuda_not_available_exception();
#endif
        }
        else if (mem_type == viennacl::MAIN_MEMORY)
        {
          elements_.switch_active_handle_id(viennacl::MAIN_MEMORY);
          elements_.ram_handle().reset(reinterpret_cast<char*>(ptr_to_mem));
          elements_.ram_handle().inc(); //prevents that the user-provided memory is deleted once the vector object is destroyed.
        }

        elements_.raw_size(sizeof(SCALARTYPE) * internal_size());
      }

#ifdef VIENNACL_WITH_OPENCL
      explicit matrix_base(cl_mem mem, size_type rows, size_type columns, viennacl::context ctx = viennacl::context())
        : size1_(rows), size2_(columns),
          start1_(0), start2_(0),
          stride1_(1), stride2_(1),
          internal_size1_(rows), internal_size2_(columns)
      {
        elements_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
        elements_.opencl_handle() = mem;
        elements_.opencl_handle().inc();  //prevents that the user-provided memory is deleted once the vector object is destroyed.
        elements_.opencl_handle().context(ctx.opencl_context());
        elements_.raw_size(sizeof(SCALARTYPE)*internal_size());
      }

      explicit matrix_base(cl_mem mem, viennacl::context ctx,
                           size_type mat_size1, size_type mat_start1, difference_type mat_stride1, size_type mat_internal_size1,
                           size_type mat_size2, size_type mat_start2, difference_type mat_stride2, size_type mat_internal_size2)
        : size1_(mat_size1), size2_(mat_size2),
          start1_(mat_start1), start2_(mat_start2),
          stride1_(mat_stride1), stride2_(mat_stride2),
          internal_size1_(mat_internal_size1), internal_size2_(mat_internal_size2)
      {
        elements_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
        elements_.opencl_handle() = mem;
        elements_.opencl_handle().inc();  //prevents that the user-provided memory is deleted once the vector object is destroyed.
        elements_.opencl_handle().context(ctx.opencl_context());
        elements_.raw_size(sizeof(SCALARTYPE)*internal_size());
      }
#endif


      self_type & operator=(const self_type & other)  //enables implicit conversions
      {
        if (internal_size() == 0)
        {
          if (other.internal_size() == 0)
            return *this;
          resize(other.size1(), other.size2(), false);
        }

        viennacl::linalg::am(*this,
                             other, cpu_value_type(1.0), 1, false, false);
        return *this;
      }

      /** @brief Creates the matrix from the supplied random matrix. */
      /*template<class DISTRIBUTION>
      matrix(rand::random_matrix_t<SCALARTYPE, DISTRIBUTION> const & m) : rows_(m.size1), columns_(m.size2)
      {
        if (internal_size() > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size());
          rand::buffer_dumper<SCALARTYPE, DISTRIBUTION>::dump(elements_,m.distribution,0,internal_size());
        }
      }*/



      /** @brief Implementation of the operation m1 = m2 @ alpha, where @ denotes either multiplication or division, and alpha is either a CPU or a GPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename LHS, typename RHS, typename OP>
      self_type & operator=(const matrix_expression<const LHS, const RHS, OP> & proxy)
      {
        assert(  (viennacl::traits::size1(proxy) == size1() || size1() == 0)
              && (viennacl::traits::size2(proxy) == size2() || size2() == 0)
              && bool("Incompatible matrix sizes!"));

        if (internal_size() == 0 && viennacl::traits::size1(proxy) > 0 && viennacl::traits::size2(proxy) > 0)
        {
          size1_ = viennacl::traits::size1(proxy);
          size2_ = viennacl::traits::size2(proxy);
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(proxy));
          if (size1_ != internal_size1_ || size2_ != internal_size2_)
            clear();
        }

        if (internal_size() > 0)
          linalg::detail::op_executor<self_type, op_assign, matrix_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }


      // A = trans(B). Currently achieved in CPU memory
      self_type & operator=(const matrix_expression< const self_type,
                                                     const self_type,
                                                     op_trans> & proxy)
      {
        assert( (handle() != proxy.lhs().handle()) && bool("Self-assignment of matrix transpose not implemented"));
        assert( ( (proxy.lhs().size1() == size2()) || (size2() == 0) ) && bool("Matrix dimensions do not match!"));
        assert( ( (proxy.lhs().size2() == size1()) || (size1() == 0) ) && bool("Matrix dimensions do not match!"));

        if (internal_size() == 0 && viennacl::traits::size1(proxy) > 0 && viennacl::traits::size2(proxy) > 0)
        {
          size1_ = viennacl::traits::size1(proxy);
          size2_ = viennacl::traits::size2(proxy);
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
        }

        std::vector<SCALARTYPE> temp(proxy.lhs().internal_size());

        viennacl::backend::memory_read(proxy.lhs().handle(), 0, sizeof(SCALARTYPE)*proxy.lhs().internal_size(), &(temp[0]));

        // now transpose it
        std::vector<SCALARTYPE> temp_trans(internal_size());

        for (vcl_size_t i=0; i<proxy.lhs().size1(); ++i)
          for (vcl_size_t j=0; j<proxy.lhs().size2(); ++j)
            temp_trans[F::mem_index(start2() + stride2() * j,
                                    start1() + stride1() * i,
                                    internal_size1(), internal_size2())]
              = temp[F::mem_index(proxy.lhs().start1() + proxy.lhs().stride1() * i,
                                  proxy.lhs().start2() + proxy.lhs().stride2() * j,
                                  proxy.lhs().internal_size1(), proxy.lhs().internal_size2())];

        // write back
        viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(proxy), &(temp_trans[0]));

        return *this;
      }

      template <typename LHS, typename RHS, typename OP>
      self_type & operator+=(const matrix_expression<const LHS, const RHS, OP> & proxy)
      {
        assert(  (viennacl::traits::size1(proxy) == size1())
              && (viennacl::traits::size2(proxy) == size2())
              && bool("Incompatible matrix sizes!"));
        assert( (size1() > 0) && bool("Vector not yet initialized!") );
        assert( (size2() > 0) && bool("Vector not yet initialized!") );

        linalg::detail::op_executor<self_type, op_inplace_add, matrix_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }

      template <typename LHS, typename RHS, typename OP>
      self_type & operator-=(const matrix_expression<const LHS, const RHS, OP> & proxy)
      {
        assert(  (viennacl::traits::size1(proxy) == size1())
              && (viennacl::traits::size2(proxy) == size2())
              && bool("Incompatible matrix sizes!"));
        assert( (size1() > 0) && bool("Vector not yet initialized!") );
        assert( (size2() > 0) && bool("Vector not yet initialized!") );

        linalg::detail::op_executor<self_type, op_inplace_sub, matrix_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }

      /** @brief Assigns the supplied identity matrix to the matrix. */
      self_type & operator = (identity_matrix<SCALARTYPE> const & m)
      {
        assert( (m.size1() == size1_ || size1_ == 0) && bool("Size mismatch!") );
        assert( (m.size2() == size2_ || size2_ == 0) && bool("Size mismatch!") );

        if (internal_size() == 0)
        {
          size1_ = m.size1();
          size2_ = m.size2();
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
          if (internal_size() > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), m.context());
            clear();
          }
        }
        else
          viennacl::linalg::matrix_assign(*this, SCALARTYPE(0));

        if (internal_size() > 0)
          viennacl::linalg::matrix_diagonal_assign(*this, m(0,0));

        return *this;
      }

      /** @brief Assigns the supplied zero matrix to the matrix. */
      self_type & operator = (zero_matrix<SCALARTYPE> const & m)
      {
        assert( (m.size1() == size1_ || size1_ == 0) && bool("Size mismatch!") );
        assert( (m.size2() == size2_ || size2_ == 0) && bool("Size mismatch!") );

        if (internal_size() == 0)
        {
          size1_ = m.size1();
          size2_ = m.size2();
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
          if (internal_size() > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), m.context());
            clear();
          }
        }
        else
          viennacl::linalg::matrix_assign(*this, SCALARTYPE(0));

        return *this;
      }

      /** @brief Assigns the supplied scalar vector to the matrix. */
      self_type & operator = (scalar_matrix<SCALARTYPE> const & m)
      {
        assert( (m.size1() == size1_ || size1_ == 0) && bool("Size mismatch!") );
        assert( (m.size2() == size2_ || size2_ == 0) && bool("Size mismatch!") );

        if (internal_size() == 0)
        {
          size1_ = m.size1();
          size2_ = m.size2();
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
          if (internal_size() > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), m.context());
            clear();
          }
        }

        if (internal_size() > 0)
        {
          viennacl::linalg::matrix_assign(*this, m(0,0));
        }

        return *this;
      }


      //read-write access to an element of the matrix/matrix_range/matrix_slice
      /** @brief Read-write access to a single element of the matrix/matrix_range/matrix_slice
      */
      entry_proxy<SCALARTYPE> operator()(size_type row_index, size_type col_index)
      {
        return entry_proxy<SCALARTYPE>(F::mem_index(start1_ + stride1_ * row_index, start2_ + stride2_ * col_index, internal_size1(), internal_size2()), elements_);
      }

      /** @brief Read access to a single element of the matrix/matrix_range/matrix_slice
      */
      const_entry_proxy<SCALARTYPE> operator()(size_type row_index, size_type col_index) const
      {
        return const_entry_proxy<SCALARTYPE>(F::mem_index(start1_ + stride1_ * row_index, start2_ + stride2_ * col_index, internal_size1(), internal_size2()), elements_);
      }

      //
      // Operator overloads for enabling implicit conversions:
      //
      self_type & operator += (const self_type & other)
      {
        viennacl::linalg::ambm(*this,
                                *this, SCALARTYPE(1.0), 1, false, false,
                                other, SCALARTYPE(1.0), 1, false, false);
        return *this;
      }

      self_type & operator -= (const self_type & other)
      {
        viennacl::linalg::ambm(*this,
                                *this, SCALARTYPE(1.0), 1, false, false,
                                other, SCALARTYPE(1.0), 1, false, true);
        return *this;
      }

      /** @brief Scales a matrix by a CPU scalar value
      */
      self_type & operator *= (SCALARTYPE val)
      {
        //viennacl::linalg::inplace_mult(*this, val);
        viennacl::linalg::am(*this,
                              *this, val, 1, false, false);
        return *this;
      }

      /** @brief Scales this matrix by a CPU scalar value
      */
      self_type & operator /= (SCALARTYPE val)
      {
        //viennacl::linalg::inplace_mult(*this, static_cast<SCALARTYPE>(1) / val);
        viennacl::linalg::am(*this,
                              *this, val, 1, true, false);
        return *this;
      }


      /** @brief Sign flip for the matrix. Emulated to be equivalent to -1.0 * matrix */
      matrix_expression<const self_type, const SCALARTYPE, op_mult> operator-() const
      {
        return matrix_expression<const self_type, const SCALARTYPE, op_mult>(*this, SCALARTYPE(-1));
      }

      /** @brief Returns the number of rows */
      size_type size1() const { return size1_;}
      /** @brief Returns the number of columns */
      size_type size2() const { return size2_; }

      /** @brief Returns the number of rows */
      size_type start1() const { return start1_;}
      /** @brief Returns the number of columns */
      size_type start2() const { return start2_; }

      /** @brief Returns the number of rows */
      size_type stride1() const { return stride1_;}
      /** @brief Returns the number of columns */
      size_type stride2() const { return stride2_; }

      /** @brief Resets all entries to zero */
      void clear()
      {
        viennacl::linalg::matrix_assign(*this, SCALARTYPE(0), true);
      }


      /** @brief Returns the internal number of rows. Usually required for launching OpenCL kernels only */
      size_type internal_size1() const { return internal_size1_; }
      /** @brief Returns the internal number of columns. Usually required for launching OpenCL kernels only */
      size_type internal_size2() const { return internal_size2_; }
      /** @brief Returns the total amount of allocated memory in multiples of sizeof(SCALARTYPE) */
      size_type internal_size() const { return internal_size1() * internal_size2(); }

      /** @brief Returns the OpenCL handle, non-const-version */
            handle_type & handle()       { return elements_; }
      /** @brief Returns the OpenCL handle, const-version */
      const handle_type & handle() const { return elements_; }


      viennacl::memory_types memory_domain() const
      {
        return elements_.get_active_handle_id();
      }

    protected:

      void set_handle(viennacl::backend::mem_handle const & h)
      {
        elements_ = h;
      }

      void switch_memory_context(viennacl::context new_ctx)
      {
        viennacl::backend::switch_memory_context<SCALARTYPE>(elements_, new_ctx);
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
        assert( (rows > 0 && columns > 0) && bool("Check failed in matrix::resize(): Number of rows and columns must be positive!"));

        if (preserve && internal_size() > 0)
        {
          //get old entries:
          std::vector< SCALARTYPE > old_entries(internal_size());
          viennacl::backend::memory_read(elements_, 0, sizeof(SCALARTYPE)*internal_size(), &(old_entries[0]));

          //set up entries of new matrix:
          std::vector< SCALARTYPE > new_entries(  viennacl::tools::align_to_multiple<vcl_size_t>(rows,    alignment)
                                                * viennacl::tools::align_to_multiple<vcl_size_t>(columns, alignment));
          for (size_type i=0; i<rows; ++i)
          {
            if (i >= size1_)
              continue;

            for (size_type j=0; j<columns; ++j)
            {
              if (j >= size2_)
                continue;
              new_entries[F::mem_index(i, j, viennacl::tools::align_to_multiple<vcl_size_t>(rows, alignment), viennacl::tools::align_to_multiple<vcl_size_t>(columns, alignment))]
                  = old_entries[F::mem_index(i, j, internal_size1(), internal_size2())];
            }
          }

          //copy new entries to GPU:
          size1_ = rows;
          size2_ = columns;
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*new_entries.size(), viennacl::traits::context(elements_), &(new_entries[0]));
        }
        else //discard old entries:
        {
          size1_ = rows;
          size2_ = columns;
          internal_size1_ = viennacl::tools::align_to_multiple<size_type>(size1_, alignment);
          internal_size2_ = viennacl::tools::align_to_multiple<size_type>(size2_, alignment);

          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(elements_));
          clear();
        }
      }

    private:
      size_type size1_;
      size_type size2_;
      size_type start1_;
      size_type start2_;
      difference_type stride1_;
      difference_type stride2_;
      size_type internal_size1_;
      size_type internal_size2_;
      handle_type elements_;
  }; //matrix



  /** @brief A dense matrix class
  *
  * @tparam SCALARTYPE   The underlying scalar type (either float or double)
  * @tparam F            Storage layout: Either row_major or column_major (at present only row_major is supported)
  * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
  */
  template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
  class matrix : public matrix_base<SCALARTYPE, F>
  {
      typedef matrix<SCALARTYPE, F, ALIGNMENT>          self_type;
      typedef matrix_base<SCALARTYPE, F>                base_type;
    public:
      typedef typename base_type::size_type             size_type;

      /** @brief The default constructor. Does not allocate any memory. */
      explicit matrix() : base_type() {}

      /** @brief Creates the matrix with the given dimensions
      *
      * @param rows     Number of rows
      * @param columns  Number of columns
      * @param ctx      Optional context in which the matrix is created (one out of multiple OpenCL contexts, CUDA, host)
      */
      explicit matrix(size_type rows, size_type columns, viennacl::context ctx = viennacl::context()) : base_type(rows, columns, ctx) {}

#ifdef VIENNACL_WITH_OPENCL
      explicit matrix(cl_mem mem, size_type rows, size_type columns) : base_type(mem, rows, columns) {}
#endif

      template <typename LHS, typename RHS, typename OP>
      matrix(matrix_expression< LHS, RHS, OP> const & proxy) : base_type(proxy) {}

      /** @brief Creates the matrix from the supplied identity matrix. */
      matrix(identity_matrix<SCALARTYPE> const & m) : base_type(m.size1(), m.size2(), m.context())
      {
        if (base_type::internal_size() > 0)
          base_type::operator=(m);
      }

      /** @brief Creates the matrix from the supplied zero matrix. */
      matrix(zero_matrix<SCALARTYPE> const & m) : base_type(m.size1(), m.size2(), m.context())
      {
        if (base_type::internal_size() > 0)
          base_type::operator=(m);
      }

      /** @brief Creates the matrix from the supplied scalar matrix. */
      matrix(scalar_matrix<SCALARTYPE> const & m) : base_type(m.size1(), m.size2(), m.context())
      {
        if (base_type::internal_size() > 0)
          base_type::operator=(m);
      }

      matrix(const base_type & other) : base_type(other.size1(), other.size2(), viennacl::traits::context(other))
      {
        base_type::operator=(other);
      }


      //copy constructor:
      matrix(const self_type & other) : base_type(other.size1(), other.size2(), viennacl::traits::context(other))
      {
        base_type::operator=(other);
      }


      /*template <typename M1>
      self_type & operator=(const matrix_expression< const M1, const M1, op_trans> & proxy)
      {
        self_type temp(proxy.lhs());
        *this = trans(temp);
        return *this;
      }*/

      using base_type::operator=;

      /** @brief Resizes the matrix.
      *   Existing entries can optionally be preserved
      *
      * @param rows       New number of rows
      * @param columns    New number of columns
      * @param preserve   If true, existing values are preserved.
      */
      void resize(size_type rows, size_type columns, bool preserve = true)
      {
        base_type::resize(rows, columns, preserve);
      }

  }; //matrix



  /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
  *
  * @param s            STL output stream
  * @param gpu_matrix   A dense ViennaCL matrix
  */
  template<class SCALARTYPE, typename F>
  std::ostream & operator<<(std::ostream & s, const matrix_base<SCALARTYPE, F> & gpu_matrix)
  {
    typedef typename matrix_base<SCALARTYPE, F>::size_type      size_type;

    std::vector<SCALARTYPE> tmp(gpu_matrix.internal_size());
    viennacl::backend::memory_read(gpu_matrix.handle(), 0, sizeof(SCALARTYPE) * gpu_matrix.internal_size(), &(tmp[0]));

    s << "[" << gpu_matrix.size1() << "," << gpu_matrix.size2() << "]";

    s << "(";
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      s << "(";
      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
      {
        s << tmp[F::mem_index(i * gpu_matrix.stride1() + gpu_matrix.start1(), j * gpu_matrix.stride2() + gpu_matrix.start2(), gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
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
  template<typename NumericT, typename F>
  matrix_expression< const matrix_base<NumericT, F>, const matrix_base<NumericT, F>, op_trans>
  trans(const matrix_base<NumericT, F> & mat)
  {
    return matrix_expression< const matrix_base<NumericT, F>, const matrix_base<NumericT, F>, op_trans>(mat, mat);
  }

  //diag():
  template<typename NumericT, typename F>
  vector_expression< const matrix_base<NumericT, F>, const int, op_matrix_diag>
  diag(const matrix_base<NumericT, F> & A, int k = 0)
  {
    return vector_expression< const matrix_base<NumericT, F>, const int, op_matrix_diag>(A, k);
  }

  template<typename NumericT>
  matrix_expression< const vector_base<NumericT>, const int, op_vector_diag>
  diag(const vector_base<NumericT> & v, int k = 0)
  {
    return matrix_expression< const vector_base<NumericT>, const int, op_vector_diag>(v, k);
  }

  // row():
  template<typename NumericT, typename F>
  vector_expression< const matrix_base<NumericT, F>, const unsigned int, op_row>
  row(const matrix_base<NumericT, F> & A, unsigned int i)
  {
    return vector_expression< const matrix_base<NumericT, F>, const unsigned int, op_row>(A, i);
  }

  // column():
  template<typename NumericT, typename F>
  vector_expression< const matrix_base<NumericT, F>, const unsigned int, op_column>
  column(const matrix_base<NumericT, F> & A, unsigned int j)
  {
    return vector_expression< const matrix_base<NumericT, F>, const unsigned int, op_column>(A, j);
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

    assert( (gpu_matrix.size1() == cpu_matrix.size1()) && (gpu_matrix.size2() == cpu_matrix.size2()) && bool("Matrix dimensions mismatch.") );

    std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
    }

    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(SCALARTYPE) * data.size(), viennacl::traits::context(gpu_matrix), &(data[0]));
    //gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
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

    assert( (gpu_matrix.size1() == cpu_matrix.size()) && bool("Matrix dimensions mismatch.") );

    std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      assert( (gpu_matrix.size2() == cpu_matrix[i].size()) && bool("Matrix dimensions mismatch.") );

      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix[i][j];
    }

    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(SCALARTYPE) * data.size(), viennacl::traits::context(gpu_matrix), &(data[0]));
    //gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
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
    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(SCALARTYPE) * (cpu_matrix_end - cpu_matrix_begin), viennacl::traits::context(gpu_matrix), cpu_matrix_begin);
    /*gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE,
                                                                          sizeof(SCALARTYPE) * (cpu_matrix_end - cpu_matrix_begin),
                                                                          cpu_matrix_begin);*/
  }


  #ifdef VIENNACL_WITH_EIGEN
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
      assert( (gpu_matrix.size1() == static_cast<vcl_size_t>(cpu_matrix.rows()))
              && (gpu_matrix.size2() == static_cast<vcl_size_t>(cpu_matrix.cols()))
              && bool("matrix size mismatch")
            );
    }

    std::vector<float> data(gpu_matrix.internal_size());
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
    }

    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(float) * data.size(), viennacl::traits::context(gpu_matrix), &(data[0]));
    //gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
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
      assert( (gpu_matrix.size1() == static_cast<vcl_size_t>(cpu_matrix.rows()))
              && (gpu_matrix.size2() == static_cast<vcl_size_t>(cpu_matrix.cols()))
              && bool("matrix size mismatch")
            );
    }

    std::vector<double> data(gpu_matrix.internal_size());
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix(i,j);
    }

    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(double) * data.size(), viennacl::traits::context(gpu_matrix), &(data[0]));
    //gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
  }
  #endif

  #ifdef VIENNACL_WITH_MTL4
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
              && bool("matrix size mismatch")
            );
    }

    std::vector<SCALARTYPE> data(gpu_matrix.internal_size());
    for (size_type i = 0; i < gpu_matrix.size1(); ++i)
    {
      for (size_type j = 0; j < gpu_matrix.size2(); ++j)
        data[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())] = cpu_matrix[i][j];
    }

    viennacl::backend::memory_create(gpu_matrix.handle(), sizeof(SCALARTYPE) * data.size(), viennacl::traits::context(gpu_matrix), &(data[0]));
    //gpu_matrix.elements_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, data);
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
      assert( viennacl::traits::size1(cpu_matrix) == gpu_matrix.size1() && bool("Matrix dimensions mismatch: rows"));

      std::vector<SCALARTYPE> temp_buffer(gpu_matrix.internal_size());
      viennacl::backend::memory_read(gpu_matrix.handle(), 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), &(temp_buffer[0]));

      //now copy entries to cpu_matrix:
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        assert( viennacl::traits::size2(cpu_matrix) == gpu_matrix.size2() && bool("Matrix dimensions mismatch: columns"));
        for (size_type j = 0; j < gpu_matrix.size2(); ++j)
          cpu_matrix(i,j) = temp_buffer[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
      }
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

    if ( (gpu_matrix.size1() > 0) && (gpu_matrix.size2() > 0) )
    {
      assert( (cpu_matrix.size() == gpu_matrix.size1()) && bool("Matrix dimensions mismatch: rows"));

      std::vector<SCALARTYPE> temp_buffer(gpu_matrix.internal_size());
      viennacl::backend::memory_read(gpu_matrix.handle(), 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), &(temp_buffer[0]));

      //now copy entries to cpu_matrix:
      for (size_type i = 0; i < gpu_matrix.size1(); ++i)
      {
        assert( (cpu_matrix[i].size() == gpu_matrix.size2()) && bool("Matrix dimensions mismatch: columns"));

        for (size_type j = 0; j < gpu_matrix.size2(); ++j)
          cpu_matrix[i][j] = temp_buffer[F::mem_index(i, j, gpu_matrix.internal_size1(), gpu_matrix.internal_size2())];
      }
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
    viennacl::backend::memory_read(gpu_matrix.handle(), 0, sizeof(SCALARTYPE)*gpu_matrix.internal_size(), cpu_matrix_begin);
  }



  /////////////////////// matrix operator overloads to follow ////////////////////////////////////////////


  // operator +
  /** @brief Generic 'catch-all' overload, which enforces a temporary if the expression tree gets too deep. */
  template <typename LHS1, typename RHS1, typename OP1,
            typename LHS2, typename RHS2, typename OP2>
  matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                     const matrix_expression<const LHS2, const RHS2, OP2>,
                     op_add>
  operator + (matrix_expression<const LHS1, const RHS1, OP1> const & proxy1,
              matrix_expression<const LHS2, const RHS2, OP2> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                              const matrix_expression<const LHS2, const RHS2, OP2>,
                              op_add>(proxy1, proxy2);
  }

  template <typename LHS1, typename RHS1, typename OP1,
            typename NumericT, typename F>
  matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                     const matrix_base<NumericT, F>,
                     op_add>
  operator + (matrix_expression<const LHS1, const RHS1, OP1> const & proxy1,
              matrix_base<NumericT, F> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                              const matrix_base<NumericT, F>,
                              op_add>(proxy1, proxy2);
  }

  template <typename NumericT, typename F,
            typename LHS2, typename RHS2, typename OP2>
  matrix_expression< const matrix_base<NumericT, F>,
                     const matrix_expression<const LHS2, const RHS2, OP2>,
                     op_add>
  operator + (matrix_base<NumericT, F> const & proxy1,
              matrix_expression<const LHS2, const RHS2, OP2> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return  matrix_expression< const matrix_base<NumericT, F>,
                               const matrix_expression<const LHS2, const RHS2, OP2>,
                               op_add>(proxy1, proxy2);
  }

  /** @brief Operator overload for m1 + m2, where m1 and m2 are either dense matrices, matrix ranges, or matrix slices. No mixing of different storage layouts allowed at the moment. */
  template <typename NumericT, typename F>
  matrix_expression< const matrix_base<NumericT, F>, const matrix_base<NumericT, F>, op_add >
  operator + (const matrix_base<NumericT, F> & m1, const matrix_base<NumericT, F> & m2)
  {
    return matrix_expression< const matrix_base<NumericT, F>,
                              const matrix_base<NumericT, F>,
                              op_add > (m1, m2);
  }


  // operator -
  template <typename LHS1, typename RHS1, typename OP1,
            typename LHS2, typename RHS2, typename OP2>
  matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                     const matrix_expression<const LHS2, const RHS2, OP2>,
                     op_sub>
  operator - (matrix_expression<const LHS1, const RHS1, OP1> const & proxy1,
              matrix_expression<const LHS2, const RHS2, OP2> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                              const matrix_expression<const LHS2, const RHS2, OP2>,
                              op_sub>(proxy1, proxy2);
  }

  template <typename LHS1, typename RHS1, typename OP1,
            typename NumericT, typename F>
  matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                     const matrix_base<NumericT, F>,
                     op_sub>
  operator - (matrix_expression<const LHS1, const RHS1, OP1> const & proxy1,
              matrix_base<NumericT, F> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return matrix_expression< const matrix_expression<const LHS1, const RHS1, OP1>,
                              const matrix_base<NumericT, F>,
                              op_sub>(proxy1, proxy2);
  }

  template <typename NumericT, typename F,
            typename LHS2, typename RHS2, typename OP2>
  matrix_expression< const matrix_base<NumericT, F>,
                     const matrix_expression<const LHS2, const RHS2, OP2>,
                     op_sub>
  operator - (matrix_base<NumericT, F> const & proxy1,
              matrix_expression<const LHS2, const RHS2, OP2> const & proxy2)
  {
    assert(    (viennacl::traits::size1(proxy1) == viennacl::traits::size1(proxy2))
            && (viennacl::traits::size2(proxy1) == viennacl::traits::size2(proxy2))
            && bool("Incompatible matrix sizes!"));
    return  matrix_expression< const matrix_base<NumericT, F>,
                               const matrix_expression<const LHS2, const RHS2, OP2>,
                               op_sub>(proxy1, proxy2);
  }

  /** @brief Operator overload for m1 - m2, where m1 and m2 are either dense matrices, matrix ranges, or matrix slices. No mixing of different storage layouts allowed at the moment. */
  template <typename NumericT, typename F>
  matrix_expression< const matrix_base<NumericT, F>, const matrix_base<NumericT, F>, op_sub >
  operator - (const matrix_base<NumericT, F> & m1, const matrix_base<NumericT, F> & m2)
  {
    return matrix_expression< const matrix_base<NumericT, F>,
                              const matrix_base<NumericT, F>,
                              op_sub > (m1, m2);
  }



  // operator *
  /** @brief Operator overload for the expression alpha * m1, where alpha is a host scalar (float or double) and m1 is a ViennaCL matrix.
  *
  * @param value   The host scalar (float or double)
  * @param m1      A ViennaCL matrix
  */
  template <typename S1, typename NumericT, typename F>
  typename viennacl::enable_if<    viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_base<NumericT, F>, const S1, op_mult>
                              >::type
  operator * (S1 const & value, matrix_base<NumericT, F> const & m1)
  {
    return matrix_expression< const matrix_base<NumericT, F>, const S1, op_mult>(m1, value);
  }


  /** @brief Operator overload for the multiplication of a matrix expression with a scalar from the right, e.g. (beta * m1) * alpha. Here, beta * m1 is wrapped into a matrix_expression and then multiplied with alpha from the right.
  *
  * @param proxy   Left hand side matrix expression
  * @param val     Right hand side scalar
  */
  template <typename LHS, typename RHS, typename OP, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_expression< LHS, RHS, OP>, const S1, op_mult> >::type
  operator * (matrix_expression< LHS, RHS, OP> const & proxy,
              S1 const & val)
  {
    return matrix_expression< const matrix_expression< LHS, RHS, OP>, const S1, op_mult>(proxy, val);
  }


  /** @brief Operator overload for the multiplication of a matrix expression with a ViennaCL scalar from the left, e.g. alpha * (beta * m1). Here, beta * m1 is wrapped into a matrix_expression and then multiplied with alpha from the left.
  *
  * @param val     Right hand side scalar
  * @param proxy   Left hand side matrix expression
  */
  template <typename S1, typename LHS, typename RHS, typename OP>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_expression< LHS, RHS, OP>, const S1, op_mult> >::type
  operator * (S1 const & val,
              matrix_expression< LHS, RHS, OP> const & proxy)
  {
    return matrix_expression< const matrix_expression< LHS, RHS, OP>, const S1, op_mult>(proxy, val);
  }

  /** @brief Scales the matrix by a GPU scalar 'alpha' and returns an expression template
  */
  template <typename NumericT, typename F, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_base<NumericT, F>, const S1, op_mult> >::type
  operator * (matrix_base<NumericT, F> const & m1, S1 const & s1)
  {
    return matrix_expression< const matrix_base<NumericT, F>, const S1, op_mult>(m1, s1);
  }


  // operator *=

  /** @brief Scales a matrix by a GPU scalar value
  */
  template <typename NumericT, typename F, typename S1>
  typename viennacl::enable_if< viennacl::is_scalar<S1>::value,
                                matrix_base<NumericT, F> &
                              >::type
  operator *= (matrix_base<NumericT, F> & m1, S1 const & gpu_val)
  {
    //viennacl::linalg::inplace_mult(*this, gpu_val);
    viennacl::linalg::am(m1,
                         m1, gpu_val, 1, false, (viennacl::is_flip_sign_scalar<S1>::value ? true : false));
    return m1;
  }


  // operator /


  /** @brief Operator overload for the division of a matrix expression by a scalar from the right, e.g. (beta * m1) / alpha. Here, beta * m1 is wrapped into a matrix_expression and then divided by alpha.
  *
  * @param proxy   Left hand side matrix expression
  * @param val     Right hand side scalar
  */
  template <typename LHS, typename RHS, typename OP, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_expression<const LHS, const RHS, OP>, const S1, op_div> >::type
  operator / (matrix_expression<const LHS, const RHS, OP> const & proxy,
              S1 const & val)
  {
    return matrix_expression< const matrix_expression<const LHS, const RHS, OP>, const S1, op_div>(proxy, val);
  }


  /** @brief Returns an expression template for scaling the matrix by a GPU scalar 'alpha'
  */
  template <typename NumericT, typename F, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                matrix_expression< const matrix_base<NumericT, F>, const S1, op_div> >::type
  operator / (matrix_base<NumericT, F> const & m1, S1 const & s1)
  {
    return matrix_expression< const matrix_base<NumericT, F>, const S1, op_div>(m1, s1);
  }


  // operator /=

  /** @brief Scales a matrix by a GPU scalar value
  */
  template <typename NumericT, typename F, typename S1>
  typename viennacl::enable_if< viennacl::is_scalar<S1>::value,
                                matrix_base<NumericT, F> &
                              >::type
  operator /= (matrix_base<NumericT, F> & m1, S1 const & gpu_val)
  {
    //viennacl::linalg::inplace_divide(*this, gpu_val);
    viennacl::linalg::am(m1,
                         m1, gpu_val, 1, true, (viennacl::is_flip_sign_scalar<S1>::value ? true : false));
    return m1;
  }





  // outer_prod(v1, v2) * val;
  template <typename NumericT, typename S1>
  typename viennacl::enable_if< viennacl::is_scalar<S1>::value,
                                viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                                             const S1,
                                                             op_mult>
                              >::type
  operator*(const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod> & proxy,
            const S1 & val)
  {
    return viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                        const S1,
                                        op_mult>(proxy, val);
  }

  template <typename NumericT, typename S1>
  typename viennacl::enable_if< viennacl::is_cpu_scalar<S1>::value,
                                viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                                              const NumericT,
                                                              op_mult>
                              >::type
  operator*(const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod> & proxy,
            const S1 & val)
  {
    return viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                        const NumericT,
                                        op_mult>(proxy, NumericT(val));
  }

  // val * outer_prod(v1, v2);
  template <typename NumericT, typename S1>
  typename viennacl::enable_if< viennacl::is_scalar<S1>::value,
                                viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                                             const S1,
                                                             op_mult>
                              >::type
  operator*(const S1 & val,
            const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod> & proxy)
  {
    return viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                        const S1,
                                        op_mult>(proxy, val);
  }

  template<typename NumericT, typename S1>
  typename viennacl::enable_if< viennacl::is_cpu_scalar<S1>::value,
                                viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                                             const NumericT,
                                                             op_mult>
                              >::type
  operator*(const S1 & val,
            const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod> & proxy)
  {
    return viennacl::matrix_expression< const viennacl::matrix_expression< const vector_base<NumericT>, const vector_base<NumericT>, op_prod>,
                                        const NumericT,
                                        op_mult>(proxy, NumericT(val));
  }



  //
  // Specify available operations:
  //

  /** \cond */

  namespace linalg
  {
    namespace detail
    {

      // x = y
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_base<T, F> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_base<T, F> const & rhs)
        {
          viennacl::linalg::am(lhs, rhs, T(1), 1, false, false);
        }
      };

      // x += y
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_base<T, F> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_base<T, F> const & rhs)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, rhs, T(1), 1, false, false);
        }
      };

      // x -= y
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_base<T, F> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_base<T, F> const & rhs)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, rhs, T(1), 1, false, true);
        }
      };

      ///////////// x  OP  y * alpha ////////////////////////


      // x = alpha * y
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> const & proxy)
        {
          viennacl::linalg::am(lhs, proxy.lhs(), proxy.rhs(), 1, false, false);
        }
      };

      // x += alpha * y
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> const & proxy)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, false);
        }
      };

      // x -= alpha * y
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult> const & proxy)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, true);
        }
      };


      ///////////// x  OP  vec_expr * alpha ////////////////////////

      // x = alpha * vec_expr
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs = temp * proxy.rhs();
          }
      };

      // x += alpha * vec_expr
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs += temp * proxy.rhs();
          }
      };

      // x -= alpha * vec_expr
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs -= temp * proxy.rhs();
          }
      };


      ///////////// x  OP  y / alpha ////////////////////////

      // x = y / alpha
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::am(lhs, proxy.lhs(), proxy.rhs(), 1, true, false);
        }
      };

      // x += y / alpha
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, true, false);
        }
      };

      // x -= y / alpha
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::ambm(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, true, true);
        }
      };


      ///////////// x  OP  vec_expr / alpha ////////////////////////

      // x = vec_expr / alpha
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs = temp / proxy.rhs();
          }
      };

      // x += vec_expr / alpha
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs += temp / proxy.rhs();
          }
      };

      // x -= vec_expr / alpha
      template <typename T, typename F, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            matrix<T, F> temp(proxy.lhs());
            lhs -= temp / proxy.rhs();
          }
      };



      // generic x = vec_expr1 + vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const LHS, const RHS, op_add> >
      {
        // generic x = vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs = temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_assign, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x = y + z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = alpha * y + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = y / alpha + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = y + beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = y + z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x = alpha * y + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = alpha * y + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x = y / alpha + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = y / alpha + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };

      // dense = sparse * dense
      template <typename T, typename F1, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F1>, op_assign, matrix_expression<const LHS, const RHS, op_prod> >
      {
        template < typename SparseMatrixType, typename F2 >
        static void apply(matrix_base<T, F1> & lhs, matrix_expression<const SparseMatrixType,
                                                                     const viennacl::matrix_base<T, F2>,
                                                                     viennacl::op_prod> const & proxy)
        {
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), lhs);
        }

        // dense = sparse * trans(dense)
        template < typename SparseMatrixType, typename F2 >
        static void apply(matrix_base<T, F1> & lhs, matrix_expression<const SparseMatrixType,
                                                                     const viennacl::matrix_expression< const viennacl::matrix_base<T, F2>,
                                                                                                        const viennacl::matrix_base<T, F2>,
                                                                                                        viennacl::op_trans >,
                                                                     viennacl::op_prod> const & proxy)
        {
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), lhs);
        }

      };

      // generic x += vec_expr1 + vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const LHS, const RHS, op_add> >
      {
        // generic x += vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs += temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_inplace_add, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x += y + z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += alpha * y + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += y / alpha + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += y + beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += y + z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x += alpha * y + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += alpha * y + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x += y / alpha + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += y / alpha + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };



      // generic x -= vec_expr1 + vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const LHS, const RHS, op_add> >
      {
        // generic x -= vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs -= temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_inplace_sub, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x -= y + z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= alpha * y + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= y / alpha + z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= y + beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= y + z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x -= alpha * y + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= alpha * y + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x -= y / alpha + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= y / alpha + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };



      ///////////////////////



      // generic x = vec_expr1 - vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const LHS, const RHS, op_sub> >
      {
        // generic x = vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs = temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_assign, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x = y - z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = alpha * y - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = y / alpha - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = y - beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = y - z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x = alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x = y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };


      // generic x += vec_expr1 - vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const LHS, const RHS, op_sub> >
      {
        // generic x += vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs += temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_inplace_add, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x += y - z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += alpha * y - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += y / alpha - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += y - beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += y - z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x += alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x += y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };



      // generic x -= vec_expr1 - vec_expr2:
      template <typename T, typename F, typename LHS, typename RHS>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const LHS, const RHS, op_sub> >
      {
        // generic x -= vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            matrix_base<T, F> temp(proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs -= temp;
          }
          else
          {
            op_executor<matrix_base<T, F>, op_inplace_sub, LHS>::apply(lhs, proxy.lhs());
            op_executor<matrix_base<T, F>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x -= y - z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= alpha * y - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= y / alpha - z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  const matrix_base<T, F>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= y - beta * z
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= y - z / beta
        template <typename ScalarType>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x -= alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_mult>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x -= y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F>, const ScalarType1, op_div>,
                                                                  const matrix_expression<const matrix_base<T, F>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::ambm_m(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };


      //////////////////// diag(), row(), column() operations ////////////////////////////////////////

      template <typename T, typename F, typename LHS>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const LHS, const int, op_vector_diag> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const vector_base<T>, const int, op_vector_diag> const & proxy)
        {
          viennacl::linalg::matrix_diag_from_vector(proxy.lhs(), proxy.rhs(), lhs);
        }
      };


      template <typename T, typename LHS>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const int, op_matrix_diag> >
      {
        template <typename F>
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const int, op_matrix_diag> const & proxy)
        {
          viennacl::linalg::matrix_diag_to_vector(proxy.lhs(), proxy.rhs(), lhs);
        }
      };

      template <typename T, typename LHS>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const unsigned int, op_row> >
      {
        template <typename F>
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const unsigned int, op_row> const & proxy)
        {
          viennacl::linalg::matrix_row(proxy.lhs(), proxy.rhs(), lhs);
        }
      };


      template <typename T, typename LHS>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const unsigned int, op_column> >
      {
        template <typename F>
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const unsigned int, op_column> const & proxy)
        {
          viennacl::linalg::matrix_column(proxy.lhs(), proxy.rhs(), lhs);
        }
      };


      //////////////////// Element-wise operations ////////////////////////////////////////

      // generic x = mat_expr1 .* mat_expr2:
      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const LHS, const RHS, op_element_binary<OP> > >
      {
        // x = y .* z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          viennacl::linalg::element_op(lhs, proxy);
        }

        // x = y .* mat_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_expression<const LHS2, const RHS2, OP2>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(proxy.lhs(), temp));
        }

        // x = mat_expr .* z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.lhs());
          viennacl::linalg::element_op(lhs, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp, proxy.rhs()));
        }

        // x = mat_expr .* mat_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>,
                                                                  const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp1(proxy.lhs());
          matrix<T, F> temp2(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp1, temp2));
        }
      };

      // generic x += mat_expr .* mat_expr:
      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const LHS, const RHS, op_element_binary<OP> > >
      {
        // x += y .* z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          viennacl::matrix<T, F> temp(proxy);
          lhs += temp;
        }

        // x += y .* mat_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_expression<const LHS2, const RHS2, OP2>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          matrix<T, F> temp2(temp.size1(), temp.size2());
          viennacl::linalg::element_op(temp2, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(proxy.lhs(), temp));
          lhs += temp2;
        }

        // x += mat_expr .* z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.lhs());
          matrix<T, F> temp2(temp.size1(), temp.size2());
          viennacl::linalg::element_op(temp2, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp, proxy.rhs()));
          lhs += temp2;
        }

        // x += mat_expr .* mat_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>,
                                                                  const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp1(proxy.lhs());
          matrix<T, F> temp2(proxy.rhs());
          matrix<T, F> temp3(temp1.size1(), temp1.size2());
          viennacl::linalg::element_op(temp3, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp1, temp2));
          lhs += temp3;
        }
      };

      // generic x -= mat_expr1 .* mat_expr2:
      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const LHS, const RHS, op_element_binary<OP> > >
      {

        // x -= y .* z
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          viennacl::matrix<T, F> temp(proxy);
          lhs -= temp;
        }

        // x -= y .* mat_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_expression<const LHS2, const RHS2, OP2>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          matrix<T, F> temp2(temp.size1(), temp.size2());
          viennacl::linalg::element_op(temp2, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(proxy.lhs(), temp));
          lhs -= temp2;
        }

        // x -= mat_expr .* z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>, const matrix_base<T, F>, op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.lhs());
          matrix<T, F> temp2(temp.size1(), temp.size2());
          viennacl::linalg::element_op(temp2, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp, proxy.rhs()));
          lhs -= temp2;
        }

        // x -= mat_expr .* mat_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS1, const RHS1, OP1>,
                                                                     const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                     op_element_binary<OP> > const & proxy)
        {
          matrix<T, F> temp1(proxy.lhs());
          matrix<T, F> temp2(proxy.rhs());
          matrix<T, F> temp3(temp1.size1(), temp1.size2());
          viennacl::linalg::element_op(temp3, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_binary<OP> >(temp1, temp2));
          lhs -= temp3;
        }
      };

      //////////////// unary expressions

      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x = OP(y)
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> > const & proxy)
        {
          viennacl::linalg::element_op(lhs, proxy);
        }

        // x = OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                     const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                     op_element_unary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> >(temp, temp));
        }
      };

      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x += OP(y)
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy);
          lhs += temp;
        }

        // x += OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                  const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_unary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          viennacl::linalg::element_op(temp, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> >(temp, temp)); // inplace operation is safe here
          lhs += temp;
        }
      };

      template <typename T, typename F, typename LHS, typename RHS, typename OP>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x -= OP(y)
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy);
          lhs -= temp;
        }

        // x -= OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                     const matrix_expression<const LHS2, const RHS2, OP2>,
                                                                     op_element_unary<OP> > const & proxy)
        {
          matrix<T, F> temp(proxy.rhs());
          viennacl::linalg::element_op(temp, viennacl::matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_element_unary<OP> >(temp, temp)); // inplace operation is safe here
          lhs -= temp;
        }
      };



      //////////////// Matrix - Matrix products ////////////////

      // C = A * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(0));
        }
      };

      // C = A * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_base<T, F1>,
                                                                         const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                         op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(0));
        }
      };

      // C = A^T * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                         const matrix_base<T, F2>,
                                                                         op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_base<T, F2>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(0));
        }
      };

      // C = A^T * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                         const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                         op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(0));
        }
      };


      // C += A * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(1.0));
        }
      };

      // C += A * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_base<T, F1>,
                                                                              const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(1.0));
        }
      };

      // C += A^T * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                              const matrix_base<T, F2>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_base<T, F2>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(1.0));
        }
      };

      // C += A^T * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                              const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(1.0), T(1.0));
        }
      };


      // C -= A * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F2>, op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(-1.0), T(1.0));
        }
      };

      // C -= A * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_base<T, F1>,
                                                                              const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_base<T, F1>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(-1.0), T(1.0));
        }
      };

      // C -= A^T * B
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                              const matrix_base<T, F2>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_base<T, F2>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(-1.0), T(1.0));
        }
      };

      // C -= A^T * B^T
      template <typename T, typename F, typename F1, typename F2>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                              const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                              op_mat_mat_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const matrix_expression<const matrix_base<T, F1>, const matrix_base<T, F1>, op_trans>,
                                                                     const matrix_expression<const matrix_base<T, F2>, const matrix_base<T, F2>, op_trans>,
                                                                     op_mat_mat_prod> const & rhs)
        {
          viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs, T(-1.0), T(1.0));
        }
      };

      ////////////////// Matrix-Vector Products ///////////////

      // y = A * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> const & rhs)
        {
          // check for x = A * x
          if (op_aliasing(lhs, rhs.rhs()))
          {
            vector_base<T> temp(rhs);
            lhs = temp;
          }
          else
            viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs);
        }
      };

      // y = A^T * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                      const vector_base<T>,
                                                                      op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                  const vector_base<T>,
                                                                  op_prod> const & rhs)
        {
          // check for x = A^T * x
          if (op_aliasing(lhs, rhs.rhs()))
          {
            vector_base<T> temp(rhs);
            lhs = temp;
          }
          else
            viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs);
        }
      };


      // y += A * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> const & rhs)
        {
          vector_base<T> temp(rhs);
          lhs += temp;
        }
      };

      // y += A^T * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                           const vector_base<T>,
                                                                           op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                  const vector_base<T>,
                                                                  op_prod> const & rhs)
        {
          vector_base<T> temp(rhs);
          lhs += temp;
        }
      };


      // y -= A * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_base<T, F>, const vector_base<T>, op_prod> const & rhs)
        {
          vector_base<T> temp(rhs);
          lhs -= temp;
        }
      };

      // y -= A^T * x
      template <typename T, typename F>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                           const vector_base<T>,
                                                                           op_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const matrix_expression<const matrix_base<T, F>, const matrix_base<T, F>, op_trans>,
                                                                  const vector_base<T>,
                                                                  op_prod> const & rhs)
        {
          vector_base<T> temp(rhs);
          lhs -= temp;
        }
      };



      ////////////////// Rank-1 Updates ///////////////

      // A = v1 * v2^T
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> const & rhs)
        {
          lhs.clear();
          viennacl::linalg::scaled_rank_1_update(lhs, T(1.0), 1, false, false, rhs.lhs(), rhs.rhs());
        }
      };

      // A = alpha * v1 * v2^T
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_assign, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                          const ScalarType,
                                                                          op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                      const ScalarType,
                                                                      op_mult> const & rhs)
        {
          lhs.clear();
          viennacl::linalg::scaled_rank_1_update(lhs, rhs.rhs(), 1, false, false, rhs.lhs().lhs(), rhs.lhs().rhs());
        }
      };

      // A += v1 * v2^T
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> const & rhs)
        {
          viennacl::linalg::scaled_rank_1_update(lhs, T(1.0), 1, false, false, rhs.lhs(), rhs.rhs());
        }
      };

      // A += alpha * v1 * v2^T
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_add, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                               const ScalarType,
                                                                               op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                      const ScalarType,
                                                                      op_mult> const & rhs)
        {
          viennacl::linalg::scaled_rank_1_update(lhs, rhs.rhs(), 1, false, false, rhs.lhs().lhs(), rhs.lhs().rhs());
        }
      };

      // A -= v1 * v2^T
      template <typename T, typename F>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression<const vector_base<T>, const vector_base<T>, op_prod> const & rhs)
        {
          viennacl::linalg::scaled_rank_1_update(lhs, T(1.0), 1, false, true, rhs.lhs(), rhs.rhs());
        }
      };

      // A -= alpha * v1 * v2^T
      template <typename T, typename F, typename ScalarType>
      struct op_executor<matrix_base<T, F>, op_inplace_sub, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                               const ScalarType,
                                                                               op_mult> >
      {
        static void apply(matrix_base<T, F> & lhs, matrix_expression< const matrix_expression<const vector_base<T>, const vector_base<T>, op_prod>,
                                                                      const ScalarType,
                                                                      op_mult> const & rhs)
        {
          viennacl::linalg::scaled_rank_1_update(lhs, rhs.rhs(), 1, false, true, rhs.lhs().lhs(), rhs.lhs().rhs());
        }
      };


    } // namespace detail

  } // namespace linalg

  /** \endcond */

} //namespace viennacl

#endif
