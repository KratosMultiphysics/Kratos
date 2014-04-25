#ifndef VIENNACL_VECTOR_HPP_
#define VIENNACL_VECTOR_HPP_

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

/** @file  viennacl/vector.hpp
    @brief The vector type with operator-overloads and proxy classes is defined here.
           Linear algebra operations such as norms and inner products are located in linalg/vector_operations.hpp
*/


#include "viennacl/forwards.h"
#include "viennacl/backend/memory.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/entry_proxy.hpp"
#include "viennacl/linalg/detail/op_executor.hpp"
#include "viennacl/linalg/vector_operations.hpp"
#include "viennacl/meta/result_of.hpp"
//#include "viennacl/rand/utils.hpp"
#include "viennacl/context.hpp"
#include "viennacl/traits/handle.hpp"

namespace viennacl
{

  /** @brief Common base class for representing vectors where the entries are not all stored explicitly.
    *
    * Typical examples are zero_vector or scalar_vector.
    */
  template<typename SCALARTYPE>
  class implicit_vector_base
  {
    protected:
      typedef vcl_size_t        size_type;
      implicit_vector_base(size_type s, vcl_size_t i, std::pair<SCALARTYPE, bool> v, viennacl::context ctx) : size_(s), index_(std::make_pair(true,i)), value_(v), ctx_(ctx){ }
      implicit_vector_base(size_type s, std::pair<SCALARTYPE, bool> v, viennacl::context ctx) : size_(s), index_(std::make_pair(false,0)), value_(v), ctx_(ctx){ }

    public:
      typedef SCALARTYPE const & const_reference;
      typedef SCALARTYPE cpu_value_type;

      viennacl::context context() const { return ctx_; }

      size_type size() const { return size_; }

      cpu_value_type  value() const { return value_.first; }

      bool is_value_static() const { return value_.second; }

      vcl_size_t index() const { return index_.second; }

      bool has_index() const { return index_.first; }

      cpu_value_type operator()(size_type i) const {
        if(index_.first)
          return (i==index_.second)?value_.first:0;
        return value_.first;
      }

      cpu_value_type operator[](size_type i) const {
        if(index_.first)
          return (i==index_.second)?value_.first:0;
        return
            value_.first;
      }

    protected:
      size_type size_;
      std::pair<bool, vcl_size_t> index_;
      std::pair<SCALARTYPE, bool> value_;
      viennacl::context ctx_;
  };

  /** @brief Represents a vector consisting of 1 at a given index and zeros otherwise.*/
  template <typename SCALARTYPE>
  class unit_vector : public implicit_vector_base<SCALARTYPE>
  {
      typedef implicit_vector_base<SCALARTYPE> base_type;
    public:
      typedef typename base_type::size_type size_type;
      unit_vector(size_type s, size_type ind, viennacl::context ctx = viennacl::context()) : base_type(s, ind, std::make_pair(SCALARTYPE(1),true), ctx)
      {
        assert( (ind < s) && bool("Provided index out of range!") );
      }
  };


  /** @brief Represents a vector consisting of zeros only. */
  template <typename SCALARTYPE>
  class zero_vector : public implicit_vector_base<SCALARTYPE>
  {
      typedef implicit_vector_base<SCALARTYPE> base_type;
    public:
      typedef typename base_type::size_type size_type;
      typedef SCALARTYPE        const_reference;
      zero_vector(size_type s, viennacl::context ctx = viennacl::context()) : base_type(s, std::make_pair(SCALARTYPE(0),true), ctx) {}
  };

  /** @brief Represents a vector consisting of ones only. */
  template <typename SCALARTYPE>
  class one_vector : public implicit_vector_base<SCALARTYPE>
  {
      typedef implicit_vector_base<SCALARTYPE> base_type;
    public:
      typedef typename base_type::size_type size_type;
      typedef SCALARTYPE        const_reference;
      one_vector(size_type s, viennacl::context ctx = viennacl::context()) : base_type(s, std::make_pair(SCALARTYPE(1),true), ctx) {}
  };


  /** @brief Represents a vector consisting of scalars 's' only, i.e. v[i] = s for all i. To be used as an initializer for viennacl::vector, vector_range, or vector_slize only. */
  template <typename SCALARTYPE>
  class scalar_vector : public implicit_vector_base<SCALARTYPE>
  {
      typedef implicit_vector_base<SCALARTYPE> base_type;
    public:
      typedef typename base_type::size_type size_type;
      typedef SCALARTYPE const & const_reference;

      scalar_vector(size_type s, SCALARTYPE val, viennacl::context ctx = viennacl::context()) : base_type(s, std::make_pair(val,false), ctx) {}
  };


//#ifdef VIENNACL_WITH_OPENCL
//  template<class SCALARTYPE, class DISTRIBUTION>
//  rand::random_vector_t<SCALARTYPE, DISTRIBUTION> random_vector(unsigned int size, DISTRIBUTION const & distribution){
//      return rand::random_vector_t<SCALARTYPE,DISTRIBUTION>(size,distribution);
//  }
//#endif


  //
  // Vector expression
  //

  /** @brief An expression template class that represents a binary operation that yields a vector
  *
  * In contrast to full expression templates as introduced by Veldhuizen, ViennaCL does not allow nested expressions.
  * The reason is that this requires automated GPU viennacl::ocl::kernel generation, which then has to be compiles just-in-time.
  * For performance-critical applications, one better writes the appropriate viennacl::ocl::kernels by hand.
  *
  * Assumption: dim(LHS) >= dim(RHS), where dim(scalar) = 0, dim(vector) = 1 and dim(matrix = 2)
  *
  * @tparam LHS   left hand side operand
  * @tparam RHS   right hand side operand
  * @tparam OP    the operator
  */
  template <typename LHS, typename RHS, typename OP>
  class vector_expression
  {
      typedef typename viennacl::result_of::reference_if_nonscalar<LHS>::type     lhs_reference_type;
      typedef typename viennacl::result_of::reference_if_nonscalar<RHS>::type     rhs_reference_type;

    public:
      enum { alignment = 1 };

      /** @brief Extracts the vector type from the two operands.
      */
      typedef vcl_size_t       size_type;

      vector_expression(LHS & l, RHS & r) : lhs_(l), rhs_(r) {}

      /** @brief Get left hand side operand
      */
      lhs_reference_type lhs() const { return lhs_; }
      /** @brief Get right hand side operand
      */
      rhs_reference_type rhs() const { return rhs_; }

      /** @brief Returns the size of the result vector */
      size_type size() const { return viennacl::traits::size(*this); }

    private:
      /** @brief The left hand side operand */
      lhs_reference_type lhs_;
      /** @brief The right hand side operand */
      rhs_reference_type rhs_;
  };

  /** @brief A STL-type const-iterator for vector elements. Elements can be accessed, but cannot be manipulated. VERY SLOW!!
  *
  * Every dereference operation initiates a transfer from the GPU to the CPU. The overhead of such a transfer is around 50us, so 20.000 dereferences take one second.
  * This is four orders of magnitude slower than similar dereferences on the CPU. However, increments and comparisons of iterators is as fast as for CPU types.
  * If you need a fast iterator, copy the whole vector to the CPU first and iterate over the CPU object, e.g.
  * std::vector<float> temp;
  * copy(gpu_vector, temp);
  * for (std::vector<float>::const_iterator iter = temp.begin();
  *      iter != temp.end();
  *      ++iter)
  * {
  *   //do something
  * }
  * Note that you may obtain inconsistent data if entries of gpu_vector are manipulated elsewhere in the meanwhile.
  *
  * @tparam SCALARTYPE  The underlying floating point type (either float or double)
  * @tparam ALIGNMENT   Alignment of the underlying vector, @see vector
  */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class const_vector_iterator
  {
      typedef const_vector_iterator<SCALARTYPE, ALIGNMENT>    self_type;
    public:
      typedef scalar<SCALARTYPE>            value_type;
      typedef vcl_ptrdiff_t                 difference_type;
      typedef viennacl::backend::mem_handle handle_type;

      //const_vector_iterator() {}

      /** @brief Constructor
      *   @param vec    The vector over which to iterate
      *   @param index  The starting index of the iterator
      *   @param start  First index of the element in the vector pointed to be the iterator (for vector_range and vector_slice)
      *   @param stride Stride for the support of vector_slice
      */
      const_vector_iterator(vector_base<SCALARTYPE> const & vec,
                            vcl_size_t index,
                            vcl_size_t start = 0,
                            vcl_ptrdiff_t stride = 1) : elements_(vec.handle()), index_(index), start_(start), stride_(stride) {}

      /** @brief Constructor for vector-like treatment of arbitrary buffers
      *   @param elements  The buffer over which to iterate
      *   @param index     The starting index of the iterator
      *   @param start     First index of the element in the vector pointed to be the iterator (for vector_range and vector_slice)
      *   @param stride    Stride for the support of vector_slice
      */
      const_vector_iterator(handle_type const & elements,
                            vcl_size_t index,
                            vcl_size_t start = 0,
                            vcl_ptrdiff_t stride = 1) : elements_(elements), index_(index), start_(start), stride_(stride) {}

      /** @brief Dereferences the iterator and returns the value of the element. For convenience only, performance is poor due to OpenCL overhead! */
      value_type operator*(void) const
      {
          value_type result;
          result = const_entry_proxy<SCALARTYPE>(start_ + index_ * stride_, elements_);
          return result;
      }
      self_type operator++(void) { index_ += stride_; return *this; }
      self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }

      bool operator==(self_type const & other) const { return index_ == other.index_; }
      bool operator!=(self_type const & other) const { return index_ != other.index_; }

//        self_type & operator=(self_type const & other)
//        {
//           index_ = other._index;
//           elements_ = other._elements;
//           return *this;
//        }

      difference_type operator-(self_type const & other) const
      {
        assert( (other.start_ == start_) && (other.stride_ == stride_) && bool("Iterators are not from the same vector (proxy)!"));
        return static_cast<difference_type>(index_) - static_cast<difference_type>(other.index_);
      }
      self_type operator+(difference_type diff) const { return self_type(elements_, index_ + diff * stride_, start_, stride_); }

      //vcl_size_t index() const { return index_; }
      /** @brief Offset of the current element index with respect to the beginning of the buffer */
      vcl_size_t offset() const { return start_ + index_ * stride_; }

      /** @brief Index increment in the underlying buffer when incrementing the iterator to the next element */
      vcl_size_t stride() const { return stride_; }
      handle_type const & handle() const { return elements_; }

    protected:
      /** @brief  The index of the entry the iterator is currently pointing to */
      handle_type const & elements_;
      vcl_size_t index_;  //offset from the beginning of elements_
      vcl_size_t start_;
      vcl_ptrdiff_t stride_;
  };


  /** @brief A STL-type iterator for vector elements. Elements can be accessed and manipulated. VERY SLOW!!
  *
  * Every dereference operation initiates a transfer from the GPU to the CPU. The overhead of such a transfer is around 50us, so 20.000 dereferences take one second.
  * This is four orders of magnitude slower than similar dereferences on the CPU. However, increments and comparisons of iterators is as fast as for CPU types.
  * If you need a fast iterator, copy the whole vector to the CPU first and iterate over the CPU object, e.g.
  * std::vector<float> temp;
  * copy(gpu_vector, temp);
  * for (std::vector<float>::const_iterator iter = temp.begin();
  *      iter != temp.end();
  *      ++iter)
  * {
  *   //do something
  * }
  * copy(temp, gpu_vector);
  * Note that you may obtain inconsistent data if you manipulate entries of gpu_vector in the meanwhile.
  *
  * @tparam SCALARTYPE  The underlying floating point type (either float or double)
  * @tparam ALIGNMENT   Alignment of the underlying vector, @see vector
  */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class vector_iterator : public const_vector_iterator<SCALARTYPE, ALIGNMENT>
  {
      typedef const_vector_iterator<SCALARTYPE, ALIGNMENT>  base_type;
      typedef vector_iterator<SCALARTYPE, ALIGNMENT>        self_type;
    public:
      typedef typename base_type::handle_type               handle_type;
      typedef typename base_type::difference_type           difference_type;

      vector_iterator() : base_type(), elements_(NULL) {}
      vector_iterator(handle_type & elements,
                      vcl_size_t index,
                      vcl_size_t start = 0,
                      vcl_ptrdiff_t stride = 1)  : base_type(elements, index, start, stride), elements_(elements) {}
      /** @brief Constructor
      *   @param vec    The vector over which to iterate
      *   @param index  The starting index of the iterator
      *   @param start  Offset from the beginning of the underlying vector (for ranges and slices)
      *   @param stride Stride for slices
      */
      vector_iterator(vector_base<SCALARTYPE> & vec,
                      vcl_size_t index,
                      vcl_size_t start = 0,
                      vcl_ptrdiff_t stride = 1) : base_type(vec, index, start, stride), elements_(vec.handle()) {}
      //vector_iterator(base_type const & b) : base_type(b) {}

      typename base_type::value_type operator*(void)
      {
          typename base_type::value_type result;
          result = entry_proxy<SCALARTYPE>(base_type::start_ + base_type::index_ * base_type::stride_, elements_);
          return result;
      }

      difference_type operator-(self_type const & other) const { difference_type result = base_type::index_; return (result - static_cast<difference_type>(other.index_)); }
      self_type operator+(difference_type diff) const { return self_type(elements_, base_type::index_ + diff * base_type::stride_, base_type::start_, base_type::stride_); }

      handle_type       & handle()       { return elements_; }
      handle_type const & handle() const { return base_type::elements_; }

      //operator base_type() const
      //{
      //  return base_type(base_type::elements_, base_type::index_, base_type::start_, base_type::stride_);
      //}
    private:
      handle_type & elements_;
  };


  /** @brief Common base class for dense vectors, vector ranges, and vector slices.
    *
    * @tparam SCALARTYPE   The floating point type, either 'float' or 'double'
    */
  template<class SCALARTYPE, typename SizeType /* see forwards.h for default type */, typename DistanceType /* see forwards.h for default type */>
  class vector_base
  {
      typedef vector_base<SCALARTYPE>         self_type;

    public:
      typedef scalar<SCALARTYPE>                                value_type;
      typedef SCALARTYPE                                        cpu_value_type;
      typedef viennacl::backend::mem_handle                     handle_type;
      typedef SizeType                                          size_type;
      typedef DistanceType                                      difference_type;
      typedef const_vector_iterator<SCALARTYPE, 1>              const_iterator;
      typedef vector_iterator<SCALARTYPE, 1>                    iterator;

      static const size_type alignment = 128;

      /** @brief Default constructor in order to be compatible with various containers.
      */
      explicit vector_base() : size_(0), start_(0), stride_(1), internal_size_(0) { /* Note: One must not call ::init() here because a vector might have been created globally before the backend has become available */ }

      /** @brief An explicit constructor for wrapping an existing vector into a vector_range or vector_slice.
       *
       *
       *
       * @param h          The existing memory handle from a vector/vector_range/vector_slice
       * @param vec_size   The length (i.e. size) of the buffer
       * @param vec_start  The offset from the beginning of the buffer identified by 'h'
       * @param vec_stride Increment between two elements in the original buffer (in multiples of SCALARTYPE)
      */
      explicit vector_base(viennacl::backend::mem_handle & h,
                           size_type vec_size, size_type vec_start, difference_type vec_stride)
        : size_(vec_size), start_(vec_start), stride_(vec_stride), internal_size_(vec_size), elements_(h) {}

      /** @brief Creates a vector and allocates the necessary memory */
      explicit vector_base(size_type vec_size, viennacl::context ctx = viennacl::context())
        : size_(vec_size), start_(0), stride_(1), internal_size_(viennacl::tools::align_to_multiple<size_type>(size_, alignment))
      {
        if (size_ > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), ctx);
          clear();
        }
      }

      // CUDA or host memory:
      explicit vector_base(SCALARTYPE * ptr_to_mem, viennacl::memory_types mem_type, size_type vec_size, vcl_size_t start = 0, difference_type stride = 1)
        : size_(vec_size), start_(start), stride_(stride), internal_size_(vec_size)
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

        elements_.raw_size(sizeof(SCALARTYPE) * vec_size);

      }

#ifdef VIENNACL_WITH_OPENCL
      /** @brief Create a vector from existing OpenCL memory
      *
      * Note: The provided memory must take an eventual ALIGNMENT into account, i.e. existing_mem must be at least of size internal_size()!
      * This is trivially the case with the default alignment, but should be considered when using vector<> with an alignment parameter not equal to 1.
      *
      * @param existing_mem   An OpenCL handle representing the memory
      * @param vec_size       The size of the vector.
      */
      explicit vector_base(cl_mem existing_mem, size_type vec_size, size_type start = 0, difference_type stride = 1, viennacl::context ctx = viennacl::context())
        : size_(vec_size), start_(start), stride_(stride), internal_size_(vec_size)
      {
        elements_.switch_active_handle_id(viennacl::OPENCL_MEMORY);
        elements_.opencl_handle() = existing_mem;
        elements_.opencl_handle().inc();  //prevents that the user-provided memory is deleted once the vector object is destroyed.
        elements_.opencl_handle().context(ctx.opencl_context());
        elements_.raw_size(sizeof(SCALARTYPE) * vec_size);
      }
#endif

      /** @brief Creates the vector from the supplied random vector. */
      /*template<class DISTRIBUTION>
      vector(rand::random_vector_t<SCALARTYPE, DISTRIBUTION> v) : size_(v.size)
      {
        if(size_ > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size());
          rand::buffer_dumper<SCALARTYPE, DISTRIBUTION>::dump(elements_,v.distribution,0,size_);
        }
      } */

      template <typename LHS, typename RHS, typename OP>
      explicit vector_base(vector_expression<const LHS, const RHS, OP> const & proxy)
        : size_(viennacl::traits::size(proxy)), start_(0), stride_(1), internal_size_(viennacl::tools::align_to_multiple<size_type>(size_, alignment))
      {
        if (size_ > 0)
        {
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(proxy));
          clear();
        }
        self_type::operator=(proxy);
      }


      //
      // operator=
      //


      /** @brief Assignment operator. Other vector needs to be of the same size, or this vector is not yet initialized.
      */
      self_type & operator=(const self_type & vec)
      {
        assert( ( (vec.size() == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        if (vec.size() > 0)
        {
          if (size_ == 0)
          {
            size_ = vec.size();
            internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
            elements_.switch_active_handle_id(vec.handle().get_active_handle_id());
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(vec));
            pad();
          }

          viennacl::linalg::av(*this,
                               vec, cpu_value_type(1.0), 1, false, false);
        }

        return *this;
      }


      /** @brief Implementation of the operation v1 = v2 @ alpha, where @ denotes either multiplication or division, and alpha is either a CPU or a GPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename LHS, typename RHS, typename OP>
      self_type & operator=(const vector_expression<const LHS, const RHS, OP> & proxy)
      {
        assert( ( (viennacl::traits::size(proxy) == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        // initialize the necessary buffer
        if (size() == 0)
        {
          size_ = viennacl::traits::size(proxy);
          internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
          viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(proxy));
          pad();
        }

        linalg::detail::op_executor<self_type, op_assign, vector_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }

      // assign vector range or vector slice
      template <typename T>
      self_type &
      operator = (const vector_base<T> & v1)
      {
        assert( ( (v1.size() == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        if (size() == 0)
        {
          size_ = v1.size();
          if (size_ > 0)
          {
            internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), viennacl::traits::context(v1));
            pad();
          }
        }

        viennacl::linalg::av(*this,
                             v1, SCALARTYPE(1.0), 1, false, false);

        return *this;
      }

      /** @brief Creates the vector from the supplied unit vector. */
      self_type & operator = (unit_vector<SCALARTYPE> const & v)
      {
        assert( ( (v.size() == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        if (size() == 0)
        {
          size_ = v.size();
          internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
          if (size_ > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), v.context());
            clear();
          }
        }
        else
          viennacl::linalg::vector_assign(*this, SCALARTYPE(0));

        if (size_ > 0)
          this->operator()(v.index()) = SCALARTYPE(1);

        return *this;
      }

      /** @brief Creates the vector from the supplied zero vector. */
      self_type & operator = (zero_vector<SCALARTYPE> const & v)
      {
        assert( ( (v.size() == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        if (size() == 0)
        {
          size_ = v.size();
          internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
          if (size_ > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), v.context());
            clear();
          }
        }
        else
          viennacl::linalg::vector_assign(*this, SCALARTYPE(0));

        return *this;
      }

      /** @brief Creates the vector from the supplied scalar vector. */
      self_type & operator = (scalar_vector<SCALARTYPE> const & v)
      {
        assert( ( (v.size() == size()) || (size() == 0) )
                && bool("Incompatible vector sizes!"));

        if (size() == 0)
        {
          size_ = v.size();
          internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
          if (size_ > 0)
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*internal_size(), v.context());
            pad();
          }
        }

        if (size_ > 0)
          viennacl::linalg::vector_assign(*this, v[0]);

        return *this;
      }



      ///////////////////////////// Matrix Vector interaction start ///////////////////////////////////

      //Note: The following operator overloads are defined in matrix_operations.hpp, compressed_matrix_operations.hpp and coordinate_matrix_operations.hpp
      //This is certainly not the nicest approach and will most likely by changed in the future, but it works :-)

      //matrix<>
      /** @brief Operator overload for v1 = A * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F>
      self_type & operator=(const viennacl::vector_expression< const matrix_base<SCALARTYPE, F>, const vector_base<SCALARTYPE>, viennacl::op_prod> & proxy)
      {
        assert(viennacl::traits::size1(proxy.lhs()) == size() && bool("Size check failed for v1 = A * v2: size1(A) != size(v1)"));

        // check for the special case x = A * x
        if (viennacl::traits::handle(proxy.rhs()) == viennacl::traits::handle(*this))
        {
          viennacl::vector<SCALARTYPE> result(viennacl::traits::size1(proxy.lhs()));
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
          *this = result;
        }
        else
        {
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        }
        return *this;
      }


      //transposed_matrix_proxy:
      /** @brief Operator overload for v1 = trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F>
      self_type & operator=(const vector_expression< const matrix_expression< const matrix_base<SCALARTYPE, F>, const matrix_base<SCALARTYPE, F>, op_trans >,
                                                     const vector_base<SCALARTYPE>,
                                                     op_prod> & proxy)
      {
        assert(viennacl::traits::size1(proxy.lhs()) == size() && bool("Size check failed in v1 = trans(A) * v2: size2(A) != size(v1)"));

        // check for the special case x = trans(A) * x
        if (viennacl::traits::handle(proxy.rhs()) == viennacl::traits::handle(*this))
        {
          viennacl::vector<SCALARTYPE> result(viennacl::traits::size1(proxy.lhs()));
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), result);
          *this = result;
        }
        else
        {
          viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
        }
        return *this;
      }

      ///////////////////////////// Matrix Vector interaction end ///////////////////////////////////


      //read-write access to an element of the vector
      /** @brief Read-write access to a single element of the vector
      */
      entry_proxy<SCALARTYPE> operator()(size_type index)
      {
        assert( (size() > 0)  && bool("Cannot apply operator() to vector of size zero!"));
        assert( index < size() && bool("Index out of bounds!") );

        return entry_proxy<SCALARTYPE>(start_ + stride_ * index, elements_);
      }

      /** @brief Read-write access to a single element of the vector
      */
      entry_proxy<SCALARTYPE> operator[](size_type index)
      {
        assert( (size() > 0)  && bool("Cannot apply operator() to vector of size zero!"));
        assert( index < size() && bool("Index out of bounds!") );

        return entry_proxy<SCALARTYPE>(start_ + stride_ * index, elements_);
      }


      /** @brief Read access to a single element of the vector
      */
      const_entry_proxy<SCALARTYPE> operator()(size_type index) const
      {
        assert( (size() > 0)  && bool("Cannot apply operator() to vector of size zero!"));
        assert( index < size() && bool("Index out of bounds!") );

        return const_entry_proxy<SCALARTYPE>(start_ + stride_ * index, elements_);
      }

      /** @brief Read access to a single element of the vector
      */
      const_entry_proxy<SCALARTYPE> operator[](size_type index) const
      {
        assert( (size() > 0)  && bool("Cannot apply operator() to vector of size zero!"));
        assert( index < size() && bool("Index out of bounds!") );

        return const_entry_proxy<SCALARTYPE>(start_ + stride_ * index, elements_);
      }

      //
      // Operator overloads with implicit conversion (thus cannot be made global without introducing additional headache)
      //
      self_type & operator += (const self_type & vec)
      {
        assert(vec.size() == size() && bool("Incompatible vector sizes!"));

        if (size() > 0)
          viennacl::linalg::avbv(*this,
                                  *this, SCALARTYPE(1.0), 1, false, false,
                                  vec,   SCALARTYPE(1.0), 1, false, false);
        return *this;
      }

      self_type & operator -= (const self_type & vec)
      {
        assert(vec.size() == size() && bool("Incompatible vector sizes!"));

        if (size() > 0)
          viennacl::linalg::avbv(*this,
                                  *this, SCALARTYPE(1.0),  1, false, false,
                                  vec,   SCALARTYPE(-1.0), 1, false, false);
        return *this;
      }

      template <typename LHS, typename RHS, typename OP>
      self_type & operator += (const vector_expression<const LHS, const RHS, OP> & proxy)
      {
        assert( (viennacl::traits::size(proxy) == size()) && bool("Incompatible vector sizes!"));
        assert( (size() > 0) && bool("Vector not yet initialized!") );

        linalg::detail::op_executor<self_type, op_inplace_add, vector_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }

      template <typename LHS, typename RHS, typename OP>
      self_type & operator -= (const vector_expression<const LHS, const RHS, OP> & proxy)
      {
        assert( (viennacl::traits::size(proxy) == size()) && bool("Incompatible vector sizes!"));
        assert( (size() > 0) && bool("Vector not yet initialized!") );

        linalg::detail::op_executor<self_type, op_inplace_sub, vector_expression<const LHS, const RHS, OP> >::apply(*this, proxy);

        return *this;
      }

      /** @brief Scales a vector (or proxy) by a CPU scalar value
      */
      self_type & operator *= (SCALARTYPE val)
      {
        if (size() > 0)
          viennacl::linalg::av(*this,
                                *this, val, 1, false, false);
        return *this;
      }

      /** @brief Scales this vector by a CPU scalar value
      */
      self_type & operator /= (SCALARTYPE val)
      {
        if (size() > 0)
          viennacl::linalg::av(*this,
                               *this, val, 1, true, false);
        return *this;
      }


      /** @brief Scales the vector by a CPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const self_type, const SCALARTYPE, op_mult>
      operator * (SCALARTYPE value) const
      {
        return vector_expression< const self_type, const SCALARTYPE, op_mult>(*this, value);
      }


      /** @brief Scales the vector by a CPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const self_type, const SCALARTYPE, op_div>
      operator / (SCALARTYPE value) const
      {
        return vector_expression< const self_type, const SCALARTYPE, op_div>(*this, value);
      }


      /** @brief Sign flip for the vector. Emulated to be equivalent to -1.0 * vector */
      vector_expression<const self_type, const SCALARTYPE, op_mult> operator-() const
      {
        return vector_expression<const self_type, const SCALARTYPE, op_mult>(*this, SCALARTYPE(-1.0));
      }

      //
      //// iterators:
      //

      /** @brief Returns an iterator pointing to the beginning of the vector  (STL like)*/
      iterator begin()
      {
        return iterator(*this, 0, start_, stride_);
      }

      /** @brief Returns an iterator pointing to the end of the vector (STL like)*/
      iterator end()
      {
        return iterator(*this, size(), start_, stride_);
      }

      /** @brief Returns a const-iterator pointing to the beginning of the vector (STL like)*/
      const_iterator begin() const
      {
        return const_iterator(*this, 0, start_, stride_);
      }

      /** @brief Returns a const-iterator pointing to the end of the vector (STL like)*/
      const_iterator end() const
      {
        return const_iterator(*this, size(), start_, stride_);
      }

      /** @brief Swaps the entries of the two vectors
      */
      self_type & swap(self_type & other)
      {
        viennacl::linalg::vector_swap(*this, other);
        return *this;
      };


      /** @brief Returns the length of the vector (cf. std::vector)
      */
      size_type size() const { return size_; }

      /** @brief Returns the internal length of the vector, which is given by size() plus the extra memory due to padding the memory with zeros up to a multiple of 'ALIGNMENT'
      */
      size_type internal_size() const { return internal_size_; }

      /** @brief Returns the offset within the buffer
      */
      size_type start() const { return start_; }

      /** @brief Returns the stride within the buffer (in multiples of sizeof(SCALARTYPE))
      */
      size_type stride() const { return stride_; }


      /** @brief Returns true is the size is zero */
      bool empty() const { return size_ == 0; }

      /** @brief Returns the memory handle. */
      const handle_type & handle() const { return elements_; }

      /** @brief Returns the memory handle. */
      handle_type & handle() { return elements_; }

      /** @brief Resets all entries to zero. Does not change the size of the vector.
      */
      void clear()
      {
        viennacl::linalg::vector_assign(*this, cpu_value_type(0.0), true);
      }

      viennacl::memory_types memory_domain() const
      {
        return elements_.get_active_handle_id();
      }

    protected:

      void set_handle(viennacl::backend::mem_handle const & h)
      {
        elements_ = h;
      }

      /** @brief Swaps the handles of two vectors by swapping the OpenCL handles only, no data copy
      */
      self_type & fast_swap(self_type & other)
      {
        assert(this->size_ == other.size_ && bool("Vector size mismatch"));
        this->elements_.swap(other.elements_);
        return *this;
      }

      /** @brief Pads vectors with alignment > 1 with trailing zeros if the internal size is larger than the visible size */
      void pad()
      {
        if (internal_size() != size())
        {
          std::vector<SCALARTYPE> pad(internal_size() - size());
          viennacl::backend::memory_write(elements_, sizeof(SCALARTYPE) * size(), sizeof(SCALARTYPE) * pad.size(), &(pad[0]));
        }
      }

      void switch_memory_context(viennacl::context new_ctx)
      {
        viennacl::backend::switch_memory_context<SCALARTYPE>(elements_, new_ctx);
      }

      //TODO: Think about implementing the following public member functions
      //void insert_element(unsigned int i, SCALARTYPE val){}
      //void erase_element(unsigned int i){}

      //enlarge or reduce allocated memory and set unused memory to zero
      /** @brief Resizes the allocated memory for the vector. Pads the memory to be a multiple of 'ALIGNMENT'
      *
      *  @param new_size  The new size of the vector
      *  @param preserve  If true, old entries of the vector are preserved, otherwise eventually discarded.
      */
      void resize(size_type new_size, bool preserve = true)
      {
        resize_impl(new_size, viennacl::traits::context(*this), preserve);
      }

      /** @brief Resizes the allocated memory for the vector. Convenience function for setting an OpenCL context in case reallocation is needed
      *
      *  @param new_size  The new size of the vector
      *  @param ctx       The context within which the new memory should be allocated
      *  @param preserve  If true, old entries of the vector are preserved, otherwise eventually discarded.
      */
      void resize(size_type new_size, viennacl::context ctx, bool preserve = true)
      {
        resize_impl(new_size, ctx, preserve);
      }

    private:

      void resize_impl(size_type new_size, viennacl::context ctx, bool preserve = true)
      {
        assert(new_size > 0 && bool("Positive size required when resizing vector!"));

        if (new_size != size_)
        {
          vcl_size_t new_internal_size = viennacl::tools::align_to_multiple<vcl_size_t>(new_size, alignment);

          std::vector<SCALARTYPE> temp(size_);
          if (preserve && size_ > 0)
            fast_copy(*this, temp);
          temp.resize(new_size);  //drop all entries above new_size
          temp.resize(new_internal_size); //enlarge to fit new internal size

          if (new_internal_size != internal_size())
          {
            viennacl::backend::memory_create(elements_, sizeof(SCALARTYPE)*new_internal_size, ctx, NULL);
          }

          fast_copy(temp, *this);
          size_ = new_size;
          internal_size_ = viennacl::tools::align_to_multiple<size_type>(size_, alignment);
          pad();
        }

      }

      size_type       size_;
      size_type       start_;
      difference_type stride_;
      size_type       internal_size_;
      handle_type elements_;
  }; //vector_base



  // forward definition in forwards.h!
  /** @brief A vector class representing a linear memory sequence on the GPU. Inspired by boost::numeric::ublas::vector
  *
  *  This is the basic vector type of ViennaCL. It is similar to std::vector and boost::numeric::ublas::vector and supports various linear algebra operations.
  * By default, the internal length of the vector is padded to a multiple of 'ALIGNMENT' in order to speed up several GPU viennacl::ocl::kernels.
  *
  * @tparam SCALARTYPE  The floating point type, either 'float' or 'double'
  * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
  */
  template<class SCALARTYPE, unsigned int ALIGNMENT>
  class vector : public vector_base<SCALARTYPE>
  {
    typedef vector<SCALARTYPE, ALIGNMENT>         self_type;
    typedef vector_base<SCALARTYPE>               base_type;

  public:
    typedef typename base_type::size_type                  size_type;
    typedef typename base_type::difference_type            difference_type;

    /** @brief Default constructor in order to be compatible with various containers.
    */
    explicit vector() : base_type() { /* Note: One must not call ::init() here because the vector might have been created globally before the backend has become available */ }

    /** @brief An explicit constructor for the vector, allocating the given amount of memory (plus a padding specified by 'ALIGNMENT')
    *
    * @param vec_size   The length (i.e. size) of the vector.
    */
    explicit vector(size_type vec_size) : base_type(vec_size) {}

    explicit vector(size_type vec_size, viennacl::context ctx) : base_type(vec_size, ctx) {}

    explicit vector(SCALARTYPE * ptr_to_mem, viennacl::memory_types mem_type, size_type vec_size, size_type start = 0, difference_type stride = 1)
        : base_type(ptr_to_mem, mem_type, vec_size, start, stride) {}

#ifdef VIENNACL_WITH_OPENCL
    /** @brief Create a vector from existing OpenCL memory
    *
    * Note: The provided memory must take an eventual ALIGNMENT into account, i.e. existing_mem must be at least of size internal_size()!
    * This is trivially the case with the default alignment, but should be considered when using vector<> with an alignment parameter not equal to 1.
    *
    * @param existing_mem   An OpenCL handle representing the memory
    * @param vec_size       The size of the vector.
    */
    explicit vector(cl_mem existing_mem, size_type vec_size, size_type start = 0, difference_type stride = 1) : base_type(existing_mem, vec_size, start, stride) {}

    /** @brief An explicit constructor for the vector, allocating the given amount of memory (plus a padding specified by 'ALIGNMENT') and the OpenCL context provided
    *
    * @param vec_size   The length (i.e. size) of the vector.
    * @param ctx        The context
    */
    explicit vector(size_type vec_size, viennacl::ocl::context const & ctx) : base_type(vec_size, ctx) {}
#endif

    template <typename LHS, typename RHS, typename OP>
    vector(vector_expression<const LHS, const RHS, OP> const & proxy) : base_type(proxy) {}

    vector(const base_type & v) : base_type(v.size(), viennacl::traits::context(v))
    {
      if (v.size() > 0)
        base_type::operator=(v);
    }

    vector(const self_type & v) : base_type(v.size(), viennacl::traits::context(v))
    {
      if (v.size() > 0)
        base_type::operator=(v);
    }

    /** @brief Creates the vector from the supplied unit vector. */
    vector(unit_vector<SCALARTYPE> const & v) : base_type(v.size())
    {
      if (v.size() > 0)
        this->operator()(v.index()) = SCALARTYPE(1);;
    }

    /** @brief Creates the vector from the supplied zero vector. */
    vector(zero_vector<SCALARTYPE> const & v) : base_type(v.size(), v.context())
    {
      if (v.size() > 0)
        viennacl::linalg::vector_assign(*this, SCALARTYPE(0.0));
    }

    /** @brief Creates the vector from the supplied scalar vector. */
    vector(scalar_vector<SCALARTYPE> const & v) : base_type(v.size(), v.context())
    {
      if (v.size() > 0)
        viennacl::linalg::vector_assign(*this, v[0]);
    }

    // the following is used to circumvent an issue with Clang 3.0 when 'using base_type::operator=;' directly
    template <typename T>
    self_type & operator=(T const & other)
    {
      base_type::operator=(other);
      return *this;
    }

    using base_type::operator+=;
    using base_type::operator-=;

    //enlarge or reduce allocated memory and set unused memory to zero
    /** @brief Resizes the allocated memory for the vector. Pads the memory to be a multiple of 'ALIGNMENT'
    *
    *  @param new_size  The new size of the vector
    *  @param preserve  If true, old entries of the vector are preserved, otherwise eventually discarded.
    */
    void resize(size_type new_size, bool preserve = true)
    {
      base_type::resize(new_size, preserve);
    }

    void resize(size_type new_size, viennacl::context ctx, bool preserve = true)
    {
      base_type::resize(new_size, ctx, preserve);
    }

    /** @brief Swaps the handles of two vectors by swapping the OpenCL handles only, no data copy
    */
    self_type & fast_swap(self_type & other)
    {
      base_type::fast_swap(other);
      return *this;
    }

    void switch_memory_context(viennacl::context new_ctx)
    {
      base_type::switch_memory_context(new_ctx);
    }

  }; //vector

  /** @brief Tuple class holding pointers to multiple vectors. Mainly used as a temporary object returned from viennacl::tie(). */
  template <typename ScalarT>
  class vector_tuple
  {
    typedef vector_base<ScalarT>   VectorType;

  public:
      // 2 vectors

      vector_tuple(VectorType const & v0, VectorType const & v1) : const_vectors_(2), non_const_vectors_()
      {
        const_vectors_[0] = &v0;
        const_vectors_[1] = &v1;
      }
      vector_tuple(VectorType       & v0, VectorType       & v1) : const_vectors_(2), non_const_vectors_(2)
      {
        const_vectors_[0] = &v0; non_const_vectors_[0] = &v0;
        const_vectors_[1] = &v1; non_const_vectors_[1] = &v1;
      }

      // 3 vectors

      vector_tuple(VectorType const & v0, VectorType const & v1, VectorType const & v2) : const_vectors_(3), non_const_vectors_()
      {
        const_vectors_[0] = &v0;
        const_vectors_[1] = &v1;
        const_vectors_[2] = &v2;
      }
      vector_tuple(VectorType       & v0, VectorType       & v1, VectorType       & v2) : const_vectors_(3), non_const_vectors_(3)
      {
        const_vectors_[0] = &v0; non_const_vectors_[0] = &v0;
        const_vectors_[1] = &v1; non_const_vectors_[1] = &v1;
        const_vectors_[2] = &v2; non_const_vectors_[2] = &v2;
      }

      // 4 vectors

      vector_tuple(VectorType const & v0, VectorType const & v1, VectorType const & v2, VectorType const & v3) : const_vectors_(4), non_const_vectors_()
      {
        const_vectors_[0] = &v0;
        const_vectors_[1] = &v1;
        const_vectors_[2] = &v2;
        const_vectors_[3] = &v3;
      }
      vector_tuple(VectorType       & v0, VectorType       & v1, VectorType       & v2, VectorType       & v3) : const_vectors_(4), non_const_vectors_(4)
      {
        const_vectors_[0] = &v0; non_const_vectors_[0] = &v0;
        const_vectors_[1] = &v1; non_const_vectors_[1] = &v1;
        const_vectors_[2] = &v2; non_const_vectors_[2] = &v2;
        const_vectors_[3] = &v3; non_const_vectors_[3] = &v3;
      }

      // add more overloads here

      // generic interface:

      vector_tuple(std::vector<VectorType const *> const & vecs) : const_vectors_(vecs.size()), non_const_vectors_()
      {
        for (vcl_size_t i=0; i<vecs.size(); ++i)
          const_vectors_[i] = vecs[i];
      }

      vector_tuple(std::vector<VectorType *> const & vecs) : const_vectors_(vecs.size()), non_const_vectors_(vecs.size())
      {
        for (vcl_size_t i=0; i<vecs.size(); ++i)
        {
              const_vectors_[i] = vecs[i];
          non_const_vectors_[i] = vecs[i];
        }
      }

      vcl_size_t size()       const { return non_const_vectors_.size(); }
      vcl_size_t const_size() const { return     const_vectors_.size(); }

      VectorType       &       at(vcl_size_t i) const { return *(non_const_vectors_.at(i)); }
      VectorType const & const_at(vcl_size_t i) const { return     *(const_vectors_.at(i)); }

  private:
    std::vector<VectorType const *>   const_vectors_;
    std::vector<VectorType *>         non_const_vectors_;
  };

  // 2 args
  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT> const & v0, vector_base<ScalarT> const & v1) { return vector_tuple<ScalarT>(v0, v1); }

  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT>       & v0, vector_base<ScalarT>       & v1) { return vector_tuple<ScalarT>(v0, v1); }

  // 3 args
  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT> const & v0, vector_base<ScalarT> const & v1, vector_base<ScalarT> const & v2) { return vector_tuple<ScalarT>(v0, v1, v2); }

  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT>       & v0, vector_base<ScalarT>       & v1, vector_base<ScalarT>       & v2) { return vector_tuple<ScalarT>(v0, v1, v2); }

  // 4 args
  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT> const & v0, vector_base<ScalarT> const & v1, vector_base<ScalarT> const & v2, vector_base<ScalarT> const & v3)
  {
    return vector_tuple<ScalarT>(v0, v1, v2, v3);
  }

  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT>       & v0, vector_base<ScalarT>       & v1, vector_base<ScalarT>       & v2, vector_base<ScalarT>       & v3)
  {
    return vector_tuple<ScalarT>(v0, v1, v2, v3);
  }

  // 5 args
  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT> const & v0,
                            vector_base<ScalarT> const & v1,
                            vector_base<ScalarT> const & v2,
                            vector_base<ScalarT> const & v3,
                            vector_base<ScalarT> const & v4)
  {
    typedef vector_base<ScalarT> const *       VectorPointerType;
    std::vector<VectorPointerType> vec(5);
    vec[0] = &v0;
    vec[1] = &v1;
    vec[2] = &v2;
    vec[3] = &v3;
    vec[4] = &v4;
    return vector_tuple<ScalarT>(vec);
  }

  template <typename ScalarT>
  vector_tuple<ScalarT> tie(vector_base<ScalarT> & v0,
                            vector_base<ScalarT> & v1,
                            vector_base<ScalarT> & v2,
                            vector_base<ScalarT> & v3,
                            vector_base<ScalarT> & v4)
  {
    typedef vector_base<ScalarT> *       VectorPointerType;
    std::vector<VectorPointerType> vec(5);
    vec[0] = &v0;
    vec[1] = &v1;
    vec[2] = &v2;
    vec[3] = &v3;
    vec[4] = &v4;
    return vector_tuple<ScalarT>(vec);
  }

  // TODO: Add more arguments to tie() here. Maybe use some preprocessor magic to accomplish this.

  //
  //////////////////// Copy from GPU to CPU //////////////////////////////////
  //


  /** @brief STL-like transfer of a GPU vector to the CPU. The cpu type is assumed to reside in a linear piece of memory, such as e.g. for std::vector.
  *
  * This method is faster than the plain copy() function, because entries are
  * directly written to the cpu vector, starting with &(*cpu.begin()) However,
  * keep in mind that the cpu type MUST represent a linear piece of
  * memory, otherwise you will run into undefined behavior.
  *
  * @param gpu_begin  GPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_end    GPU iterator pointing to the end of the vector (STL-like)
  * @param cpu_begin  Output iterator for the cpu vector. The cpu vector must be at least as long as the gpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void fast_copy(const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_begin,
                  const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_end,
                  CPU_ITERATOR cpu_begin )
  {
    if (gpu_begin != gpu_end)
    {
      if (gpu_begin.stride() == 1)
      {
        viennacl::backend::memory_read(gpu_begin.handle(),
                                      sizeof(SCALARTYPE)*gpu_begin.offset(),
                                      sizeof(SCALARTYPE)*gpu_begin.stride() * (gpu_end - gpu_begin),
                                      &(*cpu_begin));
      }
      else
      {
        vcl_size_t gpu_size = (gpu_end - gpu_begin);
        std::vector<SCALARTYPE> temp_buffer(gpu_begin.stride() * gpu_size);
        viennacl::backend::memory_read(gpu_begin.handle(), sizeof(SCALARTYPE)*gpu_begin.offset(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));

        for (vcl_size_t i=0; i<gpu_size; ++i)
        {
          (&(*cpu_begin))[i] = temp_buffer[i * gpu_begin.stride()];
        }
      }
    }
  }

  /** @brief Transfer from a gpu vector to a cpu vector. Convenience wrapper for viennacl::linalg::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  *
  * @param gpu_vec    A gpu vector.
  * @param cpu_vec    The cpu vector. Type requirements: Output iterator pointing to entries linear in memory can be obtained via member function .begin()
  */
  template <typename NumericT, typename CPUVECTOR>
  void fast_copy(vector_base<NumericT> const & gpu_vec, CPUVECTOR & cpu_vec )
  {
    viennacl::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  }


  /** @brief Asynchronous version of fast_copy(), copying data from device to host. The host iterator cpu_begin needs to reside in a linear piece of memory, such as e.g. for std::vector.
  *
  * This method allows for overlapping data transfer with host computation and returns immediately if the gpu vector has a unit-stride.
  * In order to wait for the transfer to complete, use viennacl::backend::finish().
  * Note that data pointed to by cpu_begin must not be modified prior to completion of the transfer.
  *
  * @param gpu_begin  GPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_end    GPU iterator pointing to the end of the vector (STL-like)
  * @param cpu_begin  Output iterator for the cpu vector. The cpu vector must be at least as long as the gpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void async_copy(const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_begin,
                  const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_end,
                  CPU_ITERATOR cpu_begin )
  {
    if (gpu_begin != gpu_end)
    {
      if (gpu_begin.stride() == 1)
      {
        viennacl::backend::memory_read(gpu_begin.handle(),
                                       sizeof(SCALARTYPE)*gpu_begin.offset(),
                                       sizeof(SCALARTYPE)*gpu_begin.stride() * (gpu_end - gpu_begin),
                                       &(*cpu_begin),
                                       true);
      }
      else // no async copy possible, so fall-back to fast_copy
        fast_copy(gpu_begin, gpu_end, cpu_begin);
    }
  }

  /** @brief Transfer from a gpu vector to a cpu vector. Convenience wrapper for viennacl::linalg::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  *
  * @param gpu_vec    A gpu vector.
  * @param cpu_vec    The cpu vector. Type requirements: Output iterator pointing to entries linear in memory can be obtained via member function .begin()
  */
  template <typename NumericT, typename CPUVECTOR>
  void async_copy(vector_base<NumericT> const & gpu_vec, CPUVECTOR & cpu_vec )
  {
    viennacl::async_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  }


  /** @brief STL-like transfer for the entries of a GPU vector to the CPU. The cpu type does not need to lie in a linear piece of memory.
  *
  * @param gpu_begin  GPU constant iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_end    GPU constant iterator pointing to the end of the vector (STL-like)
  * @param cpu_begin  Output iterator for the cpu vector. The cpu vector must be at least as long as the gpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void copy(const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_begin,
            const const_vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_end,
            CPU_ITERATOR cpu_begin )
  {
    assert(gpu_end - gpu_begin >= 0 && bool("Iterators incompatible"));
    if (gpu_end - gpu_begin != 0)
    {
      std::vector<SCALARTYPE> temp_buffer(gpu_end - gpu_begin);
      fast_copy(gpu_begin, gpu_end, temp_buffer.begin());

      //now copy entries to cpu_vec:
      std::copy(temp_buffer.begin(), temp_buffer.end(), cpu_begin);
    }
  }

  /** @brief STL-like transfer for the entries of a GPU vector to the CPU. The cpu type does not need to lie in a linear piece of memory.
  *
  * @param gpu_begin  GPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_end    GPU iterator pointing to the end of the vector (STL-like)
  * @param cpu_begin  Output iterator for the cpu vector. The cpu vector must be at least as long as the gpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void copy(const vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_begin,
            const vector_iterator<SCALARTYPE, ALIGNMENT> & gpu_end,
            CPU_ITERATOR cpu_begin )

  {
    viennacl::copy(const_vector_iterator<SCALARTYPE, ALIGNMENT>(gpu_begin),
                    const_vector_iterator<SCALARTYPE, ALIGNMENT>(gpu_end),
                    cpu_begin);
  }

  /** @brief Transfer from a gpu vector to a cpu vector. Convenience wrapper for viennacl::linalg::copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  *
  * @param gpu_vec    A gpu vector
  * @param cpu_vec    The cpu vector. Type requirements: Output iterator can be obtained via member function .begin()
  */
  template <typename NumericT, typename CPUVECTOR>
  void copy(vector_base<NumericT> const & gpu_vec, CPUVECTOR & cpu_vec )
  {
    viennacl::copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  }



  #ifdef VIENNACL_WITH_EIGEN
  template <unsigned int ALIGNMENT>
  void copy(vector<float, ALIGNMENT> const & gpu_vec,
            Eigen::VectorXf & eigen_vec)
  {
    viennacl::fast_copy(gpu_vec.begin(), gpu_vec.end(), &(eigen_vec[0]));
  }

  template <unsigned int ALIGNMENT>
  void copy(vector<double, ALIGNMENT> & gpu_vec,
            Eigen::VectorXd & eigen_vec)
  {
    viennacl::fast_copy(gpu_vec.begin(), gpu_vec.end(), &(eigen_vec[0]));
  }
  #endif


  //
  //////////////////// Copy from CPU to GPU //////////////////////////////////
  //

  /** @brief STL-like transfer of a CPU vector to the GPU. The cpu type is assumed to reside in a linear piece of memory, such as e.g. for std::vector.
  *
  * This method is faster than the plain copy() function, because entries are
  * directly read from the cpu vector, starting with &(*cpu.begin()). However,
  * keep in mind that the cpu type MUST represent a linear piece of
  * memory, otherwise you will run into undefined behavior.
  *
  * @param cpu_begin  CPU iterator pointing to the beginning of the cpu vector (STL-like)
  * @param cpu_end    CPU iterator pointing to the end of the vector (STL-like)
  * @param gpu_begin  Output iterator for the gpu vector. The gpu iterator must be incrementable (cpu_end - cpu_begin) times, otherwise the result is undefined.
  */
  template <typename CPU_ITERATOR, typename SCALARTYPE, unsigned int ALIGNMENT>
  void fast_copy(CPU_ITERATOR const & cpu_begin,
                  CPU_ITERATOR const & cpu_end,
                  vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin)
  {
    if (cpu_end - cpu_begin > 0)
    {
      if (gpu_begin.stride() == 1)
      {
        viennacl::backend::memory_write(gpu_begin.handle(),
                                        sizeof(SCALARTYPE)*gpu_begin.offset(),
                                        sizeof(SCALARTYPE)*gpu_begin.stride() * (cpu_end - cpu_begin), &(*cpu_begin));
      }
      else //writing to slice:
      {
        vcl_size_t cpu_size = (cpu_end - cpu_begin);
        std::vector<SCALARTYPE> temp_buffer(gpu_begin.stride() * cpu_size);

        viennacl::backend::memory_read(gpu_begin.handle(), sizeof(SCALARTYPE)*gpu_begin.offset(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));

        for (vcl_size_t i=0; i<cpu_size; ++i)
          temp_buffer[i * gpu_begin.stride()] = (&(*cpu_begin))[i];

        viennacl::backend::memory_write(gpu_begin.handle(), sizeof(SCALARTYPE)*gpu_begin.offset(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));
      }
    }
  }


  /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  *
  * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
  * @param gpu_vec    The gpu vector.
  */
  template <typename CPUVECTOR, typename NumericT>
  void fast_copy(const CPUVECTOR & cpu_vec, vector_base<NumericT> & gpu_vec)
  {
    viennacl::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  }

  /** @brief Asynchronous version of fast_copy(), copying data from host to device. The host iterator cpu_begin needs to reside in a linear piece of memory, such as e.g. for std::vector.
  *
  * This method allows for overlapping data transfer with host computation and returns immediately if the gpu vector has a unit-stride.
  * In order to wait for the transfer to complete, use viennacl::backend::finish().
  * Note that data pointed to by cpu_begin must not be modified prior to completion of the transfer.
  *
  * @param cpu_begin  CPU iterator pointing to the beginning of the cpu vector (STL-like)
  * @param cpu_end    CPU iterator pointing to the end of the vector (STL-like)
  * @param gpu_begin  Output iterator for the gpu vector. The gpu iterator must be incrementable (cpu_end - cpu_begin) times, otherwise the result is undefined.
  */
  template <typename CPU_ITERATOR, typename SCALARTYPE, unsigned int ALIGNMENT>
  void async_copy(CPU_ITERATOR const & cpu_begin,
                  CPU_ITERATOR const & cpu_end,
                  vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin)
  {
    if (cpu_end - cpu_begin > 0)
    {
      if (gpu_begin.stride() == 1)
      {
        viennacl::backend::memory_write(gpu_begin.handle(),
                                        sizeof(SCALARTYPE)*gpu_begin.offset(),
                                        sizeof(SCALARTYPE)*gpu_begin.stride() * (cpu_end - cpu_begin), &(*cpu_begin),
                                        true);
      }
      else // fallback to blocking copy. There's nothing we can do to prevent this
        fast_copy(cpu_begin, cpu_end, gpu_begin);
    }
  }


  /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  *
  * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
  * @param gpu_vec    The gpu vector.
  */
  template <typename CPUVECTOR, typename NumericT>
  void async_copy(const CPUVECTOR & cpu_vec, vector_base<NumericT> & gpu_vec)
  {
    viennacl::async_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  }

  //from cpu to gpu. Safe assumption: cpu_vector does not necessarily occupy a linear memory segment, but is not larger than the allocated memory on the GPU
  /** @brief STL-like transfer for the entries of a GPU vector to the CPU. The cpu type does not need to lie in a linear piece of memory.
  *
  * @param cpu_begin  CPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param cpu_end    CPU iterator pointing to the end of the vector (STL-like)
  * @param gpu_begin  Output iterator for the gpu vector. The gpu vector must be at least as long as the cpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
  void copy(CPU_ITERATOR const & cpu_begin,
            CPU_ITERATOR const & cpu_end,
            vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin)
  {
    assert(cpu_end - cpu_begin > 0 && bool("Iterators incompatible"));
    if (cpu_begin != cpu_end)
    {
      //we require that the size of the gpu_vector is larger or equal to the cpu-size
      std::vector<SCALARTYPE> temp_buffer(cpu_end - cpu_begin);
      std::copy(cpu_begin, cpu_end, temp_buffer.begin());
      viennacl::fast_copy(temp_buffer.begin(), temp_buffer.end(), gpu_begin);
    }
  }

  // for things like copy(std_vec.begin(), std_vec.end(), vcl_vec.begin() + 1);

  /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  *
  * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
  * @param gpu_vec    The gpu vector.
  */
  template <typename CPUVECTOR, typename T>
  void copy(const CPUVECTOR & cpu_vec, vector_base<T> & gpu_vec)
  {
    viennacl::copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  }


  #ifdef VIENNACL_WITH_EIGEN
  template <unsigned int ALIGNMENT>
  void copy(Eigen::VectorXf const & eigen_vec,
            vector<float, ALIGNMENT> & gpu_vec)
  {
    std::vector<float> entries(eigen_vec.size());
    for (vcl_size_t i = 0; i<entries.size(); ++i)
      entries[i] = eigen_vec(i);
    viennacl::fast_copy(entries.begin(), entries.end(), gpu_vec.begin());
  }

  template <unsigned int ALIGNMENT>
  void copy(Eigen::VectorXd const & eigen_vec,
            vector<double, ALIGNMENT> & gpu_vec)
  {
    std::vector<double> entries(eigen_vec.size());
    for (vcl_size_t i = 0; i<entries.size(); ++i)
      entries[i] = eigen_vec(i);
    viennacl::fast_copy(entries.begin(), entries.end(), gpu_vec.begin());
  }
  #endif



  //
  //////////////////// Copy from GPU to GPU //////////////////////////////////
  //
  /** @brief Copy (parts of a) GPU vector to another GPU vector
  *
  * @param gpu_src_begin    GPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_src_end      GPU iterator pointing to the end of the vector (STL-like)
  * @param gpu_dest_begin   Output iterator for the gpu vector. The gpu_dest vector must be at least as long as the gpu_src vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
  void copy(const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_begin,
            const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_end,
            vector_iterator<SCALARTYPE, ALIGNMENT_DEST> gpu_dest_begin)
  {
    assert(gpu_src_end - gpu_src_begin >= 0);
    assert(gpu_src_begin.stride() == 1 && bool("ViennaCL ERROR: copy() for GPU->GPU not implemented for slices! Use operator= instead for the moment."));

    if (gpu_src_begin.stride() == 1 && gpu_dest_begin.stride() == 1)
    {
      if (gpu_src_begin != gpu_src_end)
        viennacl::backend::memory_copy(gpu_src_begin.handle(), gpu_dest_begin.handle(),
                                        sizeof(SCALARTYPE) * gpu_src_begin.offset(),
                                        sizeof(SCALARTYPE) * gpu_dest_begin.offset(),
                                        sizeof(SCALARTYPE) * (gpu_src_end.offset() - gpu_src_begin.offset()));
    }
    else
    {
      assert( false && bool("not implemented yet"));
    }
  }

  /** @brief Copy (parts of a) GPU vector to another GPU vector
  *
  * @param gpu_src_begin   GPU iterator pointing to the beginning of the gpu vector (STL-like)
  * @param gpu_src_end     GPU iterator pointing to the end of the vector (STL-like)
  * @param gpu_dest_begin  Output iterator for the gpu vector. The gpu vector must be at least as long as the cpu vector!
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
  void copy(vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_begin,
            vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_end,
            vector_iterator<SCALARTYPE, ALIGNMENT_DEST> gpu_dest_begin)
  {
    viennacl::copy(static_cast<const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> >(gpu_src_begin),
                    static_cast<const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> >(gpu_src_end),
                    gpu_dest_begin);
  }

  /** @brief Transfer from a ViennaCL vector to another ViennaCL vector. Convenience wrapper for viennacl::linalg::copy(gpu_src_vec.begin(), gpu_src_vec.end(), gpu_dest_vec.begin());
  *
  * @param gpu_src_vec    A gpu vector
  * @param gpu_dest_vec    The cpu vector. Type requirements: Output iterator can be obtained via member function .begin()
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
  void copy(vector<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_vec,
            vector<SCALARTYPE, ALIGNMENT_DEST> & gpu_dest_vec )
  {
    viennacl::copy(gpu_src_vec.begin(), gpu_src_vec.end(), gpu_dest_vec.begin());
  }






  //global functions for handling vectors:
  /** @brief Output stream. Output format is ublas compatible.
  * @param os   STL output stream
  * @param val  The vector that should be printed
  */
  template <typename T>
  std::ostream & operator<<(std::ostream & os, vector_base<T> const & val)
  {
    std::vector<T> tmp(val.size());
    viennacl::copy(val.begin(), val.end(), tmp.begin());
    os << "[" << val.size() << "](";
    for (typename std::vector<T>::size_type i=0; i<val.size(); ++i)
    {
      if (i > 0)
        os << ",";
      os << tmp[i];
    }
    os << ")";
    return os;
  }

  template <typename LHS, typename RHS, typename OP>
  std::ostream & operator<<(std::ostream & os, vector_expression<LHS, RHS, OP> const & proxy)

  {
    typedef typename viennacl::result_of::cpu_value_type<typename LHS::value_type>::type ScalarType;
    viennacl::vector<ScalarType> result = proxy;
    os << result;
    return os;
  }

  /** @brief Swaps the contents of two vectors, data is copied
  *
  * @param vec1   The first vector
  * @param vec2   The second vector
  */
  template <typename T>
  void swap(vector_base<T> & vec1, vector_base<T> & vec2)
  {
    viennacl::linalg::vector_swap(vec1, vec2);
  }

  /** @brief Swaps the content of two vectors by swapping OpenCL handles only, NO data is copied
  *
  * @param v1   The first vector
  * @param v2   The second vector
  */
  template <typename SCALARTYPE, unsigned int ALIGNMENT>
  vector<SCALARTYPE, ALIGNMENT> & fast_swap(vector<SCALARTYPE, ALIGNMENT> & v1,
                                            vector<SCALARTYPE, ALIGNMENT> & v2)
  {
    return v1.fast_swap(v2);
  }





  //
  //
  ////////// operations /////////////////////////////////////////////////////////////////////////////////
  //
  //


  //
  // operator *=
  //

  /** @brief Scales this vector by a GPU scalar value
  */
  template <typename T, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                vector_base<T> &
                              >::type
  operator *= (vector_base<T> & v1, S1 const & gpu_val)
  {
    if (v1.size() > 0)
      viennacl::linalg::av(v1,
                           v1, gpu_val, 1, false, (viennacl::is_flip_sign_scalar<S1>::value ? true : false));
    return v1;
  }


  //
  // operator /=
  //


  /** @brief Scales this vector by a GPU scalar value
  */
  template <typename T, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                vector_base<T> &
                              >::type
  operator /= (vector_base<T> & v1, S1 const & gpu_val)
  {
    if (v1.size() > 0)
      viennacl::linalg::av(v1,
                           v1, gpu_val, 1, true, (viennacl::is_flip_sign_scalar<S1>::value ? true : false));
    return v1;
  }


  //
  // operator +
  //


  /** @brief Operator overload for the addition of two vector expressions.
  *
  * @param proxy1  Left hand side vector expression
  * @param proxy2  Right hand side vector expression
  */
  template <typename LHS1, typename RHS1, typename OP1,
            typename LHS2, typename RHS2, typename OP2>
  vector_expression< const vector_expression< LHS1, RHS1, OP1>,
                     const vector_expression< LHS2, RHS2, OP2>,
                     viennacl::op_add>
  operator + (vector_expression<LHS1, RHS1, OP1> const & proxy1,
              vector_expression<LHS2, RHS2, OP2> const & proxy2)
  {
    assert(proxy1.size() == proxy2.size() && bool("Incompatible vector sizes!"));
    return   vector_expression< const vector_expression<LHS1, RHS1, OP1>,
                                const vector_expression<LHS2, RHS2, OP2>,
                                viennacl::op_add>(proxy1, proxy2);
  }

  /** @brief Operator overload for the addition of a vector expression with a vector or another vector expression. This is the default implementation for all cases that are too complex in order to be covered within a single kernel, hence a temporary vector is created.
  *
  * @param proxy   Left hand side vector expression
  * @param vec     Right hand side vector (also -range and -slice is allowed)
  */
  template <typename LHS, typename RHS, typename OP, typename T>
  vector_expression< const vector_expression<LHS, RHS, OP>,
                     const vector_base<T>,
                     viennacl::op_add>
  operator + (vector_expression<LHS, RHS, OP> const & proxy,
              vector_base<T> const & vec)
  {
    assert(proxy.size() == vec.size() && bool("Incompatible vector sizes!"));
    return vector_expression< const vector_expression<LHS, RHS, OP>,
                              const vector_base<T>,
                              viennacl::op_add>(proxy, vec);
  }

  /** @brief Operator overload for the addition of a vector with a vector expression. This is the default implementation for all cases that are too complex in order to be covered within a single kernel, hence a temporary vector is created.
  *
  * @param proxy   Left hand side vector expression
  * @param vec     Right hand side vector (also -range and -slice is allowed)
  */
  template <typename T, typename LHS, typename RHS, typename OP>
  vector_expression< const vector_base<T>,
                     const vector_expression<LHS, RHS, OP>,
                     viennacl::op_add>
  operator + (vector_base<T> const & vec,
              vector_expression<LHS, RHS, OP> const & proxy)
  {
    assert(proxy.size() == vec.size() && bool("Incompatible vector sizes!"));
    return vector_expression< const vector_base<T>,
                              const vector_expression<LHS, RHS, OP>,
                              viennacl::op_add>(vec, proxy);
  }

  /** @brief Returns an expression template object for adding up two vectors, i.e. v1 + v2
  */
  template <typename T>
  vector_expression< const vector_base<T>, const vector_base<T>, op_add>
  operator + (const vector_base<T> & v1, const vector_base<T> & v2)
  {
    return vector_expression< const vector_base<T>, const vector_base<T>, op_add>(v1, v2);
  }



  //
  // operator -
  //

  /** @brief Operator overload for the subtraction of two vector expressions.
  *
  * @param proxy1  Left hand side vector expression
  * @param proxy2  Right hand side vector expression
  */
  template <typename LHS1, typename RHS1, typename OP1,
            typename LHS2, typename RHS2, typename OP2>
  vector_expression< const vector_expression< LHS1, RHS1, OP1>,
                     const vector_expression< LHS2, RHS2, OP2>,
                     viennacl::op_sub>
  operator - (vector_expression<LHS1, RHS1, OP1> const & proxy1,
              vector_expression<LHS2, RHS2, OP2> const & proxy2)
  {
    assert(proxy1.size() == proxy2.size() && bool("Incompatible vector sizes!"));
    return   vector_expression< const vector_expression<LHS1, RHS1, OP1>,
                                const vector_expression<LHS2, RHS2, OP2>,
                                viennacl::op_sub>(proxy1, proxy2);
  }


  /** @brief Operator overload for the subtraction of a vector expression with a vector or another vector expression. This is the default implementation for all cases that are too complex in order to be covered within a single kernel, hence a temporary vector is created.
  *
  * @param proxy   Left hand side vector expression
  * @param vec     Right hand side vector (also -range and -slice is allowed)
  */
  template <typename LHS, typename RHS, typename OP, typename T>
  vector_expression< const vector_expression<LHS, RHS, OP>,
                     const vector_base<T>,
                     viennacl::op_sub>
  operator - (vector_expression<LHS, RHS, OP> const & proxy,
              vector_base<T> const & vec)
  {
    assert(proxy.size() == vec.size() && bool("Incompatible vector sizes!"));
    return vector_expression< const vector_expression<LHS, RHS, OP>,
                              const vector_base<T>,
                              viennacl::op_sub>(proxy, vec);
  }

  /** @brief Operator overload for the subtraction of a vector expression with a vector or another vector expression. This is the default implementation for all cases that are too complex in order to be covered within a single kernel, hence a temporary vector is created.
  *
  * @param proxy   Left hand side vector expression
  * @param vec     Right hand side vector (also -range and -slice is allowed)
  */
  template <typename T, typename LHS, typename RHS, typename OP>
  vector_expression< const vector_base<T>,
                     const vector_expression<LHS, RHS, OP>,
                     viennacl::op_sub>
  operator - (vector_base<T> const & vec,
              vector_expression<LHS, RHS, OP> const & proxy)
  {
    assert(proxy.size() == vec.size() && bool("Incompatible vector sizes!"));
    return vector_expression< const vector_base<T>,
                              const vector_expression<LHS, RHS, OP>,
                              viennacl::op_sub>(vec, proxy);
  }

  /** @brief Returns an expression template object for subtracting two vectors, i.e. v1 - v2
  */
  template <typename T>
  vector_expression< const vector_base<T>, const vector_base<T>, op_sub>
  operator - (const vector_base<T> & v1, const vector_base<T> & v2)
  {
    return vector_expression< const vector_base<T>, const vector_base<T>, op_sub>(v1, v2);
  }


  //
  // operator *
  //


  /** @brief Operator overload for the expression alpha * v1, where alpha is a host scalar (float or double) and v1 is a ViennaCL vector.
  *
  * @param value   The host scalar (float or double)
  * @param vec     A ViennaCL vector
  */
  template <typename S1, typename T>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                vector_expression< const vector_base<T>, const S1, op_mult> >::type
  operator * (S1 const & value, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const S1, op_mult>(vec, value);
  }

  /** @brief Operator overload for the expression alpha * v1, where alpha is a char
  *
  * @param value   The host scalar (float or double)
  * @param vec     A ViennaCL vector
  */
  template <typename T>
  vector_expression< const vector_base<T>, const T, op_mult>
  operator * (char value, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const T, op_mult>(vec, value);
  }

  /** @brief Operator overload for the expression alpha * v1, where alpha is a short
  *
  * @param value   The host scalar (float or double)
  * @param vec     A ViennaCL vector
  */
  template <typename T>
  vector_expression< const vector_base<T>, const T, op_mult>
  operator * (short value, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const T, op_mult>(vec, value);
  }

  /** @brief Operator overload for the expression alpha * v1, where alpha is a int
  *
  * @param value   The host scalar (float or double)
  * @param vec     A ViennaCL vector
  */
  template <typename T>
  vector_expression< const vector_base<T>, const T, op_mult>
  operator * (int value, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const T, op_mult>(vec, value);
  }

  /** @brief Operator overload for the expression alpha * v1, where alpha is a long
  *
  * @param value   The host scalar (float or double)
  * @param vec     A ViennaCL vector
  */
  template <typename T>
  vector_expression< const vector_base<T>, const T, op_mult>
  operator * (long value, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const T, op_mult>(vec, value);
  }




  /** @brief Operator overload for the expression alpha * v1, where alpha is a scalar expression and v1 is a ViennaCL vector.
  *
  * @param expr    The scalar expression
  * @param vec     A ViennaCL vector
  */
  template <typename LHS, typename RHS, typename OP, typename T>
  vector_expression< const vector_base<T>, const scalar_expression<LHS, RHS, OP>, op_mult>
  operator * (scalar_expression<LHS, RHS, OP> const & expr, vector_base<T> const & vec)
  {
    return vector_expression< const vector_base<T>, const scalar_expression<LHS, RHS, OP>, op_mult>(vec, expr);
  }

  /** @brief Scales the vector by a scalar 'alpha' and returns an expression template
  */
  template <typename T, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                vector_expression< const vector_base<T>, const S1, op_mult> >::type
  operator * (vector_base<T> const & vec, S1 const & value)
  {
    return vector_expression< const vector_base<T>, const S1, op_mult>(vec, value);
  }

  template <typename T>
  vector_expression< const vector_base<T>, const T, op_mult>
  operator * (vector_base<T> const & vec, T const & value)
  {
    return vector_expression< const vector_base<T>, const T, op_mult>(vec, value);
  }

  /** @brief Operator overload for the multiplication of a vector expression with a scalar from the right, e.g. (beta * vec1) * alpha. Here, beta * vec1 is wrapped into a vector_expression and then multiplied with alpha from the right.
  *
  * @param proxy   Left hand side vector expression
  * @param val     Right hand side scalar
  */
  template <typename LHS, typename RHS, typename OP, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_mult>  >::type
  operator * (vector_expression< LHS, RHS, OP> const & proxy,
              S1 const & val)
  {
    return viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_mult>(proxy, val);
  }

  /** @brief Operator overload for the multiplication of a vector expression with a ViennaCL scalar from the left, e.g. alpha * (beta * vec1). Here, beta * vec1 is wrapped into a vector_expression and then multiplied with alpha from the left.
  *
  * @param val     Right hand side scalar
  * @param proxy   Left hand side vector expression
  */
  template <typename S1, typename LHS, typename RHS, typename OP>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_mult>  >::type
  operator * (S1 const & val,
              vector_expression<LHS, RHS, OP> const & proxy)
  {
    return viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_mult>(proxy, val);
  }

  //
  // operator /
  //

  /** @brief Operator overload for the division of a vector expression by a scalar from the right, e.g. (beta * vec1) / alpha. Here, beta * vec1 is wrapped into a vector_expression and then divided by alpha.
  *
  * @param proxy   Left hand side vector expression
  * @param val     Right hand side scalar
  */
  template <typename S1, typename LHS, typename RHS, typename OP>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_div>  >::type
  operator / (vector_expression< LHS, RHS, OP> const & proxy,
              S1 const & val)
  {
    return viennacl::vector_expression<const vector_expression<LHS, RHS, OP>, const S1, op_div>(proxy, val);
  }


  /** @brief Returns an expression template for scaling the vector by a GPU scalar 'alpha'
  */
  template <typename T, typename S1>
  typename viennacl::enable_if< viennacl::is_any_scalar<S1>::value,
                                vector_expression< const vector_base<T>, const S1, op_div> >::type
  operator / (vector_base<T> const & v1, S1 const & s1)
  {
    return vector_expression<const vector_base<T>, const S1, op_div>(v1, s1);
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
      template <typename T>
      struct op_executor<vector_base<T>, op_assign, vector_base<T> >
      {
        static void apply(vector_base<T> & lhs, vector_base<T> const & rhs)
        {
          viennacl::linalg::av(lhs, rhs, T(1), 1, false, false);
        }
      };

      // x = inner_prod(z, {y0, y1, ...})
      template <typename T>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const vector_base<T>, const vector_tuple<T>, op_inner_prod> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_tuple<T>, op_inner_prod> const & rhs)
        {
          viennacl::linalg::inner_prod_impl(rhs.lhs(), rhs.rhs(), lhs);
        }
      };

      // x += y
      template <typename T>
      struct op_executor<vector_base<T>, op_inplace_add, vector_base<T> >
      {
        static void apply(vector_base<T> & lhs, vector_base<T> const & rhs)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, rhs, T(1), 1, false, false);
        }
      };

      // x -= y
      template <typename T>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_base<T> >
      {
        static void apply(vector_base<T> & lhs, vector_base<T> const & rhs)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, rhs, T(1), 1, false, true);
        }
      };

      ///////////// x  OP  y * alpha ////////////////////////


      // x = alpha * y
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const vector_base<T>, const ScalarType, op_mult> >
      {
        // generic case: ScalarType is a scalar expression
        template <typename LHS, typename RHS, typename OP>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar_expression<LHS, RHS, OP>, op_mult> const & proxy)
        {
          T alpha = proxy.rhs();
          viennacl::linalg::av(lhs, proxy.lhs(), alpha, 1, false, false);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar<T>, op_mult> const & proxy)
        {
          viennacl::linalg::av(lhs, proxy.lhs(), proxy.rhs(), 1, false, false);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const T, op_mult> const & proxy)
        {
          viennacl::linalg::av(lhs, proxy.lhs(), proxy.rhs(), 1, false, false);
        }
      };

      // x += alpha * y
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const vector_base<T>, const ScalarType, op_mult> >
      {
        // generic case: ScalarType is a scalar expression
        template <typename LHS, typename RHS, typename OP>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar_expression<LHS, RHS, OP>, op_mult> const & proxy)
        {
          T alpha = proxy.rhs();
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), alpha, 1, false, false);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar<T>, op_mult> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, false);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const T, op_mult> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, false);
        }
      };

      // x -= alpha * y
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const vector_base<T>, const ScalarType, op_mult> >
      {
        // generic case: ScalarType is a scalar expression
        template <typename LHS, typename RHS, typename OP>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar_expression<LHS, RHS, OP>, op_mult> const & proxy)
        {
          T alpha = proxy.rhs();
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), alpha, 1, false, true);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const scalar<T>, op_mult> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, true);
        }

        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const T, op_mult> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, false, true);
        }
      };


      ///////////// x  OP  vec_expr * alpha ////////////////////////

      // x = alpha * vec_expr
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs = temp * proxy.rhs();
          }
      };

      // x += alpha * vec_expr
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs += temp * proxy.rhs();
          }
      };

      // x -= alpha * vec_expr
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_mult> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs -= temp * proxy.rhs();
          }
      };


      ///////////// x  OP  y / alpha ////////////////////////

      // x = y / alpha
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const vector_base<T>, const ScalarType, op_div> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::av(lhs, proxy.lhs(), proxy.rhs(), 1, true, false);
        }
      };

      // x += y / alpha
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const vector_base<T>, const ScalarType, op_div> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, true, false);
        }
      };

      // x -= y / alpha
      template <typename T, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const vector_base<T>, const ScalarType, op_div> >
      {
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const ScalarType, op_div> const & proxy)
        {
          viennacl::linalg::avbv(lhs, lhs, T(1), 1, false, false, proxy.lhs(), proxy.rhs(), 1, true, true);
        }
      };


      ///////////// x  OP  vec_expr / alpha ////////////////////////

      // x = vec_expr / alpha
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs = temp / proxy.rhs();
          }
      };

      // x += vec_expr / alpha
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs += temp / proxy.rhs();
          }
      };

      // x -= vec_expr / alpha
      template <typename T, typename LHS, typename RHS, typename OP, typename ScalarType>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> >
      {
          static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS, const RHS, OP>, const ScalarType, op_div> const & proxy)
          {
            vector<T> temp(proxy.lhs());
            lhs -= temp / proxy.rhs();
          }
      };



      // generic x = vec_expr1 + vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const RHS, op_add> >
      {
        // generic x = vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs = temp;
          }
          else
          {
            op_executor<vector_base<T>, op_assign, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x = y + z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = alpha * y + z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = y / alpha + z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs(), T(1), 1, false, false);
        }

        // x = y + beta * z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = y + z / beta
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x = alpha * y + beta * z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = alpha * y + z / beta
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x = y / alpha + beta * z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  const vector_expression<const vector_base<T>, const T, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x = y / alpha + z / beta
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  const vector_expression<const vector_base<T>, const T, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };


      // generic x += vec_expr1 + vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const LHS, const RHS, op_add> >
      {
        // generic x += vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs += temp;
          }
          else
          {
            op_executor<vector_base<T>, op_inplace_add, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x += y + z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += alpha * y + z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += y / alpha + z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x += y + beta * z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += y + z / beta
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x += alpha * y + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += alpha * y + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x += y / alpha + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x += y / alpha + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };



      // generic x -= vec_expr1 + vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const LHS, const RHS, op_add> >
      {
        // generic x -= vec_expr1 + vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_add> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(temp, proxy.rhs());
            lhs -= temp;
          }
          else
          {
            op_executor<vector_base<T>, op_inplace_sub, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x -= y + z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= alpha * y + z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= y / alpha + z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  const vector_base<T>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x -= y + beta * z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= y + z / beta
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x -= alpha * y + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= alpha * y + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x -= y / alpha + beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x -= y / alpha + z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_add> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };



      ///////////////////////



      // generic x = vec_expr1 - vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const RHS, op_sub> >
      {
        // generic x = vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs = temp;
          }
          else
          {
            op_executor<vector_base<T>, op_assign, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x = y - z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = alpha * y - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = y / alpha - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs(), T(1), 1, false, true);
        }

        // x = y - beta * z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = y - z / beta
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs(), T(1), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x = alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x = y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x = y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv(lhs,
                                 proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                 proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };


      // generic x += vec_expr1 - vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const LHS, const RHS, op_sub> >
      {
        // generic x += vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs += temp;
          }
          else
          {
            op_executor<vector_base<T>, op_inplace_add, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x += y - z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += alpha * y - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += y / alpha - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs(), T(1), 1, false, true);
        }

        // x += y - beta * z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += y - z / beta
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x += alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }

        // x += y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, true);
        }

        // x += y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, false,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, true);
        }
      };



      // generic x -= vec_expr1 - vec_expr2:
      template <typename T, typename LHS, typename RHS>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const LHS, const RHS, op_sub> >
      {
        // generic x -= vec_expr1 - vec_expr2:
        template <typename LHS1, typename RHS1>
        static void apply(vector_base<T> & lhs, vector_expression<const LHS1, const RHS1, op_sub> const & proxy)
        {
          bool op_aliasing_lhs = op_aliasing(lhs, proxy.lhs());
          bool op_aliasing_rhs = op_aliasing(lhs, proxy.rhs());

          if (op_aliasing_lhs || op_aliasing_rhs)
          {
            vector_base<T> temp(proxy.lhs());
            op_executor<vector_base<T>, op_inplace_sub, RHS>::apply(temp, proxy.rhs());
            lhs -= temp;
          }
          else
          {
            op_executor<vector_base<T>, op_inplace_sub, LHS>::apply(lhs, proxy.lhs());
            op_executor<vector_base<T>, op_inplace_add, RHS>::apply(lhs, proxy.rhs());
          }
        }

        // x -= y - z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= alpha * y - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= y / alpha - z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  const vector_base<T>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs(), T(1), 1, false, false);
        }

        // x -= y - beta * z
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= y - z / beta
        template <typename ScalarType>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs(), T(1), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x -= alpha * y - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= alpha * y - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_mult>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, false, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }

        // x -= y / alpha - beta * z
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_mult>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, false, false);
        }

        // x -= y / alpha - z / beta
        template <typename ScalarType1, typename ScalarType2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const vector_base<T>, const ScalarType1, op_div>,
                                                                  const vector_expression<const vector_base<T>, const ScalarType2, op_div>,
                                                                  op_sub> const & proxy)
        {
          viennacl::linalg::avbv_v(lhs,
                                   proxy.lhs().lhs(), proxy.lhs().rhs(), 1, true, true,
                                   proxy.rhs().lhs(), proxy.rhs().rhs(), 1, true, false);
        }
      };


















      //////////////////// Element-wise operations ////////////////////////////////////////

      // generic x = vec_expr1 .* vec_expr2:
      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const RHS, op_element_binary<OP> > >
      {
        // x = y .* z  or  x = y ./ z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          viennacl::linalg::element_op(lhs, proxy);
        }

        // x = y .* vec_expr  or  x = y ./ vec_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_expression<const LHS2, const RHS2, OP2>, op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(proxy.lhs(), temp));
        }

        // x = vec_expr .* z  or  x = vec_expr ./ z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.lhs());
          viennacl::linalg::element_op(lhs, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp, proxy.rhs()));
        }

        // x = vec_expr .* vec_expr  or  z = vec_expr .* vec_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_binary<OP> > const & proxy)
        {
          vector<T> temp1(proxy.lhs());
          vector<T> temp2(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp1, temp2));
        }
      };

      // generic x += vec_expr1 .* vec_expr2:
      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const LHS, const RHS, op_element_binary<OP> > >
      {
        // x += y .* z  or  x += y ./ z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          viennacl::vector<T> temp(proxy);
          lhs += temp;
        }

        // x += y .* vec_expr  or  x += y ./ vec_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_expression<const LHS2, const RHS2, OP2>,  op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          vector<T> temp2(temp.size());
          viennacl::linalg::element_op(temp2, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(proxy.lhs(), temp));
          lhs += temp2;
        }

        // x += vec_expr .* z  or  x += vec_expr ./ z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.lhs());
          vector<T> temp2(temp.size());
          viennacl::linalg::element_op(temp2, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp, proxy.rhs()));
          lhs += temp2;
        }

        // x += vec_expr .* vec_expr  or  x += vec_expr ./ vec_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_binary<OP> > const & proxy)
        {
          vector<T> temp1(proxy.lhs());
          vector<T> temp2(proxy.rhs());
          vector<T> temp3(temp1.size());
          viennacl::linalg::element_op(temp3, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp1, temp2));
          lhs += temp3;
        }
      };

      // generic x -= vec_expr1 .* vec_expr2:
      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const LHS, const RHS, op_element_binary<OP> > >
      {

        // x -= y .* z  or  x -= y ./ z
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          viennacl::vector<T> temp(proxy);
          lhs -= temp;
        }

        // x -= y .* vec_expr  or  x -= y ./ vec_expr
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_expression<const LHS2, const RHS2, OP2>, op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          vector<T> temp2(temp.size());
          viennacl::linalg::element_op(temp2, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(proxy.lhs(), temp));
          lhs -= temp2;
        }

        // x -= vec_expr .* z  or  x -= vec_expr ./ z
        template <typename LHS1, typename RHS1, typename OP1>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>, const vector_base<T>, op_element_binary<OP> > const & proxy)
        {
          vector<T> temp(proxy.lhs());
          vector<T> temp2(temp.size());
          viennacl::linalg::element_op(temp2, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp, proxy.rhs()));
          lhs -= temp2;
        }

        // x -= vec_expr .* vec_expr  or  x -= vec_expr ./ vec_expr
        template <typename LHS1, typename RHS1, typename OP1,
                  typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS1, const RHS1, OP1>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_binary<OP> > const & proxy)
        {
          vector<T> temp1(proxy.lhs());
          vector<T> temp2(proxy.rhs());
          vector<T> temp3(temp1.size());
          viennacl::linalg::element_op(temp3, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_binary<OP> >(temp1, temp2));
          lhs -= temp3;
        }
      };

      //////////////// unary expressions

      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_assign, vector_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x = OP(y)
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> > const & proxy)
        {
          viennacl::linalg::element_op(lhs, proxy);
        }

        // x = OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_unary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          viennacl::linalg::element_op(lhs, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> >(temp, temp));
        }
      };

      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x += OP(y)
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> > const & proxy)
        {
          vector<T> temp(proxy);
          lhs += temp;
        }

        // x += OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_unary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          viennacl::linalg::element_op(temp, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> >(temp, temp)); // inplace operation is safe here
          lhs += temp;
        }
      };

      template <typename T, typename LHS, typename RHS, typename OP>
      struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const LHS, const RHS, op_element_unary<OP> > >
      {
        // x -= OP(y)
        static void apply(vector_base<T> & lhs, vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> > const & proxy)
        {
          vector<T> temp(proxy);
          lhs -= temp;
        }

        // x -= OP(vec_expr)
        template <typename LHS2, typename RHS2, typename OP2>
        static void apply(vector_base<T> & lhs, vector_expression<const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  const vector_expression<const LHS2, const RHS2, OP2>,
                                                                  op_element_unary<OP> > const & proxy)
        {
          vector<T> temp(proxy.rhs());
          viennacl::linalg::element_op(temp, viennacl::vector_expression<const vector_base<T>, const vector_base<T>, op_element_unary<OP> >(temp, temp)); // inplace operation is safe here
          lhs -= temp;
        }
      };

    } // namespace detail

  } // namespace linalg

  /** \endcond */

} // namespace viennacl

#endif
