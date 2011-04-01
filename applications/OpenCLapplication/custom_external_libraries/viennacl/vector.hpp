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

/** @file vector.hpp
    @brief The vector type with operator-overloads and proxy classes is defined here. 
           Linear algebra operations such as norms and inner products are located in linalg/vector_operations.hpp
*/

#ifndef _VIENNACL_VECTOR_HPP_
#define _VIENNACL_VECTOR_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/scalar.hpp"
#include "viennacl/tools/tools.hpp"
#include "viennacl/tools/entry_proxy.hpp"
#include "viennacl/linalg/vector_operations.hpp"

namespace viennacl
{
    
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
      public:
        /** @brief Extracts the vector type from the two operands.
        */
        typedef typename viennacl::tools::VECTOR_EXTRACTOR<LHS, RHS>::ResultType    VectorType;
      
        vector_expression(LHS & lhs, RHS & rhs) : _lhs(lhs), _rhs(rhs) {}
        
        /** @brief Get left hand side operand
        */
        LHS & lhs() const { return _lhs; }
        /** @brief Get right hand side operand
        */
        RHS & rhs() const { return _rhs; }
        
        /** @brief Returns the size of the result vector */
        unsigned int size() const { return viennacl::tools::VECTOR_SIZE_DEDUCER<LHS, RHS, OP>::size(_lhs, _rhs); }
        
      private:
        /** @brief The left hand side operand */
        LHS & _lhs;
        /** @brief The right hand side operand */
        RHS & _rhs;
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
        typedef long                          difference_type;
        
        const_vector_iterator() {};
        /** @brief Constructor
        *   @param vec    The vector over which to iterate
        *   @param index  The starting index of the iterator
        */        
        const_vector_iterator(vector<SCALARTYPE, ALIGNMENT> const & vec,      unsigned int index)  : _elements(vec.handle()), _index(index) {};
        const_vector_iterator(viennacl::ocl::handle<cl_mem> const & elements, unsigned int index)  : _elements(elements), _index(index) {};

        
        value_type operator*(void) const 
        { 
           value_type result;
           result = entry_proxy<SCALARTYPE>(_index, _elements);
           return result;
        }
        self_type operator++(void) { ++_index; return *this; }
        self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }
        
        bool operator==(self_type const & other) const { return _index == other._index; }
        bool operator!=(self_type const & other) const { return _index != other._index; }
        
//        self_type & operator=(self_type const & other)
//        {
//           _index = other._index;
//           _elements = other._elements;
//           return *this;
//        }   

        difference_type operator-(self_type const & other) const { difference_type result = _index; return result - other._index; }
        self_type operator+(difference_type diff) const { return self_type(_elements, _index + diff); }
        
        unsigned int index() const { return _index; }
        viennacl::ocl::handle<cl_mem> const & handle() const { return _elements; }

      protected:
        /** @brief  The index of the entry the iterator is currently pointing to */
        viennacl::ocl::handle<cl_mem> _elements;
        unsigned int _index;
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
        typedef vector_iterator<SCALARTYPE, ALIGNMENT>    self_type;
      public:
        vector_iterator() : base_type(){};
        vector_iterator(viennacl::ocl::handle<cl_mem> const & elements, unsigned int index)  : base_type(elements, index) {};
        /** @brief Constructor
        *   @param vec    The vector over which to iterate
        *   @param index  The starting index of the iterator
        */        
        vector_iterator(vector<SCALARTYPE, ALIGNMENT> & vec, unsigned int index) : base_type(vec, index) {};
        vector_iterator(base_type const & b) : base_type(b) {};

        typename base_type::value_type operator*(void)  
        { 
           typename base_type::value_type result;
           result = entry_proxy<SCALARTYPE>(base_type::_index, base_type::_elements); 
           return result;
        }
        
        viennacl::ocl::handle<cl_mem> handle() { return base_type::_elements; }
        
        operator base_type() const
        {
          return base_type(base_type::_elements, base_type::_index);
        }
    };

    // forward definition in VCLForwards.h!
    /** @brief A vector class representing a linear memory sequence on the GPU. Inspired by boost::numeric::ublas::vector
    *
    *  This is the basic vector type of ViennaCL. It is similar to std::vector and boost::numeric::ublas::vector and supports various linear algebra operations.
    * By default, the internal length of the vector is padded to a multiple of 'ALIGNMENT' in order to speed up several GPU viennacl::ocl::kernels.
    *
    * @tparam SCALARTYPE  The floating point type, either 'float' or 'double'
    * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    class vector
    {
      
    public:
      typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;
      typedef const_vector_iterator<SCALARTYPE, ALIGNMENT>      const_iterator;
      typedef vector_iterator<SCALARTYPE, ALIGNMENT>            iterator;

      /** @brief Default constructor in order to be compatible with various containers.
      */
      vector() : _size(0) { viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::init();  }

      /** @brief An explicit constructor for the vector, allocating the given amount of memory (plus a padding specified by 'ALIGNMENT')
      *
      * @param vec_size   The length (i.e. size) of the vector.
      */
      explicit vector(unsigned int vec_size) : _size(vec_size)
      {
        viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::init(); 
        
        if (_size > 0)
          _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
        
        //force entries above _size to zero:
        if (_size < internal_size())
        {
          std::vector<SCALARTYPE> temp(internal_size() - _size);
          cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), _elements, CL_TRUE, sizeof(SCALARTYPE)*_size, sizeof(SCALARTYPE)*(internal_size() - _size), &(temp[0]), 0, NULL, NULL);
          //assert(err == CL_SUCCESS);
          VIENNACL_ERR_CHECK(err);
        }
      }

      /** @brief Create a vector from existing OpenCL memory
      *
      * Note: The provided memory must take an eventual ALIGNMENT into account, i.e. existing_mem must be at least of size internal_size()!
      * This is trivially the case with the default alignment, but should be considered when using vector<> with an alignment parameter not equal to 1.
      *
      * @param existing_mem   An OpenCL handle representing the memory
      * @param vec_size       The size of the vector. 
      */
      explicit vector(cl_mem existing_mem, unsigned int vec_size) : _size(vec_size),  _elements(existing_mem)
      {
        _elements.inc();  //prevents that the user-provided memory is deleted once the vector object is destroyed.
      }
      
      /** @brief The copy constructor
      *
      * Entries of 'vec' are directly copied to this vector.
      */
      vector(const vector<SCALARTYPE, ALIGNMENT> & vec) :
        _size(vec.size())
      {
        viennacl::linalg::kernels::vector<SCALARTYPE, 1>::init(); 
        
        if (size() != 0)
        {
          _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*internal_size());
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), vec.handle(), _elements, 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
          //assert(err == CL_SUCCESS);
          VIENNACL_ERR_CHECK(err);
        }
      }

      /** @brief Assignment operator. This vector is resized if 'vec' is of a different size.
      */
      vector<SCALARTYPE, ALIGNMENT> & operator=(const vector<SCALARTYPE, ALIGNMENT> & vec)
      {
        resize(vec.size());
        if (size() != 0)
        {
          cl_int err;
          err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), vec.handle(), _elements, 0, 0, sizeof(SCALARTYPE)*internal_size(), 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
        }
        return *this;
      }


      /** @brief Implementation of the operation v1 = alpha * v2, where alpha is a GPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename VectorType>   //use template to cover const/non-const of VectorType:
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< VectorType,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy)
      {
        resize(proxy.lhs().size());
        //std::cout << "vector::operator=(vec_times_scalar_proxy)" << std::endl; 
        viennacl::linalg::mult(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      /** @brief Implementation of the operation v1 = alpha * v2, where alpha is a CPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename VectorType>   //use template to cover const/non-const of VectorType:
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< VectorType,
                                                                           const SCALARTYPE,
                                                                           op_prod> & proxy)
      {
        resize(proxy.lhs().size());
        viennacl::linalg::mult(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      /** @brief Implementation of the operation v1 = v2 / alpha, where alpha is a GPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename VectorType>   //use template to cover const/non-const of VectorType:
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< VectorType,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_div> & proxy)
      {
        resize(proxy.lhs().size());
        //std::cout << "vector::operator=(vec_times_scalar_proxy)" << std::endl; 
        viennacl::linalg::divide(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      /** @brief Implementation of the operation v1 = v2 / alpha, where alpha is a CPU scalar
      *
      * @param proxy  An expression template proxy class.
      */
      template <typename VectorType>   //use template to cover const/non-const of VectorType:
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< VectorType,
                                                                           const SCALARTYPE,
                                                                           op_div> & proxy)
      {
        resize(proxy.lhs().size());
        //std::cout << "vector::operator=(vec_times_scalar_proxy)" << std::endl; 
        viennacl::linalg::mult(proxy.lhs(), static_cast<SCALARTYPE>(1.0) / proxy.rhs(), *this);
        return *this;
      }

      //v1 = v2 + v3; 
      /** @brief Implementation of the operation v1 = v2 + v3
      *
      * @param proxy  An expression template proxy class.
      */
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                           vector<SCALARTYPE, ALIGNMENT>,
                                                                           op_add> & proxy)
      {
        resize(proxy.lhs().size());
        //std::cout << "vector::operator=(vec_times_scalar_proxy)" << std::endl; 
        viennacl::linalg::add(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }
      
      //v1 = v2 - v3; 
      /** @brief Implementation of the operation v1 = v2 - v3
      *
      * @param proxy  An expression template proxy class.
      */
      vector<SCALARTYPE, ALIGNMENT> & operator = (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                           vector<SCALARTYPE, ALIGNMENT>,
                                                                           op_sub> & proxy)
      {
        resize(proxy.lhs().size());
        //std::cout << "vector::operator=(vec_times_scalar_proxy)" << std::endl; 
        viennacl::linalg::sub(proxy.lhs(), proxy.rhs(), *this);
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
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator=(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                op_prod> & proxy);

      /** @brief Operator overload for v1 += A * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator+=(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);
                                                
      /** @brief Operator overload for v1 -= A * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator-=(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);

      /** @brief Operator overload for v1 + A * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator+(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      /** @brief Operator overload for v1 - A * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator-(const vector_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      //transposed_matrix_proxy:
      /** @brief Operator overload for v1 = trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator=(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                  const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                  op_trans >,
                                                                         const vector<SCALARTYPE, ALIGNMENT>,
                                                                         op_prod> & proxy);

      /** @brief Operator overload for v1 += trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator+=(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                   const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                   op_trans >,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);
                                                
      /** @brief Operator overload for v1 -= trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator-=(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                   const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                   op_trans >,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);

      /** @brief Operator overload for v1 + trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator+(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                op_trans >,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      /** @brief Operator overload for v1 - trans(A) * v2, where v1, v2 are vectors and A is a dense matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <typename F, unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator-(const vector_expression< const matrix_expression< const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                const matrix<SCALARTYPE, F, MAT_ALIGNMENT>,
                                                                                                op_trans >,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);
                                                                       
                                                                       
                                                                       
      //////////// compressed_matrix<>
      /** @brief Operator overload for v1 = A * v2, where v1, v2 are vectors and A is a sparse matrix of type compressed_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator=(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                op_prod> & proxy);

      /** @brief Operator overload for v1 += A * v2, where v1, v2 are vectors and A is a sparse matrix of type compressed_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator+=(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);
                                                
      /** @brief Operator overload for v1 -= A * v2, where v1, v2 are vectors and A is a sparse matrix of type compressed_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator-=(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);

      /** @brief Operator overload for v1 + A * v2, where v1, v2 are vectors and A is a sparse matrix of type compressed_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator+(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      /** @brief Operator overload for v1 - A * v2, where v1, v2 are vectors and A is a sparse matrix of type compressed_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator-(const vector_expression< const compressed_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);


      //coordinate_matrix<>
      /** @brief Operator overload for v1 = A * v2, where v1, v2 are vectors and A is a sparse matrix of type coordinate_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator=(const vector_expression< const coordinate_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                const vector<SCALARTYPE, ALIGNMENT>,
                                                op_prod> & proxy);

      /** @brief Operator overload for v1 += A * v2, where v1, v2 are vectors and A is a sparse matrix of type coordinate_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator+=(const vector_expression< const coordinate_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);
                                                
      /** @brief Operator overload for v1 -= A * v2, where v1, v2 are vectors and A is a sparse matrix of type coordinate_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> & operator-=(const vector_expression< const coordinate_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                          const vector<SCALARTYPE, ALIGNMENT>,
                                                                          op_prod> & proxy);

      /** @brief Operator overload for v1 + A * v2, where v1, v2 are vectors and A is a sparse matrix of type coordinate_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator+(const vector_expression< const coordinate_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      /** @brief Operator overload for v1 - A * v2, where v1, v2 are vectors and A is a sparse matrix of type coordinate_matrix.
      *
      * @param proxy An expression template proxy class
      */
      template <unsigned int MAT_ALIGNMENT>
      vector<SCALARTYPE, ALIGNMENT> operator-(const vector_expression< const coordinate_matrix<SCALARTYPE, MAT_ALIGNMENT>,
                                                                       const vector<SCALARTYPE, ALIGNMENT>,
                                                                       op_prod> & proxy);

      ///////////////////////////// Matrix Vector interaction end ///////////////////////////////////

      //enlarge or reduce allocated memory and set unused memory to zero
      /** @brief Resizes the allocated memory for the vector. Pads the memory to be a multiple of 'ALIGNMENT'
      *
      *  @param new_size  The new size of the vector
      *  @param preserve  If true, old entries of the vector are preserved, otherwise eventually discarded.
      */
      void resize(unsigned int new_size, bool preserve = true)
      {
        assert(new_size > 0);
        
        if (new_size != _size)
        {
          unsigned int new_internal_size = viennacl::tools::roundUpToNextMultiple<unsigned int>(new_size, ALIGNMENT);
        
          std::vector<SCALARTYPE> temp(_size);
          if (preserve && _size > 0)
            fast_copy(*this, temp);
          temp.resize(new_size);  //drop all entries above new_size
          temp.resize(new_internal_size); //enlarge to fit new internal size
          
          if (new_internal_size != internal_size())
          {
            _elements = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(SCALARTYPE)*new_internal_size);
          }
          
          fast_copy(temp, *this);
          _size = new_size;
        }
        
      }
      

      //read-write access to an element of the vector
      /** @brief Read-write access to a single element of the vector
      */
      entry_proxy<SCALARTYPE> operator()(unsigned int index)
      {
        return entry_proxy<SCALARTYPE>(index, _elements);
      }

      /** @brief Read-write access to a single element of the vector
      */
      entry_proxy<SCALARTYPE> operator[](unsigned int index)
      {
        return entry_proxy<SCALARTYPE>(index, _elements);
      }


      /** @brief Read access to a single element of the vector
      */
      scalar<SCALARTYPE> operator()(unsigned int index) const
      {
        scalar<SCALARTYPE> tmp;
        cl_int err;
        err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), _elements, tmp.handle(), sizeof(SCALARTYPE)*index, 0, sizeof(SCALARTYPE), 0, NULL, NULL);
        //assert(err == CL_SUCCESS);
        VIENNACL_ERR_CHECK(err);
        return tmp;
      }
      
      /** @brief Read access to a single element of the vector
      */
      scalar<SCALARTYPE> operator[](unsigned int index) const
      {
        return operator()(index);
      }
      
      /** @brief Inplace addition of a vector
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector<SCALARTYPE, ALIGNMENT> & vec)
      {
        viennacl::linalg::inplace_add(*this, vec);
        return *this;
      }

      /** @brief Inplace addition of a scaled vector, i.e. v1 += alpha * v2, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace addition of a scaled vector, i.e. v1 += alpha * v2, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace addition of a scaled vector, i.e. v1 += alpha * v2, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                           const SCALARTYPE,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace addition of a scaled vector, i.e. v1 += alpha * v2, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                           const SCALARTYPE,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace addition of a scaled vector, i.e. v1 += alpha * v2, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator += (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_div> & proxy)
      {
        viennacl::linalg::inplace_div_add(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }



      /** @brief Inplace subtraction of a vector
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector<SCALARTYPE, ALIGNMENT> & vec)
      {
        viennacl::linalg::inplace_sub(*this, vec);
        return *this;
      }

      /** @brief Inplace subtraction of a scaled vector, i.e. v1 -= alpha * v2, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_sub(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace subtraction of a scaled vector, i.e. v1 -= alpha * v2, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                           const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_sub(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }

      /** @brief Inplace subtraction of a scaled vector, i.e. v1 -= alpha * v2, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                            const SCALARTYPE,
                                                                            op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), -proxy.rhs());
        return *this;
      }

      /** @brief Inplace subtraction of a scaled vector, i.e. v1 -= alpha * v2, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                            const SCALARTYPE,
                                                                            op_prod> & proxy)
      {
        viennacl::linalg::inplace_mul_add(*this, proxy.lhs(), -proxy.rhs());
        return *this;
      }
      
      /** @brief Inplace subtraction of a scaled vector, i.e. v1 -= alpha * v2, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> & operator -= (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                            const scalar<SCALARTYPE>,
                                                                            op_div> & proxy)
      {
        viennacl::linalg::inplace_div_sub(*this, proxy.lhs(), proxy.rhs());
        return *this;
      }
      
      
      

      /** @brief Scales this vector by a CPU scalar value
      */
      vector<SCALARTYPE, ALIGNMENT> & operator *= (SCALARTYPE val)
      {
        viennacl::linalg::inplace_mult(*this, val);
        return *this;
      }

      /** @brief Scales this vector by a GPU scalar value
      */
      vector<SCALARTYPE, ALIGNMENT> & operator *= (scalar<SCALARTYPE> const & gpu_val)
      {
        viennacl::linalg::inplace_mult(*this, gpu_val);
        return *this;
      }

      /** @brief Scales this vector by a CPU scalar value
      */
      vector<SCALARTYPE, ALIGNMENT> & operator /= (SCALARTYPE val)
      {
        viennacl::linalg::inplace_mult(*this, static_cast<SCALARTYPE>(1) / val);
        return *this;
      }
      
      /** @brief Scales this vector by a CPU scalar value
      */
      vector<SCALARTYPE, ALIGNMENT> & operator /= (scalar<SCALARTYPE> const & gpu_val)
      {
        viennacl::linalg::inplace_divide(*this, gpu_val);
        return *this;
      }
      
      
      
      // free addition
      
      /** @brief Adds up two vectors
      */
      vector<SCALARTYPE, ALIGNMENT> operator + (const vector<SCALARTYPE, ALIGNMENT> & vec) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(internal_size());
        viennacl::linalg::add(*this, vec, result);
        return result;
      }
      
      /** @brief Adds up two vectors, i.e. result = v1 + v2 * alpha, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator + (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                         const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        viennacl::linalg::mul_add(proxy.lhs(), proxy.rhs(), *this, result);
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 + v2 * alpha, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator + (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                         const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        viennacl::linalg::mul_add(proxy.lhs(), proxy.rhs(), *this, result);
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 + v2 * alpha, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator + (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                         const SCALARTYPE,
                                                                         op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        viennacl::linalg::mul_add(proxy.lhs(), proxy.rhs(), *this, result);
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 + v2 * alpha, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator + (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                         const SCALARTYPE,
                                                                         op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        viennacl::linalg::mul_add(proxy.lhs(), proxy.rhs(), *this, result);
        return result;
      }


      //free subtraction:
      /** @brief Implementation of    result = v1 - v2
      */
      vector<SCALARTYPE, ALIGNMENT> operator - (const vector<SCALARTYPE, ALIGNMENT> & vec) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        viennacl::linalg::sub(*this, vec, result);
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 - v2 * alpha, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator - (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                         const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        result = *this;
        viennacl::linalg::inplace_mul_sub(result, proxy.lhs(), proxy.rhs());
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 - v2 * alpha, where alpha is a GPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator - (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                         const scalar<SCALARTYPE>,
                                                                           op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        result = *this;
        viennacl::linalg::inplace_mul_sub(result, proxy.lhs(), proxy.rhs());
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 - v2 * alpha, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator - (const vector_expression< vector<SCALARTYPE, ALIGNMENT>,
                                                                         const SCALARTYPE,
                                                                         op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        result = *this;
        viennacl::linalg::inplace_mul_add(result, proxy.lhs(), -proxy.rhs());
        return result;
      }

      /** @brief Adds up two vectors, i.e. result = v1 - v2 * alpha, where alpha is a CPU scalar
      */
      vector<SCALARTYPE, ALIGNMENT> operator - (const vector_expression< const vector<SCALARTYPE, ALIGNMENT>,
                                                                         const SCALARTYPE,
                                                                         op_prod> & proxy) const
      {
        vector<SCALARTYPE, ALIGNMENT> result(_size);
        result = *this;
        viennacl::linalg::inplace_mul_add(result, proxy.lhs(), -proxy.rhs());
        return result;
      }

      
      //free multiplication
      /** @brief Scales the vector by a CPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const SCALARTYPE, op_prod> 
      operator * (SCALARTYPE value) const
      {
        return vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const SCALARTYPE, op_prod>(*this, value);
      }

      /** @brief Scales the vector by a GPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const scalar<SCALARTYPE>, op_prod> 
      operator * (scalar<SCALARTYPE> const & value) const
      {
        return vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const scalar<SCALARTYPE>, op_prod>(*this, value);
      }

      //free division
      /** @brief Scales the vector by a CPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const SCALARTYPE, op_div> 
      operator / (SCALARTYPE value) const
      {
        return vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const SCALARTYPE, op_div>(*this, value);
      }

      /** @brief Scales the vector by a GPU scalar 'alpha' and returns an expression template
      */
      vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const scalar<SCALARTYPE>, op_div> 
      operator / (scalar<SCALARTYPE> const & value) const
      {
        return vector_expression< const vector<SCALARTYPE, ALIGNMENT>, const scalar<SCALARTYPE>, op_div>(*this, value);
      }
      
      
      //// iterators:
      /** @brief Returns an iterator pointing to the beginning of the vector  (STL like)*/
      iterator begin()
      {
        return iterator(*this, 0);
      }

      /** @brief Returns an iterator pointing to the end of the vector (STL like)*/
      iterator end()
      {
        return iterator(*this, size());
      }

      /** @brief Returns a const-iterator pointing to the beginning of the vector (STL like)*/
      const_iterator begin() const
      {
        return const_iterator(*this, 0);
      }

      /** @brief Returns a const-iterator pointing to the end of the vector (STL like)*/
      const_iterator end() const
      {
        return const_iterator(*this, size());
      }

      /** @brief Swaps the entries of the two vectors
      */
      vector<SCALARTYPE, ALIGNMENT> & swap(vector<SCALARTYPE, ALIGNMENT> & other)
      {
        swap(*this, other);
        return *this;
      };
      
      /** @brief Swaps the handles of two vectors by swapping the OpenCL handles only, no data copy
      */ 
      vector<SCALARTYPE, ALIGNMENT> & fast_swap(vector<SCALARTYPE, ALIGNMENT> & other) 
      { 
        assert(this->_size == other._size); 
        this->_elements.swap(other._elements); 
        return *this; 
      };       
      
      /** @brief Returns the length of the vector (cf. std::vector)
      */
      unsigned int size() const { return _size; }
      
      /** @brief Returns the maximum possible size of the vector, which is given by 128 MByte due to limitations by OpenCL.
      */
      unsigned int max_size() const
      {
        return (128*1024*1024) / sizeof(SCALARTYPE);  //128 MB is maximum size of memory chunks in OpenCL!
      }
      /** @brief Returns the internal length of the vector, which is given by size() plus the extra memory due to padding the memory with zeros up to a multiple of 'ALIGNMENT'
      */
      unsigned int internal_size() const { return viennacl::tools::roundUpToNextMultiple<unsigned int>(_size, ALIGNMENT); }
      
      /** @brief Returns true is the size is zero */
      bool empty() { return _size == 0; }
      
      /** @brief Returns the OpenCL memory viennacl::ocl::handle. Typically used for launching compute viennacl::ocl::kernels */
      const viennacl::ocl::handle<cl_mem> & handle() const { return _elements; }

      /** @brief Resets all entries to zero. Does not change the size of the vector.
      */
      void clear()
      {
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "clear");
        
        viennacl::ocl::enqueue(k(_elements, internal_size()));
      }
      //void swap(vector & other){}
      

      //TODO: Think about implementing the following public member functions
      //void insert_element(unsigned int i, SCALARTYPE val){}
      //void erase_element(unsigned int i){}
      
    private:
      unsigned int _size;
      viennacl::ocl::handle<cl_mem> _elements;
    }; //vector
    

    //
    //////////////////// Copy from GPU to CPU //////////////////////////////////
    //
    
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
      assert(gpu_end - gpu_begin >= 0);
      if (gpu_end - gpu_begin != 0)
      {
        std::vector<SCALARTYPE> temp_buffer(gpu_end - gpu_begin);
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(),
                                         gpu_begin.handle(), CL_TRUE, 0, 
                                         sizeof(SCALARTYPE)*(gpu_end - gpu_begin),
                                         &(temp_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
        
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
      copy(const_vector_iterator<SCALARTYPE, ALIGNMENT>(gpu_begin),
           const_vector_iterator<SCALARTYPE, ALIGNMENT>(gpu_end),
           cpu_begin);
    }
    
    /** @brief Transfer from a gpu vector to a cpu vector. Convenience wrapper for viennacl::linalg::copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
    *
    * @param gpu_vec    A gpu vector
    * @param cpu_vec    The cpu vector. Type requirements: Output iterator can be obtained via member function .begin()
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPUVECTOR>
    void copy(vector<SCALARTYPE, ALIGNMENT> const & gpu_vec,
              CPUVECTOR & cpu_vec )
    {
      viennacl::copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
    }

    //from gpu to cpu. Type assumption: cpu_vec lies in a linear memory chunk
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
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(),
                                         gpu_begin.handle(), CL_TRUE, 0,
                                         sizeof(SCALARTYPE)*(gpu_end - gpu_begin),
                                         &(*cpu_begin), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        viennacl::ocl::get_queue().finish();
      }
    }

    /** @brief Transfer from a gpu vector to a cpu vector. Convenience wrapper for viennacl::linalg::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
    *
    * @param gpu_vec    A gpu vector.
    * @param cpu_vec    The cpu vector. Type requirements: Output iterator can be obtained via member function .begin()
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPUVECTOR>
    void fast_copy(vector<SCALARTYPE, ALIGNMENT> const & gpu_vec,
                   CPUVECTOR & cpu_vec )
    {
      viennacl::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
    }



    #ifdef VIENNACL_HAVE_EIGEN
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
      if (cpu_begin != cpu_end)
      {
        //we require that the size of the gpu_vector is larger or equal to the cpu-size
        std::vector<SCALARTYPE> temp_buffer(cpu_end - cpu_begin);
        std::copy(cpu_begin, cpu_end, temp_buffer.begin());
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(),
                                          gpu_begin.handle(), CL_TRUE, 0,
                                          sizeof(SCALARTYPE)*(cpu_end - cpu_begin),
                                          &(temp_buffer[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
      }
    }

    // for things like copy(std_vec.begin(), std_vec.end(), vcl_vec.begin() + 1);
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
    void copy(CPU_ITERATOR const & cpu_begin,
              CPU_ITERATOR const & cpu_end,
              const_vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin)
    {
      copy(cpu_begin, cpu_end, vector_iterator<SCALARTYPE, ALIGNMENT>(gpu_begin));
    }

    /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
    *
    * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
    * @param gpu_vec    The gpu vector.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPUVECTOR>
    void copy(const CPUVECTOR & cpu_vec, vector<SCALARTYPE, ALIGNMENT> & gpu_vec)
    {
      viennacl::copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
    }

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
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPU_ITERATOR>
    void fast_copy(CPU_ITERATOR const & cpu_begin,
                   CPU_ITERATOR const & cpu_end,
                   vector_iterator<SCALARTYPE, ALIGNMENT> gpu_begin)
    {
      if (cpu_begin != cpu_end)
      {
        //we require that the size of the gpu_vector is larger or equal to the cpu-size
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), 
                                          gpu_begin.handle(), CL_TRUE, 0, 
                                          sizeof(SCALARTYPE)*(cpu_end - cpu_begin), &(*cpu_begin), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
      }
    }


    /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
    *
    * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
    * @param gpu_vec    The gpu vector.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename CPUVECTOR>
    void fast_copy(const CPUVECTOR & cpu_vec, vector<SCALARTYPE, ALIGNMENT> & gpu_vec)
    {
      viennacl::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
    }

    #ifdef VIENNACL_HAVE_EIGEN
    template <unsigned int ALIGNMENT>
    void copy(Eigen::VectorXf const & eigen_vec,
              vector<float, ALIGNMENT> & gpu_vec)
    {
      std::vector<float> entries(eigen_vec.size());
      for (size_t i = 0; i<entries.size(); ++i)
        entries[i] = eigen_vec(i);
      viennacl::fast_copy(entries.begin(), entries.end(), gpu_vec.begin());
    }
    
    template <unsigned int ALIGNMENT>
    void copy(Eigen::VectorXd const & eigen_vec,
              vector<double, ALIGNMENT> & gpu_vec)
    {
      std::vector<double> entries(eigen_vec.size());
      for (size_t i = 0; i<entries.size(); ++i)
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
      if (gpu_src_begin != gpu_src_end)
      {
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(),
                                          gpu_src_begin.handle(),  //src handle
                                          gpu_dest_begin.handle(), //dest handle
                                          sizeof(SCALARTYPE) * gpu_src_begin.index(), //src offset
                                          sizeof(SCALARTYPE) * gpu_dest_begin.index(), //dest offset
                                          sizeof(SCALARTYPE) * (gpu_src_end.index() - gpu_src_begin.index()), //data length
                                          0, //don't know -> check!! (something related to increment?)
                                          NULL, NULL);
        VIENNACL_ERR_CHECK(err);
      }
    }

    /** @brief Copy (parts of a) GPU vector to another GPU vector
    *
    * @param gpu_src_begin   GPU iterator pointing to the beginning of the gpu vector (STL-like)
    * @param gpu_src_end     GPU iterator pointing to the end of the vector (STL-like)
    * @param gpu_dest_begin  Output iterator for the gpu vector. The gpu vector must be at least as long as the cpu vector!
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT_SRC, unsigned int ALIGNMENT_DEST>
    void copy(const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_begin,
              const_vector_iterator<SCALARTYPE, ALIGNMENT_SRC> const & gpu_src_end,
              const_vector_iterator<SCALARTYPE, ALIGNMENT_DEST> gpu_dest_begin)
    {
      copy(gpu_src_begin, gpu_src_end, vector_iterator<SCALARTYPE, ALIGNMENT_DEST>(gpu_dest_begin));
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
    * @param s    STL output stream
    * @param val  The vector that should be printed
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream & s, vector<SCALARTYPE,ALIGNMENT> const & val)
    {
      viennacl::ocl::get_queue().finish();
      std::vector<SCALARTYPE> tmp(val.size());
      copy(val.begin(), val.end(), tmp.begin());
      std::cout << "[" << val.size() << "](";
      for (unsigned int i=0; i<val.size(); ++i)
      {
        if (i > 0)
          s << ",";
        s << tmp[i];
      }
      std::cout << ")" << std::endl;
      return s;
    }

    /** @brief Swaps the contents of two vectors, data is copied
    *
    * @param vec1   The first vector
    * @param vec2   The second vector
    * @param NUM_THREADS  The number of threads per work group
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    void swap(viennacl::vector<SCALARTYPE, ALIGNMENT> & vec1,
              viennacl::vector<SCALARTYPE, ALIGNMENT> & vec2, 
              size_t NUM_THREADS = 0)
    {
      assert(vec1.size() == vec2.size());

      viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::vector<SCALARTYPE, ALIGNMENT>::program_name(), "swap");

      viennacl::ocl::enqueue(k(vec1, vec2, vec1.size()));
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
    
    
    
    ////////// operations /////////////
    /** @brief Operator overload for the expression alpha * v1, where alpha is a host scalar (float or double) and v1 is a ViennaCL vector.
    *
    * @param value   The host scalar (float or double)
    * @param vec     A ViennaCL vector
    */
    template <typename SCALARTYPE, unsigned int A>
    vector_expression< const vector<SCALARTYPE, A>, const SCALARTYPE, op_prod> operator * (SCALARTYPE const & value, vector<SCALARTYPE, A> const & vec)
    {
      return vector_expression< const vector<SCALARTYPE, A>, const SCALARTYPE, op_prod>(vec, value);
    }

    /** @brief Operator overload for the expression alpha * v1, where alpha is a ViennaCL scalar (float or double) and v1 is a ViennaCL vector.
    *
    * @param value   The ViennaCL scalar
    * @param vec     A ViennaCL vector
    */
    template <typename SCALARTYPE, unsigned int A>
    vector_expression< const vector<SCALARTYPE, A>, const scalar<SCALARTYPE>, op_prod> operator * (scalar<SCALARTYPE> const & value, vector<SCALARTYPE, A> const & vec)
    {
        return vector_expression< const vector<SCALARTYPE, A>, const scalar<SCALARTYPE>, op_prod>(vec, value);
    }


    //addition and subtraction of two vector_expressions:
    /** @brief Operator overload for the addition of two vector expressions.
    *
    * @param proxy1  Left hand side vector expression
    * @param proxy2  Right hand side vector expression
    */
    template <typename LHS1, typename RHS1, typename OP1,
              typename LHS2, typename RHS2, typename OP2>
    typename vector_expression< LHS1, RHS1, OP1>::VectorType
    operator + (vector_expression< LHS1, RHS1, OP1> const & proxy1,
                vector_expression< LHS2, RHS2, OP2> const & proxy2)
    {
      assert(proxy1.size() == proxy2.size());
      typename vector_expression< LHS1, RHS1, OP1>::VectorType result(proxy1.size());
      result = proxy1;
      result += proxy2;
      return result;
    }

    /** @brief Operator overload for the subtraction of two vector expressions.
    *
    * @param proxy1  Left hand side vector expression
    * @param proxy2  Right hand side vector expression
    */
    template <typename LHS1, typename RHS1, typename OP1,
              typename LHS2, typename RHS2, typename OP2>
    typename vector_expression< LHS1, RHS1, OP1>::VectorType
    operator - (vector_expression< LHS1, RHS1, OP1> const & proxy1,
                vector_expression< LHS2, RHS2, OP2> const & proxy2)
    {
      assert(proxy1.size() == proxy2.size());
      typename vector_expression< LHS1, RHS1, OP1>::VectorType result(proxy1.size());
      result = proxy1;
      result -= proxy2;
      return result;
    }
    
    //////////// one vector expression from left /////////////////////////////////////////
    
    /** @brief Operator overload for the addition of a vector expression from the left, e.g. alpha * vec1 + vec2. Here, alpha * vec1 is wrapped into a vector_expression and then added to vec2.
    *
    * @param proxy   Left hand side vector expression
    * @param vec     Right hand side vector
    */
    template <typename SCALARTYPE, unsigned int A, typename LHS, typename RHS, typename OP>
    vector<SCALARTYPE, A> operator + (vector_expression< LHS, RHS, OP> const & proxy,
                                      vector<SCALARTYPE, A> const & vec)
    {
      assert(proxy.size() == vec.size());
      vector<SCALARTYPE, A> result(vec.size());
      result = proxy;
      result += vec;
      return result;
    }

    /** @brief Operator overload for the subtraction of a vector expression from the left, e.g. alpha * vec1 + vec2. Here, alpha * vec1 is wrapped into a vector_expression and then added to vec2.
    *
    * @param proxy   Left hand side vector expression
    * @param vec     Right hand side vector
    */
    template <typename SCALARTYPE, unsigned int A, typename LHS, typename RHS, typename OP>
    vector<SCALARTYPE, A> operator - (vector_expression< LHS, RHS, OP> const & proxy,
                                      vector<SCALARTYPE, A> const & vec)
    {
      assert(proxy.size() == vec.size());
      vector<SCALARTYPE, A> result(vec.size());
      result = proxy;
      result -= vec;
      return result;
    }


    /** @brief Operator overload for the multiplication of a vector expression with a scalar from the right, e.g. (beta * vec1) * alpha. Here, beta * vec1 is wrapped into a vector_expression and then multiplied with alpha from the right.
    *
    * @param proxy   Left hand side vector expression
    * @param val     Right hand side scalar
    */
    template <typename SCALARTYPE, typename LHS, typename RHS, typename OP>
    vector<SCALARTYPE> operator * (vector_expression< LHS, RHS, OP> const & proxy,
                                   scalar<SCALARTYPE> const & val)
    {
      vector<SCALARTYPE> result(proxy.size());
      result = proxy;
      result *= val;
      return result;
    }

    /** @brief Operator overload for the division of a vector expression by a scalar from the right, e.g. (beta * vec1) / alpha. Here, beta * vec1 is wrapped into a vector_expression and then divided by alpha.
    *
    * @param proxy   Left hand side vector expression
    * @param val     Right hand side scalar
    */
    template <typename SCALARTYPE, typename LHS, typename RHS, typename OP>
    vector<SCALARTYPE> operator / (vector_expression< LHS, RHS, OP> const & proxy,
                                      scalar<SCALARTYPE> const & val)
    {
      vector<SCALARTYPE> result(proxy.size());
      result = proxy;
      result /= val;
      return result;
    }


    //////////// one vector expression from right (on scalar) ///////////////////////
    
    /** @brief Operator overload for the multiplication of a vector expression with a ViennaCL scalar from the left, e.g. alpha * (beta * vec1). Here, beta * vec1 is wrapped into a vector_expression and then multiplied with alpha from the left.
    *
    * @param val     Right hand side scalar
    * @param proxy   Left hand side vector expression
    */
    template <typename SCALARTYPE, typename LHS, typename RHS, typename OP>
    vector<SCALARTYPE> operator * (scalar<SCALARTYPE> const & val,
                                   vector_expression< LHS, RHS, OP> const & proxy)
    {
      vector<SCALARTYPE> result(proxy.size());
      result = proxy;
      result *= val;
      return result;
    }
    
    /** @brief Operator overload for the multiplication of a vector expression with a host scalar (float or double) from the left, e.g. alpha * (beta * vec1). Here, beta * vec1 is wrapped into a vector_expression and then multiplied with alpha from the left.
    *
    * @param val     Right hand side scalar
    * @param proxy   Left hand side vector expression
    */
    template <typename SCALARTYPE, typename LHS, typename RHS, typename OP>
    viennacl::vector<SCALARTYPE> operator * (SCALARTYPE val,
                                   viennacl::vector_expression< LHS, RHS, OP> const & proxy)
    {
      viennacl::vector<SCALARTYPE> result(proxy.size());
      result = proxy;
      result *= val;
      return result;
    }

}

#endif
