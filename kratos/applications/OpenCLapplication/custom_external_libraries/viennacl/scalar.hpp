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

#ifndef _VIENNACL_SCALAR_HPP_
#define _VIENNACL_SCALAR_HPP_

/** @file scalar.hpp
    @brief Implementation of the ViennaCL scalar class
*/

#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/linalg/kernels/scalar_kernels.h"

#include <iostream>

namespace viennacl
{
    /** @brief A proxy for scalar expressions (e.g. from inner vector products)
    * 
    * assumption: dim(LHS) >= dim(RHS), where dim(scalar) = 0, dim(vector) = 1 and dim(matrix = 2)
    * @tparam LHS   The left hand side operand
    * @tparam RHS   The right hand side operand
    * @tparam OP    The operation tag
    */
    template <typename LHS, typename RHS, typename OP>
    class scalar_expression
    {
        typedef typename LHS::value_type          DummyType; //Visual C++ 2005 does not allow to write LHS::value_type::value_type
      public:
        typedef typename DummyType::value_type    ScalarType;
        
        scalar_expression(LHS & lhs, RHS & rhs) : _lhs(lhs), _rhs(rhs) {}
        
        /** @brief Returns the left hand side operand */
        LHS & get_lhs() const { return _lhs; }
        /** @brief Returns the left hand side operand */
        RHS & get_rhs() const { return _rhs; }

        /** @brief Conversion operator to a ViennaCL scalar */
        operator ScalarType () const
        {
          viennacl::scalar<ScalarType> temp;
          temp = *this;
          return temp;
        }

      private:
        LHS & _lhs;
        RHS & _rhs;
    };
    
    /** @brief This class represents a single scalar value on the GPU and behaves mostly like a built-in scalar type like float or double.
    *
    * Since every read and write operation requires a CPU->GPU or GPU->CPU transfer, this type should be used with care.
    * The advantage of this type is that the GPU command queue can be filled without blocking read operations.
    *
    * @tparam TYPE  Either float or double. Checked at compile time.
    */
    template<class TYPE>
    class scalar
    {
    public:
      /** @brief Returns the underlying host scalar type. */
      typedef typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<TYPE>::ResultType   value_type;
      
      /** @brief Allocates the memory for the scalar, but does not set it to zero. */
      scalar()
      {
        viennacl::linalg::kernels::scalar<TYPE, 1>::init(); 
        val_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(TYPE));
      }
      /** @brief Allocates the memory for the scalar and sets it to the supplied value. */
      scalar(TYPE val)
      {
        viennacl::linalg::kernels::scalar<TYPE, 1>::init(); 
        val_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(TYPE), &val);
      }
      
      /** @brief Wraps an existing memory entry into a scalar 
      *
      * @param mem    The OpenCL memory handle
      * @param size   Ignored - Only necessary to avoid ambiguities. Users are advised to set this parameter to '1'.
      */
      explicit scalar(cl_mem mem, size_t size) : val_(mem) { val_.inc(); }

      /** @brief Allocates memory for the scalar and sets it to the result of supplied expression. */
      template <typename T1, typename T2, typename OP>
      scalar(scalar_expression<T1, T2, OP> const & proxy)
      {
        viennacl::linalg::kernels::scalar<TYPE, 1>::init(); 
        val_ = viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(TYPE));
        *this = proxy;
      }

      //copy constructor
      /** @brief Copy constructor. Allocates new memory for the scalar and copies the value of the supplied scalar */
      scalar(const scalar & other) : val_(viennacl::ocl::current_context().create_memory(CL_MEM_READ_WRITE, sizeof(TYPE)))
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), other.handle(), handle(), 0, 0, sizeof(TYPE), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
      }

      /** @brief Reads the value of the scalar from the GPU and returns the float or double value. */
      operator TYPE() const
      {
        TYPE tmp;
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle(), val_, CL_TRUE, 0, sizeof(TYPE), &tmp, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        return tmp;
      } 
      
      /** @brief Assigns a vector entry. */
      scalar<TYPE> & operator= (entry_proxy<TYPE> const & other)
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), other.handle(), handle(), other.index() * sizeof(TYPE), 0, sizeof(TYPE), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        return *this;
      }

      /** @brief Assigns the value from another scalar. */
      scalar<TYPE> & operator= (scalar<TYPE> const & other)
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle(), other.handle(), handle(), 0, 0, sizeof(TYPE), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (float cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (double cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (long cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (unsigned long cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (int cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }

      scalar<TYPE> & operator= (unsigned int cpu_other)
      {
        //copy value:
        TYPE other = cpu_other;
        cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle(), handle(), CL_TRUE, 0, sizeof(TYPE), &other, 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        return *this;
      }
      /** @brief Sets the scalar to the result of supplied inner product expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_inner_prod> const & proxy)
      {
        viennacl::linalg::inner_prod_impl(proxy.get_lhs(), proxy.get_rhs(), *this);
        return *this;
      }

      /** @brief Sets the scalar to the result of supplied norm_1 expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_norm_1> const & proxy)
      {
        viennacl::linalg::norm_1_impl(proxy.get_lhs(), *this);
        return *this;
      }

      /** @brief Sets the scalar to the result of supplied norm_2 expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_norm_2> const & proxy)
      {
        viennacl::linalg::norm_2_impl(proxy.get_lhs(), *this);
        return *this;
      }

      /** @brief Sets the scalar to the result of supplied norm_inf expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_norm_inf> const & proxy)
      {
        viennacl::linalg::norm_inf_impl(proxy.get_lhs(), *this);
        return *this;
      }

      /** @brief Inplace addition of a ViennaCL scalar */
      scalar<TYPE> & operator += (scalar<TYPE> const & other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "inplace_add");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        
        viennacl::ocl::enqueue(k(val_, other.val_));
        return *this;
      }
      /** @brief Inplace addition of a host scalar (float or double) */
      scalar<TYPE> & operator += (TYPE other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_inplace_add");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_));        
        return *this;
      }


      /** @brief Inplace subtraction of a ViennaCL scalar */
      scalar<TYPE> & operator -= (scalar<TYPE> const & other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "inplace_sub");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        
        viennacl::ocl::enqueue(k(val_, other.val_));
        return *this;
      }
      /** @brief Inplace subtraction of a host scalar (float or double) */
      scalar<TYPE> & operator -= (TYPE other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_inplace_sub");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_));        
        return *this;
      }


      /** @brief Inplace multiplication with a ViennaCL scalar */
      scalar<TYPE> & operator *= (scalar<TYPE> const & other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "inplace_mul");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        
        viennacl::ocl::enqueue(k(val_, other.val_));
        return *this;
      }
      /** @brief Inplace  multiplication with a host scalar (float or double) */
      scalar<TYPE> & operator *= (TYPE other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_inplace_mul");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_));        
        return *this;
      }


      //////////////// operator /=    ////////////////////////////
      /** @brief Inplace division with a ViennaCL scalar */
      scalar<TYPE> & operator /= (scalar<TYPE> const & other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "inplace_div");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);
        
        viennacl::ocl::enqueue(k(val_, other.val_));
        return *this;
      }
      /** @brief Inplace division with a host scalar (float or double) */
      scalar<TYPE> & operator /= (TYPE other)
      {
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_inplace_div");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_));        
        return *this;
      }
      
      
      //////////////// operator + ////////////////////////////
      /** @brief Addition of two ViennaCL scalars */
      scalar<TYPE> operator + (scalar<TYPE> const & other)
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "add");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_, result));        
        return result;
      }
      /** @brief Addition of a ViennaCL scalar with a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator + (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "add");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, result, result));        
        return result;
      }
      /** @brief Addition of a ViennaCL scalar with a host scalar (float, double) */
      scalar<TYPE> operator + (TYPE other)
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_add");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other, result));        
        return result;
      }


      //////////////// operator - ////////////////////////////
      /** @brief Subtraction of two ViennaCL scalars */
      scalar<TYPE> operator - (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "sub");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_, result));        
        return result;
      }
      /** @brief Subtraction of a ViennaCL scalar from a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator - (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = *this;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "sub");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, result, result));        
        return result;
      }
      /** @brief Subtraction of a host scalar (float, double) from a ViennaCL scalar */
      scalar<TYPE> operator - (TYPE other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_sub");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other, result));        
        return result;
        
        return result;
      }

      //////////////// operator * ////////////////////////////
      /** @brief Multiplication of two ViennaCL scalars */
      scalar<TYPE> operator * (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "mul");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_, result));        
        return result;
      }
      /** @brief Multiplication of a ViennaCL scalar with a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator * (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "mul");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, result, result));        
        return result;
      }
      /** @brief Multiplication of a host scalar (float, double) with a ViennaCL scalar */
      scalar<TYPE> operator * (TYPE other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_mul");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other, result));        
        return result;
      }
      
      //////////////// operator /    ////////////////////////////
      /** @brief Division of two ViennaCL scalars */
      scalar<TYPE> operator / (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "divide");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other.val_, result));        
        return result;
      }
      /** @brief Division of a ViennaCL scalar by a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator / (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "divide");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, result, result));        
        return result;
      }
      /** @brief Division of a ViennaCL scalar by a host scalar (float, double)*/
      scalar<TYPE> operator / (TYPE other) const
      {
        scalar<TYPE> result;
        //get kernel:
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(viennacl::linalg::kernels::scalar<TYPE, 1>::program_name(), "cpu_div");
        k.local_work_size(0, 1);
        k.global_work_size(0, 1);

        viennacl::ocl::enqueue(k(val_, other, result));        
        return result;
      }

      /** @brief Returns the OpenCL handle */
      const viennacl::ocl::handle<cl_mem> & handle() const { return val_; }
      
    private:
      viennacl::ocl::handle<cl_mem> val_;
    };
    
    
    //stream operators:
    /** @brief Allows to directly print the value of a scalar to an output stream */
    template<class SCALARTYPE>
    std::ostream & operator<<(std::ostream & s, const scalar<SCALARTYPE> & val)
    {
      SCALARTYPE temp = val;
      s << temp;
      return s;
    }

    /** @brief Allows to directly read a value of a scalar from an input stream */
    template<class SCALARTYPE>
    std::istream & operator>>(std::istream & s, const scalar<SCALARTYPE> & val)
    {
      SCALARTYPE temp;
      s >> temp;
      val = temp;
      return s;
    }

} //namespace viennacl

#endif
