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

#ifndef _VIENNACL_SCALAR_HPP_
#define _VIENNACL_SCALAR_HPP_

#include "viennacl/forwards.h"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/kernel.hpp"
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
        _val = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(TYPE));
      }
      /** @brief Allocates the memory for the scalar and sets it to the supplied value. */
      scalar(TYPE val)
      {
        viennacl::linalg::kernels::scalar<TYPE, 1>::init(); 
        _val = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(TYPE), &val);
      }

      /** @brief Allocates memory for the scalar and sets it to the result of supplied expression. */
      template <typename T1, typename T2, typename OP>
      scalar(scalar_expression<T1, T2, OP> const & proxy)
      {
        viennacl::linalg::kernels::scalar<TYPE, 1>::init(); 
        _val = viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(TYPE));
        *this = proxy;
      }

      //copy constructor
      /** @brief Copy constructor. Allocates new memory for the scalar and copies the value of the supplied scalar */
      scalar(const scalar & other) : _val(viennacl::ocl::device().createMemory(CL_MEM_READ_WRITE, sizeof(TYPE)))
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), other.handle().get(), handle().get(), 0, 0, sizeof(TYPE), 0, NULL, NULL);
        CL_ERR_CHECK(err);
      }

      /** @brief Reads the value of the scalar from the GPU and returns the float or double value. */
      operator TYPE() const
      {
        TYPE tmp;
        cl_int err;
        err = clEnqueueReadBuffer(viennacl::ocl::device().queue().get(), _val.get(), CL_TRUE, 0, sizeof(TYPE), &tmp, 0, NULL, NULL);
        assert(err == CL_SUCCESS);
        viennacl::ocl::finish();
        return tmp;
      } 
      
      /** @brief Assigns a vector entry. */
      scalar<TYPE> & operator= (vector_entry_proxy<TYPE> const & other)
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), other.handle().get(), handle().get(), other.index() * sizeof(TYPE), 0, sizeof(TYPE), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        return *this;
      }

      /** @brief Assigns the value from another scalar. */
      scalar<TYPE> & operator= (scalar<TYPE> const & other)
      {
        //copy value:
        cl_int err = clEnqueueCopyBuffer(viennacl::ocl::device().queue().get(), other.handle().get(), handle().get(), 0, 0, sizeof(TYPE), 0, NULL, NULL);
        CL_ERR_CHECK(err);
        
        return *this;
      }

      /** @brief Sets the scalar to the result of supplied inner product expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_inner_prod> const & proxy)
      {
        viennacl::linalg::inner_prod_impl(proxy.get_lhs(), proxy.get_rhs(), *this);
        return *this;
      }

      /** @brief Sets the scalar to the result of supplied norm_2 expression. */
      template <typename T1, typename T2>
      scalar<TYPE> & operator= (scalar_expression<T1, T2, op_norm_2> const & proxy)
      {
        viennacl::linalg::norm_2_impl(proxy.get_lhs(), *this);
        return *this;
      }

      /** @brief Inplace addition of a ViennaCL scalar */
      scalar<TYPE> & operator += (scalar<TYPE> const & other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_add.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_add.setArgument(pos++, other._val);

        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_add.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }
      /** @brief Inplace addition of a host scalar (float or double) */
      scalar<TYPE> & operator += (TYPE other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_add.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_add.setArgument(pos++, other);

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_add.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }


      /** @brief Inplace subtraction of a ViennaCL scalar */
      scalar<TYPE> & operator -= (scalar<TYPE> const & other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_sub.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_sub.setArgument(pos++, other._val);

        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_sub.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }
      /** @brief Inplace subtraction of a host scalar (float or double) */
      scalar<TYPE> & operator -= (TYPE other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_sub.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_sub.setArgument(pos++, other);

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_sub.start1D(1,1);  //one thread is enough for such a trivial operation
        return *this;
      }


      /** @brief Inplace multiplication with a ViennaCL scalar */
      scalar<TYPE> & operator *= (scalar<TYPE> const & other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_mul.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_mul.setArgument(pos++, other._val);

        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_mul.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }
      /** @brief Inplace  multiplication with a host scalar (float or double) */
      scalar<TYPE> & operator *= (TYPE other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_mul.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_mul.setArgument(pos++, other);

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_mul.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }


      //////////////// operator /=    ////////////////////////////
      /** @brief Inplace division with a ViennaCL scalar */
      scalar<TYPE> & operator /= (scalar<TYPE> const & other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_div.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_div.setArgument(pos++, other._val);

        viennacl::linalg::kernels::scalar<TYPE, 1>::inplace_div.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }
      /** @brief Inplace division with a host scalar (float or double) */
      scalar<TYPE> & operator /= (TYPE other)
      {
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_div.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_div.setArgument(pos++, other);

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_inplace_div.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return *this;
      }
      
      
      //////////////// operator + ////////////////////////////
      /** @brief Addition of two ViennaCL scalars */
      scalar<TYPE> operator + (scalar<TYPE> const & other)
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, other._val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::add.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      /** @brief Addition of a ViennaCL scalar with a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator + (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::add.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::add.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      /** @brief Addition of a ViennaCL scalar with a host scalar (float, double) */
      scalar<TYPE> operator + (TYPE other)
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_add.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_add.setArgument(pos++, other);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_add.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_add.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }


      //////////////// operator - ////////////////////////////
      /** @brief Subtraction of two ViennaCL scalars */
      scalar<TYPE> operator - (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, other._val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }
      /** @brief Subtraction of a ViennaCL scalar from a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator - (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = *this;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      /** @brief Subtraction of a host scalar (float, double) from a ViennaCL scalar */
      scalar<TYPE> operator - (TYPE other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::sub.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }

      //////////////// operator * ////////////////////////////
      /** @brief Multiplication of two ViennaCL scalars */
      scalar<TYPE> operator * (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, other._val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      /** @brief Multiplication of a ViennaCL scalar with a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator * (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::mul.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      /** @brief Multiplication of a host scalar (float, double) with a ViennaCL scalar */
      scalar<TYPE> operator * (TYPE other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_mul.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_mul.setArgument(pos++, other);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_mul.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_mul.start1D(1,1);  //one thread is enough for such a trivial operation
        return result;
      }
      
      //////////////// operator /    ////////////////////////////
      /** @brief Division of two ViennaCL scalars */
      scalar<TYPE> operator / (scalar<TYPE> const & other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, other._val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }
      /** @brief Division of a ViennaCL scalar by a scalar expression */
      template <typename T1, typename T2, typename OP>
      scalar<TYPE> operator / (scalar_expression<T1, T2, OP> const & proxy) const
      {
        scalar<TYPE> result = proxy;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, result.handle());
        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::divide.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }
      /** @brief Division of a ViennaCL scalar by a host scalar (float, double)*/
      scalar<TYPE> operator / (TYPE other) const
      {
        scalar<TYPE> result;
        unsigned int pos = 0;
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_div.setArgument(pos++, _val);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_div.setArgument(pos++, other);
        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_div.setArgument(pos++, result.handle());

        viennacl::linalg::kernels::scalar<TYPE, 1>::cpu_div.start1D(1,1);  //one thread is enough for such a trivial operation
        
        return result;
      }

      /** @brief Returns the OpenCL handle */
      const viennacl::ocl::handle<cl_mem> & handle() const { return _val; }

    private:
      viennacl::ocl::handle<cl_mem> _val;
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
    std::istream & operator>>(std::ostream & s, const scalar<SCALARTYPE> & val)
    {
      SCALARTYPE temp;
      s >> temp;
      val = temp;
      return s;
    }

} //namespace viennacl

#endif
