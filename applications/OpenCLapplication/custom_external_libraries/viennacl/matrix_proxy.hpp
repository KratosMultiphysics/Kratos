#ifndef VIENNACL_MATRIX_PROXY_HPP_
#define VIENNACL_MATRIX_PROXY_HPP_

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

/** @file matrix_proxy.hpp
    @brief Proxy classes for matrices.
*/

#include "viennacl/forwards.h"
#include "viennacl/range.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/matrix_operations.hpp"

namespace viennacl
{

  template <typename MatrixType>
  class matrix_range
  {
    public:
      typedef typename MatrixType::value_type     value_type;
      typedef typename viennacl::result_of::cpu_value_type<value_type>::type    cpu_value_type;
      typedef range::size_type                    size_type;
      typedef range::difference_type              difference_type;
      typedef value_type                          reference;
      typedef const value_type &                  const_reference;
      
      matrix_range(MatrixType & A, 
                   range const & row_range,
                   range const & col_range) : A_(&A), row_range_(row_range), col_range_(col_range) {}
                   
      size_type start1() const { return row_range_.start(); }
      size_type size1() const { return row_range_.size(); }

      size_type start2() const { return col_range_.start(); }
      size_type size2() const { return col_range_.size(); }
      
      ////////// operator= //////////////////////////
      
      matrix_range<MatrixType> & operator = (const matrix_range<MatrixType> & other) 
      {
        viennacl::linalg::assign(*this, other);
        return *this;
      }
      
      template <typename MatrixType2>
      matrix_range<MatrixType> & operator = (const MatrixType2 & other) 
      {
        viennacl::linalg::assign(*this, other);
        return *this;
      }

      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator = (const matrix_expression< MatrixType1,
                                                                      MatrixType2,
                                                                      op_prod > & proxy) 
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, 1.0, 0.0);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & 
      operator = (const matrix_expression< MatrixType1,
                                           MatrixType2,
                                           op_add > & proxy) 
      {
        viennacl::linalg::add(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & 
      operator = (const matrix_expression< MatrixType1,
                                           MatrixType2,
                                           op_sub > & proxy) 
      {
        viennacl::linalg::sub(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }


      ////////// operator+= //////////////////////////

      matrix_range<MatrixType> & operator += (MatrixType const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }
      
      matrix_range<MatrixType> & operator += (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }

      matrix_range<MatrixType> & operator += (matrix_slice<MatrixType> const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator += (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, 1.0, 1.0);
        return *this;
      }
      
      
      ////////// operator-= //////////////////////////

      matrix_range<MatrixType> & operator -= (MatrixType const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }
      
      matrix_range<MatrixType> & operator -= (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }

      matrix_range<MatrixType> & operator -= (matrix_slice<MatrixType> const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }

      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator -= (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, -1.0, 1.0);
        return *this;
      }


      ////////// operator*= //////////////////////////

      template <typename T>
      matrix_range<MatrixType> & operator *= (T const & val)
      {
        viennacl::linalg::inplace_mult(*this, val);
        return *this;
      }
      
      ////////// operator/= //////////////////////////

      template <typename T>
      matrix_range<MatrixType> & operator /= (T const & val)
      {
        viennacl::linalg::inplace_divide(*this, val);
        return *this;
      }

      matrix_range<MatrixType> & operator /= (cpu_value_type val)
      {
        viennacl::linalg::inplace_mult(*this, cpu_value_type(1.0) / val);
        return *this;
      }


      ////////// operator+ //////////////////////////
      
      template <typename MatrixType2>
      typename viennacl::enable_if< viennacl::is_matrix<MatrixType2>::value,
                                    matrix_expression< const matrix_range<MatrixType>,
                                                       const MatrixType2,
                                                       op_add > >::type
      operator + (const MatrixType2 & other) 
      {
        return matrix_expression< const matrix_range<MatrixType>,
                                  const MatrixType2,
                                  op_add > (*this, other);
      }
      
      ////////// operator- //////////////////////////
      
      template <typename MatrixType2>
      typename viennacl::enable_if< viennacl::is_matrix<MatrixType2>::value,
                                    matrix_expression< const matrix_range<MatrixType>,
                                                       const MatrixType2,
                                                       op_sub > >::type
      operator - (const MatrixType2 & other) 
      {
        return matrix_expression< const matrix_range<MatrixType>,
                                  const MatrixType2,
                                  op_sub > (*this, other);
      }
      
      
      

      //const_reference operator()(size_type i, size_type j) const { return A_(start1() + i, start2() + i); }
      //reference operator()(size_type i, size_type j) { return A_(start1() + i, start2() + i); }

      MatrixType & get() { return *A_; }
      const MatrixType & get() const { return *A_; }

    private:
      MatrixType * A_;
      range row_range_;
      range col_range_;
  };

  
  /** @brief Returns an expression template class representing a transposed matrix */
  template <typename MatrixType>
  matrix_expression< const matrix_range<MatrixType>,
                     const matrix_range<MatrixType>,
                     op_trans> trans(const matrix_range<MatrixType> & mat)
  {
    return matrix_expression< const matrix_range<MatrixType>,
                              const matrix_range<MatrixType>,
                              op_trans>(mat, mat);
  }
  
  
  
  
  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////
  
  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_range<matrix<SCALARTYPE, row_major, 1> > & gpu_matrix_range )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2()) );
    
     if ( gpu_matrix_range.start2() != 0 ||  gpu_matrix_range.size2() !=  gpu_matrix_range.get().size2())
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size2());
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
       {
         for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[j] = cpu_matrix(i,j);
         
         std::size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.get().internal_size2() + gpu_matrix_range.start2();
         std::size_t num_entries = gpu_matrix_range.size2();
         cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[i*gpu_matrix_range.get().internal_size2() + j] = cpu_matrix(i,j);
       
       std::size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.get().internal_size2();
       std::size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;
     }
  }
  
  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_range<matrix<SCALARTYPE, column_major, 1> > & gpu_matrix_range )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2()) );
    
     if ( gpu_matrix_range.start1() != 0 ||  gpu_matrix_range.size1() != gpu_matrix_range.get().size1())
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1());
       
       //copy each stride separately:
       for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
           entries[i] = cpu_matrix(i,j);
         
         std::size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.get().internal_size1() + gpu_matrix_range.start1();
         std::size_t num_entries = gpu_matrix_range.size1();
         cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[i + j*gpu_matrix_range.get().internal_size1()] = cpu_matrix(i,j);
       
       std::size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.get().internal_size1();
       std::size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;
     }
    
  }


  /////////////////////////////////////////////////////////////
  ///////////////////////// GPU to CPU ////////////////////////
  /////////////////////////////////////////////////////////////
  
  
  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_range<matrix<SCALARTYPE, row_major, 1> > const & gpu_matrix_range,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2()) );
    
     if ( gpu_matrix_range.start2() != 0 ||  gpu_matrix_range.size2() !=  gpu_matrix_range.get().size2())
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size2());
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
       {
         std::size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.get().internal_size2() + gpu_matrix_range.start2();
         std::size_t num_entries = gpu_matrix_range.size2();
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
        
        for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
          cpu_matrix(i,j) = entries[j];
         
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       std::size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.get().internal_size2();
       std::size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;

       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
           cpu_matrix(i,j) = entries[i*gpu_matrix_range.get().internal_size2() + j];
    }
    
  }
  
  
  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_range<matrix<SCALARTYPE, column_major, 1> > const & gpu_matrix_range,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2()) );
    
     if ( gpu_matrix_range.start1() != 0 ||  gpu_matrix_range.size1() !=  gpu_matrix_range.get().size1())
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1());
       
       //copy each stride separately:
       for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         std::size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.get().internal_size1() + gpu_matrix_range.start1();
         std::size_t num_entries = gpu_matrix_range.size1();
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
        
        for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
          cpu_matrix(i,j) = entries[i];
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       //copy each stride separately:
       std::size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.get().internal_size1();
       std::size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;
       
       for (std::size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (std::size_t j=0; j < gpu_matrix_range.size2(); ++j)
           cpu_matrix(i,j) = entries[i + j*gpu_matrix_range.get().internal_size1()];
     }
    
  }


  template<typename MatrixType>
  std::ostream & operator<<(std::ostream & s, matrix_range<MatrixType> const & proxy)
  {
    MatrixType temp = proxy;
    s << temp;
    return s;
  }

  template<typename MatrixType>
  std::ostream & operator<<(std::ostream & s, matrix_range<const MatrixType> const & proxy)
  {
    MatrixType temp = proxy;
    s << temp;
    return s;
  }


  //
  // Convenience function
  //
  template <typename MatrixType>
  matrix_range<MatrixType> project(MatrixType & A, viennacl::range const & r1, viennacl::range const & r2)
  {
    return matrix_range<MatrixType>(A, r1, r2);
  }






//
//
//
/////////////////////////////// Slice /////////////////////////////////////////////
//
//
//









  template <typename MatrixType>
  class matrix_slice
  {
    public:
      typedef typename MatrixType::value_type     value_type;
      typedef typename viennacl::result_of::cpu_value_type<value_type>::type    cpu_value_type;
      typedef slice::size_type                    size_type;
      typedef slice::difference_type              difference_type;
      typedef value_type                          reference;
      typedef const value_type &                  const_reference;
      
      matrix_slice(MatrixType & A, 
                   slice const & row_slice,
                   slice const & col_slice) : A_(&A), row_slice_(row_slice), col_slice_(col_slice) {}
                   
      size_type start1() const { return row_slice_.start(); }
      size_type stride1() const { return row_slice_.stride(); }
      size_type size1() const { return row_slice_.size(); }

      size_type start2() const { return col_slice_.start(); }
      size_type stride2() const { return col_slice_.stride(); }
      size_type size2() const { return col_slice_.size(); }
      
      ////////// operator= //////////////////////////

      matrix_slice<MatrixType> & operator = (const matrix_slice<MatrixType> & other) 
      {
        viennacl::linalg::assign(*this, other);
        return *this;
      }
      
      template <typename MatrixType2>
      matrix_slice<MatrixType> & operator = (const MatrixType2 & other) 
      {
        viennacl::linalg::assign(*this, other);
        return *this;
      }

      
      template <typename MatrixType1, typename MatrixType2>
      matrix_slice<MatrixType> & operator = (const matrix_expression< MatrixType1,
                                                                      MatrixType2,
                                                                      op_prod > & proxy) 
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, 1.0, 0.0);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_slice<MatrixType> & 
      operator = (const matrix_expression< MatrixType1,
                                           MatrixType2,
                                           op_add > & proxy) 
      {
        viennacl::linalg::add(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }

      template <typename MatrixType1, typename MatrixType2>
      matrix_slice<MatrixType> & 
      operator = (const matrix_expression< MatrixType1,
                                           MatrixType2,
                                           op_sub > & proxy) 
      {
        viennacl::linalg::sub(proxy.lhs(), proxy.rhs(), *this);
        return *this;
      }


      ////////// operator+= //////////////////////////

      matrix_slice<MatrixType> & operator += (MatrixType const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }
      
      matrix_slice<MatrixType> & operator += (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }

      matrix_slice<MatrixType> & operator += (matrix_slice<MatrixType> const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_slice<MatrixType> & operator += (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, 1.0, 1.0);
        return *this;
      }
      
      
      ////////// operator-= //////////////////////////
      matrix_slice<MatrixType> & operator -= (MatrixType const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }
      
      matrix_slice<MatrixType> & operator -= (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }

      matrix_slice<MatrixType> & operator -= (matrix_slice<MatrixType> const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_slice<MatrixType> & operator -= (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this, -1.0, 1.0);
        return *this;
      }


      ////////// operator*= //////////////////////////

      template <typename T>
      matrix_slice<MatrixType> & operator *= (T const & val)
      {
        viennacl::linalg::inplace_mult(*this, val);
        return *this;
      }
      
      ////////// operator/= //////////////////////////

      template <typename T>
      matrix_slice<MatrixType> & operator /= (T const & val)
      {
        viennacl::linalg::inplace_divide(*this, val);
        return *this;
      }

      matrix_slice<MatrixType> & operator /= (cpu_value_type val)
      {
        viennacl::linalg::inplace_mult(*this, cpu_value_type(1.0) / val);
        return *this;
      }


      ////////// operator+ //////////////////////////
      
      template <typename MatrixType2>
      typename viennacl::enable_if< viennacl::is_matrix<MatrixType2>::value,
                                    matrix_expression< const matrix_slice<MatrixType>,
                                                       const MatrixType2,
                                                       op_add > >::type
      operator + (const MatrixType2 & other) 
      {
        return matrix_expression< const matrix_slice<MatrixType>,
                                  const MatrixType2,
                                  op_add > (*this, other);
      }
      
      ////////// operator- //////////////////////////
      
      template <typename MatrixType2>
      typename viennacl::enable_if< viennacl::is_matrix<MatrixType2>::value,
                                    matrix_expression< const matrix_slice<MatrixType>,
                                                       const MatrixType2,
                                                       op_sub > >::type
      operator - (const MatrixType2 & other) 
      {
        return matrix_expression< const matrix_slice<MatrixType>,
                                  const MatrixType2,
                                  op_sub > (*this, other);
      }
      
      
      

      //const_reference operator()(size_type i, size_type j) const { return A_(start1() + i, start2() + i); }
      //reference operator()(size_type i, size_type j) { return A_(start1() + i, start2() + i); }

      MatrixType & get() { return *A_; }
      const MatrixType & get() const { return *A_; }

    private:
      MatrixType * A_;
      slice row_slice_;
      slice col_slice_;
  };

  
  /** @brief Returns an expression template class representing a transposed matrix */
  template <typename MatrixType>
  matrix_expression< const matrix_slice<MatrixType>,
                     const matrix_slice<MatrixType>,
                     op_trans> trans(const matrix_slice<MatrixType> & mat)
  {
    return matrix_expression< const matrix_slice<MatrixType>,
                              const matrix_slice<MatrixType>,
                              op_trans>(mat, mat);
  }
  
  
  
  
  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////
  
  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_slice<matrix<SCALARTYPE, row_major, 1> > & gpu_matrix_slice )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2()) );
    
     if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
     {
       std::size_t num_entries = gpu_matrix_slice.size2() * gpu_matrix_slice.stride2(); //no. of entries per stride
       
       std::vector<SCALARTYPE> entries(num_entries);
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_slice.size1(); ++i)
       {
         std::size_t start_offset = (gpu_matrix_slice.start1() + i * gpu_matrix_slice.stride1()) * gpu_matrix_slice.get().internal_size2() + gpu_matrix_slice.start2();
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
         VIENNACL_ERR_CHECK(err);
         
         for (std::size_t j=0; j < gpu_matrix_slice.size2(); ++j)
           entries[j * gpu_matrix_slice.stride2()] = cpu_matrix(i,j);
         
         err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                    gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                    sizeof(SCALARTYPE)*start_offset,
                                    sizeof(SCALARTYPE)*num_entries,
                                    &(entries[0]), 0, NULL, NULL);
         
         VIENNACL_ERR_CHECK(err);
       }
     }
  }
  
  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_slice<matrix<SCALARTYPE, column_major, 1> > & gpu_matrix_slice )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2()) );
    
    if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
    {
      std::size_t num_entries = gpu_matrix_slice.size1() * gpu_matrix_slice.stride1(); //no. of entries per stride
      
      std::vector<SCALARTYPE> entries(num_entries);
      
      //copy each column stride separately:
      for (std::size_t j=0; j < gpu_matrix_slice.size2(); ++j)
      {
        std::size_t start_offset = gpu_matrix_slice.start1() + (gpu_matrix_slice.start2() + j * gpu_matrix_slice.stride2()) * gpu_matrix_slice.get().internal_size1();
        
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                        gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                        sizeof(SCALARTYPE)*start_offset,
                                        sizeof(SCALARTYPE)*num_entries,
                                        &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        
        for (std::size_t i=0; i < gpu_matrix_slice.size1(); ++i)
          entries[i * gpu_matrix_slice.stride1()] = cpu_matrix(i,j);
        
        err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(),
                                   gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                   sizeof(SCALARTYPE)*start_offset,
                                   sizeof(SCALARTYPE)*num_entries,
                                   &(entries[0]), 0, NULL, NULL);
        
        VIENNACL_ERR_CHECK(err);
      }
    }
    
  }


  /////////////////////////////////////////////////////////////
  ///////////////////////// GPU to CPU ////////////////////////
  /////////////////////////////////////////////////////////////
  
  
  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_slice<matrix<SCALARTYPE, row_major, 1> > const & gpu_matrix_slice,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2()) );
    
     if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
     {
       std::size_t num_entries = gpu_matrix_slice.size2() * gpu_matrix_slice.stride2(); //no. of entries per stride
       
       std::vector<SCALARTYPE> entries(num_entries);
       
       //copy each stride separately:
       for (std::size_t i=0; i < gpu_matrix_slice.size1(); ++i)
       {
         std::size_t start_offset = (gpu_matrix_slice.start1() + i * gpu_matrix_slice.stride1()) * gpu_matrix_slice.get().internal_size2() + gpu_matrix_slice.start2();
         
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
         
         for (std::size_t j=0; j < gpu_matrix_slice.size2(); ++j)
           cpu_matrix(i,j) = entries[j * gpu_matrix_slice.stride2()];
         
        VIENNACL_ERR_CHECK(err);
       }
     }
    
  }
  
  
  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_slice<matrix<SCALARTYPE, column_major, 1> > const & gpu_matrix_slice,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2()) );
    
    if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
    {
      std::size_t num_entries = gpu_matrix_slice.size1() * gpu_matrix_slice.stride1(); //no. of entries per stride
      
      std::vector<SCALARTYPE> entries(num_entries);
      
      //copy each column stride separately:
      for (std::size_t j=0; j < gpu_matrix_slice.size2(); ++j)
      {
        std::size_t start_offset = gpu_matrix_slice.start1() + (gpu_matrix_slice.start2() + j * gpu_matrix_slice.stride2()) * gpu_matrix_slice.get().internal_size1();
        
        cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                        gpu_matrix_slice.get().handle().get(), CL_TRUE, 
                                        sizeof(SCALARTYPE)*start_offset,
                                        sizeof(SCALARTYPE)*num_entries,
                                        &(entries[0]), 0, NULL, NULL);
        
        for (std::size_t i=0; i < gpu_matrix_slice.size1(); ++i)
          cpu_matrix(i,j) = entries[i * gpu_matrix_slice.stride1()];
        
        VIENNACL_ERR_CHECK(err);
      }
    }
    
  }


  template<typename MatrixType>
  std::ostream & operator<<(std::ostream & s, matrix_slice<MatrixType> const & proxy)
  {
    MatrixType temp = proxy;
    s << temp;
    return s;
  }

  template<typename MatrixType>
  std::ostream & operator<<(std::ostream & s, matrix_slice<const MatrixType> const & proxy)
  {
    MatrixType temp = proxy;
    s << temp;
    return s;
  }


  //
  // Convenience function
  //
  template <typename MatrixType>
  matrix_slice<MatrixType> project(MatrixType & A, viennacl::slice const & r1, viennacl::slice const & r2)
  {
    return matrix_slice<MatrixType>(A, r1, r2);
  }



}

#endif