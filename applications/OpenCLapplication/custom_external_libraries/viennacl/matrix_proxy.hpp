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
      
      /** @brief Copy-constructor: Writes the entries from the matrix_range to the wrapped matrix.
       * 
       * Note: A generic overload of operator=() is insufficient, because then the compiler generates the copy-CTOR!
       * 
       * @param other    The submatrix to be assigned
       */
      matrix_range<MatrixType> & operator = (const matrix_range<MatrixType> & other) 
      {
        assert(size1() == other.size1());
        assert(size2() == other.size2());

        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< MatrixType >::ResultType    KernelClass;
        
        std::size_t block_size = 16;
        
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "assign");
        k.global_work_size(0, block_size*block_size);
        k.global_work_size(1, block_size*block_size);
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        
        viennacl::ocl::enqueue(k(viennacl::traits::handle(*A_),
                                        cl_uint(start1()),             cl_uint(start2()), 
                                        cl_uint(size1()),              cl_uint(size2()),
                                        cl_uint(A_->internal_size1()), cl_uint(A_->internal_size2()),
                                viennacl::traits::handle(other), 
                                        cl_uint(viennacl::traits::start1(other)),            cl_uint(viennacl::traits::start2(other)), 
                                        cl_uint(viennacl::traits::size1(other)),             cl_uint(viennacl::traits::size2(other)),
                                        cl_uint(viennacl::traits::internal_size1(other)),    cl_uint(viennacl::traits::internal_size2(other))
                                )
                              );

        return *this;
      }

      template <typename MatrixType2>
      matrix_range<MatrixType> & operator = (const MatrixType2 & other) 
      {
        assert(size1() == other.size1());
        assert(size2() == other.size2());

        typedef typename viennacl::tools::MATRIX_KERNEL_CLASS_DEDUCER< MatrixType >::ResultType    KernelClass;
        
        std::size_t block_size = 16;
        
        viennacl::ocl::kernel & k = viennacl::ocl::get_kernel(KernelClass::program_name(), "assign");
        k.global_work_size(0, block_size*block_size);
        k.global_work_size(1, block_size*block_size);
        k.local_work_size(0, block_size);
        k.local_work_size(1, block_size);
        
        viennacl::ocl::enqueue(k(viennacl::traits::handle(*A_),
                                        cl_uint(start1()),             cl_uint(start2()), 
                                        cl_uint(size1()),              cl_uint(size2()),
                                        cl_uint(A_->internal_size1()), cl_uint(A_->internal_size2()),
                                viennacl::traits::handle(other), 
                                        cl_uint(viennacl::traits::start1(other)),            cl_uint(viennacl::traits::start2(other)), 
                                        cl_uint(viennacl::traits::size1(other)),             cl_uint(viennacl::traits::size2(other)),
                                        cl_uint(viennacl::traits::internal_size1(other)),    cl_uint(viennacl::traits::internal_size2(other))
                                )
                              );

        return *this;
      }

      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator = (const matrix_expression< MatrixType1,
                                                                      MatrixType2,
                                                                      op_prod > & proxy) 
      {
        viennacl::linalg::prod_impl(proxy.lhs(), proxy.rhs(), *this);
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

      matrix_range<MatrixType> & operator += (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_add(*this, other);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator += (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        MatrixType temp = proxy;
        viennacl::linalg::inplace_add(*this, temp);
        return *this;
      }
      
      
      ////////// operator-= //////////////////////////
      matrix_range<MatrixType> & operator -= (matrix_range<MatrixType> const & other)
      {
        viennacl::linalg::inplace_sub(*this, other);
        return *this;
      }
      
      template <typename MatrixType1, typename MatrixType2>
      matrix_range<MatrixType> & operator -= (const matrix_expression< MatrixType1,
                                                                       MatrixType2,
                                                                       op_prod > & proxy)
      {
        MatrixType temp = proxy;
        viennacl::linalg::inplace_sub(*this, temp);
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
       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
       {
         for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[j] = cpu_matrix(i,j);
         
         size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.get().internal_size2() + gpu_matrix_range.start2();
         size_t num_entries = gpu_matrix_range.size2();
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
       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[i*gpu_matrix_range.get().internal_size2() + j] = cpu_matrix(i,j);
       
       size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.get().internal_size2();
       size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
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
       for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
           entries[i] = cpu_matrix(i,j);
         
         size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.get().internal_size1() + gpu_matrix_range.start1();
         size_t num_entries = gpu_matrix_range.size1();
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
       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[i + j*gpu_matrix_range.get().internal_size1()] = cpu_matrix(i,j);
       
       size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.get().internal_size1();
       size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
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
       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
       {
         size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.get().internal_size2() + gpu_matrix_range.start2();
         size_t num_entries = gpu_matrix_range.size2();
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
        
        for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
          cpu_matrix(i,j) = entries[j];
         
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.get().internal_size2();
       size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;

       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
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
       for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.get().internal_size1() + gpu_matrix_range.start1();
         size_t num_entries = gpu_matrix_range.size1();
         cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                          gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                          sizeof(SCALARTYPE)*start_offset,
                                          sizeof(SCALARTYPE)*num_entries,
                                          &(entries[0]), 0, NULL, NULL);
        VIENNACL_ERR_CHECK(err);
        //std::cout << "Strided copy worked!" << std::endl;
        
        for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
          cpu_matrix(i,j) = entries[i];
       }
     }
     else
     {
       //full block can be copied: 
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.size2());
       
       //copy each stride separately:
       size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.get().internal_size1();
       size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
       //std::cout << "start_offset: " << start_offset << std::endl;
       cl_int err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(),
                                         gpu_matrix_range.get().handle().get(), CL_TRUE, 
                                         sizeof(SCALARTYPE)*start_offset,
                                         sizeof(SCALARTYPE)*num_entries,
                                         &(entries[0]), 0, NULL, NULL);
       VIENNACL_ERR_CHECK(err);
       //std::cout << "Block copy worked!" << std::endl;
       
       for (size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (size_t j=0; j < gpu_matrix_range.size2(); ++j)
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

  /*template <typename MatrixType>
  matrix_range<MatrixType> project(MatrixType const & A, viennacl::range const & r1, viennacl::range const & r2)
  {
    return matrix_range<MatrixType>(A, r1, r2);
  }*/

  //TODO: Think about const-matrix...
  /*template <typename MatrixType>
  matrix_range<const MatrixType> project(MatrixType const & A, viennacl::range const & r1, viennacl::range const & r2)
  {
    return matrix_range<MatrixType>(A, r1, r2);
  }*/


}

#endif