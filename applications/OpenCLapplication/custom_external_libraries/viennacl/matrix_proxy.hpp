#ifndef VIENNACL_MATRIX_PROXY_HPP_
#define VIENNACL_MATRIX_PROXY_HPP_

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

/** @file matrix_proxy.hpp
    @brief Proxy classes for matrices.
*/

#include "viennacl/forwards.h"
#include "viennacl/range.hpp"
#include "viennacl/matrix.hpp"
#include "viennacl/linalg/matrix_operations.hpp"

namespace viennacl
{

  /** @brief Class for representing non-strided submatrices of a bigger matrix A.
    *
    * In MATLAB notation, this could for example refer to the submatrix A(3:8, 6:10) of a matrix A.
    */
  template <typename MatrixType>
  class matrix_range : public matrix_base<typename MatrixType::cpu_value_type, typename MatrixType::orientation_functor>
  {
      typedef matrix_base<typename MatrixType::cpu_value_type,
                          typename MatrixType::orientation_functor>    base_type;
      typedef matrix_range<MatrixType>                                 self_type;

    public:
      typedef typename MatrixType::orientation_category       orientation_category;

      typedef typename MatrixType::value_type     value_type;
      typedef typename viennacl::result_of::cpu_value_type<value_type>::type    cpu_value_type;
      typedef range::size_type                    size_type;
      typedef range::difference_type              difference_type;
      typedef value_type                          reference;
      typedef const value_type &                  const_reference;

      matrix_range(MatrixType & A,
                   range const & row_range,
                   range const & col_range) : base_type(A.handle(),
                                                        row_range.size(), row_range.start(), 1, A.internal_size1(),
                                                        col_range.size(), col_range.start(), 1, A.internal_size2()) {}

      using base_type::operator=;

  };


  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////

  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_range<matrix<SCALARTYPE, row_major, 1> > & gpu_matrix_range )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2())
           && bool("Matrix size mismatch!"));

    if ( gpu_matrix_range.start2() != 0)
    {
      std::vector<SCALARTYPE> entries(gpu_matrix_range.size2());

      //copy each stride separately:
      for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
      {
        for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
          entries[j] = cpu_matrix(i,j);

        vcl_size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.internal_size2() + gpu_matrix_range.start2();
        vcl_size_t num_entries = gpu_matrix_range.size2();
        viennacl::backend::memory_write(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
      //std::cout << "Strided copy worked!" << std::endl;
      }
    }
    else
    {
      //full block can be copied:
      std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.internal_size2());

      //copy each stride separately:
      for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
        for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
          entries[i*gpu_matrix_range.internal_size2() + j] = cpu_matrix(i,j);

      vcl_size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.internal_size2();
      vcl_size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.internal_size2();
      viennacl::backend::memory_write(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
      //std::cout << "Block copy worked!" << std::endl;
    }
  }

  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_range<matrix<SCALARTYPE, column_major, 1> > & gpu_matrix_range )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2())
           && bool("Matrix size mismatch!"));

     if ( gpu_matrix_range.start1() != 0 ||  gpu_matrix_range.size1() != gpu_matrix_range.size1())
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1());

       //copy each stride separately:
       for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
           entries[i] = cpu_matrix(i,j);

         vcl_size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.internal_size1() + gpu_matrix_range.start1();
         vcl_size_t num_entries = gpu_matrix_range.size1();
         viennacl::backend::memory_write(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
        //std::cout << "Strided copy worked!" << std::endl;
       }
     }
     else
     {
       //full block can be copied:
       std::vector<SCALARTYPE> entries(gpu_matrix_range.internal_size1()*gpu_matrix_range.size2());

       //copy each stride separately:
       for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
           entries[i + j*gpu_matrix_range.internal_size1()] = cpu_matrix(i,j);

       vcl_size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.internal_size1();
       vcl_size_t num_entries = gpu_matrix_range.internal_size1() * gpu_matrix_range.size2();
       viennacl::backend::memory_write(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
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
           && (cpu_matrix.size2() == gpu_matrix_range.size2())
           && bool("Matrix size mismatch!"));

     if ( gpu_matrix_range.start2() != 0)
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size2());

       //copy each stride separately:
       for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
       {
         vcl_size_t start_offset = (gpu_matrix_range.start1() + i) * gpu_matrix_range.internal_size2() + gpu_matrix_range.start2();
         vcl_size_t num_entries = gpu_matrix_range.size2();
         viennacl::backend::memory_read(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
        //std::cout << "Strided copy worked!" << std::endl;

        for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
          cpu_matrix(i,j) = entries[j];
       }
     }
     else
     {
       //full block can be copied:
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1()*gpu_matrix_range.internal_size2());

       vcl_size_t start_offset = gpu_matrix_range.start1() * gpu_matrix_range.internal_size2();
       vcl_size_t num_entries = gpu_matrix_range.size1() * gpu_matrix_range.size2();
         viennacl::backend::memory_read(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
       //std::cout << "Block copy worked!" << std::endl;

       for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
           cpu_matrix(i,j) = entries[i*gpu_matrix_range.internal_size2() + j];
    }

  }


  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_range<matrix<SCALARTYPE, column_major, 1> > const & gpu_matrix_range,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_range.size1())
           && (cpu_matrix.size2() == gpu_matrix_range.size2())
           && bool("Matrix size mismatch!"));

     if ( gpu_matrix_range.start1() != 0)
     {
       std::vector<SCALARTYPE> entries(gpu_matrix_range.size1());

       //copy each stride separately:
       for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
       {
         vcl_size_t start_offset = (gpu_matrix_range.start2() + j) * gpu_matrix_range.internal_size1() + gpu_matrix_range.start1();
         vcl_size_t num_entries = gpu_matrix_range.size1();
         viennacl::backend::memory_read(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
        //std::cout << "Strided copy worked!" << std::endl;

        for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
          cpu_matrix(i,j) = entries[i];
       }
     }
     else
     {
       //full block can be copied:
       std::vector<SCALARTYPE> entries(gpu_matrix_range.internal_size1()*gpu_matrix_range.size2());

       //copy each stride separately:
       vcl_size_t start_offset = gpu_matrix_range.start2() * gpu_matrix_range.internal_size1();
       vcl_size_t num_entries = gpu_matrix_range.internal_size1() * gpu_matrix_range.size2();
       viennacl::backend::memory_read(gpu_matrix_range.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
       //std::cout << "Block copy worked!" << std::endl;

       for (vcl_size_t i=0; i < gpu_matrix_range.size1(); ++i)
         for (vcl_size_t j=0; j < gpu_matrix_range.size2(); ++j)
           cpu_matrix(i,j) = entries[i + j*gpu_matrix_range.internal_size1()];
     }

  }


  //
  // Convenience function
  //
  template <typename MatrixType>
  matrix_range<MatrixType> project(MatrixType & A, viennacl::range const & r1, viennacl::range const & r2)
  {
    assert(r1.size() <= A.size1() && r2.size() <= A.size2() && bool("Size of range invalid!"));

    return matrix_range<MatrixType>(A, r1, r2);
  }


  template <typename MatrixType>
  matrix_range<MatrixType> project(matrix_range<MatrixType> & A, viennacl::range const & r1, viennacl::range const & r2)
  {
    assert(r1.size() <= A.size1() && r2.size() <= A.size2() && bool("Size of range invalid!"));

    return matrix_range<MatrixType>(A,
                                    viennacl::range(A.start1() + r1.start(), A.start1() + r1.start() + r1.size()),
                                    viennacl::range(A.start2() + r2.start(), A.start2() + r2.start() + r2.size())
                                   );
  }




//
//
//
/////////////////////////////// Slice /////////////////////////////////////////////
//
//
//





  /** @brief Class for representing strided submatrices of a bigger matrix A.
    *
    * In MATLAB notation, this could for example refer to the submatrix A(3:2:8, 6:3:16) of a matrix A.
    */
  template <typename MatrixType>
  class matrix_slice : public matrix_base<typename MatrixType::cpu_value_type, typename MatrixType::orientation_functor>
  {
      typedef matrix_base<typename MatrixType::cpu_value_type,
                          typename MatrixType::orientation_functor>    base_type;
      typedef matrix_slice<MatrixType>                                 self_type;

    public:
      typedef typename MatrixType::orientation_category       orientation_category;

      typedef typename MatrixType::value_type     value_type;
      typedef typename viennacl::result_of::cpu_value_type<value_type>::type    cpu_value_type;
      typedef range::size_type                    size_type;
      typedef range::difference_type              difference_type;
      typedef value_type                          reference;
      typedef const value_type &                  const_reference;

      matrix_slice(MatrixType & A,
                   slice const & row_slice,
                   slice const & col_slice) : base_type(A.handle(),
                                                        row_slice.size(), row_slice.start(), row_slice.stride(), A.internal_size1(),
                                                        col_slice.size(), col_slice.start(), col_slice.stride(), A.internal_size2()) {}

      using base_type::operator=;

  };



  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////

  //row_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_slice<matrix<SCALARTYPE, row_major, 1> > & gpu_matrix_slice )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2())
           && bool("Matrix size mismatch!"));

     if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
     {
       vcl_size_t num_entries = gpu_matrix_slice.size2() * gpu_matrix_slice.stride2(); //no. of entries per stride

       std::vector<SCALARTYPE> entries(num_entries);

       //copy each stride separately:
       for (vcl_size_t i=0; i < gpu_matrix_slice.size1(); ++i)
       {
         vcl_size_t start_offset = (gpu_matrix_slice.start1() + i * gpu_matrix_slice.stride1()) * gpu_matrix_slice.internal_size2() + gpu_matrix_slice.start2();
         viennacl::backend::memory_read(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));

         for (vcl_size_t j=0; j < gpu_matrix_slice.size2(); ++j)
           entries[j * gpu_matrix_slice.stride2()] = cpu_matrix(i,j);

         viennacl::backend::memory_write(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
       }
     }
  }

  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(const CPU_MATRIX & cpu_matrix,
            matrix_slice<matrix<SCALARTYPE, column_major, 1> > & gpu_matrix_slice )
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2())
           && bool("Matrix size mismatch!"));


    if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
    {
      vcl_size_t num_entries = gpu_matrix_slice.size1() * gpu_matrix_slice.stride1(); //no. of entries per stride

      std::vector<SCALARTYPE> entries(num_entries);

      //copy each column stride separately:
      for (vcl_size_t j=0; j < gpu_matrix_slice.size2(); ++j)
      {
        vcl_size_t start_offset = gpu_matrix_slice.start1() + (gpu_matrix_slice.start2() + j * gpu_matrix_slice.stride2()) * gpu_matrix_slice.internal_size1();

        viennacl::backend::memory_read(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));

        for (vcl_size_t i=0; i < gpu_matrix_slice.size1(); ++i)
          entries[i * gpu_matrix_slice.stride1()] = cpu_matrix(i,j);

        viennacl::backend::memory_write(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));
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
           && (cpu_matrix.size2() == gpu_matrix_slice.size2())
           && bool("Matrix size mismatch!"));

     if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
     {
       vcl_size_t num_entries = gpu_matrix_slice.size2() * gpu_matrix_slice.stride2(); //no. of entries per stride

       std::vector<SCALARTYPE> entries(num_entries);

       //copy each stride separately:
       for (vcl_size_t i=0; i < gpu_matrix_slice.size1(); ++i)
       {
         vcl_size_t start_offset = (gpu_matrix_slice.start1() + i * gpu_matrix_slice.stride1()) * gpu_matrix_slice.internal_size2() + gpu_matrix_slice.start2();

         viennacl::backend::memory_read(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));

         for (vcl_size_t j=0; j < gpu_matrix_slice.size2(); ++j)
           cpu_matrix(i,j) = entries[j * gpu_matrix_slice.stride2()];
       }
     }

  }


  //column_major:
  template <typename CPU_MATRIX, typename SCALARTYPE>
  void copy(matrix_slice<matrix<SCALARTYPE, column_major, 1> > const & gpu_matrix_slice,
            CPU_MATRIX & cpu_matrix)
  {
    assert( (cpu_matrix.size1() == gpu_matrix_slice.size1())
           && (cpu_matrix.size2() == gpu_matrix_slice.size2())
           && bool("Matrix size mismatch!"));

    if ( (gpu_matrix_slice.size1() > 0) && (gpu_matrix_slice.size1() > 0) )
    {
      vcl_size_t num_entries = gpu_matrix_slice.size1() * gpu_matrix_slice.stride1(); //no. of entries per stride

      std::vector<SCALARTYPE> entries(num_entries);

      //copy each column stride separately:
      for (vcl_size_t j=0; j < gpu_matrix_slice.size2(); ++j)
      {
        vcl_size_t start_offset = gpu_matrix_slice.start1() + (gpu_matrix_slice.start2() + j * gpu_matrix_slice.stride2()) * gpu_matrix_slice.internal_size1();

        viennacl::backend::memory_read(gpu_matrix_slice.handle(), sizeof(SCALARTYPE)*start_offset, sizeof(SCALARTYPE)*num_entries, &(entries[0]));

        for (vcl_size_t i=0; i < gpu_matrix_slice.size1(); ++i)
          cpu_matrix(i,j) = entries[i * gpu_matrix_slice.stride1()];
      }
    }

  }


  //
  // Convenience function
  //
  template <typename MatrixType>
  matrix_slice<MatrixType> project(MatrixType & A, viennacl::slice const & r1, viennacl::slice const & r2)
  {
    assert(r1.size() <= A.size1() && r2.size() <= A.size2() && bool("Size of slice invalid!"));

    return matrix_slice<MatrixType>(A, r1, r2);
  }

  template <typename MatrixType>
  matrix_slice<MatrixType> project(matrix_range<MatrixType> & A, viennacl::slice const & r1, viennacl::slice const & r2)
  {
    assert(r1.size() <= A.size1() && r2.size() <= A.size2() && bool("Size of slice invalid!"));

    return matrix_slice<MatrixType>(A,
                                    viennacl::slice(A.start1() + r1.start(), r1.stride(), r1.size()),
                                    viennacl::slice(A.start2() + r2.start(), r2.stride(), r2.size())
                                   );
  }

  template <typename MatrixType>
  matrix_slice<MatrixType> project(matrix_slice<MatrixType> & A, viennacl::slice const & r1, viennacl::slice const & r2)
  {
    assert(r1.size() <= A.size1() && r2.size() <= A.size2() && bool("Size of slice invalid!"));

    return matrix_slice<MatrixType>(A,
                                    viennacl::slice(A.start1() + r1.start(), A.stride1() * r1.stride(), r1.size()),
                                    viennacl::slice(A.start2() + r2.start(), A.stride2() * r2.stride(), r2.size())
                                   );
  }

  // TODO: Allow mix of range/slice

}

#endif
