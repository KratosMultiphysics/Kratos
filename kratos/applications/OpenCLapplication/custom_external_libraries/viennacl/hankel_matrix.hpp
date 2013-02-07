#ifndef VIENNACL_HANKEL_MATRIX_HPP
#define VIENNACL_HANKEL_MATRIX_HPP

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


/** @file hankel_matrix.hpp
    @brief Implementation of the hankel_matrix class for efficient manipulation of Hankel matrices.  Experimental in 1.2.x.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/context.hpp"

#include "viennacl/toeplitz_matrix.hpp"
#include "viennacl/fft.hpp"

#include "viennacl/linalg/hankel_matrix_operations.hpp"

namespace viennacl {
    /** @brief A Hankel matrix class
    *
    * @tparam SCALARTYPE   The underlying scalar type (either float or double)
    * @tparam ALIGNMENT    The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    class hankel_matrix 
    {
    public:
        /**
         * @brief The default constructor. Does not allocate any memory.
         *
         */
        explicit hankel_matrix()
        {
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /**
         * @brief         Creates the matrix with the given size
         *
         * @param rows      Number of rows of the matrix
         * @param cols      Number of columns of the matrix
         */
        explicit hankel_matrix(std::size_t rows, std::size_t cols) : elements_(rows, cols)
        {
          assert(rows == cols && "Hankel matrix must be square!");
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /** @brief Resizes the matrix.
        *   Existing entries can be preserved
        *
        * @param sz         New size of matrix
        * @param preserve   If true, existing values are preserved.
        */
        void resize(size_t sz, bool preserve = true)
        {
            elements_.resize(sz, preserve);
        }

        /** @brief Returns the OpenCL handle
        *
        *   @return OpenCL handle
        */
        viennacl::ocl::handle<cl_mem> handle() const { return elements_.handle(); }

        /**
         * @brief Returns an internal viennacl::toeplitz_matrix, which represents a Hankel matrix elements
         *
         */
        toeplitz_matrix<SCALARTYPE, ALIGNMENT> & elements() { return elements_; }
        toeplitz_matrix<SCALARTYPE, ALIGNMENT> const & elements() const { return elements_; }

        /**
         * @brief Returns the number of rows of the matrix
         */
        std::size_t size1() const { return elements_.size1(); }
        
        /**
         * @brief Returns the number of columns of the matrix
         */
        std::size_t size2() const { return elements_.size2(); }

        /** @brief Returns the internal size of matrix representtion.
        *   Usually required for launching OpenCL kernels only
        *
        *   @return Internal size of matrix representation
        */
        std::size_t internal_size() const { return elements_.internal_size(); }

        /**
         * @brief Read-write access to a element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @param col_index  Column index of accessed element
         * @return Proxy for matrix entry
         */
        entry_proxy<SCALARTYPE> operator()(unsigned int row_index, unsigned int col_index)
        {
            assert(row_index < size1() && col_index < size2() && "Invalid access");
            
            return elements_(size1() - row_index - 1, col_index);
        }

        /**
         * @brief += operation for Hankel matrices
         *
         * @param that Matrix which will be added
         * @return Result of addition
         */
        hankel_matrix<SCALARTYPE, ALIGNMENT>& operator +=(hankel_matrix<SCALARTYPE, ALIGNMENT>& that)
        {
            elements_ += that.elements();
            return *this;
        }

    private:
        hankel_matrix(hankel_matrix const & t) {}
        hankel_matrix & operator=(hankel_matrix const & t) {}
      
        toeplitz_matrix<SCALARTYPE, ALIGNMENT> elements_;
    };

    /** @brief Copies a Hankel matrix from the std::vector to the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param cpu_vec   A std::vector on the host.
    * @param gpu_mat   A hankel_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(std::vector<SCALARTYPE> const & cpu_vec, hankel_matrix<SCALARTYPE, ALIGNMENT> & gpu_mat)
    {
        assert((gpu_mat.size1() * 2 - 1)  == cpu_vec.size() && "Size mismatch");

        copy(cpu_vec, gpu_mat.elements());
    }

    /** @brief Copies a Hankel matrix from the OpenCL device (either GPU or multi-core CPU) to the std::vector
    *
    *
    * @param gpu_mat   A hankel_matrix from ViennaCL
    * @param cpu_vec   A std::vector on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(hankel_matrix<SCALARTYPE, ALIGNMENT> const & gpu_mat, std::vector<SCALARTYPE> & cpu_vec)
    {
        assert((gpu_mat.size1() * 2 - 1)  == cpu_vec.size() && "Size mismatch");

        copy(gpu_mat.elements(), cpu_vec);
    }

    /** @brief Copies a Hankel matrix from the OpenCL device (either GPU or multi-core CPU) to the matrix-like object
    *
    *
    * @param han_src   A hankel_matrix from ViennaCL
    * @param com_dst   A matrix-like object
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(hankel_matrix<SCALARTYPE, ALIGNMENT> const & han_src, MATRIXTYPE& com_dst)
    {
        std::size_t size = han_src.size1();
        assert(size == com_dst.size1() && "Size mismatch");
        assert(size == com_dst.size2() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(size * 2 - 1);
        copy(han_src, tmp);

        for (std::size_t i = 0; i < size; i++)
            for (std::size_t j = 0; j < size; j++)
                com_dst(i, j) = tmp[i + j];
    }
    
    /** @brief Copies a the matrix-like object to the Hankel matrix from the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param com_src   A std::vector on the host
    * @param han_dst   A hankel_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(MATRIXTYPE const & com_src, hankel_matrix<SCALARTYPE, ALIGNMENT>& han_dst)
    {
        std::size_t size = han_dst.size1();
        assert(size == com_src.size1() && "Size mismatch");
        assert(size == com_src.size2() && "Size mismatch");

        std::vector<SCALARTYPE> tmp(2*size - 1);

        for (std::size_t i = 0; i < size; i++)
            tmp[i] = com_src(0, i);

        for (std::size_t i = 1; i < size; i++)
            tmp[size + i - 1] = com_src(size - 1, i);

        viennacl::copy(tmp, han_dst);
    }

    /*template <typename SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(hankel_matrix<SCALARTYPE, ALIGNMENT>& mat,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& vec,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& result)
    {
        prod_impl(mat.elements(), vec, result);
        fft::reverse(result);
    }*/

    template<class SCALARTYPE, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream & s, hankel_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix)
    {
        std::size_t size = gpu_matrix.size1();
        std::vector<SCALARTYPE> tmp(2*size - 1);
        copy(gpu_matrix, tmp);
        s << "[" << size << "," << size << "](";

        for(std::size_t i = 0; i < size; i++) {
            s << "(";
            for(std::size_t j = 0; j < size; j++) {
                s << tmp[i + j];
                //s << (int)i - (int)j;
                if(j < (size - 1)) s << ",";
            }
            s << ")";
        }
        s << ")";
        return s;
    }
}
#endif // _VIENNACL_HANKEL_MATRIX_HPP
