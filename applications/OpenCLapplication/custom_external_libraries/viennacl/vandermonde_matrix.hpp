#ifndef VIENNACL_VANDERMONDE_MATRIX_HPP
#define VIENNACL_VANDERMONDE_MATRIX_HPP

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

#include <cmath>

/** @file vandermonde_matrix.hpp
    @brief Implementation of the vandermonde_matrix class for efficient manipulation of Vandermonde matrices.  Experimental in 1.2.x.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/context.hpp"

#include "viennacl/fft.hpp"

#include "viennacl/linalg/vandermonde_matrix_operations.hpp"

namespace viennacl {
    /** @brief A Vandermonde matrix class
    *
    * @tparam SCALARTYPE   The underlying scalar type (either float or double)
    * @tparam ALIGNMENT    The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    class vandermonde_matrix
    {
      public:
        /**
         * @brief The default constructor. Does not allocate any memory.
         *
         */
        explicit vandermonde_matrix()
        {
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /**
         * @brief         Creates the matrix with the given size
         *
         * @param rows      Number of rows of the matrix
         * @param cols      Number of columns of the matrix
         */
        explicit vandermonde_matrix(std::size_t rows, std::size_t cols) : elements_(rows)
        {
          assert(rows == cols && "Vandermonde matrix must be square in this release!");
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /** @brief Resizes the matrix.
        *   Existing entries can be preserved
        *
        * @param sz         New size of matrix
        * @param preserve   If true, existing values are preserved.
        */
        void resize(std::size_t sz, bool preserve = true) {
            elements_.resize(sz, preserve);
        }

        /** @brief Returns the OpenCL handle
        *
        *   @return OpenCL handle
        */
        viennacl::ocl::handle<cl_mem> handle() const { return elements_.handle(); }

        /**
         * @brief Returns an internal viennacl::vector, which represents a Vandermonde matrix elements
         *
         */
        viennacl::vector<SCALARTYPE, ALIGNMENT> & elements() { return elements_; }
        viennacl::vector<SCALARTYPE, ALIGNMENT> const & elements() const { return elements_; }

        /**
         * @brief Returns the number of rows of the matrix
         */
        std::size_t size1() const { return elements_.size(); }
        
        /**
         * @brief Returns the number of columns of the matrix
         */
        std::size_t size2() const { return elements_.size(); }

        /** @brief Returns the internal size of matrix representtion.
        *   Usually required for launching OpenCL kernels only
        *
        *   @return Internal size of matrix representation
        */
        std::size_t internal_size() const { return elements_.internal_size(); }

        /**
         * @brief Read-write access to a base element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @return Proxy for matrix entry
         */
        entry_proxy<SCALARTYPE> operator()(std::size_t row_index)
        {
            return elements_[row_index];
        }

        /**
         * @brief Read access to a element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @param col_index  Column index of accessed element
         * @return Proxy for matrix entry
         */
        SCALARTYPE operator()(std::size_t row_index, std::size_t col_index) const
        {
            assert(row_index < size1() && col_index < size2() && "Invalid access");
            
            return pow(elements_[row_index], static_cast<int>(col_index));
        }

    private:
        vandermonde_matrix(vandermonde_matrix const & t) {}
        vandermonde_matrix & operator=(vandermonde_matrix const & t) {}
        
        viennacl::vector<SCALARTYPE, ALIGNMENT> elements_;
    };

    /** @brief Copies a Vandermonde matrix from the std::vector to the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param cpu_vec   A std::vector on the host.
    * @param gpu_mat   A vandermonde_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(std::vector<SCALARTYPE>& cpu_vec, vandermonde_matrix<SCALARTYPE, ALIGNMENT>& gpu_mat)
    {
        assert(cpu_vec.size() == gpu_mat.size1()  && "Size mismatch");
        copy(cpu_vec, gpu_mat.elements());
    }

    /** @brief Copies a Vandermonde matrix from the OpenCL device (either GPU or multi-core CPU) to the std::vector
    *
    *
    * @param gpu_mat   A vandermonde_matrix from ViennaCL
    * @param cpu_vec   A std::vector on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(vandermonde_matrix<SCALARTYPE, ALIGNMENT>& gpu_mat, std::vector<SCALARTYPE>& cpu_vec)
    {
        assert(cpu_vec.size() == gpu_mat.size1() && "Size mismatch");
        copy(gpu_mat.elements(), cpu_vec);
    }

    /** @brief Copies a Vandermonde matrix from the OpenCL device (either GPU or multi-core CPU) to the matrix-like object
    *
    *
    * @param vander_src   A vandermonde_matrix from ViennaCL
    * @param com_dst   A matrix-like object
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(vandermonde_matrix<SCALARTYPE, ALIGNMENT>& vander_src, MATRIXTYPE& com_dst)
    {
        std::size_t size = vander_src.size1();
        assert(size == com_dst.size1() && "Size mismatch");
        assert(size == com_dst.size2() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(size);
        copy(vander_src, tmp);

        for(std::size_t i = 0; i < size; i++) {
            for(std::size_t j = 0; j < size; j++) {
                com_dst(i, j) = pow(tmp[i], static_cast<int>(j));
            }
        }
    }
    
    /** @brief Copies a the matrix-like object to the Vandermonde matrix from the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param com_src   A std::vector on the host
    * @param vander_dst   A vandermonde_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(MATRIXTYPE& com_src, vandermonde_matrix<SCALARTYPE, ALIGNMENT>& vander_dst) 
    {
        std::size_t size = vander_dst.size1();
        assert(size == com_src.size1() && "Size mismatch");
        assert(size == com_src.size2() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(size);

        for(std::size_t i = 0; i < size; i++)
            tmp[i] = com_src(i, 1);

        copy(tmp, vander_dst);
    }

    /*template <typename SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(vandermonde_matrix<SCALARTYPE, ALIGNMENT>& mat,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& vec,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& result) {
        assert(mat.size1() == vec.size());

        fft::vandermonde_prod<SCALARTYPE>(mat.handle(), vec.handle(), result.handle(), mat.size1());
    } */

    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param gpu_matrix   A ViennaCL Vandermonde matrix
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream& s, vandermonde_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix)
    {
        std::size_t size = gpu_matrix.size1();
        std::vector<SCALARTYPE> tmp(size);
        copy(gpu_matrix, tmp);
        s << "[" << size << "," << size << "](\n";

        for(std::size_t i = 0; i < size; i++) {
            s << "(";
            for(std::size_t j = 0; j < size; j++) {
                s << pow(tmp[i], j);
                if(j < (size - 1)) s << ",";
            }
            s << ")";
        }
        s << ")";
        return s;
    }

}

#endif // _VIENNACL_VANDERMONDE_MATRIX_HPP
