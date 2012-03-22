#ifndef VIENNACL_TOEPLITZ_MATRIX_HPP
#define VIENNACL_TOEPLITZ_MATRIX_HPP

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

/** @file toeplitz_matrix.hpp
    @brief Implementation of the toeplitz_matrix class for efficient manipulation of Toeplitz matrices.  Experimental in 1.2.x.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/context.hpp"

#include "viennacl/fft.hpp"

#include "viennacl/linalg/toeplitz_matrix_operations.hpp"


namespace viennacl {

    /** @brief A Toeplitz matrix class
    *
    * @tparam SCALARTYPE   The underlying scalar type (either float or double)
    * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    class toeplitz_matrix
    {
      public:

        /**
         * @brief The default constructor. Does not allocate any memory.
         *
         */
        explicit toeplitz_matrix()
        {
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /** @brief         Creates the matrix with the given size
        *
        * @param rows      Number of rows of the matrix
        * @param cols      Number of columns of the matrix
        */
        explicit toeplitz_matrix(std::size_t rows, std::size_t cols) : elements_(rows * 2)
        {
          assert(rows == cols && "Toeplitz matrix must be square!");
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }
        

        /** @brief Resizes the matrix.
        *   Existing entries can be preserved
        *
        * @param sz         New size of matrix
        * @param preserve   If true, existing values are preserved.
        */
        void resize(size_t sz, bool preserve = true) {
            elements_.resize(sz * 2, preserve);
        }

        /** @brief Returns the OpenCL handle
        *
        *   @return OpenCL handle
        */
        viennacl::ocl::handle<cl_mem> handle() const { return elements_.handle(); }

        /**
         * @brief Returns an internal viennacl::vector, which represents a Toeplitz matrix elements
         *
         */
        viennacl::vector<SCALARTYPE, ALIGNMENT> & elements() { return elements_; }
        viennacl::vector<SCALARTYPE, ALIGNMENT> const & elements() const { return elements_; }


        /**
         * @brief Returns the number of rows of the matrix
         */
        std::size_t size1() const { return elements_.size() / 2; }
        
        /**
         * @brief Returns the number of columns of the matrix
         */
        std::size_t size2() const { return elements_.size() / 2; }

        /** @brief Returns the internal size of matrix representtion.
        *   Usually required for launching OpenCL kernels only
        *
        *   @return Internal size of matrix representation
        */
        std::size_t internal_size() const { return elements_.internal_size(); }


        /**
         * @brief Read-write access to a single element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @param col_index  Column index of accessed element
         * @return Proxy for matrix entry
         */
        entry_proxy<SCALARTYPE> operator()(std::size_t row_index, std::size_t col_index) 
        {
            assert(row_index < size1() && col_index < size2() && "Invalid access");
            
            int index = static_cast<int>(col_index) - static_cast<int>(row_index);
            
            if (index < 0)
              index = -index;
            else if
              (index > 0) index = 2 * size1() - index;
            return elements_[index];
        }


        /**
         * @brief += operation for Toeplitz matrices
         *
         * @param that Matrix which will be added
         * @return Result of addition
         */
        toeplitz_matrix<SCALARTYPE, ALIGNMENT>& operator +=(toeplitz_matrix<SCALARTYPE, ALIGNMENT>& that) {
            elements_ += that.elements();
            return *this;
        }

    private:
        toeplitz_matrix(toeplitz_matrix const & t) {}
        toeplitz_matrix & operator=(toeplitz_matrix const & t) {}
        
      
        viennacl::vector<SCALARTYPE, ALIGNMENT> elements_;
    };

    /** @brief Copies a Toeplitz matrix from the std::vector to the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param cpu_vec   A std::vector on the host.
    * @param gpu_mat   A toeplitz_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(std::vector<SCALARTYPE> const & cpu_vec, toeplitz_matrix<SCALARTYPE, ALIGNMENT>& gpu_mat)
    {
        std::size_t size = gpu_mat.size1();
        assert((size * 2 - 1)  == cpu_vec.size() && "Size mismatch");
        std::vector<SCALARTYPE> rvrs(cpu_vec.size());
        std::copy(cpu_vec.begin(), cpu_vec.end(), rvrs.begin());
        std::reverse(rvrs.begin(), rvrs.end());

        std::vector<SCALARTYPE> tmp(size * 2);
        std::copy(rvrs.begin() + size - 1, rvrs.end(), tmp.begin());
        std::copy(rvrs.begin(), rvrs.begin() + size - 1, tmp.begin() + size + 1);
        tmp[size] = 0.0;
        copy(tmp, gpu_mat.elements());
    }

    /** @brief Copies a Toeplitz matrix from the OpenCL device (either GPU or multi-core CPU) to the std::vector
    *
    *
    * @param gpu_mat   A toeplitz_matrix from ViennaCL
    * @param cpu_vec   A std::vector on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(toeplitz_matrix<SCALARTYPE, ALIGNMENT> const & gpu_mat, std::vector<SCALARTYPE> & cpu_vec)
    {
        std::size_t size = gpu_mat.size1();
        assert((size * 2 - 1)  == cpu_vec.size() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(size * 2);
        copy(gpu_mat.elements(), tmp);
        std::reverse(tmp.begin(), tmp.end());

        std::copy(tmp.begin(), tmp.begin() + size - 1, cpu_vec.begin() + size);
        std::copy(tmp.begin() + size, tmp.end(), cpu_vec.begin());

    }

    /** @brief Copies a Toeplitz matrix from the OpenCL device (either GPU or multi-core CPU) to the matrix-like object
    *
    *
    * @param tep_src   A toeplitz_matrix from ViennaCL
    * @param com_dst   A matrix-like object
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(toeplitz_matrix<SCALARTYPE, ALIGNMENT> const & tep_src, MATRIXTYPE & com_dst)
    {
        std::size_t size = tep_src.size1();
        assert(size == com_dst.size1() && "Size mismatch");
        assert(size == com_dst.size2() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(tep_src.size1() * 2 - 1);
        copy(tep_src, tmp);

        for(std::size_t i = 0; i < size; i++)
            for(std::size_t j = 0; j < size; j++)
                com_dst(i, j) = tmp[static_cast<int>(j) - static_cast<int>(i) + static_cast<int>(size) - 1];
    }

    /** @brief Copies a the matrix-like object to the Toeplitz matrix from the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param com_src   A std::vector on the host
    * @param tep_dst   A toeplitz_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(MATRIXTYPE const & com_src, toeplitz_matrix<SCALARTYPE, ALIGNMENT>& tep_dst)
    {
        std::size_t size = tep_dst.size1();
        assert(size == com_src.size1() && "Size mismatch");
        assert(size == com_src.size2() && "Size mismatch");

        std::vector<SCALARTYPE> tmp(2*size - 1);

        for(int i = size - 1; i >= 0; i--)
            tmp[size - i - 1] = com_src(i, 0);

        for(std::size_t i = 1; i < size; i++)
            tmp[size + i - 1] = com_src(0, i);

        copy(tmp, tep_dst);
    }

    /*template <typename SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
    void prod_impl(toeplitz_matrix<SCALARTYPE, ALIGNMENT>& mat,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& vec,
                   vector<SCALARTYPE, VECTOR_ALIGNMENT>& result) {
        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tep(mat.elements().size() * 2);
        fft::real_to_complex(mat.elements(), tep, mat.elements().size());

        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp(vec.size() * 4);
        viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp2(vec.size() * 4);

        tmp.clear();
        copy(vec, tmp);
        fft::real_to_complex(tmp, tmp2, vec.size() * 2);
        fft::convolve(tep, tmp2, tmp);
        fft::complex_to_real(tmp, tmp2, vec.size() * 2);
        copy(tmp2.begin(), tmp2.begin() + vec.size(), result.begin());
    }*/

    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param gpu_matrix   A ViennaCL Toeplitz matrix
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream & s, toeplitz_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix)
    {
        std::size_t size = gpu_matrix.size1();
        std::vector<SCALARTYPE> tmp(2*size - 1);
        copy(gpu_matrix, tmp);
        s << "[" << size << "," << size << "](";

        for(std::size_t i = 0; i < size; i++) {
            s << "(";
            for(std::size_t j = 0; j < size; j++) {
                s << tmp[(int)j - (int)i + (int)size - 1];
                //s << (int)i - (int)j;
                if(j < (size - 1)) s << ",";
            }
            s << ")";
        }
        s << ")";
        return s;
    }

}

#endif // _VIENNACL_TOEPLITZ_MATRIX_HPP
