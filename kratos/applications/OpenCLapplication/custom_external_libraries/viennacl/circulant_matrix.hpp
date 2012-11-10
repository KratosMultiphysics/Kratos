#ifndef _VIENNACL_CIRCULANT_MATRIX_HPP
#define _VIENNACL_CIRCULANT_MATRIX_HPP

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

/** @file circulant_matrix.hpp
    @brief Implementation of the circulant_matrix class for efficient manipulation of circulant matrices.  Experimental in 1.2.x.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/context.hpp"

#include "viennacl/linalg/circulant_matrix_operations.hpp"

#include "viennacl/fft.hpp"

namespace viennacl 
{
    /** @brief A Circulant matrix class
    *
    * @tparam SCALARTYPE  The underlying scalar type (either float or double)
    * @tparam ALIGNMENT   The internal memory size is given by (size()/ALIGNMENT + 1) * ALIGNMENT. ALIGNMENT must be a power of two. Best values or usually 4, 8 or 16, higher values are usually a waste of memory.
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    class circulant_matrix
    {
      public:     
        /**
         * @brief The default constructor. Does not allocate any memory.
         *
         */
        explicit circulant_matrix()
        {
          viennacl::linalg::kernels::fft<SCALARTYPE, 1>::init();
        }

        /**
         * @brief         Creates the matrix with the given size
         *
         * @param rows      Number of rows of the matrix
         * @param cols      Number of columns of the matrix
         */
        explicit circulant_matrix(std::size_t rows, std::size_t cols) : elements_(rows)
        {
          assert(rows == cols && "Circulant matrix must be square!");
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
         * @brief Returns an internal viennacl::vector, which represents a circulant matrix elements
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
         * @brief Read-write access to a single element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @param col_index  Column index of accessed element
         * @return Proxy for matrix entry
         */
        entry_proxy<SCALARTYPE> operator()(std::size_t row_index, std::size_t col_index)
        {
            int index = static_cast<int>(row_index) - static_cast<int>(col_index);

            assert(row_index < size1() && col_index < size2() && "Invalid access");
            
            while (index < 0)
              index += size1();
            return elements_[index];
        }

        /**
         * @brief += operation for circulant matrices
         *
         * @param that Matrix which will be added
         * @return Result of addition
         */
        circulant_matrix<SCALARTYPE, ALIGNMENT>& operator +=(circulant_matrix<SCALARTYPE, ALIGNMENT>& that)
        {
            elements_ += that.elements();
            return *this;
        }

    private:
        circulant_matrix(circulant_matrix const & t) {}
        circulant_matrix & operator=(circulant_matrix const & t) {}
      
        viennacl::vector<SCALARTYPE, ALIGNMENT> elements_;
    };

    /** @brief Copies a circulant matrix from the std::vector to the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param cpu_vec   A std::vector on the host.
    * @param gpu_mat   A circulant_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(std::vector<SCALARTYPE>& cpu_vec, circulant_matrix<SCALARTYPE, ALIGNMENT>& gpu_mat)
    {
        assert(cpu_vec.size() == gpu_mat.size1() && "Size mismatch");
        copy(cpu_vec, gpu_mat.elements());
    }

    /** @brief Copies a circulant matrix from the OpenCL device (either GPU or multi-core CPU) to the std::vector
    *
    *
    * @param gpu_mat   A circulant_matrix from ViennaCL
    * @param cpu_vec   A std::vector on the host.
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT>
    void copy(circulant_matrix<SCALARTYPE, ALIGNMENT>& gpu_mat, std::vector<SCALARTYPE>& cpu_vec)
    {
        assert(cpu_vec.size() == gpu_mat.size1() && "Size mismatch");
        copy(gpu_mat.elements(), cpu_vec);
    }

    /** @brief Copies a circulant matrix from the OpenCL device (either GPU or multi-core CPU) to the matrix-like object
    *
    *
    * @param circ_src   A circulant_matrix from ViennaCL
    * @param com_dst   A matrix-like object
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(circulant_matrix<SCALARTYPE, ALIGNMENT>& circ_src, MATRIXTYPE& com_dst) {
        std::size_t size = circ_src.size1();
        assert(size == com_dst.size1() && "Size mismatch");
        assert(size == com_dst.size2() && "Size mismatch");
        std::vector<SCALARTYPE> tmp(size);
        copy(circ_src, tmp);

        for (std::size_t i = 0; i < size; i++) {
            for (std::size_t j = 0; j < size; j++) {
                int index = static_cast<int>(i) - static_cast<int>(j);
                if (index < 0)
                  index = size + index;
                com_dst(i, j) = tmp[index];
            }
        }
    }

    /** @brief Copies a the matrix-like object to the circulant matrix from the OpenCL device (either GPU or multi-core CPU)
    *
    *
    * @param com_src   A std::vector on the host
    * @param circ_dst   A circulant_matrix from ViennaCL
    */
    template <typename SCALARTYPE, unsigned int ALIGNMENT, typename MATRIXTYPE>
    void copy(MATRIXTYPE& com_src, circulant_matrix<SCALARTYPE, ALIGNMENT>& circ_dst) {
        std::size_t size = circ_dst.size1();
        assert(size == com_src.size1() && "Size mismatch");
        assert(size == com_src.size2() && "Size mismatch");

        std::vector<SCALARTYPE> tmp(size);

        for(std::size_t i = 0; i < size; i++) tmp[i] = com_src(i, 0);

        copy(tmp, circ_dst);
    }

    /*namespace linalg
    {
      template <typename SCALARTYPE, unsigned int ALIGNMENT, unsigned int VECTOR_ALIGNMENT>
      void prod_impl(circulant_matrix<SCALARTYPE, ALIGNMENT> const & mat,
                      vector<SCALARTYPE, VECTOR_ALIGNMENT> const & vec,
                      vector<SCALARTYPE, VECTOR_ALIGNMENT>& result) {
          viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> circ(mat.elements().size() * 2);
          fft::real_to_complex(mat.elements(), circ, mat.elements().size());

          viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp(vec.size() * 2);
          viennacl::vector<SCALARTYPE, VECTOR_ALIGNMENT> tmp2(vec.size() * 2);

          fft::real_to_complex(vec, tmp, vec.size());
          fft::convolve(circ, tmp, tmp2);
          fft::complex_to_real(tmp2, result, vec.size());
      }
    }*/

    /** @brief Prints the matrix. Output is compatible to boost::numeric::ublas
    *
    * @param s            STL output stream
    * @param gpu_matrix   A ViennaCL circulant matrix
    */
    template<class SCALARTYPE, unsigned int ALIGNMENT>
    std::ostream & operator<<(std::ostream& s, circulant_matrix<SCALARTYPE, ALIGNMENT>& gpu_matrix)
    {
        std::size_t size = gpu_matrix.size1();
        std::vector<SCALARTYPE> tmp(size);
        copy(gpu_matrix, tmp);
        s << "[" << size << "," << size << "](";

        for(std::size_t i = 0; i < size; i++) {
            s << "(";
            for(std::size_t j = 0; j < size; j++) {
                int index = (int)i - (int)j;
                if(index < 0) index = size + index;
                s << tmp[index];
                //s << index;
                if(j < (size - 1)) s << ",";
            }
            s << ")";
        }
        s << ")";
        return s;
    }
}

#endif // _VIENNACL_CIRCULANT_MATRIX_HPP
