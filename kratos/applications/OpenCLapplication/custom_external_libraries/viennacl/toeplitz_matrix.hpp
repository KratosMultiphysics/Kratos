#ifndef VIENNACL_TOEPLITZ_MATRIX_HPP
#define VIENNACL_TOEPLITZ_MATRIX_HPP

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

/** @file toeplitz_matrix.hpp
    @brief Implementation of the toeplitz_matrix class for efficient manipulation of Toeplitz matrices.  Experimental.
*/

#include "viennacl/forwards.h"
#include "viennacl/vector.hpp"
#include "viennacl/ocl/backend.hpp"

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
        typedef viennacl::backend::mem_handle                                                              handle_type;
        typedef scalar<typename viennacl::tools::CHECK_SCALAR_TEMPLATE_ARGUMENT<SCALARTYPE>::ResultType>   value_type;

        /**
         * @brief The default constructor. Does not allocate any memory.
         *
         */
        explicit toeplitz_matrix() {}

        /** @brief         Creates the matrix with the given size
        *
        * @param rows      Number of rows of the matrix
        * @param cols      Number of columns of the matrix
        */
        explicit toeplitz_matrix(vcl_size_t rows, vcl_size_t cols) : elements_(rows * 2)
        {
          assert(rows == cols && bool("Toeplitz matrix must be square!"));
          (void)cols;  // avoid 'unused parameter' warning in optimized builds
        }


        /** @brief Resizes the matrix.
        *   Existing entries can be preserved
        *
        * @param sz         New size of matrix
        * @param preserve   If true, existing values are preserved.
        */
        void resize(vcl_size_t sz, bool preserve = true) {
            elements_.resize(sz * 2, preserve);
        }

        /** @brief Returns the OpenCL handle
        *
        *   @return OpenCL handle
        */
        handle_type const & handle() const { return elements_.handle(); }

        /**
         * @brief Returns an internal viennacl::vector, which represents a Toeplitz matrix elements
         *
         */
        viennacl::vector<SCALARTYPE, ALIGNMENT> & elements() { return elements_; }
        viennacl::vector<SCALARTYPE, ALIGNMENT> const & elements() const { return elements_; }


        /**
         * @brief Returns the number of rows of the matrix
         */
        vcl_size_t size1() const { return elements_.size() / 2; }

        /**
         * @brief Returns the number of columns of the matrix
         */
        vcl_size_t size2() const { return elements_.size() / 2; }

        /** @brief Returns the internal size of matrix representtion.
        *   Usually required for launching OpenCL kernels only
        *
        *   @return Internal size of matrix representation
        */
        vcl_size_t internal_size() const { return elements_.internal_size(); }


        /**
         * @brief Read-write access to a single element of the matrix
         *
         * @param row_index  Row index of accessed element
         * @param col_index  Column index of accessed element
         * @return Proxy for matrix entry
         */
        entry_proxy<SCALARTYPE> operator()(vcl_size_t row_index, vcl_size_t col_index)
        {
            assert(row_index < size1() && col_index < size2() && bool("Invalid access"));

            long index = static_cast<long>(col_index) - static_cast<long>(row_index);

            if (index < 0)
              index = -index;
            else if
              (index > 0) index = 2 * static_cast<long>(size1()) - index;
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
        toeplitz_matrix(toeplitz_matrix const &) {}
        toeplitz_matrix & operator=(toeplitz_matrix const & t);


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
        assert( (gpu_mat.size1() == 0 || (gpu_mat.size1() * 2 - 1)  == cpu_vec.size()) && bool("Size mismatch"));

        vcl_size_t size = gpu_mat.size1();
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
        assert((gpu_mat.size1() * 2 - 1)  == cpu_vec.size() && bool("Size mismatch"));

        vcl_size_t size = gpu_mat.size1();
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
        assert(tep_src.size1() == viennacl::traits::size1(com_dst) && bool("Size mismatch"));
        assert(tep_src.size2() == viennacl::traits::size2(com_dst) && bool("Size mismatch"));

        vcl_size_t size = tep_src.size1();
        std::vector<SCALARTYPE> tmp(tep_src.size1() * 2 - 1);
        copy(tep_src, tmp);

        for(vcl_size_t i = 0; i < size; i++)
            for(vcl_size_t j = 0; j < size; j++)
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
        assert( (tep_dst.size1() == 0 || tep_dst.size1() == viennacl::traits::size1(com_src)) && bool("Size mismatch"));
        assert( (tep_dst.size2() == 0 || tep_dst.size2() == viennacl::traits::size2(com_src)) && bool("Size mismatch"));

        vcl_size_t size = tep_dst.size1();

        std::vector<SCALARTYPE> tmp(2*size - 1);

        for(long i = static_cast<long>(size) - 1; i >= 0; i--)
            tmp[size - i - 1] = com_src(i, 0);

        for(vcl_size_t i = 1; i < size; i++)
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
        vcl_size_t size = gpu_matrix.size1();
        std::vector<SCALARTYPE> tmp(2*size - 1);
        copy(gpu_matrix, tmp);
        s << "[" << size << "," << size << "](";

        for(vcl_size_t i = 0; i < size; i++) {
            s << "(";
            for(vcl_size_t j = 0; j < size; j++) {
                s << tmp[static_cast<int>(j) - static_cast<int>(i) + static_cast<int>(size - 1)];
                //s << (int)i - (int)j;
                if(j < (size - 1)) s << ",";
            }
            s << ")";
        }
        s << ")";
        return s;
    }

    //
    // Specify available operations:
    //

    /** \cond */

    namespace linalg
    {
      namespace detail
      {
        // x = A * y
        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              // check for the special case x = A * x
              if (viennacl::traits::handle(lhs) == viennacl::traits::handle(rhs.rhs()))
              {
                viennacl::vector<T> temp(lhs);
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
                lhs = temp;
              }
              else
                viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), lhs);
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs += temp;
            }
        };

        template <typename T, unsigned int A>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, const vector_base<T>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), rhs.rhs(), temp);
              lhs -= temp;
            }
        };


        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_assign, vector_expression<const toeplitz_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs());
              viennacl::linalg::prod_impl(rhs.lhs(), temp, lhs);
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_add, vector_expression<const toeplitz_matrix<T, A>, vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs());
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs += temp_result;
            }
        };

        // x = A * vec_op
        template <typename T, unsigned int A, typename LHS, typename RHS, typename OP>
        struct op_executor<vector_base<T>, op_inplace_sub, vector_expression<const toeplitz_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> >
        {
            static void apply(vector_base<T> & lhs, vector_expression<const toeplitz_matrix<T, A>, const vector_expression<const LHS, const RHS, OP>, op_prod> const & rhs)
            {
              viennacl::vector<T> temp(rhs.rhs());
              viennacl::vector<T> temp_result(lhs);
              viennacl::linalg::prod_impl(rhs.lhs(), temp, temp_result);
              lhs -= temp_result;
            }
        };

      } // namespace detail
    } // namespace linalg

    /** \endcond */

}

#endif // VIENNACL_TOEPLITZ_MATRIX_HPP
