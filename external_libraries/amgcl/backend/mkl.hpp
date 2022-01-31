#ifndef AMGCL_BACKEND_MKL_HPP
#define AMGCL_BACKEND_MKL_HPP

/*
The MIT License

Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/backend/mkl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Intel Math Kernel library backend.
 */

#include <vector>
#include <amgcl/backend/builtin.hpp>

#include <mkl.h>

namespace amgcl {
namespace backend {

/// Sparse matrix stored in CRS format.
struct mkl_crs : public crs<double, int> {
    typedef crs<double, int> Base;

    template <
        class PtrRange,
        class ColRange,
        class ValRange
        >
    mkl_crs(size_t nrows, size_t ncols,
        const PtrRange &ptr_range,
        const ColRange &col_range,
        const ValRange &val_range
        )
    : Base(nrows, ncols, ptr_range, col_range, val_range)
    {}

    template <class Matrix>
    mkl_crs(const Matrix &A) : Base(A) {}
};

/// BLAS vector.
class mkl_vec {
    public:
        mkl_vec() {}
        mkl_vec(size_t n) : buf(n) {}

        template <class Iterator>
        mkl_vec(Iterator beg, Iterator end) : buf(beg, end) {}

        size_t size() const {
            return buf.size();
        }

        double operator[](size_t i) const {
            return buf[i];
        }

        double& operator[](size_t i) {
            return buf[i];
        }

        const double* data() const {
            return buf.data();
        }

        double* data() {
            return buf.data();
        }
    private:
        std::vector<double> buf;
};

/// Intel Math Kernel library backend.
struct mkl {
    typedef double value_type;
    typedef int    index_type;

    typedef mkl_crs matrix;
    typedef mkl_vec vector;
    typedef vector  matrix_diagonal;
    typedef solver::skyline_lu<value_type> direct_solver;

    /// Backend parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
    };

    /// Copy matrix from builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(
            std::shared_ptr< typename builtin<value_type>::matrix > A,
            const params&
            )
    {
        return std::make_shared<matrix>(*A);
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(typename builtin<value_type>::vector const &x, const params&)
    {
        return std::make_shared<vector>(x.data(), x.data() + x.size());
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(
            std::shared_ptr< typename builtin<value_type>::vector > x,
            const params &prm
            )
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return std::make_shared<vector>(size);
    }

    /// Create direct solver for coarse level
    static std::shared_ptr<direct_solver>
    create_solver(std::shared_ptr< typename builtin<value_type>::matrix > A, const params&)
    {
        return std::make_shared<direct_solver>(*A);
    }
};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template <>
struct value_type < mkl_crs > {
    typedef double type;
};

template <>
struct value_type < mkl_vec > {
    typedef double type;
};

template <>
struct rows_impl< mkl_crs > {
    static size_t get(const mkl_crs &A) {
        return A.nrows;
    }
};

template <>
struct cols_impl< mkl_crs > {
    static size_t get(const mkl_crs &A) {
        return A.ncols;
    }
};

template <>
struct nonzeros_impl< mkl_crs > {
    static size_t get(const mkl_crs &A) {
        return A.ptr[A.nrows];
    }
};

template <typename Alpha, typename Beta>
struct spmv_impl<
  Alpha, mkl_crs, mkl_vec,
  Beta, mkl_vec
  >
{
    static void apply(Alpha alpha, const mkl_crs &A, const mkl_vec &x,
            Beta beta, mkl_vec &y)
    {
        MKL_INT m = A.nrows;
        MKL_INT k = A.ncols;
        mkl_dcsrmv("N", &m, &k, &alpha, "G__C",
                const_cast<double*>(&A.val[0]),
                const_cast<int*   >(&A.col[0]),
                const_cast<int*   >(&A.ptr[0]),
                const_cast<int*   >(&A.ptr[1]),
                const_cast<double*>(x.data()),
                &beta,
                y.data());
    }
};

template <>
struct residual_impl< mkl_crs, mkl_vec, mkl_vec, mkl_vec >
{
    static void apply(const mkl_vec &rhs, const mkl_crs &A, const mkl_vec &x,
            mkl_vec &r)
    {
        cblas_dcopy(rhs.size(), rhs.data(), 1, r.data(), 1);
        spmv_impl<double, mkl_crs, mkl_vec, double, mkl_vec>::apply(-1, A, x, 1, r);
    }
};

template <>
struct clear_impl< mkl_vec >
{
    static void apply(mkl_vec &x)
    {
        std::fill_n(x.data(), x.size(), 0.0);
    }
};

template <>
struct copy_impl< mkl_vec, mkl_vec >
{
    static void apply(const mkl_vec &x, mkl_vec &y)
    {
        cblas_dcopy(x.size(), x.data(), 1, y.data(), 1);
    }
};

template <>
struct inner_product_impl< mkl_vec, mkl_vec >
{
    static double get(const mkl_vec &x, const mkl_vec &y)
    {
        return cblas_ddot(x.size(), x.data(), 1, y.data(), 1);
    }
};

template <typename A, typename B>
struct axpby_impl<
  A, mkl_vec,
  B, mkl_vec
  >
{
    static void apply(A a, const mkl_vec &x, B b, mkl_vec &y)
    {
        cblas_dscal(y.size(), b, y.data(), 1);
        cblas_daxpy(y.size(), a, x.data(), 1, y.data(), 1);
    }
};

template <typename A, typename B, typename C>
struct axpbypcz_impl<
  A, mkl_vec,
  B, mkl_vec,
  C, mkl_vec
  >
{
    static void apply(
            A a, const mkl_vec &x,
            B b, const mkl_vec &y,
            C c,       mkl_vec &z
            )
    {
        cblas_dscal(z.size(), c, z.data(), 1);
        cblas_daxpy(z.size(), a, x.data(), 1, z.data(), 1);
        cblas_daxpy(z.size(), b, y.data(), 1, z.data(), 1);
    }
};

template <typename A, typename B>
struct vmul_impl<
  A, mkl_vec, mkl_vec,
  B, mkl_vec >
{
    static void apply(A a, const mkl_vec &x, const mkl_vec &y, B b, mkl_vec &z)
    {
        cblas_dsbmv(CblasRowMajor, CblasLower, z.size(), 0, a, x.data(), 1, y.data(), 1, b, z.data(), 1);
    }
};
} // namespace backend
} // namespace amgcl

#endif
