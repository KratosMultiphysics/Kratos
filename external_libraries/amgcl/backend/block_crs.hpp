#ifndef AMGCL_BACKEND_BLOCK_CRS_HPP
#define AMGCL_BACKEND_BLOCK_CRS_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/backend/block_crs.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse matrix in block-CRS format.
 */

#include <algorithm>
#include <numeric>

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/skyline_lu.hpp>

namespace amgcl {
namespace backend {

/// Sparse matrix in Block CRS format.
/**
 * \param V Value type.
 * \param C Column number type.
 * \param P Index type.
 */
template < typename V, typename C, typename P >
struct bcrs {
    typedef V val_type;
    typedef C col_type;
    typedef P ptr_type;

    size_t block_size;
    size_t nrows, ncols;
    size_t brows, bcols;

    std::vector<ptr_type> ptr;
    std::vector<col_type> col;
    std::vector<val_type> val;

    /// Converts matrix in CRS format to Block CRS format.
    /**
     * \param A          Input matrix.
     * \param block_size Block size.
     *
     * \note Input matrix dimensions are *not* required to be divisible by
     * block_size.
     */
    template < class Matrix >
    bcrs(const Matrix &A, size_t block_size)
        : block_size(block_size), nrows( rows(A) ), ncols( cols(A) ),
          brows((nrows + block_size - 1) / block_size),
          bcols((ncols + block_size - 1) / block_size),
          ptr(brows + 1, 0)
    {
#pragma omp parallel
        {
            std::vector<ptrdiff_t> marker(bcols, -1);

            // Count number of nonzeros in block matrix.
#pragma omp for
            for(ptr_type ib = 0; ib < static_cast<ptr_type>(brows); ++ib) {
                ptr_type ia = ib * block_size;

                for(size_t k = 0; k < block_size && ia < static_cast<ptr_type>(nrows); ++k, ++ia) {
                    for(auto a = backend::row_begin(A, ia); a; ++a) {
                        col_type cb = a.col() / block_size;

                        if (marker[cb] != static_cast<col_type>(ib)) {
                            marker[cb]  = static_cast<col_type>(ib);
                            ++ptr[ib + 1];
                        }
                    }
                }
            }

#pragma omp single
            {
                std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());
                col.resize(ptr.back());
                val.resize(ptr.back() * block_size * block_size, 0);
            }

            std::fill(marker.begin(), marker.end(), -1);

            // Fill the block matrix.
#pragma omp for
            for(ptr_type ib = 0; ib < static_cast<ptr_type>(brows); ++ib) {
                ptr_type ia = ib * block_size;
                ptr_type row_beg = ptr[ib];
                ptr_type row_end = row_beg;

                for(size_t k = 0; k < block_size && ia < static_cast<ptr_type>(nrows); ++k, ++ia) {
                    for(auto a = backend::row_begin(A, ia); a; ++a) {
                        col_type cb = a.col() / block_size;
                        col_type cc = a.col() % block_size;
                        val_type va = a.value();

                        if (marker[cb] < row_beg) {
                            marker[cb] = row_end;
                            col[row_end] = cb;
                            val[block_size * (block_size * row_end + k) + cc] = va;
                            ++row_end;
                        } else {
                            val[block_size * (block_size * marker[cb] + k) + cc] = va;
                        }
                    }
                }
            }
        }
    }
};

/// block_crs backend definition.
/**
 * \param real Value type.
 * \ingroup backends
 */
template <typename real>
struct block_crs {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    typedef bcrs<real, index_type, index_type> matrix;
    typedef typename builtin<real>::vector     vector;
    typedef typename builtin<real>::vector     matrix_diagonal;
    typedef solver::skyline_lu<value_type>     direct_solver;

    struct provides_row_iterator : std::false_type {};

    /// Backend parameters.
    struct params {
        /// Block size to use with the created matrices.
        size_t block_size;

        params(size_t block_size = 4) : block_size(block_size) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : AMGCL_PARAMS_IMPORT_VALUE(p, block_size)
        {
            check_params(p, {"block_size"});
        }
        void get(boost::property_tree::ptree &p, const std::string &path) const {
            AMGCL_PARAMS_EXPORT_VALUE(p, path, block_size);
        }
#endif
    };

    static std::string name() { return "block_crs"; }

    /// Copy matrix from builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr< typename backend::builtin<real>::matrix > A,
            const params &prm)
    {
        return std::make_shared<matrix>(*A, prm.block_size);
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(const vector &x, const params&)
    {
        return std::make_shared<vector>(x);
    }

    static std::shared_ptr< vector >
    copy_vector(const std::vector<value_type> &x, const params&)
    {
        return std::make_shared<vector>(x);
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(std::shared_ptr< vector > x, const params&)
    {
        return x;
    }

    /// Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return std::make_shared<vector>(size);
    }

    static std::shared_ptr<direct_solver>
    create_solver(
            std::shared_ptr< typename backend::builtin<real>::matrix > A,
            const params&)
    {
        return std::make_shared<direct_solver>(*A);
    }
};

//---------------------------------------------------------------------------
// Specialization of backend interface
//---------------------------------------------------------------------------
template < typename V, typename C, typename P >
struct value_type< bcrs<V, C, P> > {
    typedef V type;
};

template < typename V, typename C, typename P >
struct rows_impl< bcrs<V, C, P> > {
    static size_t get(const bcrs<V, C, P> &A) {
        return A.nrows;
    }
};

template < typename V, typename C, typename P >
struct cols_impl< bcrs<V, C, P> > {
    static size_t get(const bcrs<V, C, P> &A) {
        return A.ncols;
    }
};

template < typename V, typename C, typename P >
struct nonzeros_impl< bcrs<V, C, P> > {
    static size_t get(const bcrs<V, C, P> &A) {
        return A.ptr.back() * A.block_size * A.block_size;
    }
};

template < typename Alpha, typename Beta, typename V, typename C, typename P, class Vec1, class Vec2 >
struct spmv_impl< Alpha, bcrs<V, C, P>, Vec1, Beta, Vec2 >
{
    typedef bcrs<V, C, P>  matrix;

    static void apply(Alpha alpha, const matrix &A, const Vec1 &x, Beta beta, Vec2 &y)
    {
        const size_t nb  = A.brows;
        const size_t na  = A.nrows;
        const size_t ma  = A.ncols;
        const size_t b1 = A.block_size;
        const size_t b2 = b1 * b1;

        if (!math::is_zero(beta)) {
            if (beta != 1) {
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(na); ++i) {
                    y[i] *= beta;
                }
            }
        } else {
            backend::clear(y);
        }

#pragma omp parallel for
        for(ptrdiff_t ib = 0; ib < static_cast<ptrdiff_t>(nb); ++ib) {
            for(P jb = A.ptr[ib], eb = A.ptr[ib + 1]; jb < eb; ++jb) {
                size_t x0 = A.col[jb] * b1;
                size_t y0 = ib * b1;
                block_prod(b1, std::min(b1, ma - x0), std::min(b1, na - y0),
                        alpha, &A.val[jb * b2], &x[x0], &y[y0]
                        );
            }
        }
    }

    static void block_prod(size_t dim, size_t nx, size_t ny,
            Alpha alpha, const V *A, const V *x, V *y)
    {
        for(size_t i = 0; i < ny; ++i, ++y) {
            const V * xx = x;
            V sum = 0;
            for(size_t j = 0; j < dim; ++j, ++A, ++xx)
                if (j < nx) sum += (*A) * (*xx);
            *y += alpha * sum;
        }
    }
};

template < typename V, typename C, typename P, class Vec1, class Vec2, class Vec3 >
struct residual_impl< bcrs<V, C, P>, Vec1, Vec2, Vec3 >
{
    typedef bcrs<V, C, P>  matrix;

    static void apply(const Vec1 &rhs, const matrix &A, const Vec2 &x, Vec3 &r)
    {
        backend::copy(rhs, r);
        backend::spmv(-1, A, x, 1, r);
    }
};

} // namespace backend
} // namespace amgcl

#endif
