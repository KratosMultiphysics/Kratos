#ifndef AMGCL_SPMAT_HPP
#define AMGCL_SPMAT_HPP

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   spmat.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  A set of routines to work with sparse matrices in CRS format.
 */

#include <vector>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <cassert>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace amgcl {

/// A set of routines to work with sparse matrices in CRS format.
namespace sparse {

/// Type of matrix indices.
template <class spmat>
struct matrix_index {
    typedef typename spmat::index_type type;
};

/// Type of matrix values.
template <class spmat>
struct matrix_value {
    typedef typename spmat::value_type type;
};

/// Number of rows in a matrix.
template <class spmat>
typename matrix_index<spmat>::type matrix_rows(const spmat &A) {
    return A.rows;
}

/// Number of columns in a matrix.
template <class spmat>
typename matrix_index<spmat>::type matrix_cols(const spmat &A) {
    return A.cols;
}

/// Pointer to an array of outer indices.
template <class spmat>
const typename matrix_index<spmat>::type * matrix_outer_index(const spmat &A) {
    return &(A.row[0]);
}

/// Pointer to an array of outer indices.
template <class spmat>
typename matrix_index<spmat>::type * matrix_outer_index(spmat &A) {
    return &(A.row[0]);
}

/// Pointer to an array of inner indices.
template <class spmat>
const typename matrix_index<spmat>::type * matrix_inner_index(const spmat &A) {
    return &(A.col[0]);
}

/// Pointer to an array of inner indices.
template <class spmat>
typename matrix_index<spmat>::type * matrix_inner_index(spmat &A) {
    return &(A.col[0]);
}

/// Pointer to an array of values.
template <class spmat>
const typename matrix_value<spmat>::type * matrix_values(const spmat &A) {
    return &(A.val[0]);
}

/// Pointer to an array of values.
template <class spmat>
typename matrix_value<spmat>::type * matrix_values(spmat &A) {
    return &(A.val[0]);
}

/// Number of nonzero entries.
template <class spmat>
typename matrix_index<spmat>::type matrix_nonzeros(const spmat &A) {
    return matrix_outer_index(A)[matrix_rows(A)];
}

/// Sparse matrix in CRS format.
/**
 * \param value_t  Type of matrix values.
 * \param indexe_t Type of matrix indices.
 */
template <typename value_t = double, class index_t = long long>
struct matrix {
    typedef index_t index_type;
    typedef value_t value_type;

    /// Empty matrix construction.
    matrix() : rows(0), cols(0) {}

    /// Constructs matrix with a given number of rows, columns, and nonzero entries.
    matrix(index_t rows, index_t cols, index_t nnz = 0) :
        rows(rows), cols(cols), row(rows + 1), col(nnz), val(nnz)
    {}

    /// Resizes matrix.
    void resize(index_t rows, index_t cols, index_t nnz = 0) {
        matrix(rows, cols, nnz).swap(*this);
    }

    /// Copy constructor.
    matrix(const matrix &A) :
        rows(A.rows),
        cols(A.cols),
        row( A.row ),
        col( A.col ),
        val( A.val )
    { }

    /// Swap contents with other matrix.
    void swap(matrix &A) {
        if (this != &A) {
            std::swap(rows, A.rows);
            std::swap(cols, A.cols);

            std::swap(row, A.row);
            std::swap(col, A.col);
            std::swap(val, A.val);
        }
    }

    /// Copy constructor from a compatible type.
    template <class spmat>
    matrix(const spmat &A) :
        rows(matrix_rows(A)),
        cols(matrix_cols(A)),
        row(matrix_outer_index(A), matrix_outer_index(A) + rows + 1),
        col(matrix_inner_index(A), matrix_inner_index(A) + row[rows]),
        val(matrix_values(A), matrix_values(A) + row[rows])
    { }

    /// Allocates memory for a given number of nonzero entries.
    void reserve(index_t nnz) {
        col.resize(nnz);
        val.resize(nnz);
    }

    /// Deallocates any heap memory held by the matrix.
    void clear() {
        matrix().swap(*this);
    }

    index_t rows;               ///< Number of rows.
    index_t cols;               ///< Number of columns.

    std::vector<index_t> row;   ///< Outer indices.
    std::vector<index_t> col;   ///< Inner indices.
    std::vector<value_t> val;   ///< Values.
};

/// Proxy for sparse matrix in CRS format that does not own its data.
/**
 * \param value_t  Type of matrix values.
 * \param indexe_t Type of matrix indices.
 */
template <typename value_t, class index_t>
struct matrix_map {
    typedef index_t index_type;
    typedef value_t value_type;

    /// Stores matrix dimensions and pointers to data arrays.
    matrix_map(
            index_t rows, index_t cols,
            const index_t *row, const index_t *col, const value_t *val
            )
        : rows(rows), cols(cols), row(row), col(col), val(val)
    {}

    const index_t rows;     ///< Number of rows.
    const index_t cols;     ///< Number of columns.

    const index_t *row;     ///< Outer indices.
    const index_t *col;     ///< Inner indices.
    const value_t *val;     ///< Values.
};

/// Constructs proxy matrix class from a raw data.
template <typename value_t, class index_t>
matrix_map<value_t, index_t> map(
        index_t rows, index_t cols,
        const index_t *row, const index_t *col, const value_t *val
        )
{
    return matrix_map<value_t, index_t>(rows, cols, row, col, val);
}

/// Transpose of a sparse matrix.
template <class spmat>
matrix<
    typename matrix_value<spmat>::type,
    typename matrix_index<spmat>::type
    >
transpose(const spmat &A) {
    typedef typename matrix_index<spmat>::type index_t;
    typedef typename matrix_value<spmat>::type value_t;

    const index_t n   = matrix_rows(A);
    const index_t m   = matrix_cols(A);
    const index_t nnz = matrix_nonzeros(A);

    const index_t *Arow = matrix_outer_index(A);
    const index_t *Acol = matrix_inner_index(A);
    const value_t *Aval = matrix_values(A);

    matrix<value_t, index_t> T(m, n, nnz);

    std::fill(T.row.begin(), T.row.end(), static_cast<index_t>(0));

    for(index_t j = 0; j < nnz; ++j)
        ++( T.row[Acol[j] + 1] );

    std::partial_sum(T.row.begin(), T.row.end(), T.row.begin());

    for(index_t i = 0; i < n; i++) {
        for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
            index_t head = T.row[Acol[j]]++;

            T.col[head] = i;
            T.val[head] = Aval[j];
        }
    }

    std::rotate(T.row.begin(), T.row.end() - 1, T.row.end());
    T.row[0] = 0;

    return T;
}

/// Matrix-matrix product.
template <class spmat1, class spmat2>
matrix<
    typename matrix_value<spmat1>::type,
    typename matrix_index<spmat1>::type
    >
prod(const spmat1 &A, const spmat2 &B) {
    typedef typename matrix_index<spmat1>::type index_t;
    typedef typename matrix_value<spmat1>::type value_t;

    const index_t n   = matrix_rows(A);
    const index_t m   = matrix_cols(B);

    const index_t *Arow = matrix_outer_index(A);
    const index_t *Acol = matrix_inner_index(A);
    const value_t *Aval = matrix_values(A);

    const index_t *Brow = matrix_outer_index(B);
    const index_t *Bcol = matrix_inner_index(B);
    const value_t *Bval = matrix_values(B);

    matrix<value_t, index_t> C(n, m);

    std::fill(C.row.begin(), C.row.end(), static_cast<index_t>(0));

#pragma omp parallel
    {
        std::vector<index_t> marker(m, static_cast<index_t>(-1));

#ifdef _OPENMP
	int nt  = omp_get_num_threads();
	int tid = omp_get_thread_num();

	index_t chunk_size  = (n + nt - 1) / nt;
	index_t chunk_start = tid * chunk_size;
	index_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
	index_t chunk_start = 0;
	index_t chunk_end   = n;
#endif

        for(index_t ia = chunk_start; ia < chunk_end; ++ia) {
            for(index_t ja = Arow[ia], ea = Arow[ia + 1]; ja < ea; ++ja) {
                index_t ca = Acol[ja];
                for(index_t jb = Brow[ca], eb = Brow[ca + 1]; jb < eb; ++jb) {
                    index_t cb = Bcol[jb];

                    if (marker[cb] != ia) {
                        marker[cb] = ia;
                        ++( C.row[ia + 1] );
                    }
                }
            }
        }

        std::fill(marker.begin(), marker.end(), static_cast<index_t>(-1));

#pragma omp barrier
#pragma omp single
        {
            std::partial_sum(C.row.begin(), C.row.end(), C.row.begin());
            C.reserve(C.row.back());
        }

        for(index_t ia = chunk_start; ia < chunk_end; ++ia) {
            index_t row_beg = C.row[ia];
            index_t row_end = row_beg;

            for(index_t ja = Arow[ia], ea = Arow[ia + 1]; ja < ea; ++ja) {
                index_t ca = Acol[ja];
                value_t va = Aval[ja];

                for(index_t jb = Brow[ca], eb = Brow[ca + 1]; jb < eb; ++jb) {
                    index_t cb = Bcol[jb];
                    value_t vb = Bval[jb];

                    if (marker[cb] < row_beg) {
                        marker[cb] = row_end;
                        C.col[row_end] = cb;
                        C.val[row_end] = va * vb;
                        ++row_end;
                    } else {
                        C.val[marker[cb]] += va * vb;
                    }
                }
            }
        }
    }

    return C;
}

/// Gauss-Jordan elimination.
template <typename index_t, class value_t>
void gaussj(index_t n, value_t *a) {
    const static value_t one = static_cast<value_t>(1);
    const static value_t zero = static_cast<value_t>(0);

    std::vector<index_t> idxc(n);
    std::vector<index_t> idxr(n);
    std::vector<char>    ipiv(n, false);

    for(index_t i = 0; i < n; ++i) {
        index_t irow = 0, icol = 0;

        value_t big = zero;
        for(index_t j = 0; j < n; ++j) {
            if (ipiv[j]) continue;

            for(index_t k = 0; k < n; ++k) {
                if (!ipiv[k] && fabs(a[j * n + k]) > big) {
                    big  = fabs(a[j * n + k]);
                    irow = j;
                    icol = k;
                }
            }
        }

        ipiv[icol] = true;

        if (irow != icol)
            std::swap_ranges(
                    a + n * irow, a + n * (irow + 1),
                    a + n * icol
                    );

        idxr[i] = irow;
        idxc[i] = icol;

        if (a[icol * n + icol] == zero)
            throw std::logic_error("Singular matrix in gaussj");

        value_t pivinv = one / a[icol * n + icol];
        a[icol * n + icol] = one;

        for(value_t *v = a + icol * n, *e = a + (icol + 1) * n; v != e; ++v)
            *v *= pivinv;

        for(index_t k = 0; k < n; ++k) {
            if (k != icol) {
                value_t dum = a[k * n + icol];
                a[k * n + icol] = zero;
                for(value_t *v1 = a + n * k, *v2 = a + n * icol, *e = a + n * (k + 1); v1 != e; ++v1, ++v2)
                    *v1 -= *v2 * dum;
            }
        }
    }

    for(index_t i = n - 1; i >= 0; --i) {
        if (idxr[i] != idxc[i]) {
            for(index_t j = 0; j < n; ++j)
                std::swap(a[j * n + idxr[i]], a[j * n + idxc[i]]);
        }
    }
}

/// Inversion of a sparse matrix.
/**
 * Gauss-Jordan elimination routine is used.
 */
template <class spmat>
matrix<
    typename matrix_value<spmat>::type,
    typename matrix_index<spmat>::type
    >
inverse(const spmat &A) {
    typedef typename matrix_index<spmat>::type index_t;
    typedef typename matrix_value<spmat>::type value_t;

    const index_t n = sparse::matrix_rows(A);

    assert(n == sparse::matrix_cols(A)
            && "Inverse of a non-square matrix does not make sense");

    const index_t *Arow = matrix_outer_index(A);
    const index_t *Acol = matrix_inner_index(A);
    const value_t *Aval = matrix_values(A);

    matrix<
        typename matrix_value<spmat>::type,
        typename matrix_index<spmat>::type
    > Ainv(n, n, n * n);

    std::fill(Ainv.val.begin(), Ainv.val.end(), static_cast<value_t>(0));

    for(index_t i = 0; i < n; ++i)
        for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j)
            Ainv.val[i * n + Acol[j]] = Aval[j];

    gaussj(n, Ainv.val.data());

    Ainv.row[0] = 0;
    for(index_t i = 0, idx = 0; i < n; ) {
        for(index_t j = 0; j < n; ++j, ++idx) Ainv.col[idx] = j;

        Ainv.row[++i] = idx;
    }

    return Ainv;
}

/// Returns vector containing diagonal entrix of a sparse matrix.
template <class spmat>
std::vector< typename matrix_value<spmat>::type >
diagonal(const spmat &A) {
    typedef typename matrix_index<spmat>::type index_t;
    typedef typename matrix_value<spmat>::type value_t;

    const index_t n = sparse::matrix_rows(A);

    assert(n == sparse::matrix_cols(A)
            && "Diagonal of a non-square matrix is not well-defined");

    const index_t *Arow = matrix_outer_index(A);
    const index_t *Acol = matrix_inner_index(A);
    const value_t *Aval = matrix_values(A);

    std::vector<value_t> dia(n);

#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        value_t d = 0;
        for(index_t j = Arow[i], e = Arow[i + 1]; j < e; ++j) {
            if (Acol[j] == i) {
                d = Aval[j];
                break;
            }
        }
        dia[i] = d;
    }

    return dia;
}

//---------------------------------------------------------------------------
template <typename I, typename V>
void insertion_sort(I *col, V *val, I n) {
    for(I j = 1; j < n; ++j) {
        I c = col[j];
        V v = val[j];
        I i = j - 1;
        while(i >= 0 && col[i] > c) {
            col[i + 1] = col[i];
            val[i + 1] = val[i];
            i--;
        }
        col[i + 1] = c;
        val[i + 1] = v;
    }
}

/// Sort rows of the matrix column-wise.
template <class spmat>
void sort_rows(spmat &A) {
    typedef typename matrix_index<spmat>::type index_t;
    typedef typename matrix_value<spmat>::type value_t;

    const index_t n = sparse::matrix_rows(A);

    index_t *Arow = matrix_outer_index(A);
    index_t *Acol = matrix_inner_index(A);
    value_t *Aval = matrix_values(A);

#pragma omp parallel for schedule(dynamic, 1024)
    for(index_t i = 0; i < n; ++i) {
        index_t beg = Arow[i];
        index_t end = Arow[i + 1];
        insertion_sort(Acol + beg, Aval + beg, end - beg);
    }
}

} // namespace sparse
} // namespace amgcl

#endif
