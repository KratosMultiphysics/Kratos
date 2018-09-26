#ifndef AMGCL_BACKEND_BUILTIN_HPP
#define AMGCL_BACKEND_BUILTIN_HPP

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
 * \file   amgcl/backend/builtin.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Builtin backend.
 */

#include <vector>
#include <numeric>
#include <memory>
#include <random>
#include <type_traits>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/detail/inverse.hpp>
#include <amgcl/detail/sort_row.hpp>
#include <amgcl/detail/spgemm.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>

namespace amgcl {
namespace backend {

/// Sparse matrix stored in CRS format.
template <
    typename val_t = double,
    typename col_t = ptrdiff_t,
    typename ptr_t = col_t
    >
struct crs {
    typedef val_t val_type;
    typedef col_t col_type;
    typedef ptr_t ptr_type;

    size_t nrows, ncols, nnz;
    ptr_type * ptr;
    col_type * col;
    val_type * val;
    bool own_data;

    crs() : nrows(0), ncols(0), nnz(0), ptr(0), col(0), val(0), own_data(true)
    {}

    template <
        class PtrRange,
        class ColRange,
        class ValRange
        >
    crs(size_t nrows, size_t ncols,
        const PtrRange &ptr_range,
        const ColRange &col_range,
        const ValRange &val_range
        ) :
        nrows(nrows), ncols(ncols), nnz(0),
        ptr(0), col(0), val(0), own_data(true)
    {
        precondition(static_cast<ptrdiff_t>(nrows + 1) == std::distance(
                    std::begin(ptr_range), std::end(ptr_range)),
                "ptr_range has wrong size in crs constructor");

        nnz = ptr_range[nrows];

        precondition(static_cast<ptrdiff_t>(nnz) == std::distance(
                    std::begin(col_range), std::end(col_range)),
                "col_range has wrong size in crs constructor");

        precondition(static_cast<ptrdiff_t>(nnz) == std::distance(
                    std::begin(val_range), std::end(val_range)),
                "val_range has wrong size in crs constructor");

        ptr = new ptr_type[nrows + 1];
        col = new col_type[nnz];
        val = new val_type[nnz];

        ptr[0] = ptr_range[0];
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            ptr[i+1] = ptr_range[i+1];
            for(ptr_type j = ptr_range[i]; j < ptr_range[i+1]; ++j) {
                col[j] = col_range[j];
                val[j] = val_range[j];
            }
        }
    }

    template <class Matrix>
    crs(const Matrix &A) :
        nrows(backend::rows(A)), ncols(backend::cols(A)),
        nnz(0), ptr(0), col(0), val(0), own_data(true)
    {
        ptr = new ptr_type[nrows + 1];
        ptr[0] = 0;

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            int row_width = 0;
            for(auto a = backend::row_begin(A, i); a; ++a) ++row_width;
            ptr[i+1] = row_width;
        }

        nnz = scan_row_sizes();
        col = new col_type[nnz];
        val = new val_type[nnz];

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            ptr_type row_head = ptr[i];
            for(auto a = backend::row_begin(A, i); a; ++a) {
                col[row_head] = a.col();
                val[row_head] = a.value();

                ++row_head;
            }
        }
    }

    crs(const crs &other) :
        nrows(other.nrows), ncols(other.ncols), nnz(other.nnz),
        ptr(0), col(0), val(0), own_data(true)
    {
        if (other.ptr && other.col && other.val) {
            ptr = new ptr_type[nrows + 1];
            col = new col_type[nnz];
            val = new val_type[nnz];

            ptr[0] = other.ptr[0];
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
                ptr[i+1] = other.ptr[i+1];
                for(ptr_type j = other.ptr[i]; j < other.ptr[i+1]; ++j) {
                    col[j] = other.col[j];
                    val[j] = other.val[j];
                }
            }
        }
    }

    crs(crs &&other) :
        nrows(other.nrows), ncols(other.ncols), nnz(other.nnz),
        ptr(other.ptr), col(other.col), val(other.val),
        own_data(other.own_data)
    {
        other.nrows = 0;
        other.ncols = 0;
        other.nnz   = 0;
        other.ptr   = 0;
        other.col   = 0;
        other.val   = 0;
    }

    const crs& operator=(const crs &other) {
        free_data();

        nrows = other.nrows;
        ncols = other.ncols;
        nnz   = other.nnz;

        if (other.ptr && other.col && other.val) {
            ptr = new ptr_type[nrows + 1];
            col = new col_type[nnz];
            val = new val_type[nnz];

            ptr[0] = other.ptr[0];
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
                ptr[i+1] = other.ptr[i+1];
                for(ptr_type j = other.ptr[i]; j < other.ptr[i+1]; ++j) {
                    col[j] = other.col[j];
                    val[j] = other.val[j];
                }
            }
        }

        return *this;
    }

    const crs& operator=(crs &&other) {
        std::swap(nrows,    other.nrows);
        std::swap(ncols,    other.ncols);
        std::swap(nnz,      other.nnz);
        std::swap(ptr,      other.ptr);
        std::swap(col,      other.col);
        std::swap(val,      other.val);
        std::swap(own_data, other.own_data);

        return *this;
    }

    void free_data() {
        if (own_data) {
            delete[] ptr; ptr = 0;
            delete[] col; col = 0;
            delete[] val; val = 0;
        }
    }

    void set_size(size_t n, size_t m, bool clean_ptr = false) {
        precondition(!ptr, "matrix data has already been allocated!");

        nrows = n;
        ncols = m;

        ptr = new ptr_type[nrows + 1];

        if (clean_ptr) {
            ptr[0] = 0;
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i)
                ptr[i+1] = 0;
        }
    }

    ptr_type scan_row_sizes() {
        std::partial_sum(ptr, ptr + nrows + 1, ptr);
        return ptr[nrows];
    }

    void set_nonzeros() {
        set_nonzeros(ptr[nrows]);

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(nrows); ++i) {
            ptrdiff_t row_beg = ptr[i];
            ptrdiff_t row_end = ptr[i+1];
            for(ptrdiff_t j = row_beg; j < row_end; ++j) {
                col[j] = 0;
                val[j] = math::zero<val_type>();
            }
        }
    }

    void set_nonzeros(size_t n, bool need_values = true) {
        precondition(!col && !val, "matrix data has already been allocated!");

        nnz = n;

        col = new col_type[nnz];

        if (need_values)
            val = new val_type[nnz];
    }

    ~crs() {
        free_data();
    }

    class row_iterator {
        public:
            row_iterator(
                    const col_type * col,
                    const col_type * end,
                    const val_type * val
                    ) : m_col(col), m_end(end), m_val(val)
            {}

            operator bool() const {
                return m_col < m_end;
            }

            row_iterator& operator++() {
                ++m_col;
                ++m_val;
                return *this;
            }

            col_type col() const {
                return *m_col;
            }

            val_type value() const {
                return *m_val;
            }

        private:
            const col_type * m_col;
            const col_type * m_end;
            const val_type * m_val;
    };

    row_iterator row_begin(size_t row) const {
        ptr_type p = ptr[row];
        ptr_type e = ptr[row + 1];
        return row_iterator(col + p, col + e, val + p);
    }

    size_t bytes() const {
        if (own_data) {
            return sizeof(ptr_type) * (nrows + 1)
                 + sizeof(col_type) * nnz
                 + sizeof(val_type) * nnz;
        } else {
            return 0;
        }
    }
};

/// Sort rows of the matrix column-wise.
template < typename V, typename C, typename P >
void sort_rows(crs<V, C, P> &A) {
    const size_t n = rows(A);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        P beg = A.ptr[i];
        P end = A.ptr[i + 1];
        amgcl::detail::sort_row(A.col + beg, A.val + beg, end - beg);
    }
}

/// Transpose of a sparse matrix.
template < typename V, typename C, typename P >
std::shared_ptr< crs<V,C,P> > transpose(const crs<V, C, P> &A)
{
    const size_t n   = rows(A);
    const size_t m   = cols(A);
    const size_t nnz = nonzeros(A);

    auto T = std::make_shared< crs<V,C,P> >();
    T->set_size(m, n, true);

    for(size_t j = 0; j < nnz; ++j)
        ++( T->ptr[A.col[j] + 1] );

    T->scan_row_sizes();
    T->set_nonzeros();

    for(size_t i = 0; i < n; i++) {
        for(P j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
            P head = T->ptr[A.col[j]]++;

            T->col[head] = static_cast<C>(i);
            T->val[head] = A.val[j];
        }
    }

    std::rotate(T->ptr, T->ptr + m, T->ptr + m + 1);
    T->ptr[0] = 0;

    return T;
}

/// Matrix-matrix product.
template <class Val, class Col, class Ptr>
std::shared_ptr< crs<Val, Col, Ptr> >
product(const crs<Val,Col,Ptr> &A, const crs<Val,Col,Ptr> &B, bool sort = false) {
    auto C = std::make_shared< crs<Val,Col,Ptr> >();

#ifdef _OPENMP
    int nt = omp_get_max_threads();
#else
    int nt = 1;
#endif

    if (nt > 16) {
        spgemm_rmerge(A, B, *C);
    } else {
        spgemm_saad(A, B, *C, sort);
    }

    return C;
}


/// Scale matrix values.
template<class Val, class Col, class Ptr, class T>
void scale(crs<Val, Col, Ptr> &A, T s) {
    ptrdiff_t n = backend::rows(A);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < n; ++i) {
        for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j)
            A.val[j] *= s;
    }
}

// Reduce matrix to a pointwise one
template <class value_type>
std::shared_ptr< crs<typename math::scalar_of<value_type>::type> >
pointwise_matrix(const crs<value_type> &A, unsigned block_size) {
    typedef value_type V;
    typedef typename math::scalar_of<V>::type S;

    AMGCL_TIC("pointwise_matrix");
    const ptrdiff_t n  = A.nrows;
    const ptrdiff_t m  = A.ncols;
    const ptrdiff_t np = n / block_size;
    const ptrdiff_t mp = m / block_size;

    precondition(np * block_size == n,
            "Matrix size should be divisible by block_size");

    auto ap = std::make_shared< crs<S> >();
    crs<S> &Ap = *ap;

    Ap.set_size(np, mp, true);

#pragma omp parallel
    {
        std::vector<ptrdiff_t> j(block_size);
        std::vector<ptrdiff_t> e(block_size);

        // Count number of nonzeros in block matrix.
#pragma omp for
        for(ptrdiff_t ip = 0; ip < np; ++ip) {
            ptrdiff_t ia = ip * block_size;
            ptrdiff_t cur_col = 0;
            bool done = true;

            for(unsigned k = 0; k < block_size; ++k) {
                ptrdiff_t beg = j[k] = A.ptr[ia + k];
                ptrdiff_t end = e[k] = A.ptr[ia + k + 1];

                if (beg == end) continue;

                ptrdiff_t c = A.col[beg];

                if (done) {
                    done = false;
                    cur_col = c;
                } else {
                    cur_col = std::min(cur_col, c);
                }
            }

            while(!done) {
                cur_col /= block_size;
                ++Ap.ptr[ip + 1];

                done = true;
                ptrdiff_t col_end = (cur_col + 1) * block_size;
                for(unsigned k = 0; k < block_size; ++k) {
                    ptrdiff_t beg = j[k];
                    ptrdiff_t end = e[k];

                    while(beg < end) {
                        ptrdiff_t c = A.col[beg++];

                        if (c >= col_end) {
                            if (done) {
                                done = false;
                                cur_col = c;
                            } else {
                                cur_col = std::min(cur_col, c);
                            }

                            break;
                        }
                    }

                    j[k] = beg;
                }
            }
        }
    }

    Ap.set_nonzeros(Ap.scan_row_sizes());

#pragma omp parallel
    {
        std::vector<ptrdiff_t> j(block_size);
        std::vector<ptrdiff_t> e(block_size);

#pragma omp for
        for(ptrdiff_t ip = 0; ip < np; ++ip) {
            ptrdiff_t ia = ip * block_size;
            ptrdiff_t cur_col = 0;
            ptrdiff_t head = Ap.ptr[ip];
            bool done = true;

            for(unsigned k = 0; k < block_size; ++k) {
                ptrdiff_t beg = j[k] = A.ptr[ia + k];
                ptrdiff_t end = e[k] = A.ptr[ia + k + 1];

                if (beg == end) continue;

                ptrdiff_t c = A.col[beg];

                if (done) {
                    done = false;
                    cur_col = c;
                } else {
                    cur_col = std::min(cur_col, c);
                }
            }

            while(!done) {
                cur_col /= block_size;

                Ap.col[head] = cur_col;

                done = true;
                bool first = true;
                S cur_val = math::zero<S>();

                ptrdiff_t col_end = (cur_col + 1) * block_size;
                for(unsigned k = 0; k < block_size; ++k) {
                    ptrdiff_t beg = j[k];
                    ptrdiff_t end = e[k];

                    while(beg < end) {
                        ptrdiff_t c = A.col[beg];
                        S v = math::norm(A.val[beg]);
                        ++beg;

                        if (c >= col_end) {
                            if (done) {
                                done = false;
                                cur_col = c;
                            } else {
                                cur_col = std::min(cur_col, c);
                            }

                            break;
                        }


                        if (first) {
                            first = false;
                            cur_val = v;
                        } else {
                            cur_val = std::max(cur_val, v);
                        }
                    }

                    j[k] = beg;
                }

                Ap.val[head++] = cur_val;
            }
        }
    }

    AMGCL_TOC("pointwise_matrix");
    return ap;
}

/** NUMA-aware vector container. */
template <class T>
class numa_vector {
    public:
        typedef T value_type;

        numa_vector() : n(0), p(0) {}

        numa_vector(size_t n, bool init = true) : n(n), p(new T[n]) {
            if (init) {
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                    p[i] = math::zero<T>();
            }
        }

        void resize(size_t size, bool init = true) {
            delete[] p; p = 0;

            n = size;
            p = new T[n];

            if (init) {
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                    p[i] = math::zero<T>();
            }
        }

        template <class Vector>
        numa_vector(const Vector &other,
                typename std::enable_if<!std::is_integral<Vector>::value, int>::type = 0
                ) : n(other.size()), p(new T[n])
        {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                p[i] = other[i];
        }

        template <class Iterator>
        numa_vector(Iterator beg, Iterator end)
            : n(std::distance(beg, end)), p(new T[n])
        {
            static_assert(std::is_same<
                    std::random_access_iterator_tag,
                    typename std::iterator_traits<Iterator>::iterator_category
                    >::value,
                    "Iterator has to be random access");

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                p[i] = beg[i];
        }

        ~numa_vector() {
            delete[] p; p = 0;
        }

        inline size_t size() const {
            return n;
        }

        inline const T& operator[](size_t i) const {
            return p[i];
        }

        inline T& operator[](size_t i) {
            return p[i];
        }

        inline const T* data() const {
            return p;
        }

        inline T* data() {
            return p;
        }

        void swap(numa_vector &other) {
            std::swap(n, other.n);
            std::swap(p, other.p);
        }

    private:
        size_t n;
        T *p;
};

/// Diagonal of a matrix
template < typename V, typename C, typename P >
std::shared_ptr< numa_vector<V> > diagonal(const crs<V, C, P> &A, bool invert = false)
{
    const size_t n = rows(A);
    auto dia = std::make_shared< numa_vector<V> >(n, false);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        for(auto a = A.row_begin(i); a; ++a) {
            if (a.col() == i) {
                (*dia)[i] = invert ? math::inverse(a.value()) : a.value();
                break;
            }
        }
    }

    return dia;
}

// Estimate spectral radius of the matrix.
// Use Gershgorin disk theorem when power_iters == 0,
// Use Power method when power_iters > 0.
// When scale = true, scale the matrix by its inverse diagonal.
template <bool scale, class Matrix>
static typename math::scalar_of<typename backend::value_type<Matrix>::type>::type
spectral_radius(const Matrix &A, int power_iters = 0) {
    AMGCL_TIC("spectral radius");
    typedef typename backend::value_type<Matrix>::type   value_type;
    typedef typename math::rhs_of<value_type>::type      rhs_type;
    typedef typename math::scalar_of<value_type>::type   scalar_type;

    const ptrdiff_t n = backend::rows(A);
    scalar_type radius;

    if (power_iters <= 0) {
        // Use Gershgorin disk theorem.
        radius = 0;

#pragma omp parallel
        {
            scalar_type emax = 0;
            value_type  dia = math::identity<value_type>();

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < n; ++i) {
                scalar_type s = 0;

                for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                    ptrdiff_t  c = A.col[j];
                    value_type v = A.val[j];

                    s += math::norm(v);

                    if (scale && c == i) dia = v;
                }

                if (scale) s *= math::norm(math::inverse(dia));

                emax = std::max(emax, s);
            }

#pragma omp critical
            radius = std::max(radius, emax);
        }
    } else {
        // Power method.
        backend::numa_vector<rhs_type> b0(n, false), b1(n, false);

        // Fill the initial vector with random values.
        // Also extract the inverted matrix diagonal values.
        scalar_type b0_norm = 0;
#pragma omp parallel
        {
#ifdef _OPENMP
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            std::mt19937 rng(tid);
            std::uniform_real_distribution<scalar_type> rnd(-1, 1);

            scalar_type loc_norm = 0;

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < n; ++i) {
                rhs_type v = math::constant<rhs_type>(rnd(rng));

                b0[i] = v;
                loc_norm += math::norm(math::inner_product(v,v));
            }

#pragma omp critical
            b0_norm += loc_norm;
        }

        // Normalize b0
        b0_norm = 1 / sqrt(b0_norm);
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            b0[i] = b0_norm * b0[i];
        }

        for(int iter = 0; iter < power_iters;) {
            // b1 = scale ? (D^1 * A) * b0 : A * b0
            // b1_norm = ||b1||
            // radius = <b1,b0>
            scalar_type b1_norm = 0;
            radius = 0;
#pragma omp parallel
            {
                scalar_type loc_norm = 0;
                scalar_type loc_radi = 0;
                value_type  dia = math::identity<value_type>();

#pragma omp for nowait
                for(ptrdiff_t i = 0; i < n; ++i) {
                    rhs_type s = math::zero<rhs_type>();

                    for(ptrdiff_t j = A.ptr[i], e = A.ptr[i+1]; j < e; ++j) {
                        ptrdiff_t  c = A.col[j];
                        value_type v = A.val[j];
                        if (scale && c == i) dia = v;
                        s += v * b0[c];
                    }

                    if (scale) s = math::inverse(dia) * s;

                    loc_norm += math::norm(math::inner_product(s, s));
                    loc_radi += math::norm(math::inner_product(s, b0[i]));

                    b1[i] = s;
                }

#pragma omp critical
                {
                    b1_norm += loc_norm;
                    radius  += loc_radi;
                }
            }

            if (++iter < power_iters) {
                // b0 = b1 / b1_norm
                b1_norm = 1 / sqrt(b1_norm);
#pragma omp parallel for
                for(ptrdiff_t i = 0; i < n; ++i) {
                    b0[i] = b1_norm * b1[i];
                }
            }
        }
    }
    AMGCL_TOC("spectral radius");

    return radius < 0 ? static_cast<scalar_type>(2) : radius;
}

/**
 * The builtin backend does not have any dependencies, and uses OpenMP for
 * parallelization. Matrices are stored in the CRS format, and vectors are
 * instances of ``std::vector<value_type>``. There is no usual overhead of
 * moving the constructed hierarchy to the builtin backend, since the backend
 * is used internally during setup.
 */
template <typename ValueType>
struct builtin {
    typedef ValueType      value_type;
    typedef ptrdiff_t      index_type;

    typedef typename math::rhs_of<value_type>::type rhs_type;

    struct provides_row_iterator : std::true_type {};

    typedef crs<value_type, index_type>    matrix;
    typedef numa_vector<rhs_type>          vector;
    typedef numa_vector<value_type>        matrix_diagonal;
    typedef solver::skyline_lu<value_type> direct_solver;

    /// The backend has no parameters.
    typedef amgcl::detail::empty_params params;

    static std::string name() { return "builtin"; }

    // Copy matrix. This is a noop for builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr<matrix> A, const params&)
    {
        return A;
    }

    // Copy vector to builtin backend.
    template <class T>
    static std::shared_ptr< numa_vector<T> >
    copy_vector(const std::vector<T> &x, const params&)
    {
        return std::make_shared< numa_vector<T> >(x);
    }

    // Copy vector to builtin backend. This is a noop for builtin backend.
    template <class T>
    static std::shared_ptr< numa_vector<T> >
    copy_vector(std::shared_ptr< numa_vector<T> > x, const params&)
    {
        return x;
    }

    // Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return std::make_shared<vector>(size);
    }

    struct gather {
        std::vector<ptrdiff_t> I;

        gather(size_t /*size*/, const std::vector<ptrdiff_t> &I, const params&)
            : I(I) { }

        template <class InVec, class OutVec>
        void operator()(const InVec &vec, OutVec &vals) const {
            for(size_t i = 0; i < I.size(); ++i)
                vals[i] = vec[I[i]];
        }
    };

    struct scatter {
        std::vector<ptrdiff_t> I;

        scatter(size_t /*size*/, const std::vector<ptrdiff_t> &I, const params&)
            : I(I) { }

        template <class InVec, class OutVec>
        void operator()(const InVec &vals, OutVec &vec) const {
            for(size_t i = 0; i < I.size(); ++i)
                vec[I[i]] = vals[i];
        }
    };

    // Create direct solver for coarse level
    static std::shared_ptr<direct_solver>
    create_solver(std::shared_ptr<matrix> A, const params&) {
        return std::make_shared<direct_solver>(*A);
    }
};

template <class T>
struct is_builtin_vector : std::false_type {};

template <class V>
struct is_builtin_vector< std::vector<V> > : std::is_arithmetic<V> {};

template <class V>
struct is_builtin_vector< numa_vector<V> > : std::true_type {};

//---------------------------------------------------------------------------
// Specialization of backend interface
//---------------------------------------------------------------------------
template <typename T1, typename T2>
struct backends_compatible< builtin<T1>, builtin<T2> > : std::true_type {};

template < typename V, typename C, typename P >
struct value_type< crs<V, C, P> > {
    typedef V type;
};

template < typename V, typename C, typename P >
struct rows_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A) {
        return A.nrows;
    }
};

template < typename V, typename C, typename P >
struct cols_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A) {
        return A.ncols;
    }
};

template < class Vec >
struct bytes_impl<
    Vec,
    typename std::enable_if< is_builtin_vector<Vec>::value >::type
    >
{
    static size_t get(const Vec &x) {
        typedef typename backend::value_type<Vec>::type V;
        return x.size() * sizeof(V);
    }
};

template < typename V, typename C, typename P >
struct ptr_data_impl< crs<V, C, P> > {
    typedef const P* type;
    static type get(const crs<V, C, P> &A) {
        return &A.ptr[0];
    }
};

template < typename V, typename C, typename P >
struct col_data_impl< crs<V, C, P> > {
    typedef const C* type;
    static type get(const crs<V, C, P> &A) {
        return &A.col[0];
    }
};

template < typename V, typename C, typename P >
struct val_data_impl< crs<V, C, P> > {
    typedef const V* type;
    static type get(const crs<V, C, P> &A) {
        return &A.val[0];
    }
};

template < typename V, typename C, typename P >
struct nonzeros_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A) {
        return A.nrows == 0 ? 0 : A.ptr[A.nrows];
    }
};

template < typename V, typename C, typename P >
struct row_iterator< crs<V, C, P> > {
    typedef
        typename crs<V, C, P>::row_iterator
        type;
};

template < typename V, typename C, typename P >
struct row_begin_impl< crs<V, C, P> > {
    typedef crs<V, C, P> Matrix;
    static typename row_iterator<Matrix>::type
    get(const Matrix &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

template < typename V, typename C, typename P >
struct row_nonzeros_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A, size_t row) {
        return A.ptr[row + 1] - A.ptr[row];
    }
};

template < class Vec >
struct clear_impl<
    Vec,
    typename std::enable_if< is_builtin_vector<Vec>::value >::type
    >
{
    static void apply(Vec &x)
    {
        typedef typename backend::value_type<Vec>::type V;

        const size_t n = x.size();
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            x[i] = math::zero<V>();
        }
    }
};

template < class Vec1, class Vec2 >
struct inner_product_impl<
    Vec1, Vec2,
    typename std::enable_if<
        is_builtin_vector<Vec1>::value &&
        is_builtin_vector<Vec2>::value
        >::type
    >
{
    typedef typename value_type<Vec1>::type V;

    typedef typename math::inner_product_impl<V>::return_type return_type;

    static return_type get(const Vec1 &x, const Vec2 &y)
    {
        const size_t n = x.size();
#ifdef _OPENMP
        return_type              _sum_stat[64];
        std::vector<return_type> _sum_dyna;
        return_type              *sum;

        const int nt = omp_get_max_threads();

        if (nt < 64) {
            sum = _sum_stat;
        } else {
            _sum_dyna.resize(nt);
            sum = _sum_dyna.data();
        }
#else
        const int nt = 1;
        return_type sum[1];
#endif


#pragma omp parallel
        {
#ifdef _OPENMP
            const int tid = omp_get_thread_num();
#else
            const int tid = 0;
#endif

            return_type s = math::zero<return_type>();
            return_type c = math::zero<return_type>();

#pragma omp for nowait
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                return_type d = math::inner_product(x[i], y[i]) - c;
                return_type t = s + d;
                c = (t - s) - d;
                s = t;
            }

            sum[tid] = s;
        }

        return std::accumulate(sum, sum + nt, math::zero<return_type>());
    }
};

template <class A, class Vec1, class B, class Vec2 >
struct axpby_impl<
    A, Vec1, B, Vec2,
    typename std::enable_if<
        is_builtin_vector<Vec1>::value &&
        is_builtin_vector<Vec2>::value
        >::type
    >
{
    static void apply(A a, const Vec1 &x, B b, Vec2 &y)
    {
        const size_t n = x.size();
        if (!math::is_zero(b)) {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                y[i] = a * x[i] + b * y[i];
            }
        } else {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                y[i] = a * x[i];
            }
        }
    }
};

template < class A, class Vec1, class B, class Vec2, class C, class Vec3 >
struct axpbypcz_impl<
    A, Vec1, B, Vec2, C, Vec3,
    typename std::enable_if<
        is_builtin_vector<Vec1>::value &&
        is_builtin_vector<Vec2>::value &&
        is_builtin_vector<Vec3>::value
        >::type
    >
{
    static void apply(A a, const Vec1 &x, B b, const Vec2 &y, C c, Vec3 &z)
    {
        const size_t n = x.size();
        if (!math::is_zero(c)) {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                z[i] = a * x[i] + b * y[i] + c * z[i];
            }
        } else {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                z[i] = a * x[i] + b * y[i];
            }
        }
    }
};

template < class Alpha, class Vec1, class Vec2, class Beta, class Vec3 >
struct vmul_impl<
    Alpha, Vec1, Vec2, Beta, Vec3,
    typename std::enable_if<
        is_builtin_vector<Vec1>::value &&
        is_builtin_vector<Vec2>::value &&
        is_builtin_vector<Vec3>::value
        >::type
    >
{
    static void apply(Alpha a, const Vec1 &x, const Vec2 &y, Beta b, Vec3 &z)
    {
        const size_t n = x.size();
        if (!math::is_zero(b)) {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                z[i] = a * x[i] * y[i] + b * z[i];
            }
        } else {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
                z[i] = a * x[i] * y[i];
            }
        }
    }
};

template < class Vec1, class Vec2 >
struct copy_impl<
    Vec1, Vec2,
    typename std::enable_if<
        is_builtin_vector<Vec1>::value &&
        is_builtin_vector<Vec2>::value
        >::type
    >
{
    static void apply(const Vec1 &x, Vec2 &y)
    {
        const size_t n = x.size();
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            y[i] = x[i];
        }
    }
};

namespace detail {

template <typename V, typename C, typename P>
struct use_builtin_matrix_ops< amgcl::backend::crs<V, C, P> >
    : std::true_type
{};

} // namespace detail

} // namespace backend
} // namespace amgcl


// Allow to use boost::iterator_range as vector in builtin backend:
namespace boost { template <class Iterator> class iterator_range; }

namespace amgcl {
namespace backend {

template <class Iterator>
struct is_builtin_vector< boost::iterator_range<Iterator> > : std::true_type {};

} // namespace backend
} // namespace amgcl

#endif
