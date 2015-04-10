#ifndef AMGCL_BACKEND_BUILTIN_HPP
#define AMGCL_BACKEND_BUILTIN_HPP

/*
The MIT License

Copyright (c) 2012-2015 Denis Demidov <dennis.demidov@gmail.com>

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

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/skyline_lu.hpp>
#include <amgcl/detail/inverse.hpp>
#include <amgcl/detail/sort_row.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>

namespace amgcl {
namespace backend {

/// Sparse matrix stored in CRS format.
template <
    typename val_t = double,
    typename col_t = int,
    typename ptr_t = col_t
    >
struct crs {
    typedef val_t val_type;
    typedef col_t col_type;
    typedef ptr_t ptr_type;

    size_t nrows, ncols;
    std::vector<ptr_type> ptr;
    std::vector<col_type> col;
    std::vector<val_type> val;

    crs() : nrows(0), ncols(0) {}

    template <
        class PtrRange,
        class ColRange,
        class ValRange
        >
    crs(size_t nrows, size_t ncols,
        const PtrRange &ptr_range,
        const ColRange &col_range,
        const ValRange &val_range
        )
    : nrows(nrows), ncols(ncols),
      ptr(boost::begin(ptr_range), boost::end(ptr_range)),
      col(boost::begin(col_range), boost::end(col_range)),
      val(boost::begin(val_range), boost::end(val_range))
    {
        precondition(
                ptr.size() == nrows + 1                       &&
                static_cast<size_t>(ptr.back()) == col.size() &&
                static_cast<size_t>(ptr.back()) == val.size(),
                "Inconsistent sizes in crs constructor"
                );
    }

    template <class Matrix>
    crs(const Matrix &A) : nrows(backend::rows(A)), ncols(backend::cols(A))
    {
        ptr.reserve(nrows + 1);
        ptr.push_back(0);

        col.reserve(backend::nonzeros(A));
        val.reserve(backend::nonzeros(A));

        typedef typename backend::row_iterator<Matrix>::type row_iterator;
        for(size_t i = 0; i < nrows; ++i) {
            for(row_iterator a = backend::row_begin(A, i); a; ++a) {
                col.push_back(a.col());
                val.push_back(a.value());
            }
            ptr.push_back( static_cast<ptr_type>(col.size()) );
        }
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
        return row_iterator(&col[0] + p, &col[0] + e, &val[0] + p);
    }

};

/// Transpose of a sparse matrix.
template < typename V, typename C, typename P >
crs<V, C, P> transpose(const crs<V, C, P> &A)
{
    const size_t n   = rows(A);
    const size_t m   = cols(A);
    const size_t nnz = nonzeros(A);

    crs<V, C, P> T;
    T.nrows = m;
    T.ncols = n;
    T.ptr.resize(m+1);
    T.col.resize(nnz);
    T.val.resize(nnz);

    boost::fill(T.ptr, P());

    for(size_t j = 0; j < nnz; ++j)
        ++( T.ptr[A.col[j] + 1] );

    boost::partial_sum(T.ptr, T.ptr.begin());

    for(size_t i = 0; i < n; i++) {
        for(P j = A.ptr[i], e = A.ptr[i + 1]; j < e; ++j) {
            P head = T.ptr[A.col[j]]++;

            T.col[head] = static_cast<C>(i);
            T.val[head] = A.val[j];
        }
    }

    std::rotate(T.ptr.begin(), T.ptr.end() - 1, T.ptr.end());
    T.ptr.front() = 0;

    return T;
}

/// Matrix-matrix product.
template <class MatrixA, class MatrixB>
crs< typename value_type<MatrixA>::type >
product(const MatrixA &A, const MatrixB &B, bool sort = false) {
    typedef typename value_type<MatrixA>::type  V;
    typedef ptrdiff_t C;
    typedef ptrdiff_t P;

    typedef typename row_iterator<MatrixA>::type Aiterator;
    typedef typename row_iterator<MatrixB>::type Biterator;

    const size_t n = rows(A);
    const size_t m = cols(B);

    crs<V, C, P> c;
    c.nrows = n;
    c.ncols = m;
    c.ptr.resize(n + 1);
    boost::fill(c.ptr, P());

#pragma omp parallel
    {
        std::vector<ptrdiff_t> marker(m, -1);

#ifdef _OPENMP
        int nt  = omp_get_num_threads();
        int tid = omp_get_thread_num();

        size_t chunk_size  = (n + nt - 1) / nt;
        size_t chunk_start = tid * chunk_size;
        size_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
        size_t chunk_start = 0;
        size_t chunk_end   = n;
#endif

        for(size_t ia = chunk_start; ia < chunk_end; ++ia) {
            for(Aiterator a = A.row_begin(ia); a; ++a) {
                for(Biterator b = B.row_begin(a.col()); b; ++b) {
                    if (marker[b.col()] != static_cast<C>(ia)) {
                        marker[b.col()]  = static_cast<C>(ia);
                        ++( c.ptr[ia + 1] );
                    }
                }
            }
        }

        boost::fill(marker, -1);

#pragma omp barrier
#pragma omp single
        {
            boost::partial_sum(c.ptr, c.ptr.begin());
            c.col.resize(c.ptr.back());
            c.val.resize(c.ptr.back());
        }

        for(size_t ia = chunk_start; ia < chunk_end; ++ia) {
            P row_beg = c.ptr[ia];
            P row_end = row_beg;

            for(Aiterator a = A.row_begin(ia); a; ++a) {
                C ca = a.col();
                V va = a.value();

                for(Biterator b = B.row_begin(ca); b; ++b) {
                    C cb = b.col();
                    V vb = b.value();

                    if (marker[cb] < row_beg) {
                        marker[cb] = row_end;
                        c.col[row_end] = cb;
                        c.val[row_end] = va * vb;
                        ++row_end;
                    } else {
                        c.val[marker[cb]] += va * vb;
                    }
                }
            }

            if (sort) {
                amgcl::detail::sort_row(
                        &c.col[row_beg], &c.val[row_beg], row_end - row_beg
                        );
            }
        }
    }

    return c;
}

/// Diagonal of a matrix
template < typename V, typename C, typename P >
std::vector<V> diagonal(const crs<V, C, P> &A, bool invert = false)
{
    typedef typename crs<V, C, P>::row_iterator row_iterator;
    const size_t n = rows(A);
    std::vector<V> dia(n);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        for(row_iterator a = A.row_begin(i); a; ++a) {
            if (a.col() == i) {
                dia[i] = invert ? 1 / a.value() : a.value();
                break;
            }
        }
    }

    return dia;
}

/// Sort rows of the matrix column-wise.
template < typename V, typename C, typename P >
void sort_rows(crs<V, C, P> &A) {
    const size_t n = rows(A);

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        P beg = A.ptr[i];
        P end = A.ptr[i + 1];
        amgcl::detail::sort_row(A.col.data() + beg, A.val.data() + beg, end - beg);
    }
}

/// Invert matrix.
template < typename V, typename C, typename P >
crs<V, C, P> inverse(const crs<V, C, P> &A) {
    typedef typename crs<V, C, P>::row_iterator row_iterator;
    const size_t n = rows(A);

    crs<V, C, P> Ainv;
    Ainv.nrows = n;
    Ainv.ncols = n;
    Ainv.ptr.resize(n + 1);
    Ainv.col.resize(n * n);
    Ainv.val.resize(n * n);

    boost::fill(Ainv.val, V());

    for(size_t i = 0; i < n; ++i)
        for(row_iterator a = A.row_begin(i); a; ++a)
            Ainv.val[i * n + a.col()] = a.value();

    amgcl::detail::inverse(n, Ainv.val.data());

    Ainv.ptr[0] = 0;
    for(size_t i = 0, idx = 0; i < n; ) {
        for(size_t j = 0; j < n; ++j, ++idx) Ainv.col[idx] = static_cast<C>(j);

        Ainv.ptr[++i] = static_cast<P>(idx);
    }

    return Ainv;
}

/// Builtin backend.
/**
 * \param real Value type.
 * \ingroup backends
 */
template <typename real>
struct builtin {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    struct provides_row_iterator : boost::true_type {};

    typedef crs<value_type, index_type>    matrix;
    typedef std::vector<value_type>        vector;
    typedef solver::skyline_lu<value_type> direct_solver;

    /// Backend parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
        void get(boost::property_tree::ptree&, const std::string&) const {}
    };

    static std::string name() { return "builtin"; }

    /// Copy matrix.
    /** This is a noop for builtin backend. */
    static boost::shared_ptr<matrix>
    copy_matrix(boost::shared_ptr<matrix> A, const params&)
    {
        return A;
    }

    /// Copy vector to builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(const vector &x, const params&)
    {
        return boost::make_shared<vector>(x);
    }

    /// Copy vector to builtin backend.
    /** This is a noop for builtin backend. */
    static boost::shared_ptr<vector>
    copy_vector(boost::shared_ptr< vector > x, const params&)
    {
        return x;
    }

    /// Create vector of the specified size.
    static boost::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return boost::make_shared<vector>(size);
    }

    struct gather {
        std::vector<ptrdiff_t> I;

        gather(size_t /*src_size*/, const std::vector<ptrdiff_t> &I, const params&)
            : I(I) { }

        template <class InVec, class OutVec>
        void operator()(const InVec &vec, OutVec &vals) const {
            for(size_t i = 0; i < I.size(); ++i)
                vals[i] = vec[I[i]];
        }
    };

    /// Create direct solver for coarse level
    static boost::shared_ptr<direct_solver>
    create_solver(boost::shared_ptr<matrix> A, const params&) {
        return boost::make_shared<direct_solver>(*A);
    }
};

template <class T>
struct is_builtin_vector : boost::false_type {};

template <class V>
struct is_builtin_vector< std::vector<V> > : boost::is_arithmetic<V> {};

template <class Iterator>
struct is_builtin_vector< boost::iterator_range<Iterator> > : boost::true_type {};

//---------------------------------------------------------------------------
// Specialization of backend interface
//---------------------------------------------------------------------------
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

template < typename V, typename C, typename P >
struct ptr_data_impl< crs<V, C, P> > {
    typedef const P* type;
    static type get(const crs<V, C, P> &A) {
        return A.ptr.data();
    }
};

template < typename V, typename C, typename P >
struct col_data_impl< crs<V, C, P> > {
    typedef const C* type;
    static type get(const crs<V, C, P> &A) {
        return A.col.data();
    }
};

template < typename V, typename C, typename P >
struct val_data_impl< crs<V, C, P> > {
    typedef const V* type;
    static type get(const crs<V, C, P> &A) {
        return A.val.data();
    }
};

template < typename V, typename C, typename P >
struct nonzeros_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A) {
        return A.ptr.empty() ? 0 : A.ptr.back();
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
    typename boost::enable_if< typename is_builtin_vector<Vec>::type >::type
    >
{
    static void apply(Vec &x)
    {
        const size_t n = x.size();
#pragma omp parallel for
        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
            x[i] = 0;
        }
    }
};

template < class Vec1, class Vec2 >
struct inner_product_impl<
    Vec1, Vec2,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type
                >::type
        >::type
    >
{
    typedef typename value_type<Vec1>::type V;

    static V get(const Vec1 &x, const Vec2 &y)
    {
        const size_t n = x.size();
        V sum = 0;

#pragma omp parallel reduction(+:sum)
        {
#ifdef _OPENMP
            int nt  = omp_get_num_threads();
            int tid = omp_get_thread_num();

            size_t chunk_size  = (n + nt - 1) / nt;
            size_t chunk_start = tid * chunk_size;
            size_t chunk_end   = std::min(n, chunk_start + chunk_size);
#else
            size_t chunk_start = 0;
            size_t chunk_end   = n;
#endif

            V s = 0;
            V c = 0;
            for(size_t i = chunk_start; i < chunk_end; ++i) {
                V d = x[i] * y[i] - c;
                V t = s + d;
                c = (t - s) - d;
                s = t;
            }

            sum += s;
        }
        return sum;
    }
};

template < class Vec1, class Vec2 >
struct axpby_impl<
    Vec1, Vec2,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type
                >::type
        >::type
    >
{
    typedef typename value_type<Vec2>::type V;

    static void apply(V a, const Vec1 &x, V b, Vec2 &y)
    {
        const size_t n = x.size();
        if (b) {
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

template < class Vec1, class Vec2, class Vec3 >
struct axpbypcz_impl<
    Vec1, Vec2, Vec3,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type,
                typename is_builtin_vector<Vec3>::type
                >::type
        >::type
    >
{
    typedef typename value_type<Vec3>::type V;

    static void apply(V a, const Vec1 &x, V b, const Vec2 &y, V c, Vec3 &z)
    {
        const size_t n = x.size();
        if (c) {
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

template < class Vec1, class Vec2, class Vec3 >
struct vmul_impl<
    Vec1, Vec2, Vec3,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type,
                typename is_builtin_vector<Vec3>::type
                >::type
        >::type
    >
{
    typedef typename value_type<Vec3>::type V;

    static void apply(V a, const Vec1 &x, const Vec2 &y, V b, Vec3 &z)
    {
        const size_t n = x.size();
        if (b) {
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type
                >::type
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

template < class Vec >
struct copy_to_backend_impl<
    Vec,
    typename boost::enable_if<
            typename is_builtin_vector<Vec>::type
        >::type
    > : copy_impl< std::vector<typename value_type<Vec>::type>, Vec > {};

namespace detail {

template <typename V, typename C, typename P>
struct use_builtin_matrix_ops< amgcl::backend::crs<V, C, P> >
    : boost::true_type
{};

} // namespace detail

} // namespace backend
} // namespace amgcl

#endif
