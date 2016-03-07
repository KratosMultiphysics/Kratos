#ifndef AMGCL_BACKEND_BUILTIN_HPP
#define AMGCL_BACKEND_BUILTIN_HPP

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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

#include <boost/typeof/typeof.hpp>
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

    virtual ~crs() {}

    virtual const ptr_type* ptr_data() const { return &ptr[0]; }
    virtual const col_type* col_data() const { return &col[0]; }
    virtual const val_type* val_data() const { return &val[0]; }

    virtual ptr_type* ptr_data() { return &ptr[0]; }
    virtual col_type* col_data() { return &col[0]; }
    virtual val_type* val_data() { return &val[0]; }

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
        ptr_type p = ptr_data()[row];
        ptr_type e = ptr_data()[row + 1];
        return row_iterator(col_data() + p, col_data() + e, val_data() + p);
    }

};

/// Sort rows of the matrix column-wise.
template < typename V, typename C, typename P >
void sort_rows(crs<V, C, P> &A) {
    const size_t n = rows(A);
    BOOST_AUTO(Aptr, A.ptr_data());
    BOOST_AUTO(Acol, A.col_data());
    BOOST_AUTO(Aval, A.val_data());

#pragma omp parallel for
    for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i) {
        P beg = Aptr[i];
        P end = Aptr[i + 1];
        amgcl::detail::sort_row(Acol + beg, Aval + beg, end - beg);
    }
}

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

    const P* Aptr = A.ptr_data();
    const C* Acol = A.col_data();
    const V* Aval = A.val_data();

    for(size_t j = 0; j < nnz; ++j)
        ++( T.ptr[Acol[j] + 1] );

    boost::partial_sum(T.ptr, T.ptr.begin());

    for(size_t i = 0; i < n; i++) {
        for(P j = Aptr[i], e = Aptr[i + 1]; j < e; ++j) {
            P head = T.ptr[Acol[j]]++;

            T.col[head] = static_cast<C>(i);
            T.val[head] = Aval[j];
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

    crs<V, ptrdiff_t> C;
    C.nrows = rows(A);
    C.ncols = cols(B);

#ifdef _OPENMP
    int nt = omp_get_max_threads();
#else
    int nt = 1;
#endif

    if (nt > 4) {
        spgemm_rmerge(
                static_cast<ptrdiff_t>(C.nrows),
                ptr_data(A), col_data(A), val_data(A),
                ptr_data(B), col_data(B), val_data(B),
                C.ptr, C.col, C.val
                );
    } else {
        spgemm_saad(
                static_cast<ptrdiff_t>(C.nrows),
                static_cast<ptrdiff_t>(C.ncols),
                ptr_data(A), col_data(A), val_data(A),
                ptr_data(B), col_data(B), val_data(B),
                C.ptr, C.col, C.val, sort
                );
    }

    return C;
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
                dia[i] = invert ? math::inverse(a.value()) : a.value();
                break;
            }
        }
    }

    return dia;
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

    amgcl::detail::inverse(n, &Ainv.val[0]);

    Ainv.ptr[0] = 0;
    for(size_t i = 0, idx = 0; i < n; ) {
        for(size_t j = 0; j < n; ++j, ++idx) Ainv.col[idx] = static_cast<C>(j);

        Ainv.ptr[++i] = static_cast<P>(idx);
    }

    return Ainv;
}

/**
 * The builtin backend does not have any dependencies except for the
 * <a href="http://www.boost.org">Boost</a> libraries, and uses OpenMP for
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

    struct provides_row_iterator : boost::true_type {};

    typedef crs<value_type, index_type>    matrix;
    typedef std::vector<rhs_type>          vector;
    typedef std::vector<value_type>        matrix_diagonal;
    typedef solver::skyline_lu<value_type> direct_solver;

    /// The backend has no parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
        void get(boost::property_tree::ptree&, const std::string&) const {}
    };

    static std::string name() { return "builtin"; }

    // Copy matrix. This is a noop for builtin backend.
    static boost::shared_ptr<matrix>
    copy_matrix(boost::shared_ptr<matrix> A, const params&)
    {
        return A;
    }

    // Copy vector to builtin backend.
    template <class T>
    static boost::shared_ptr< std::vector<T> >
    copy_vector(const std::vector<T> &x, const params&)
    {
        return boost::make_shared< std::vector<T> >(x);
    }

    // Copy vector to builtin backend. This is a noop for builtin backend.
    template <class T>
    static boost::shared_ptr< std::vector<T> >
    copy_vector(boost::shared_ptr< std::vector<T> > x, const params&)
    {
        return x;
    }

    // Create vector of the specified size.
    static boost::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return boost::make_shared<vector>(size);
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
        return A.ptr_data();
    }
};

template < typename V, typename C, typename P >
struct col_data_impl< crs<V, C, P> > {
    typedef const C* type;
    static type get(const crs<V, C, P> &A) {
        return A.col_data();
    }
};

template < typename V, typename C, typename P >
struct val_data_impl< crs<V, C, P> > {
    typedef const V* type;
    static type get(const crs<V, C, P> &A) {
        return A.val_data();
    }
};

template < typename V, typename C, typename P >
struct nonzeros_impl< crs<V, C, P> > {
    static size_t get(const crs<V, C, P> &A) {
        return A.nrows == 0 ? 0 : A.ptr_data()[A.nrows];
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
        const P *Aptr = A.ptr_data();
        return Aptr[row + 1] - Aptr[row];
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type
                >::type
        >::type
    >
{
    typedef typename value_type<Vec1>::type V;

    typedef typename math::inner_product_impl<V>::return_type return_type;

    static return_type get(const Vec1 &x, const Vec2 &y)
    {
        const size_t n = x.size();
        return_type sum = math::zero<return_type>();

#pragma omp parallel
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

            return_type s = math::zero<return_type>();
            return_type c = math::zero<return_type>();
            for(size_t i = chunk_start; i < chunk_end; ++i) {
                return_type d = math::inner_product(x[i], y[i]) - c;
                return_type t = s + d;
                c = (t - s) - d;
                s = t;
            }
#pragma omp critical
            sum += s;
        }
        return sum;
    }
};

template <class A, class Vec1, class B, class Vec2 >
struct axpby_impl<
    A, Vec1, B, Vec2,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type
                >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type,
                typename is_builtin_vector<Vec3>::type
                >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_builtin_vector<Vec1>::type,
                typename is_builtin_vector<Vec2>::type,
                typename is_builtin_vector<Vec3>::type
                >::type
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
