#ifndef AMGCL_ADAPTER_REORDER_HPP
#define AMGCL_ADAPTER_REORDER_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
\file    amgcl/adapter/reorder.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   On-the-fly reodering of matrix and vectors.
\ingroup adapters
*/

#include <type_traits>
#include <boost/range/size.hpp>
#include <boost/iterator/permutation_iterator.hpp>


#include <amgcl/reorder/cuthill_mckee.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>

namespace amgcl {
namespace adapter {

template <class Matrix>
struct reordered_matrix {
    typedef typename backend::value_type<Matrix>::type value_type;
    typedef typename backend::row_iterator<Matrix>::type base_iterator;

    const Matrix &A;
    const ptrdiff_t * perm;
    const ptrdiff_t * iperm;

    reordered_matrix(const Matrix &A, const ptrdiff_t *perm, const ptrdiff_t * iperm)
        : A(A), perm(perm), iperm(iperm)
    {}

    size_t rows() const {
        return backend::rows(A);
    }

    size_t cols() const {
        return backend::cols(A);
    }

    size_t nonzeros() const {
        return backend::nonzeros(A);
    }

    struct row_iterator {
        base_iterator base;
        const ptrdiff_t * iperm;

        row_iterator(const base_iterator &base, const ptrdiff_t *iperm)
            : base(base), iperm(iperm)
        {}

        operator bool() const {
            return base;
        }

        row_iterator& operator++() {
            ++base;
            return *this;
        }

        ptrdiff_t col() const {
            return iperm[base.col()];
        }

        value_type value() const {
            return base.value();
        }
    };

    row_iterator row_begin(size_t i) const {
        return row_iterator(backend::row_begin(A, perm[i]), iperm);
    }
};

template <class Vector>
struct reordered_vector {
    typedef typename backend::value_type<typename std::decay<Vector>::type>::type raw_value_type;
    typedef typename std::conditional<
        std::is_const<Vector>::value,
        const raw_value_type,
        raw_value_type
        >::type value_type;

    Vector &x;
    const ptrdiff_t *perm;

    reordered_vector(Vector &x, const ptrdiff_t *perm) : x(x), perm(perm) {}

    size_t size() const {
        return boost::size(x);
    }

    value_type& operator[](size_t i) const {
        return x[perm[i]];
    }

    boost::permutation_iterator<
        typename std::decay<Vector>::type::iterator,
        const ptrdiff_t*
        >
    begin() {
        return boost::make_permutation_iterator(boost::begin(x), perm);
    }

    boost::permutation_iterator<
        typename std::decay<Vector>::type::const_iterator,
        const ptrdiff_t*
        >
    begin() const {
        return boost::make_permutation_iterator(boost::begin(x), perm);
    }

    boost::permutation_iterator<
        typename std::decay<Vector>::type::iterator,
        const ptrdiff_t*
        >
    end() {
        return boost::make_permutation_iterator(boost::end(x), perm + size());
    }

    boost::permutation_iterator<
        typename std::decay<Vector>::type::const_iterator,
        const ptrdiff_t*
        >
    end() const {
        return boost::make_permutation_iterator(boost::end(x), perm + size());
    }
};

} // namespace adapter

namespace backend {

//---------------------------------------------------------------------------
// Specialization of matrix interface
//---------------------------------------------------------------------------
template <class Matrix>
struct value_type< adapter::reordered_matrix<Matrix> >
    : value_type<Matrix>
{};

template <class Matrix>
struct rows_impl< adapter::reordered_matrix<Matrix> >
{
    static size_t get(const adapter::reordered_matrix<Matrix> &A) {
        return A.rows();
    }
};

template <class Matrix>
struct cols_impl< adapter::reordered_matrix<Matrix> >
{
    static size_t get(const adapter::reordered_matrix<Matrix> &A) {
        return A.cols();
    }
};

template <class Matrix>
struct nonzeros_impl< adapter::reordered_matrix<Matrix> >
{
    static size_t get(const adapter::reordered_matrix<Matrix> &A) {
        return A.nonzeros();
    }
};

template <class Matrix>
struct row_iterator< adapter::reordered_matrix<Matrix> >
{
    typedef
        typename adapter::reordered_matrix<Matrix>::row_iterator
        type;
};

template <class Matrix>
struct row_begin_impl< adapter::reordered_matrix<Matrix> >
{
    typedef adapter::reordered_matrix<Matrix> M;
    static typename row_iterator<M>::type
    get(const M &matrix, size_t row) {
        return matrix.row_begin(row);
    }
};

namespace detail {

template <class Matrix>
struct use_builtin_matrix_ops< adapter::reordered_matrix<Matrix> >
    : std::true_type
{};

} // namespace detail


template <class Vector>
struct is_builtin_vector< adapter::reordered_vector<Vector> >
    : is_builtin_vector<typename std::decay<Vector>::type>
{};

} // namespace backend

namespace adapter {

template <class ordering = amgcl::reorder::cuthill_mckee<false> >
class reorder {
    public:
        template <class Matrix>
        reorder(const Matrix &A) : n(backend::rows(A)), perm(n), iperm(n)
        {
            ordering::get(A, perm);
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) iperm[perm[i]] = i;
        }

        template <class Matrix>
        typename std::enable_if<
            !backend::is_builtin_vector<Matrix>::value,
            reordered_matrix<Matrix>
        >::type
        operator()(const Matrix &A) const {
            return reordered_matrix<Matrix>(A, perm.data(), iperm.data());
        }

        template <class Vector>
        typename std::enable_if<
            backend::is_builtin_vector<Vector>::value,
            reordered_vector<Vector>
        >::type
        operator()(Vector &x) const {
            return reordered_vector<Vector>(x, perm.data());
        }

        template <class Vector>
        typename std::enable_if<
            backend::is_builtin_vector<Vector>::value,
            reordered_vector<const Vector>
        >::type
        operator()(const Vector &x) const {
            return reordered_vector<const Vector>(x, perm.data());
        }

        template <class Vector1, class Vector2>
        void forward(const Vector1 &x, Vector2 &y) const {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) y[i] = x[perm[i]];
        }

        template <class Vector1, class Vector2>
        void inverse(const Vector1 &x, Vector2 &y) const {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) y[perm[i]] = x[i];
        }

    private:
        ptrdiff_t n;
        backend::numa_vector<ptrdiff_t> perm, iperm;
};

} // namespace adapter
} // namespace amgcl

#endif
