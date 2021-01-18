#ifndef AMGCL_BACKEND_DETAIL_MATRIX_OPS_HPP
#define AMGCL_BACKEND_DETAIL_MATRIX_OPS_HPP

/*
The MIT License

Copyright (c) 2012-2020 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file    amgcl/adapter/detail/matrix_ops.hpp
 * \author  Denis Demidov <dennis.demidov@gmail.com>
 * \brief   Sparse matrix operations for matrices that provide row_iterator.
 */

#include <type_traits>
#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace backend {
namespace detail {

template <class Matrix, class Enable = void>
struct use_builtin_matrix_ops : std::false_type {};

} // namespace detail

template <class Alpha, class Matrix, class Vector1, class Beta, class Vector2>
struct spmv_impl<
    Alpha, Matrix, Vector1, Beta, Vector2,
    typename std::enable_if<
        detail::use_builtin_matrix_ops<Matrix>::value &&
        math::static_rows<typename value_type<Matrix>::type>::value == math::static_rows<typename value_type<Vector1>::type>::value &&
        math::static_rows<typename value_type<Matrix>::type>::value == math::static_rows<typename value_type<Vector2>::type>::value
        >::type
    >
{
    static void apply(
            Alpha alpha, const Matrix &A, const Vector1 &x, Beta beta, Vector2 &y
            )
    {
        typedef typename value_type<Vector2>::type V;

        const ptrdiff_t n = static_cast<ptrdiff_t>( rows(A) );

        if (!math::is_zero(beta)) {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                V sum = math::zero<V>();
                for(typename row_iterator<Matrix>::type a = row_begin(A, i); a; ++a)
                    sum += a.value() * x[ a.col() ];
                y[i] = alpha * sum + beta * y[i];
            }
        } else {
#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                V sum = math::zero<V>();
                for(typename row_iterator<Matrix>::type a = row_begin(A, i); a; ++a)
                    sum += a.value() * x[ a.col() ];
                y[i] = alpha * sum;
            }
        }
    }
};

template <class Matrix, class Vector1, class Vector2, class Vector3>
struct residual_impl<
    Matrix, Vector1, Vector2, Vector3,
    typename std::enable_if<
        detail::use_builtin_matrix_ops<Matrix>::value &&
        math::static_rows<typename value_type<Matrix>::type>::value == math::static_rows<typename value_type<Vector1>::type>::value &&
        math::static_rows<typename value_type<Matrix>::type>::value == math::static_rows<typename value_type<Vector2>::type>::value &&
        math::static_rows<typename value_type<Matrix>::type>::value == math::static_rows<typename value_type<Vector3>::type>::value
        >::type
    >
{
    static void apply(
            Vector1 const &rhs,
            Matrix  const &A,
            Vector2 const &x,
            Vector3       &res
            )
    {
        typedef typename value_type<Vector3>::type V;

        const ptrdiff_t n = static_cast<ptrdiff_t>( rows(A) );

#pragma omp parallel for
        for(ptrdiff_t i = 0; i < n; ++i) {
            V sum = math::zero<V>();
            for(typename row_iterator<Matrix>::type a = row_begin(A, i); a; ++a)
                sum += a.value() * x[ a.col() ];
            res[i] = rhs[i] - sum;
        }
    }
};

/* Allows to do matrix-vector products with mixed scalar/nonscalar types.
 * Reinterprets pointers to the vectors data into appropriate types.
 */
template <class Alpha, class Matrix, class Vector1, class Beta, class Vector2>
struct spmv_impl<
    Alpha, Matrix, Vector1, Beta, Vector2,
    typename std::enable_if<
            detail::use_builtin_matrix_ops<Matrix>::value && (
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector1>::type>::value ||
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector2>::type>::value)
        >::type
    >
{
    static void apply(
            Alpha alpha, const Matrix &A, const Vector1 &x, Beta beta, Vector2 &y
            )
    {
        typedef typename value_type<Matrix>::type     val_type;
        typedef typename math::rhs_of<val_type>::type rhs_type;
        typedef typename math::replace_scalar<rhs_type, typename math::scalar_of<typename value_type<Vector1>::type>::type>::type x_type;
        typedef typename math::replace_scalar<rhs_type, typename math::scalar_of<typename value_type<Vector2>::type>::type>::type y_type;

        const size_t n = backend::rows(A);
        const size_t m = backend::cols(A);

        x_type const * xptr = reinterpret_cast<x_type const *>(&x[0]);
        y_type       * yptr = reinterpret_cast<y_type       *>(&y[0]);

        iterator_range<x_type const *> xrng(xptr, xptr + m);
        iterator_range<y_type       *> yrng(yptr, yptr + n);

        spmv(alpha, A, xrng, beta, yrng);
    }
};

template <class Matrix, class Vector1, class Vector2, class Vector3>
struct residual_impl<
    Matrix, Vector1, Vector2, Vector3,
    typename std::enable_if<
            detail::use_builtin_matrix_ops<Matrix>::value && (
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector1>::type>::value ||
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector2>::type>::value ||
            math::static_rows<typename value_type<Matrix>::type>::value != math::static_rows<typename value_type<Vector3>::type>::value)
        >::type
    >
{
    static void apply(
            Vector1 const &f,
            Matrix  const &A,
            Vector2 const &x,
            Vector3       &r
            )
    {
        typedef typename value_type<Matrix>::type     val_type;
        typedef typename math::rhs_of<val_type>::type rhs_type;

        typedef typename math::replace_scalar<rhs_type, typename math::scalar_of<typename value_type<Vector1>::type>::type>::type f_type;
        typedef typename math::replace_scalar<rhs_type, typename math::scalar_of<typename value_type<Vector2>::type>::type>::type x_type;
        typedef typename math::replace_scalar<rhs_type, typename math::scalar_of<typename value_type<Vector3>::type>::type>::type r_type;

        const size_t n = backend::rows(A);
        const size_t m = backend::cols(A);

        x_type const * xptr = reinterpret_cast<x_type const *>(&x[0]);
        f_type const * fptr = reinterpret_cast<f_type const *>(&f[0]);
        r_type       * rptr = reinterpret_cast<r_type       *>(&r[0]);

        iterator_range<x_type const *> xrng(xptr, xptr + m);
        iterator_range<f_type const *> frng(fptr, fptr + n);
        iterator_range<r_type       *> rrng(rptr, rptr + n);

        residual(frng, A, xrng, rrng);
    }
};

} // namespace backend
} // namespace amgcl

#endif
