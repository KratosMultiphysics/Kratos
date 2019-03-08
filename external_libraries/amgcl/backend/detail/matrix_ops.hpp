#ifndef AMGCL_BACKEND_DETAIL_MATRIX_OPS_HPP
#define AMGCL_BACKEND_DETAIL_MATRIX_OPS_HPP

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
        detail::use_builtin_matrix_ops<Matrix>::value
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

} // namespace backend
} // namespace amgcl

#endif
