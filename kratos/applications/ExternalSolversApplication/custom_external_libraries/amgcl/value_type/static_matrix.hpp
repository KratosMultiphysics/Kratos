#ifndef AMGCL_VALUE_TYPE_STATIC_MATRIX_HPP
#define AMGCL_VALUE_TYPE_STATIC_MATRIX_HPP

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
 * \file   amgcl/value_type/static_matrix.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Enable statically sized matrices as value types.
 */

#include <boost/array.hpp>
#include <boost/type_traits.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {

template <typename T, int N, int M>
struct static_matrix {
    boost::array<T, N * M> buf;

    T operator()(int i, int j) const {
        return buf[i * M + j];
    }

    T& operator()(int i, int j) {
        return buf[i * M + j];
    }

    T operator()(int i) const {
        return buf[i];
    }

    T& operator()(int i) {
        return buf[i];
    }

    const T* data() const {
        return buf.data();
    }

    T* data() {
        return buf.data();
    }

    const static_matrix& operator+=(const static_matrix &y) {
        for(int i = 0; i < N * M; ++i)
            buf[i] += y.buf[i];
        return *this;
    }

    const static_matrix& operator-=(const static_matrix &y) {
        for(int i = 0; i < N * M; ++i)
            buf[i] -= y.buf[i];
        return *this;
    }

    const static_matrix& operator*=(T c) {
        for(int i = 0; i < N * M; ++i)
            buf[i] *= c;
        return *this;
    }


    friend static_matrix operator+(static_matrix x, const static_matrix &y)
    {
        return x += y;
    }

    friend static_matrix operator-(static_matrix x, const static_matrix &y)
    {
        return x -= y;
    }

    friend static_matrix operator*(T a, static_matrix x)
    {
        return x *= a;
    }

    friend static_matrix operator-(static_matrix x)
    {
        for(int i = 0; i < N * M; ++i)
            x.buf[i] = -x.buf[i];
        return x;
    }

    friend bool operator<(const static_matrix &x, const static_matrix &y)
    {
        T xtrace = math::zero<T>();
        T ytrace = math::zero<T>();

        const int K = N < M ? N : M;

        for(int i = 0; i < K; ++i) {
            xtrace += x(i,i);
            ytrace += y(i,i);
        }

        return xtrace < ytrace;
    }
};

template <typename T, int N, int K, int M>
static_matrix<T, N, M> operator*(
        const static_matrix<T, N, K> &a,
        const static_matrix<T, K, M> &b
        )
{
    static_matrix<T, N, M> c;
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < M; ++j)
            c(i,j) = math::zero<T>();
        for(int k = 0; k < K; ++k) {
            T aik = a(i,k);
            for(int j = 0; j < M; ++j)
                c(i,j) += aik * b(k,j);
        }
    }
    return c;
}

template <class T> struct is_static_matrix : boost::false_type {};

template <class T, int N, int M>
struct is_static_matrix< static_matrix<T, N, M> > : boost::true_type {};

namespace backend {

/// Enable static matrix as a value-type.
template <typename T, int N, int M>
struct is_builtin_vector< std::vector<static_matrix<T, N, M> > > : boost::true_type {};

} // namespace backend

namespace math {

/// Scalar type of a non-scalar type.
template <class T, int N, int M>
struct scalar_of< static_matrix<T, N, M> > {
    typedef T type;
};

/// RHS type corresponding to a non-scalar type.
template <class T, int N>
struct rhs_of< static_matrix<T, N, N> > {
    typedef static_matrix<T, N, 1> type;
};

/// Specialization of conjugate transpose for static matrices.
template <typename T, int N, int M>
struct adjoint_impl< static_matrix<T, N, M> >
{
    typedef static_matrix<T, M, N> return_type;

    static static_matrix<T, M, N> get(const static_matrix<T, N, M> &x) {
        static_matrix<T, M, N> y;
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < M; ++j)
                y(j,i) = x(i,j);
        return y;
    }
};

/// Inner-product result of two static vectors.
template <class T, int N>
struct inner_product_impl< static_matrix<T, N, 1> >
{
    typedef T return_type;
    static T get(const static_matrix<T, N, 1> &x, const static_matrix<T, N, 1> &y) {
        T sum = math::zero<T>();
        for(int i = 0; i < N; ++i)
            sum += x(i) * y(i);
        return sum;
    }
};

/// Inner-product result of two static matrices.
template <class T, int N, int M>
struct inner_product_impl< static_matrix<T, N, M> >
{
    typedef static_matrix<T, M, M> return_type;

    static return_type get(const static_matrix<T, N, M> &x, const static_matrix<T, N, M> &y) {
        static_matrix<T, M, M> p;
        for(int i = 0; i < M; ++i) {
            for(int j = 0; j < M; ++j) {
                T sum = math::zero<T>();
                for(int k = 0; k < N; ++k)
                    sum += x(k,i) * y(k,j);
                p(i,j) = sum;
            }
        }
        return p;
    }
};

/// Implementation of Frobenius norm for static matrices.
template <typename T, int N, int M>
struct norm_impl< static_matrix<T, N, M> >
{
    static T get(const static_matrix<T, N, M> &x) {
        T s = math::zero<T>();
        for(int i = 0; i < N * M; ++i)
            s += x(i) * x(i);
        return sqrt(s);
    }
};

/// Specialization of zero element for static matrices.
template <typename T, int N, int M>
struct zero_impl< static_matrix<T, N, M> >
{
    static static_matrix<T, N, M> get() {
        static_matrix<T, N, M> z;
        for(int i = 0; i < N * M; ++i)
            z(i) = math::zero<T>();
        return z;
    }
};

/// Specialization of zero element for static matrices.
template <typename T, int N, int M>
struct is_zero_impl< static_matrix<T, N, M> >
{
    static bool get(const static_matrix<T, N, M> &x) {
        for(int i = 0; i < N * M; ++i)
            if (!math::is_zero(x(i))) return false;
        return true;
    }
};

/// Specialization of identity for static matrices.
template <typename T, int N>
struct identity_impl< static_matrix<T, N, N> >
{
    static static_matrix<T, N, N> get() {
        static_matrix<T, N, N> I;
        for(int i = 0; i < N; ++i)
            for(int j = 0; j < N; ++j)
                I(i,j) = static_cast<T>(i == j);
        return I;
    }
};

/// Specialization of constant for static matrices.
template <typename T, int N, int M>
struct constant_impl< static_matrix<T, N, M> >
{
    static static_matrix<T, N, M> get(T c) {
        static_matrix<T, N, M> C;
        for(int i = 0; i < N * M; ++i)
            C(i) = c;
        return C;
    }
};

/// Specialization of inversion for static matrices.
template <typename T, int N>
struct inverse_impl< static_matrix<T, N, N> >
{
    static static_matrix<T, N, N> get(static_matrix<T, N, N> A) {
        // Perform LU-factorization of A in-place
        for(int k = 0; k < N; ++k) {
            T d = A(k,k);
            assert(!math::is_zero(d));
            for(int i = k+1; i < N; ++i) {
                A(i,k) /= d;
                for(int j = k+1; j < N; ++j)
                    A(i,j) -= A(i,k) * A(k,j);
            }
        }

        // Invert identity matrix in-place to get the solution.
        static_matrix<T, N, N> y;
        for(int k = 0; k < N; ++k) {
            // Lower triangular solve:
            for(int i = 0; i < N; ++i) {
                T b = static_cast<T>(i == k);
                for(int j = 0; j < i; ++j)
                    b -= A(i,j) * y(k,j);
                y(k,i) = b;
            }

            // Upper triangular solve:
            for(int i = N; i --> 0; ) {
                for(int j = i+1; j < N; ++j)
                    y(k,i) -= A(i,j) * y(k,j);
                y(k,i) /= A(i,i);
            }
        }

        return y;
    }
};


} // namespace math

namespace relaxation {
template <class Backend> struct spai1;
} //namespace relaxation

namespace backend {

template <class Backend>
struct relaxation_is_supported<
    Backend, relaxation::spai1,
    typename boost::enable_if<
        typename amgcl::is_static_matrix<
            typename Backend::value_type
            >::type
        >::type
    > : boost::false_type
{};

} // namespace backend
} // namespace amgcl

#endif
