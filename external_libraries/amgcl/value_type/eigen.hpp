#ifndef AMGCL_VALUE_TYPE_EIGEN_HPP
#define AMGCL_VALUE_TYPE_EIGEN_HPP

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
 * \file   amgcl/value_type/eigen.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Enable statically sized eigen matrices as value types.
 */

#include <Eigen/Dense>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace backend {

/// Enable Eigen matrix as a value-type.
template <typename T, int N, int M>
struct is_builtin_vector< std::vector<Eigen::Matrix<T, N, M> > > : std::true_type {};

} // namespace backend

namespace math {

/// Scalar type of a non-scalar type.
template <class T, int N, int M>
struct scalar_of< Eigen::Matrix<T, N, M> > {
    typedef typename math::scalar_of<T>::type type;
};

/// RHS type corresponding to a non-scalar type.
template <class T, int N>
struct rhs_of< Eigen::Matrix<T, N, N> > {
    typedef Eigen::Matrix<T, N, 1> type;
};

/// Whether the value type is a statically sized matrix.
template <class T, int N, int M>
struct is_static_matrix< Eigen::Matrix<T, N, M> > : std::true_type {};

/// Number of rows for statically sized matrix types.
template <class T, int N, int M>
struct static_rows< Eigen::Matrix<T, N, M> > : std::integral_constant<int, N> {};

/// Number of columns for statically sized matrix types.
template <class T, int N, int M>
struct static_cols< Eigen::Matrix<T, N, M> > : std::integral_constant<int, M> {};

/// Specialization of conjugate transpose for eigen matrices.
template <typename T, int N, int M>
struct adjoint_impl< Eigen::Matrix<T, N, M> >
{
    typedef typename Eigen::Matrix<T, N, M>::AdjointReturnType return_type;

    static return_type get(const Eigen::Matrix<T, N, M> &x) {
        return x.adjoint();
    }
};

/// Inner-product result of two Eigen vectors.
template <class T, int N>
struct inner_product_impl< Eigen::Matrix<T, N, 1> >
{
    typedef T return_type;
    static T get(const Eigen::Matrix<T, N, 1> &x, const Eigen::Matrix<T, N, 1> &y) {
        return x.adjoint() * y;
    }
};

/// Inner-product result of two Eigen matrices.
template <class T, int N, int M>
struct inner_product_impl< Eigen::Matrix<T, N, M> >
{
    typedef Eigen::Matrix<T, M, M> return_type;

    static return_type get(const Eigen::Matrix<T, N, M> &x, const Eigen::Matrix<T, N, M> &y) {
        return x.adjoint() * y;
    }
};

/// Specialization of element norm for eigen matrices.
template <typename T, int N, int M>
struct norm_impl< Eigen::Matrix<T, N, M> >
{
    static typename math::scalar_of<T>::type get(const Eigen::Matrix<T, N, M> &x) {
        return x.norm();
    }
};

/// Specialization of zero element for eigen matrices.
template <typename T, int N, int M>
struct zero_impl< Eigen::Matrix<T, N, M> >
{
    static Eigen::Matrix<T, N, M> get() {
        return Eigen::Matrix<T, N, M>::Zero();
    }
};

/// Specialization of zero element for eigen matrices.
template <typename T, int N, int M>
struct is_zero_impl< Eigen::Matrix<T, N, M> >
{
    static bool get(const Eigen::Matrix<T, N, M> &x) {
        return x.isZero();
    }
};

/// Specialization of identity for eigen matrices.
template <typename T, int N>
struct identity_impl< Eigen::Matrix<T, N, N> >
{
    static Eigen::Matrix<T, N, N> get() {
        return Eigen::Matrix<T, N, N>::Identity();
    }
};

/// Specialization of constant for eigen matrices.
template <typename T, int N, int M>
struct constant_impl< Eigen::Matrix<T, N, M> >
{
    static Eigen::Matrix<T, N, M> get(T c) {
        return Eigen::Matrix<T, N, M>::Constant(c);
    }
};

/// Specialization of inversion for eigen matrices.
template <typename T, int N>
struct inverse_impl< Eigen::Matrix<T, N, N> >
{
    static Eigen::Matrix<T, N, N> get(const Eigen::Matrix<T, N, N> &x) {
        return x.inverse();
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
    typename std::enable_if<
        std::is_base_of<
            Eigen::MatrixBase<typename Backend::value_type>,
            typename Backend::value_type
            >::value
        >::type
    > : std::false_type
{};

} // namespace backend
} // namespace amgcl

namespace Eigen {

template <class A, class B>
bool operator<(const MatrixBase<A> &a, const MatrixBase<B> &b) {
    return a.trace() < b.trace();
}

} // namespace Eigen

#endif
