#ifndef AMGCL_BACKEND_EIGEN_HPP
#define AMGCL_BACKEND_EIGEN_HPP

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
 * \file   amgcl/backend/eigen.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse matrix in CRS format.
 */

#include <memory>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/solver/skyline_lu.hpp>

namespace amgcl {
namespace backend {

/// Eigen backend.
/**
 * This is a backend that uses types defined in the Eigen library
 * (http://eigen.tuxfamily.org).
 *
 * \param real Value type.
 * \ingroup backends
 */
template <typename real>
struct eigen {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    typedef
        Eigen::MappedSparseMatrix<value_type, Eigen::RowMajor, index_type>
        matrix;

    typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> vector;
    typedef Eigen::Matrix<value_type, Eigen::Dynamic, 1> matrix_diagonal;

    typedef solver::skyline_lu<real> direct_solver;

    struct provides_row_iterator : std::true_type {};

    /// Backend parameters.
    typedef amgcl::detail::empty_params params;

    static std::string name() { return "eigen"; }

    /// Copy matrix from builtin backend.
    static std::shared_ptr<matrix>
    copy_matrix(std::shared_ptr< typename builtin<real>::matrix > A, const params&)
    {
        const typename builtin<real>::matrix &a = *A;

        return std::shared_ptr<matrix>(
                new matrix(
                    rows(*A), cols(*A), nonzeros(*A),
                    const_cast<index_type*>(a.ptr),
                    const_cast<index_type*>(a.col),
                    const_cast<value_type*>(a.val)
                    ),
                hold_host(A)
                );
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params&)
    {
        return std::make_shared<vector>(
                Eigen::Map<const vector>(x.data(), x.size())
                );
    }

    /// Copy vector from builtin backend.
    static std::shared_ptr<vector>
    copy_vector(std::shared_ptr< typename builtin<real>::vector > x, const params &prm)
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static std::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return std::make_shared<vector>(size);
    }

    /// Create direct solver for coarse level
    static std::shared_ptr<direct_solver>
    create_solver(std::shared_ptr< typename builtin<real>::matrix > A, const params&)
    {
        return std::make_shared<direct_solver>(*A);
    }

    private:
        struct hold_host {
            typedef std::shared_ptr< crs<real, ptrdiff_t, ptrdiff_t> > host_matrix;
            host_matrix host;

            hold_host( host_matrix host ) : host(host) {}

            void operator()(matrix *ptr) const {
                delete ptr;
            }
        };

};

template < class Alpha, class M, class V1, class Beta, class V2 >
struct spmv_impl<
    Alpha, M, V1, Beta, V2,
    typename std::enable_if<
        is_eigen_sparse_matrix<M>::value &&
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value
        >::type
    >
{
    static void apply(Alpha alpha, const M &A, const V1 &x, Beta beta, V2 &y)
    {
        if (!math::is_zero(beta))
            y = alpha * A * x + beta * y;
        else
            y = alpha * A * x;
    }
};

template < class M, class V1, class V2, class V3 >
struct residual_impl<
    M, V1, V2, V3,
    typename std::enable_if<
        is_eigen_sparse_matrix<M>::value &&
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value &&
        is_eigen_type<V3>::value
        >::type
    >
{
    static void apply(const V1 &rhs, const M &A, const V2 &x, V3 &r)
    {
        r = rhs - A * x;
    }
};

template < typename V >
struct clear_impl<
    V,
    typename std::enable_if<is_eigen_type<V>::value>::type
    >
{
    static void apply(V &x)
    {
        x.setZero();
    }
};

template < class V1, class V2 >
struct inner_product_impl<
    V1, V2,
    typename std::enable_if<
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value
        >::type
    >
{
    typedef typename value_type<V1>::type real;
    static real get(const V1 &x, const V2 &y)
    {
        return x.dot(y);
    }
};

template < class A, class V1, class B, class V2 >
struct axpby_impl<
    A, V1, B, V2,
    typename std::enable_if<
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value
        >::type
    >
{
    static void apply(A a, const V1 &x, B b, V2 &y)
    {
        if (!math::is_zero(b))
            y = a * x + b * y;
        else
            y = a * x;
    }
};

template < class A, class V1, class B, class V2, class C, class V3 >
struct axpbypcz_impl<
    A, V1, B, V2, C, V3,
    typename std::enable_if<
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value &&
        is_eigen_type<V3>::value
        >::type
    >
{
    typedef typename value_type<V1>::type real;

    static void apply(
            real a, const V1 &x,
            real b, const V2 &y,
            real c,       V3 &z
            )
    {
        if (!math::is_zero(c))
            z = a * x + b * y + c * z;
        else
            z = a * x + b * y;
    }
};

template < class Alpha, class V1, class V2, class Beta, class V3 >
struct vmul_impl<
    Alpha, V1, V2, Beta, V3,
    typename std::enable_if<
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value &&
        is_eigen_type<V3>::value
        >::type
    >
{
    static void apply(Alpha a, const V1 &x, const V2 &y, Beta b, V3 &z)
    {
        if (!math::is_zero(b))
            z.array() = a * x.array() * y.array() + b * z.array();
        else
            z.array() = a * x.array() * y.array();
    }
};

template < class V1, class V2 >
struct copy_impl<
    V1, V2,
    typename std::enable_if<
        is_eigen_type<V1>::value &&
        is_eigen_type<V2>::value
        >::type
    >
{
    static void apply(const V1 &x, V2 &y)
    {
        y = x;
    }
};

} // namespace backend
} // namespace amgcl

#endif
