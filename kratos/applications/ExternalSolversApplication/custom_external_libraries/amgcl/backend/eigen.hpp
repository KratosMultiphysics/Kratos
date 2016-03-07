#ifndef AMGCL_BACKEND_EIGEN_HPP
#define AMGCL_BACKEND_EIGEN_HPP

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
 * \file   amgcl/backend/eigen.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Sparse matrix in CRS format.
 */

#include <boost/typeof/typeof.hpp>
#include <boost/type_traits.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <Eigen/Sparse>
#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>
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

    struct provides_row_iterator : boost::true_type {};

    /// Backend parameters.
    struct params {
        params() {}
        params(const boost::property_tree::ptree&) {}
        void get(boost::property_tree::ptree&, const std::string&) const {}
    };

    static std::string name() { return "eigen"; }

    /// Copy matrix from builtin backend.
    static boost::shared_ptr<matrix>
    copy_matrix(boost::shared_ptr< typename builtin<real>::matrix > A, const params&)
    {
        const typename builtin<real>::matrix &a = *A;

        BOOST_AUTO(Aptr, a.ptr_data());
        BOOST_AUTO(Acol, a.col_data());
        BOOST_AUTO(Aval, a.val_data());

        return boost::shared_ptr<matrix>(
                new matrix(
                    rows(*A), cols(*A), nonzeros(*A),
                    const_cast<index_type*>(Aptr),
                    const_cast<index_type*>(Acol),
                    const_cast<value_type*>(Aval)
                    ),
                hold_host(A)
                );
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params&)
    {
        return boost::make_shared<vector>(
                Eigen::Map<const vector>(&x[0], x.size())
                );
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(boost::shared_ptr< typename builtin<real>::vector > x, const params &prm)
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static boost::shared_ptr<vector>
    create_vector(size_t size, const params&)
    {
        return boost::make_shared<vector>(size);
    }

    /// Create direct solver for coarse level
    static boost::shared_ptr<direct_solver>
    create_solver(boost::shared_ptr< typename builtin<real>::matrix > A, const params&)
    {
        return boost::make_shared<direct_solver>(*A);
    }

    private:
        struct hold_host {
            typedef boost::shared_ptr< crs<real, ptrdiff_t, ptrdiff_t> > host_matrix;
            host_matrix host;

            hold_host( host_matrix host ) : host(host) {}

            void operator()(matrix *ptr) const {
                delete ptr;
            }
        };

};

//---------------------------------------------------------------------------
// Backend interface specialization for Eigen types
//---------------------------------------------------------------------------
template <class T, class Enable = void>
struct is_eigen_sparse_matrix : boost::false_type {};

template <class T, class Enable = void>
struct is_eigen_type : boost::false_type {};

template <class T>
struct is_eigen_sparse_matrix<
    T,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename boost::is_arithmetic<typename T::Scalar>::type,
                typename boost::is_base_of<Eigen::SparseMatrixBase<T>, T>::type
            >::type
        >::type
    > : boost::true_type
{};

template <class T>
struct is_eigen_type<
    T,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename boost::is_arithmetic<typename T::Scalar>::type,
                typename boost::is_base_of<Eigen::EigenBase<T>, T>::type
            >::type
        >::type
    > : boost::true_type
{};

template <class T>
struct value_type<
    T,
    typename boost::enable_if<
        typename is_eigen_type<T>::type>::type
    >
{
    typedef typename T::Scalar type;
};

template <class T>
struct rows_impl<
    T,
    typename boost::enable_if<typename is_eigen_sparse_matrix<T>::type>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.rows();
    }
};

template <class T>
struct cols_impl<
    T,
    typename boost::enable_if<typename is_eigen_sparse_matrix<T>::type>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.cols();
    }
};

template <class T>
struct nonzeros_impl<
    T,
    typename boost::enable_if<typename is_eigen_type<T>::type>::type
    >
{
    static size_t get(const T &matrix) {
        return matrix.nonZeros();
    }
};

template <class T>
struct row_iterator <
    T,
    typename boost::enable_if<typename is_eigen_sparse_matrix<T>::type>::type
    >
{
    typedef typename T::InnerIterator type;
};

template <class T>
struct row_begin_impl <
    T,
    typename boost::enable_if<typename is_eigen_sparse_matrix<T>::type>::type
    >
{
    typedef typename row_iterator<T>::type iterator;
    static iterator get(const T &matrix, size_t row) {
        return iterator(matrix, row);
    }
};

template < class Alpha, class M, class V1, class Beta, class V2 >
struct spmv_impl<
    Alpha, M, V1, Beta, V2,
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_sparse_matrix<M>::type,
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type
            >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_sparse_matrix<M>::type,
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type,
                typename is_eigen_type<V3>::type
            >::type
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
    typename boost::enable_if< typename is_eigen_type<V>::type >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type
            >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type
            >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type,
                typename is_eigen_type<V3>::type
            >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type,
                typename is_eigen_type<V3>::type
            >::type
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
    typename boost::enable_if<
            typename boost::mpl::and_<
                typename is_eigen_type<V1>::type,
                typename is_eigen_type<V2>::type
            >::type
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
