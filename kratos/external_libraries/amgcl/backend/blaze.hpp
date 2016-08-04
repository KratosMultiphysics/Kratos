#ifndef AMGCL_BACKEND_BLAZE_HPP
#define AMGCL_BACKEND_BLAZE_HPP

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
 * \file   amgcl/backend/viennacl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Blaze backend.
 *
 * Uses Blaze (https://code.google.com/p/blaze-lib) types and operations.
 */

#include <blaze/Math.h>

#include <amgcl/backend/builtin.hpp>
#include <amgcl/solver/skyline_lu.hpp>

namespace amgcl {
namespace backend {

/// Blaze backend
/**
 * This is a backend that uses types defined in the Blaze library
 * (https://code.google.com/p/blaze-lib).
 *
 * \param real Value type.
 * \ingroup backends
 */
template <class real>
struct blaze {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    struct provides_row_iterator : boost::false_type {};

    typedef ::blaze::CompressedMatrix<real> matrix;
    typedef ::blaze::DynamicVector<real>    vector;
    typedef ::blaze::DynamicVector<real>    matrix_diagonal;
    typedef solver::skyline_lu<real>        direct_solver;

    /// Backend parameters.
    typedef amgcl::detail::empty_params params;

    static std::string name() { return "blaze"; }

    /// Copy matrix from builtin backend.
    static boost::shared_ptr<matrix>
    copy_matrix(
            boost::shared_ptr< typename builtin<real>::matrix > A,
            const params&
            )
    {
        typedef
            typename row_iterator<typename builtin<real>::matrix>::type
            row_iterator;

        const size_t n = rows(*A);
        const size_t m = cols(*A);

        boost::shared_ptr<matrix> B = boost::make_shared<matrix>(n, m);

        B->reserve(nonzeros(*A));
        for(size_t i = 0; i < n; ++i) {
            for(row_iterator a = A->row_begin(i); a; ++a) {
                B->append(i, a.col(), a.value());
            }
            B->finalize(i);
        }

        return B;
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params&)
    {
        boost::shared_ptr<vector> v = boost::make_shared<vector>(
                x.size(), &x[0]);
        return v;
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(
            boost::shared_ptr< typename builtin<real>::vector > x,
            const params &prm
            )
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

};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template < typename V >
struct value_type < ::blaze::CompressedMatrix<V> > {
    typedef V type;
};

template < typename V >
struct value_type < ::blaze::DynamicVector<V> > {
    typedef V type;
};

template < typename V >
struct rows_impl< ::blaze::CompressedMatrix<V> > {
    typedef ::blaze::CompressedMatrix<V> matrix;

    static size_t get(const matrix &A) {
        return A.rows();
    }
};

template < typename V >
struct cols_impl< ::blaze::CompressedMatrix<V> > {
    typedef ::blaze::CompressedMatrix<V> matrix;

    static size_t get(const matrix &A) {
        return A.columns();
    }
};

template < typename V >
struct nonzeros_impl< ::blaze::CompressedMatrix<V> > {
    typedef ::blaze::CompressedMatrix<V> matrix;

    static size_t get(const matrix &A) {
        return A.nonZeros();
    }
};

template < class A, class B, typename V >
struct spmv_impl<
    A, ::blaze::CompressedMatrix<V>, ::blaze::DynamicVector<V>,
    B, ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::CompressedMatrix<V> matrix;
    typedef ::blaze::DynamicVector<V>    vector;

    static void apply(A alpha, const matrix &K, const vector &x, B beta, vector &y)
    {
        if (!math::is_zero(beta))
            y = alpha * (K * x) + beta * y;
        else
            y = alpha * (K * x);
    }
};

template < typename V >
struct residual_impl<
    ::blaze::CompressedMatrix<V>,
    ::blaze::DynamicVector<V>,
    ::blaze::DynamicVector<V>,
    ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::CompressedMatrix<V> matrix;
    typedef ::blaze::DynamicVector<V>    vector;

    static void apply(const vector &rhs, const matrix &A, const vector &x,
            vector &r)
    {
        r = rhs - A * x;
    }
};

template < typename V >
struct clear_impl< ::blaze::DynamicVector<V> >
{
    typedef ::blaze::DynamicVector<V> vector;

    static void apply(vector &x)
    {
        x = 0;
    }
};

template < typename V >
struct copy_impl<
    ::blaze::DynamicVector<V>,
    ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::DynamicVector<V> vector;

    static void apply(const vector &x, vector &y)
    {
        y = x;
    }
};

template < typename V >
struct inner_product_impl<
    ::blaze::DynamicVector<V>,
    ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::DynamicVector<V> vector;

    static V get(const vector &x, const vector &y)
    {
        return (x, y);
    }
};

template < typename A, typename B, typename V >
struct axpby_impl<
    A, ::blaze::DynamicVector<V>,
    B, ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::DynamicVector<V> vector;

    static void apply(A a, const vector &x, B b, vector &y)
    {
        if (!math::is_zero(b))
            y = a * x + b * y;
        else
            y = a * x;
    }
};

template < typename A, typename B, typename C, typename V >
struct axpbypcz_impl<
    A, ::blaze::DynamicVector<V>,
    B, ::blaze::DynamicVector<V>,
    C, ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::DynamicVector<V> vector;

    static void apply(
            V a, const vector &x,
            V b, const vector &y,
            V c,       vector &z
            )
    {
        if (!math::is_zero(c))
            z = a * x + b * y + c * z;
        else
            z = a * x + b * y;
    }
};

template < typename A, typename B, typename V >
struct vmul_impl<
    A, ::blaze::DynamicVector<V>, ::blaze::DynamicVector<V>,
    B, ::blaze::DynamicVector<V>
    >
{
    typedef ::blaze::DynamicVector<V> vector;

    static void apply(A a, const vector &x, const vector &y, B b, vector &z)
    {
        if (!math::is_zero(b))
            z = a * x * y + b * z;
        else
            z = a * x * y;
    }
};

} // namespace backend
} // namespace amgcl

#endif
