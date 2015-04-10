#ifndef AMGCL_BACKEND_VEXCL_HPP
#define AMGCL_BACKEND_VEXCL_HPP

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
 * \file   amgcl/backend/vexcl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  VexCL backend.
 */

#include <iostream>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <vexcl/vexcl.hpp>

#include <amgcl/util.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/backend/detail/default_direct_solver.hpp>

namespace amgcl {
namespace backend {

/// VexCL backend
/**
 * This is a backend that uses types defined in the VexCL GPGPU library
 * (https://github.com/ddemidov/vexcl).
 *
 * \param real Value type.
 * \ingroup backends
 */
template <typename real>
struct vexcl {
    typedef real      value_type;
    typedef ptrdiff_t index_type;

    typedef vex::SpMat<value_type, index_type, index_type> matrix;
    typedef vex::vector<value_type>                        vector;
    typedef detail::default_direct_solver<vexcl>           direct_solver;

    struct provides_row_iterator : boost::false_type {};

    /// Backend parameters.
    struct params {
        /// Command queues that identify compute devices to use with VexCL.
        std::vector< vex::backend::command_queue > q;

        params() : q(vex::current_context().queue()) {}

        params(const boost::property_tree::ptree &p) {
            std::vector<vex::backend::command_queue> *ptr = 0;
            ptr = p.get("q", ptr);
            q = ptr ? *ptr : vex::current_context().queue();
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            p.put(path + "q", &q);
        }
    };

    static std::string name() { return "vexcl"; }

    /// Copy matrix from builtin backend.
    static boost::shared_ptr<matrix>
    copy_matrix(boost::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        precondition(!prm.q.empty(), "Empty VexCL context!");

        return boost::make_shared<matrix>(prm.q, rows(*A), cols(*A),
                A->ptr.data(), A->col.data(), A->val.data()
                );
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(typename builtin<real>::vector const &x, const params &prm)
    {
        precondition(!prm.q.empty(), "Empty VexCL context!");

        return boost::make_shared<vector>(prm.q, x);
    }

    /// Copy vector from builtin backend.
    static boost::shared_ptr<vector>
    copy_vector(boost::shared_ptr< typename builtin<real>::vector > x, const params &prm)
    {
        return copy_vector(*x, prm);
    }

    /// Create vector of the specified size.
    static boost::shared_ptr<vector>
    create_vector(size_t size, const params &prm)
    {
        precondition(!prm.q.empty(), "Empty VexCL context!");

        return boost::make_shared<vector>(prm.q, size);
    }

    struct gather {
        mutable vex::gather<value_type> G;

        gather(size_t src_size, const std::vector<ptrdiff_t> &I, const params &prm)
            : G(prm.q, src_size, std::vector<size_t>(I.begin(), I.end())) { }

        void operator()(const vector &vec, std::vector<value_type> &vals) const {
            G(vec, vals);
        }
    };


    /// Create direct solver for coarse level
    static boost::shared_ptr<direct_solver>
    create_solver(boost::shared_ptr< typename builtin<real>::matrix > A, const params &prm)
    {
        return boost::make_shared<direct_solver>(A, prm);
    }
};

//---------------------------------------------------------------------------
// Backend interface implementation
//---------------------------------------------------------------------------
template < typename V, typename C, typename P >
struct rows_impl< vex::SpMat<V, C, P> > {
    static size_t get(const vex::SpMat<V, C, P> &A) {
        return A.rows();
    }
};

template < typename V, typename C, typename P >
struct cols_impl< vex::SpMat<V, C, P> > {
    static size_t get(const vex::SpMat<V, C, P> &A) {
        return A.cols();
    }
};

template < typename V, typename C, typename P >
struct nonzeros_impl< vex::SpMat<V, C, P> > {
    static size_t get(const vex::SpMat<V, C, P> &A) {
        return A.nonzeros();
    }
};

template < typename V, typename C, typename P >
struct spmv_impl<
    vex::SpMat<V, C, P>,
    vex::vector<V>,
    vex::vector<V>
    >
{
    typedef vex::SpMat<V, C, P> matrix;
    typedef vex::vector<V>      vector;

    static void apply(V alpha, const matrix &A, const vector &x,
            V beta, vector &y)
    {
        if (beta)
            y = alpha * (A * x) + beta * y;
        else
            y = alpha * (A * x);
    }
};

template < typename V, typename C, typename P >
struct residual_impl<
    vex::SpMat<V, C, P>,
    vex::vector<V>,
    vex::vector<V>,
    vex::vector<V>
    >
{
    typedef vex::SpMat<V, C, P> matrix;
    typedef vex::vector<V>      vector;

    static void apply(const vector &rhs, const matrix &A, const vector &x,
            vector &r)
    {
        r = rhs - A * x;
    }
};

template < typename V >
struct clear_impl< vex::vector<V> >
{
    static void apply(vex::vector<V> &x)
    {
        x = 0;
    }
};

template < typename V >
struct copy_impl<
    vex::vector<V>,
    vex::vector<V>
    >
{
    static void apply(const vex::vector<V> &x, vex::vector<V> &y)
    {
        y = x;
    }
};

template < typename V >
struct copy_to_backend_impl<
    vex::vector<V>
    >
{
    static void apply(const std::vector<V> &data, vex::vector<V> &x)
    {
        vex::copy(data, x);
    }
};

template < typename V >
struct inner_product_impl<
    vex::vector<V>,
    vex::vector<V>
    >
{
    static V get(const vex::vector<V> &x, const vex::vector<V> &y)
    {
        vex::Reductor<V, vex::SUM_Kahan> sum( x.queue_list() );
        return sum(x * y);
    }
};

template < typename V >
struct axpby_impl<
    vex::vector<V>,
    vex::vector<V>
    > {
    static void apply(V a, const vex::vector<V> &x, V b, vex::vector<V> &y)
    {
        if (b)
            y = a * x + b * y;
        else
            y = a * x;
    }
};

template < typename V >
struct axpbypcz_impl<
    vex::vector<V>,
    vex::vector<V>,
    vex::vector<V>
    >
{
    static void apply(
            V a, const vex::vector<V> &x,
            V b, const vex::vector<V> &y,
            V c,       vex::vector<V> &z
            )
    {
        if (c)
            z = a * x + b * y + c * z;
        else
            z = a * x + b * y;
    }
};

template < typename V >
struct vmul_impl<
    vex::vector<V>,
    vex::vector<V>,
    vex::vector<V>
    >
{
    static void apply(V a, const vex::vector<V> &x, const vex::vector<V> &y,
            V b, vex::vector<V> &z)
    {
        if (b)
            z = a * x * y + b * z;
        else
            z = a * x * y;
    }
};

} // namespace backend
} // namespace amgcl

#endif
