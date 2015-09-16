#ifndef AMGCL_SOLVERS_CG_HPP
#define AMGCL_SOLVERS_CG_HPP

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
 * \file   amgcl/solver/cg.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Conjugate Gradient method.
 */

#include <boost/tuple/tuple.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

/// Iterative solvers
namespace solver {

/**
 * \defgroup solvers
 * \brief Iterative solvers
 *
 * AMGCL provides several iterative solvers, but it should be easy to use it as
 * a preconditioner with a user-provided solver.  Each solver in AMGCL is a
 * class template. Its single template parameter specifies the backend to use.
 * This allows to preallocate necessary resources at class construction.
 * Obviously, the solver backend has to coincide with the AMG backend.
 */


/// Conjugate Gradients iterative solver.
/**
 * \param Backend Backend for temporary structures allocation.
 * \ingroup solvers
 * \sa \cite Barrett1994
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class cg {
    public:
        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        /// Solver parameters.
        struct params {
            /// Maximum number of iterations.
            size_t maxiter;

            /// Target residual error.
            value_type tol;

            params(size_t maxiter = 100, value_type tol = 1e-8)
                : maxiter(maxiter), tol(tol)
            {}

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol)
            {}

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
            }
        };

        /// Preallocates necessary data structures
        /**
         * \param n           The system size.
         * \param prm         Solver parameters.
         * \param backend_prm Backend parameters.
         */
        cg(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
          ) : prm(prm), n(n),
              r(Backend::create_vector(n, backend_prm)),
              s(Backend::create_vector(n, backend_prm)),
              p(Backend::create_vector(n, backend_prm)),
              q(Backend::create_vector(n, backend_prm)),
              inner_product(inner_product)
        { }

        /// Solves the linear system for the given system matrix.
        /**
         * \param A   System matrix.
         * \param P   Preconditioner.
         * \param rhs Right-hand side.
         * \param x   Solution vector.
         *
         * The system matrix may differ from the matrix used for the AMG
         * preconditioner construction. This may be used for the solution of
         * non-stationary problems with slowly changing coefficients. There is
         * a strong chance that AMG built for one time step will act as a
         * reasonably good preconditioner for several subsequent time steps
         * \cite Demidov2012.
         */
        template <class Matrix, class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, value_type> operator()(
                Matrix  const &A,
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            backend::residual(rhs, A, x, *r);
            value_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<value_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            value_type eps  = prm.tol * norm_rhs;
            value_type eps2 = eps * eps;
            value_type rho1 = 2 * eps2, rho2 = 0;
            value_type res_norm = norm(*r);

            size_t iter = 0;
            for(; iter < prm.maxiter && fabs(rho1) > eps2; ++iter) {
                P.apply(*r, *s);

                rho2 = rho1;
                rho1 = inner_product(*r, *s);

                if (iter)
                    backend::axpby(1, *s, rho1 / rho2, *p);
                else
                    backend::copy(*s, *p);

                backend::spmv(1, A, *p, 0, *q);

                value_type alpha = rho1 / inner_product(*q, *p);

                backend::axpby( alpha, *p, 1,  x);
                backend::axpby(-alpha, *q, 1, *r);
            }

            backend::residual(rhs, A, x, *r);
            res_norm = norm(*r);

            return boost::make_tuple(iter, res_norm / norm_rhs);
        }

        /// Solves the linear system for the same matrix that was used for the AMG preconditioner construction.
        /**
         * \param P   AMG preconditioner.
         * \param rhs Right-hand side.
         * \param x   Solution vector.
         */
        template <class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, value_type> operator()(
                Precond const &P,
                Vec1    const &rhs,
#ifdef BOOST_NO_CXX11_RVALUE_REFERENCES
                Vec2          &x
#else
                Vec2          &&x
#endif
                ) const
        {
            return (*this)(P.top_matrix(), P, rhs, x);
        }

    public:
        params prm;

    private:
        size_t n;

        boost::shared_ptr<vector> r;
        boost::shared_ptr<vector> s;
        boost::shared_ptr<vector> p;
        boost::shared_ptr<vector> q;

        InnerProduct inner_product;

        template <class Vec>
        value_type norm(const Vec &x) const {
            return sqrt(inner_product(x, x));
        }
};

} // namespace solver
} // namespace amgcl


#endif
