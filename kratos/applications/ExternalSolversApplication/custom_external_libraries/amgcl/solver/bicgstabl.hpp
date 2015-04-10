#ifndef AMGCL_SOLVER_BICGSTABL_HPP
#define AMGCL_SOLVER_BICGSTABL_HPP

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
 * \file   amgcl/solver/bicgstabl.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  BiCGStab(L) iterative method.
 */

#include <boost/tuple/tuple.hpp>
#include <boost/multi_array.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace solver {

/// BiCGStab(L) iterative solver.
/**
 * \param Backend Backend for temporary structures allocation.
 * \ingroup solvers
 * \sa \cite Sleijpen1993
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class bicgstabl {
    public:
        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        /// Solver parameters.
        struct params {
            /// Order of the method.
            int L;

            /// Maximum number of iterations.
            size_t maxiter;

            /// Target residual error.
            value_type tol;

            params(int L = 2, size_t maxiter = 100, value_type tol = 1e-8)
                : L(L), maxiter(maxiter), tol(tol)
            {
                precondition(L > 0, "L in BiCGStab(L) should be >=1");
            }

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, L),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol)
            {}

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, L);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
            }
        };

        /// \copydoc amgcl::solver::cg::cg
        bicgstabl(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
                )
            : prm(prm), n(n),
              r0( Backend::create_vector(n, backend_prm) ),
              q ( Backend::create_vector(n, backend_prm) ),
              r(prm.L + 1), u(prm.L + 1),
              tau(boost::extents[prm.L][prm.L]),
              sigma(prm.L), gamma(prm.L), gamma1(prm.L), gamma2(prm.L),
              inner_product(inner_product)
        {
            for(int i = 0; i <= prm.L; ++i) {
                r[i] = Backend::create_vector(n, backend_prm);
                u[i] = Backend::create_vector(n, backend_prm);
            }
        }

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
            const int L = prm.L;

            backend::residual(rhs, A, x, *r0);

            value_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<value_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            value_type res_norm = norm(*r0);
            value_type eps      = prm.tol * norm_rhs;

            if(res_norm < eps)
                return boost::make_tuple(0, res_norm / norm_rhs);

            backend::copy(*r0, *r[0]);
            backend::clear( *u[0] );
            value_type rho0 = 1, alpha = 0, omega = 1;

            size_t iter = 0;

            for(; res_norm > eps && iter < prm.maxiter; iter += prm.L) {
                rho0 = -omega * rho0;

                // Bi-CG part
                for(int j = 0; j < L; ++j) {
                    precondition(rho0, "Zero rho in BiCGStab(L)");

                    double rho1 = inner_product(*r[j], *r0);
                    double beta = alpha * rho1 / rho0;
                    rho0 = rho1;

                    for(int i = 0; i <= j; ++i)
                        backend::axpby(1, *r[i], -beta, *u[i]);

                    P.apply(*u[j], *q);
                    backend::spmv(1, A, *q, 0, *u[j+1]);

                    alpha = inner_product(*u[j+1], *r0);

                    if (alpha == 0) break;

                    alpha = rho0 / alpha;

                    for(int i = 0; i <= j; ++i)
                        backend::axpby(-alpha, *u[i+1], 1, *r[i]);

                    P.apply(*r[j], *q);
                    backend::spmv(1, A, *q, 0, *r[j+1]);
                    backend::axpby(alpha, *u[0], 1, x);
                }

                // MR part
                for(int j = 0; j < L; ++j) {
                    for(int i = 0; i < j; ++i) {
                        tau[i][j] = inner_product(*r[j+1], *r[i+1]) / sigma[i];
                        backend::axpby(-tau[i][j], *r[i+1], 1, *r[j+1]);
                    }
                    sigma[j] = inner_product(*r[j+1], *r[j+1]);
                    gamma1[j] = inner_product(*r[0], *r[j+1]) / sigma[j];
                }

                omega = gamma[L-1] = gamma1[L-1];
                for(int j = L-2; j >= 0; --j) {
                    gamma[j] = gamma1[j];
                    for(int i = j+1; i < L; ++i)
                        gamma[j] -= tau[j][i] * gamma[i];
                }

                for(int j = 0; j < L-1; ++j) {
                    gamma2[j] = gamma[j+1];
                    for(int i = j+1; i < L-1; ++i)
                        gamma2[j] += tau[j][i] * gamma[i+1];
                }

                // Update
                backend::axpby(gamma[0], *r[0], 1, x);
                backend::axpby(-gamma1[L-1], *r[L], 1, *r[0]);
                backend::axpby(-gamma[L-1], *u[L], 1, *u[0]);

                for(int j = 1; j < L; ++j) {
                    backend::axpby(-gamma[j-1], *u[j], 1, *u[0]);
                    backend::axpby(gamma2[j-1], *r[j], 1, x);
                    backend::axpby(-gamma1[j-1], *r[j], 1, *r[0]);
                }

                res_norm = norm(*r[0]);
            }

            P.apply(x, *q);
            backend::copy(*q, x);
            backend::residual(rhs, A, x, *r0);
            res_norm = norm(*r0);

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

        mutable boost::shared_ptr< vector > r0;
        mutable boost::shared_ptr< vector > q;

        mutable std::vector< boost::shared_ptr< vector > > r;
        mutable std::vector< boost::shared_ptr< vector > > u;

        mutable boost::multi_array<value_type, 2> tau;
        mutable std::vector<value_type> sigma;
        mutable std::vector<value_type> gamma, gamma1, gamma2;

        InnerProduct inner_product;

        template <class Vec>
        value_type norm(const Vec &x) const {
            return sqrt(inner_product(x, x));
        }
};

} // namespace solver
} // namespace amgcl


#endif
