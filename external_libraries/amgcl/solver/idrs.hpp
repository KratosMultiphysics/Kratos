#ifndef AMGCL_SOLVER_IDRS_HPP
#define AMGCL_SOLVER_IDRS_HPP

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
\file   idrs.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  IDR(s) method.

The code is ported from Matlab code published at
http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html.

This is a very stable and efficient IDR(s) variant (implemented in the MATLAB
code idrs.m given above) as described in: Martin B. van Gijzen and Peter
Sonneveld, Algorithm 913: An Elegant IDR(s) Variant that Efficiently Exploits
Bi-orthogonality Properties. ACM Transactions on Mathematical Software, Vol.
38, No. 1, pp. 5:1-5:19, 2011 (copyright ACM).
*/

#include <vector>
#include <algorithm>

#include <tuple>
#include <random>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

#ifdef MPI_VERSION
#  include <amgcl/mpi/util.hpp>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace amgcl {
namespace solver {

/// IDR(s) method (Induced Dimension Reduction)
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class idrs {
    public:
        typedef Backend backend_type;

        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        typedef typename math::scalar_of<value_type>::type scalar_type;
        typedef typename math::rhs_of<value_type>::type rhs_type;

        typedef typename math::inner_product_impl<
            typename math::rhs_of<value_type>::type
            >::return_type coef_type;

        /// Solver parameters.
        struct params {
            /// Dimension of the shadow space in IDR(s).
            unsigned s;

            /// Computation of omega.
            /**
             * If omega = 0: a standard minimum residual step is performed
             * If omega > 0: omega is increased if
             * the cosine of the angle between Ar and r < omega
             * Default: omega = 0.7;
             */
            scalar_type omega;

            /// Specifies if residual smoothing must be applied.
            bool smoothing;

            /// Residual replacement.
            /**
             * Determines the residual replacement strategy.
             *    If |r| > 1E3 |b| TOL/EPS) (EPS is the machine precision)
             *    the recursively computed residual is replaced by the true residual
             *    once |r| < |b| (to reduce the effect of large intermediate residuals
             *    on the final accuracy).
             * Default: No residual replacement.
             */
            bool replacement;

            /// Maximum number of iterations.
            unsigned maxiter;

            /// Target relative residual error.
            scalar_type tol;

            /// Target absolute residual error.
            scalar_type abstol;

            params()
                : s(4), omega(0.7), smoothing(false),
                  replacement(false), maxiter(100), tol(1e-8),
                  abstol(std::numeric_limits<scalar_type>::min())
            { }

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, s),
                  AMGCL_PARAMS_IMPORT_VALUE(p, omega),
                  AMGCL_PARAMS_IMPORT_VALUE(p, smoothing),
                  AMGCL_PARAMS_IMPORT_VALUE(p, replacement),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol),
                  AMGCL_PARAMS_IMPORT_VALUE(p, abstol)
            {
                check_params(p, {"s", "omega", "smoothing", "replacement", "maxiter", "tol", "abstol"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, s);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, omega);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, smoothing);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, replacement);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, abstol);
            }
#endif
        } prm;

        /// Preallocates necessary data structures for the system of size \p n.
        idrs(
                size_t n,
                const params &prm = params(),
                const backend_params &bprm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
             )
            : prm(prm), n(n), inner_product(inner_product),
              M(prm.s, prm.s),
              f(prm.s), c(prm.s),
              r(Backend::create_vector(n, bprm)),
              v(Backend::create_vector(n, bprm)),
              t(Backend::create_vector(n, bprm))
        {
            static const scalar_type one = math::identity<scalar_type>();
            static const scalar_type zero = math::zero<scalar_type>();

            if (prm.smoothing) {
              x_s = Backend::create_vector(n, bprm);
              r_s = Backend::create_vector(n, bprm);
            }

            G.reserve(prm.s);
            U.reserve(prm.s);
            for(unsigned i = 0; i < prm.s; ++i) {
                G.push_back(Backend::create_vector(n, bprm));
                U.push_back(Backend::create_vector(n, bprm));
            }

            // Initialize P.
            P.reserve(prm.s);
            {
                std::vector<rhs_type> p(n);

#ifdef MPI_VERSION
                int pid = amgcl::mpi::communicator(MPI_COMM_WORLD).rank;
#else
                int pid = 0;
#endif

#pragma omp parallel
                {
#ifdef _OPENMP
                    int tid = omp_get_thread_num();
                    int nt = omp_get_max_threads();
#else
                    int tid = 0;
                    int nt = 1;
#endif

                    std::mt19937 rng(pid * nt + tid);
                    std::uniform_real_distribution<scalar_type> rnd(-1, 1);

                    for(unsigned j = 0; j < prm.s; ++j) {
#pragma omp for
                        for(ptrdiff_t i = 0; i < static_cast<ptrdiff_t>(n); ++i)
                            p[i] = math::constant<rhs_type>(rnd(rng));

#pragma omp single
                        {
                            P.push_back(Backend::copy_vector(p, bprm));
                        }
                    }
                }

                for(unsigned j = 0; j < prm.s; ++j) {
                    for(unsigned k = 0; k < j; ++k) {
                        coef_type alpha = inner_product(*P[k], *P[j]);
                        backend::axpby(-alpha, *P[k], one, *P[j]);
                    }
                    scalar_type norm_pj = norm(*P[j]);
                    backend::axpby(math::inverse(norm_pj), *P[j], zero, *P[j]);
                }
            }
        }

        /* Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``std::tuple``. The solution vector
         * \p x provides initial approximation in input and holds the computed
         * solution on output.
         *
         * The system matrix may differ from the matrix used during
         * initialization. This may be used for the solution of non-stationary
         * problems with slowly changing coefficients. There is a strong chance
         * that a preconditioner built for a time step will act as a reasonably
         * good preconditioner for several subsequent time steps [DeSh12]_.
         */
        template <class Matrix, class Precond, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                Matrix  const &A,
                Precond const &Prec,
                Vec1    const &rhs,
                Vec2          &x
                ) const
        {
            static const scalar_type one = math::identity<scalar_type>();
            static const scalar_type zero = math::zero<scalar_type>();

            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(1)) {
                backend::clear(x);
                return std::make_tuple(0, norm_rhs);
            }

            scalar_type eps = std::max(prm.tol * norm_rhs, prm.abstol);

            // Compute initial residual:
            backend::residual(rhs, A, x, *r);

            scalar_type res_norm = norm(*r);
            if (res_norm <= eps) {
                // Initial guess is a good enough solution.
                return std::make_tuple(0, res_norm / norm_rhs);
            }

            if (prm.smoothing) {
                backend::copy( x, *x_s);
                backend::copy(*r, *r_s);
            }

            // Initialization.
            coef_type om = math::identity<coef_type>();

            for(unsigned i = 0; i < prm.s; ++i) {
                backend::clear(*G[i]);
                backend::clear(*U[i]);

                for(unsigned j = 0; j < prm.s; ++j)
                    M(i, j) = (i == j);
            }

            scalar_type eps_replace = norm_rhs / (
                    // Number close to machine precision:
                    1e3 * std::numeric_limits<scalar_type>::epsilon());

            // Main iteration loop, build G-spaces:
            size_t iter = 0;
            bool trueres = false;
            while(iter < prm.maxiter && res_norm > eps) {
                // New righ-hand size for small system:
                for(unsigned i = 0; i < prm.s; ++i)
                    f[i] = inner_product(*r, *P[i]);

                for(unsigned k = 0; k < prm.s; ++k) {
                    // Compute new v
                    backend::copy(*r, *v);

                    // Solve small system (Note: M is lower triangular)
                    // and make v orthogonal to P:
                    for(unsigned i = k; i < prm.s; ++i) {
                        c[i] = f[i];
                        for(unsigned j = k; j < i; ++j)
                            c[i] -= M(i, j) * c[j];
                        c[i] = math::inverse(M(i, i)) * c[i];

                        backend::axpby(-c[i], *G[i], one, *v);
                    }

                    Prec.apply(*v, *t);

                    // Compute new U[k]
                    backend::axpby(om, *t, c[k], *U[k]);
                    for(unsigned i = k+1; i < prm.s; ++i)
                        backend::axpby(c[i], *U[i], one, *U[k]);

                    // Compute new G[k], G[k] is in space G_j
                    backend::spmv(one, A, *U[k], zero, *G[k]);

                    // Bi-Orthogonalise the new basis vectors:
                    for(unsigned i = 0; i < k; ++i) {
                        coef_type alpha = inner_product(*G[k], *P[i]) / M(i, i);

                        backend::axpby(-alpha, *G[i], one, *G[k]);
                        backend::axpby(-alpha, *U[i], one, *U[k]);
                    }

                    // New column of M = P'*G  (first k-1 entries are zero)
                    for(unsigned i = k; i < prm.s; ++i)
                        M(i, k) = inner_product(*G[k], *P[i]);

                    precondition(!math::is_zero(M(k, k)), "IDR(s) breakdown: zero M[k,k]");

                    // Make r orthogonal to q_i, i = [0..k)
                    coef_type beta = math::inverse(M(k, k)) * f[k];
                    backend::axpby(-beta, *G[k], one, *r);
                    backend::axpby( beta, *U[k], one,  x);

                    res_norm = norm(*r);

                    if (prm.replacement && res_norm > eps_replace)
                        trueres = true;

                    // Smoothing
                    if (prm.smoothing) {
                        backend::axpbypcz(one, *r_s, -one, *r, zero, *t);
                        coef_type gamma = inner_product(*t, *r_s) / inner_product(*t, *t);
                        backend::axpby(-gamma, *t, one, *r_s);
                        backend::axpbypcz(-gamma, *x_s, gamma, x, one, *x_s);
                        res_norm = norm(*r_s);
                    }

                    if (res_norm <= eps || ++iter >= prm.maxiter) break;

                    // New f = P'*r (first k  components are zero)
                    for(unsigned i = k + 1; i < prm.s; ++i)
                        f[i] -= beta * M(i, k);
                }

                if (res_norm <= eps || iter >= prm.maxiter) break;

                // Now we have sufficient vectors in G_j to compute residual in G_j+1
                // Note: r is already perpendicular to P so v = r

                Prec.apply(*r, *v);
                backend::spmv(one, A, *v, zero, *t);

                // Computation of a new omega
                om = omega(*t, *r);
                precondition(!math::is_zero(om), "IDR(s) breakdown: zero omega");

                backend::axpby(-om, *t, one, *r);
                backend::axpby( om, *v, one,  x);

                res_norm = norm(*r);
                if (prm.replacement && res_norm > eps_replace)
                    trueres = true;

                // Residual replacement?
                if (trueres && res_norm < norm_rhs) {
                    trueres = 0;
                    backend::residual(rhs, A, x, *r);
                }

                // Smoothing.
                if (prm.smoothing) {
                    backend::axpbypcz(one, *r_s, -one, *r, zero, *t);
                    coef_type gamma = inner_product(*t, *r_s) / inner_product(*t, *t);
                    backend::axpby(-gamma, *t, one, *r_s);
                    backend::axpbypcz(-gamma, *x_s, gamma, x, one, *x_s);
                    res_norm = norm(*r_s);
                }

                ++iter;
            }

            if (prm.smoothing)
                backend::copy(*x_s, x);

            return std::make_tuple(iter, res_norm / norm_rhs);
        }

        /* Computes the solution for the given right-hand side \p rhs. The
         * system matrix is the same that was used for the setup of the
         * preconditioner \p P.  Returns the number of iterations made and the
         * achieved residual as a ``std::tuple``. The solution vector \p x
         * provides initial approximation in input and holds the computed
         * solution on output.
         */
        template <class Precond, class Vec1, class Vec2>
        std::tuple<size_t, scalar_type> operator()(
                Precond const &P,
                Vec1    const &rhs,
                Vec2          &x
                ) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }

        size_t bytes() const {
            size_t b = 0;

            b += M.size() * sizeof(coef_type);

            b += backend::bytes(f);
            b += backend::bytes(c);

            b += backend::bytes(*r);
            b += backend::bytes(*v);
            b += backend::bytes(*t);

            if (x_s) b += backend::bytes(*x_s);
            if (r_s) b += backend::bytes(*r_s);

            for(const auto &v : P) b += backend::bytes(*v);
            for(const auto &v : G) b += backend::bytes(*v);
            for(const auto &v : U) b += backend::bytes(*v);

            return b;
        }

        friend std::ostream& operator<<(std::ostream &os, const idrs &s) {
            return os
                << "Type:             IDR(" << s.prm.s << ")"
                << "\nUnknowns:         " << s.n
                << "\nMemory footprint: " << human_readable_memory(s.bytes())
                << std::endl;
        }

    private:
        size_t n;

        InnerProduct inner_product;

        mutable multi_array<coef_type,2> M;
        mutable std::vector<coef_type> f, c;

        std::shared_ptr<vector> r, v, t;
        std::shared_ptr<vector> x_s;
        std::shared_ptr<vector> r_s;

        std::vector< std::shared_ptr<vector> > P, G, U;


        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return std::abs(sqrt(inner_product(x, x)));
        }

        template <class Vector1, class Vector2>
        coef_type omega(const Vector1 &t, const Vector2 &s) const {
            scalar_type norm_t = norm(t);
            scalar_type norm_s = norm(s);

            coef_type   ts  = inner_product(t, s);
            scalar_type rho = math::norm(ts / (norm_t * norm_s));
            coef_type   om  = ts / (norm_t * norm_t);

            if (rho < prm.omega)
                om *= prm.omega/rho;

            return om;
        }
};

} // namespace solver
} // namespace amgcl

#endif
