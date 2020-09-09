#ifndef AMGCL_SOLVER_BICGSTABL_HPP
#define AMGCL_SOLVER_BICGSTABL_HPP

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
\file   amgcl/solver/bicgstabl.hpp
\author Denis Demidov <dennis.demidov@gmail.com>
\brief  BiCGStab(L) iterative method.

The code is ported from PETSC BCGSL [1] and is based on [2].

[1] http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPBCGSL.html
[2] Fokkema, Diederik R. Enhanced implementation of BiCGstab (l) for solving
    linear systems of equations. Universiteit Utrecht. Mathematisch Instituut,
    1996.

The original code came with the following license:

\verbatim
Copyright (c) 1991-2014, UChicago Argonne, LLC and the PETSc Development Team
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 \endverbatim
 */

#include <tuple>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/solver/precond_side.hpp>
#include <amgcl/detail/qr.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace solver {

/** BiCGStab(L) method.
 * \rst
 * Generalization of BiCGStab method [SlDi93]_.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class bicgstabl {
    public:
        typedef Backend backend_type;

        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        typedef typename math::inner_product_impl<
            typename math::rhs_of<value_type>::type
            >::return_type coef_type;


        /// Solver parameters.
        struct params {
            // Order of the method.
            int L;

            // Threshold used to decide when to refresh computed residuals.
            scalar_type delta;

            // Use a convex function of the MinRes and OR polynomials
            // after the BiCG step instead of default MinRes
            bool convex;

            // Preconditioning kind (left/right).
            preconditioner::side::type pside;

            // Maximum number of iterations.
            size_t maxiter;

            // Target relative residual error.
            scalar_type tol;

            // Target absolute residual error.
            scalar_type abstol;

            params()
                : L(2), delta(0), convex(true),
                  pside(preconditioner::side::right), maxiter(100), tol(1e-8),
                  abstol(std::numeric_limits<scalar_type>::min())
            {
            }

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, L),
                  AMGCL_PARAMS_IMPORT_VALUE(p, delta),
                  AMGCL_PARAMS_IMPORT_VALUE(p, convex),
                  AMGCL_PARAMS_IMPORT_VALUE(p, pside),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol),
                  AMGCL_PARAMS_IMPORT_VALUE(p, abstol)
            {
                check_params(p, {"L", "delta", "convex", "pside", "maxiter", "tol", "abstol"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, L);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, delta);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, convex);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, pside);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, abstol);
            }
#endif
        };

        /// Preallocates necessary data structures for the system of size \p n.
        bicgstabl(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
                )
            : prm(prm), n(n),
              Rt( Backend::create_vector(n, backend_prm) ),
              X ( Backend::create_vector(n, backend_prm) ),
              B ( Backend::create_vector(n, backend_prm) ),
              T ( Backend::create_vector(n, backend_prm) ),
              R(prm.L + 1), U(prm.L + 1),
              MZa(prm.L + 1, prm.L + 1),
              MZb(prm.L + 1, prm.L + 1),
              Y0(prm.L + 1), YL(prm.L + 1),
              inner_product(inner_product)
        {
            precondition(prm.L > 0, "L in BiCGStab(L) should be >=1");

            for(int i = 0; i <= prm.L; ++i) {
                R[i] = Backend::create_vector(n, backend_prm);
                U[i] = Backend::create_vector(n, backend_prm);
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
                const Matrix &A, const Precond &P, const Vec1 &rhs, Vec2 &&x) const
        {
            namespace side = preconditioner::side;

            static const coef_type one  = math::identity<coef_type>();
            static const coef_type zero = math::zero<coef_type>();

            const int L = prm.L;

            scalar_type norm_rhs = norm(rhs);

            // Check if there is a trivial solution
            if (norm_rhs < amgcl::detail::eps<scalar_type>(1)) {
                backend::clear(x);
                return std::make_tuple(0, norm_rhs);
            }

            if (prm.pside == side::left) {
                backend::residual(rhs, A, x, *T);
                P.apply(*T, *B);
            } else {
                backend::residual(rhs, A, x, *B);
            }

            scalar_type zeta0 = norm(*B);
            scalar_type eps = std::max(prm.tol * norm_rhs, prm.abstol);

            coef_type alpha = zero;
            coef_type rho0  = one;
            coef_type omega = one;

            // Go
            backend::copy(*B, *R[0]);
            backend::copy(*B, *Rt);
            backend::clear(*X);
            backend::clear(*U[0]);

            scalar_type zeta           = zeta0;
            scalar_type rnmax_computed = zeta0;
            scalar_type rnmax_true     = zeta0;

            size_t iter = 0;
            for(; iter < prm.maxiter && zeta >= eps; iter += L) {
                // BiCG part
                rho0 = -omega * rho0;

                for(int j = 0; j < L; ++j) {
                    coef_type rho1 = inner_product(*R[j], *Rt);
                    precondition(!math::is_zero(rho1),
                            "BiCGStab(L) breakdown: diverged (zero rho)");

                    coef_type beta = alpha * (rho1 / rho0);
                    rho0 = rho1;

                    for(int i = 0; i <= j; ++i)
                        backend::axpby(one, *R[i], -beta, *U[i]);

                    preconditioner::spmv(prm.pside, P, A, *U[j], *U[j+1], *T);

                    coef_type sigma = inner_product(*U[j+1], *Rt);
                    precondition(!math::is_zero(sigma),
                            "BiCGStab(L) breakdown: diverged (zero sigma)");
                    alpha = rho1 / sigma;

                    backend::axpby(alpha, *U[0], one, *X);

                    for(int i = 0; i <= j; ++i)
                        backend::axpby(-alpha, *U[i+1], one, *R[i]);

                    preconditioner::spmv(prm.pside, P, A, *R[j], *R[j+1], *T);

                    zeta = norm(*R[0]);

                    rnmax_computed = std::max(zeta, rnmax_computed);
                    rnmax_true     = std::max(zeta, rnmax_true);

                    // Check for early exit
                    if (zeta < eps) {
                        iter += j+1;
                        goto done;
                    }
                }

                // Polynomial part
                for(int i = 0; i <= L; ++i) {
                    for(int j = 0; j <= i; ++j) {
                        MZa(i, j) = inner_product(*R[i], *R[j]);
                    }
                }

                // Symmetrize MZa
                for (int i = 0; i <= L; ++i) {
                    for (int j = i+1; j <= L; ++j) {
                        MZa(i, j) = MZa(j, i) = math::adjoint(MZa(j, i));
                    }
                }

                std::copy(MZa.data(), MZa.data() + MZa.size(), MZb.data());

                if (prm.convex || L == 1) {
                    Y0[0] = -one;

                    qr.solve(L, L, MZa.stride(0), MZa.stride(1),
                            &MZa(1, 1), &MZb(0, 1), &Y0[1]);
                } else {
                    Y0[0] = -one;
                    Y0[L] = zero;
                    qr.solve(L-1, L-1, MZa.stride(0), MZa.stride(1),
                            &MZa(1, 1), &MZb(0, 1), &Y0[1]);

                    YL[0] = zero;
                    YL[L] = -one;
                    qr.solve(L-1, L-1, MZa.stride(0), MZa.stride(1),
                            &MZa(1, 1), &MZb(L, 1), &YL[1], /*computed=*/true);

                    coef_type dot0 = zero;
                    coef_type dot1 = zero;
                    coef_type dotA = zero;
                    for(int i = 0; i <= L; ++i) {
                        coef_type s0 = zero;
                        coef_type sL = zero;

                        for(int j = 0; j <= L; ++j) {
                            coef_type M = MZb(i, j);
                            s0 += M * Y0[j];
                            sL += M * YL[j];
                        }

                        dot0 += Y0[i] * s0;
                        dotA += YL[i] * s0;
                        dot1 += YL[i] * sL;
                    }

                    scalar_type kappa0 = sqrt(std::abs(std::real(dot0)));
                    scalar_type kappa1 = sqrt(std::abs(std::real(dot1)));
                    scalar_type kappaA = std::real(dotA);

                    if (!math::is_zero(kappa0) && !math::is_zero(kappa1)) {
                        scalar_type ghat;
                        if (kappaA < 0.7 * kappa0 * kappa1) {
                            ghat = (kappaA < 0) ? -0.7 * kappa0 / kappa1 : 0.7 * kappa0 / kappa1;
                        } else {
                            ghat = kappaA / (kappa1 * kappa1);
                        }

                        for (int i = 0; i <= L; ++i)
                            Y0[i] -= ghat * YL[i];
                    }
                }

                omega = Y0[L];
                for(int h = L; h > 0 && math::is_zero(omega); --h)
                    omega = Y0[h];
                precondition(!math::is_zero(omega),
                        "BiCGStab(L) breakdown: diverged (zero omega)");

                backend::lin_comb(L, &Y0[1], &R[0], one, *X);

                for(int i = 1; i <= L; ++i) Y0[i] = -one * Y0[i];

                backend::lin_comb(L, &Y0[1], &U[1], one, *U[0]);
                backend::lin_comb(L, &Y0[1], &R[1], one, *R[0]);

                for(int i = 1; i <= L; ++i) Y0[i] = -one * Y0[i];

                zeta = norm(*R[0]);

                // Accurate update
                if (prm.delta > 0) {
                    rnmax_computed = std::max(zeta, rnmax_computed);
                    rnmax_true     = std::max(zeta, rnmax_true);

                    bool update_x = zeta < prm.delta * zeta0 && zeta0 <= rnmax_computed;

                    if ((zeta < prm.delta * rnmax_true && zeta <= rnmax_true) || update_x) {
                        preconditioner::spmv(prm.pside, P, A, *X, *R[0], *T);
                        backend::axpby(one, *B, -one, *R[0]);
                        rnmax_true = zeta;

                        if (update_x) {
                            if (prm.pside == side::left) {
                                backend::axpby(one, *X, one, x);
                            } else {
                                backend::axpby(one, *T, one, x);
                            }
                            backend::clear(*X);
                            backend::copy(*R[0], *B);

                            rnmax_computed = zeta;
                        }
                    }
                }
            }

done:
            if (prm.pside == side::left) {
                backend::axpby(one, *X, one, x);
            } else {
                P.apply(*X, *T);
                backend::axpby(one, *T, one, x);
            }

            return std::make_tuple(iter, zeta / norm_rhs);
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
                const Precond &P, const Vec1 &rhs, Vec2 &&x) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }

        size_t bytes() const {
            size_t b = 0;

            b += backend::bytes(*Rt);
            b += backend::bytes(*X);
            b += backend::bytes(*B);
            b += backend::bytes(*T);

            for(const auto &v : R) b += backend::bytes(*v);
            for(const auto &v : U) b += backend::bytes(*v);

            b += MZa.size() * sizeof(coef_type);
            b += MZb.size() * sizeof(coef_type);

            b += backend::bytes(Y0);
            b += backend::bytes(YL);

            b += qr.bytes();

            return b;
        }

        friend std::ostream& operator<<(std::ostream &os, const bicgstabl &s) {
            return os
                << "Type:             BiCGStab(" << s.prm.L << ")"
                << "\nUnknowns:         " << s.n
                << "\nMemory footprint: " << human_readable_memory(s.bytes())
                << std::endl;
        }
    public:
        params prm;

    private:
        size_t n;

        mutable std::shared_ptr< vector > Rt;
        mutable std::shared_ptr< vector > X;
        mutable std::shared_ptr< vector > B;
        mutable std::shared_ptr< vector > T;

        mutable std::vector< std::shared_ptr< vector > > R;
        mutable std::vector< std::shared_ptr< vector > > U;

        mutable multi_array<coef_type, 2> MZa, MZb;
        mutable std::vector<coef_type> Y0, YL;
        mutable amgcl::detail::QR<coef_type> qr;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return sqrt(math::norm(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl


#endif
