#ifndef AMGCL_SOLVER_FGMRES_HPP
#define AMGCL_SOLVER_FGMRES_HPP

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
 * \file   fgmres.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Flexible GMRES method.
 */

#include <vector>
#include <cmath>

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/range/algorithm.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace solver {

/** Flexible GMRES method.
 * \rst
 * Flexible version of the GMRES method [Saad03]_.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class fgmres {
    public:
        typedef Backend backend_type;

        typedef typename Backend::vector     vector;
        typedef typename Backend::value_type value_type;
        typedef typename Backend::params     backend_params;

        typedef typename math::scalar_of<value_type>::type scalar_type;
        typedef typename math::rhs_of<value_type>::type rhs_type;
        typedef typename math::inner_product_impl<rhs_type>::return_type coef_type;

        /// Solver parameters.
        struct params {
            /// Number of inner GMRES iterations per each outer iteration.
            unsigned M;

            /// Maximum number of iterations.
            unsigned maxiter;

            /// Target residual error.
            scalar_type tol;

            params() : M(30), maxiter(100), tol(1e-8) { }

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, M),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol)
            {
                AMGCL_PARAMS_CHECK(p, (M)(maxiter)(tol));
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, M);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
            }
        } prm;

        /// Preallocates necessary data structures for the system of size \p n.
        fgmres(
                size_t n,
                const params &prm = params(),
                const backend_params &bprm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
             )
            : prm(prm), n(n), inner_product(inner_product),
              H(boost::extents[prm.M + 1][prm.M + 1], boost::fortran_storage_order()),
              y(prm.M)
        {
            vs.reserve(prm.M + 1);
            for(unsigned i = 0; i <= prm.M; ++i)
                vs.push_back(Backend::create_vector(n, bprm));

            zs.reserve(prm.M);
            for(unsigned i = 0; i < prm.M; ++i)
                zs.push_back(Backend::create_vector(n, bprm));
        }

        /* Computes the solution for the given system matrix \p A and the
         * right-hand side \p rhs.  Returns the number of iterations made and
         * the achieved residual as a ``boost::tuple``. The solution vector
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
        boost::tuple<size_t, scalar_type> operator()(
                Matrix  const &A,
                Precond const &P,
                Vec1    const &rhs,
                Vec2          &x
                ) const
        {
            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            scalar_type norm_r = math::zero<scalar_type>();

            unsigned iter = 0;
            while(true) {
                backend::residual(rhs, A, x, *vs[0]);

                // -- Check stopping condition
                if ((norm_r = norm(*vs[0])) < prm.tol * norm_rhs || iter >= prm.maxiter)
                    break;

                // -- Inner GMRES iteration
                backend::axpby(math::inverse(norm_r), *vs[0],
                        math::zero<scalar_type>(), *vs[0]);

                // H is stored in QR factorized form
                for(unsigned j = 0; j <= prm.M; ++j)
                    for(unsigned i = 0; i <= prm.M; ++i)
                        H[i][j] = math::zero<coef_type>();

                qr.compute(prm.M+1, 0, H.data(), prm.M+1);

                const scalar_type eps = std::numeric_limits<scalar_type>::epsilon();
                bool breakdown = false;

                unsigned j = 0;
                while(true) {
                    // -- Arnoldi process
                    //
                    // Build an orthonormal basis V and matrix H such that
                    //     A V_{i-1} = V_{i} H

                    vector &v_new = *vs[j+1];

                    P.apply(*vs[j], *zs[j]);
                    backend::spmv(math::identity<scalar_type>(), A, *zs[j],
                            math::zero<scalar_type>(), v_new);

                    scalar_type v_new_norm = norm(v_new);

                    for(unsigned i = 0; i <= j; ++i) {
                        H[i][j] = inner_product(*vs[i], v_new);
                        backend::axpby(-H[i][j], *vs[i], math::identity<coef_type>(), v_new);
                    }
                    H[j+1][j] = norm(v_new);

                    // Careful with denormals:
                    coef_type alpha = math::inverse(H[j+1][j]);
                    if (boost::math::isfinite(alpha))
                        backend::axpby(alpha, v_new, math::zero<coef_type>(), v_new);

                    if (!(math::norm(H[j+1][j]) > eps * v_new_norm)) {
                        // v_new essentially in the span of previous vectors,
                        // or we have nans. Bail out after updating the QR
                        // solution.
                        breakdown = true;
                    }

                    // -- GMRES optimization problem
                    //
                    // Add new column to H = Q*R
                    qr.append_cols(1);

                    // Transformed least squares problem
                    // || Q R y - norm_r * e_1 ||_2 = min!
                    // Since R = [R'; 0], solution is y = norm_r (R')^{-1} (Q^H)[:j,0]
                    //
                    // Residual is immediately known
                    qr.compute_q(j+2);
                    scalar_type inner_res = std::abs(qr.Q(0,j+1)) * norm_r;

                    // Check for termination
                    ++j, ++iter;
                    if (iter >= prm.maxiter || j >= prm.M)
                        break;

                    if (inner_res <= prm.tol * norm_rhs || breakdown)
                        break;
                }

                precondition(boost::math::isfinite(qr.R(j-1,j-1)),
                        "NaNs encountered in FGMRES");

                // The problem is triangular, but the condition number may be
                // bad (or in case of breakdown the last diagonal entry may be
                // zero), so use lstsq instead of triangular solve.
                //
                // TODO: This is triangular solve for now.
                for(unsigned i = 0; i < j; ++i) y[i] = math::adjoint(qr.Q(0, i));
                for(unsigned i = j; i --> 0; ) {
                    coef_type rii = qr.R(i,i);
                    if (math::is_zero(rii)) continue;
                    y[i] = math::inverse(rii) * y[i];
                    for(unsigned k = 0; k < i; ++k)
                        y[k] -= qr.R(k, i) * y[i];

                    y[i] *= norm_r;

                    precondition(boost::math::isfinite(y[i]),
                            "NaNs encountered in FGMRES");
                }

                // -- GMRES terminated: eval solution
                unsigned k = 0;
                for(; k + 1 < j; k += 2)
                    backend::axpbypcz(y[k], *zs[k], y[k+1], *zs[k+1],
                            math::identity<coef_type>(), x);

                for(; k < j; ++k)
                    backend::axpby(y[k], *zs[k], math::identity<coef_type>(), x);
            }

            return boost::make_tuple(iter, norm_r / norm_rhs);
        }

        /* Computes the solution for the given right-hand side \p rhs. The
         * system matrix is the same that was used for the setup of the
         * preconditioner \p P.  Returns the number of iterations made and the
         * achieved residual as a ``boost::tuple``. The solution vector \p x
         * provides initial approximation in input and holds the computed
         * solution on output.
         */
        template <class Precond, class Vec1, class Vec2>
        boost::tuple<size_t, scalar_type> operator()(
                Precond const &P,
                Vec1    const &rhs,
                Vec2          &x
                ) const
        {
            return (*this)(P.system_matrix(), P, rhs, x);
        }

    private:
        size_t n;
        InnerProduct inner_product;

        mutable boost::multi_array<coef_type, 2> H;
        mutable amgcl::detail::QR<coef_type, amgcl::detail::col_major> qr;
        mutable std::vector<coef_type> y;
        mutable std::vector< boost::shared_ptr<vector> > vs, zs;


        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return std::abs(sqrt(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl

#endif
