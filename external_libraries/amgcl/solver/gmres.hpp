#ifndef AMGCL_SOLVER_GMRES_HPP
#define AMGCL_SOLVER_GMRES_HPP

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
 * \file   gmres.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  GMRES method.
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

/** Generalized Minimal Residual (GMRES) method.
 * \rst
 * The Generalized Minimal Residual method is an extension of MINRES (which is
 * only applicable to symmetric systems) to unsymmetric systems [Barr94]_.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class gmres {
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
            /// Number of iterations before restart.
            int M;

            /// Maximum number of iterations.
            size_t maxiter;

            /// Target residual error.
            scalar_type tol;

            params(int M = 50, size_t maxiter = 100, scalar_type tol = 1e-8)
                : M(M), maxiter(maxiter), tol(tol)
            {
                precondition(M > 0, "M in GMRES(M) should be >=1");
            }

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
        };

        /// Preallocates necessary data structures for the system of size \p n.
        gmres(
                size_t n,
                const params &prm = params(),
                const backend_params &backend_prm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
             )
            : prm(prm), n(n),
              H(boost::extents[prm.M + 1][prm.M]),
              s(prm.M + 1), cs(prm.M + 1), sn(prm.M + 1), y(prm.M + 1),
              r( Backend::create_vector(n, backend_prm) ),
              w( Backend::create_vector(n, backend_prm) ),
              inner_product(inner_product)
        {
            v.reserve(prm.M + 1);
            for(int i = 0; i <= prm.M; ++i)
                v.push_back( Backend::create_vector(n, backend_prm) );
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
            size_t iter = 0;

            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            scalar_type eps = prm.tol * norm_rhs;

            scalar_type res_norm = restart(A, rhs, P, x);
            if (res_norm < eps)
                return boost::make_tuple(0, res_norm / norm_rhs);

            do {
                for(int i = 0; i < prm.M && iter < prm.maxiter; ++i, ++iter) {
                    res_norm = iteration(A, P, i);

                    if (res_norm < eps) {
                        update(x, i);
                        return boost::make_tuple(iter + 1, res_norm / norm_rhs);
                    };
                }

                update(x, prm.M-1);
                res_norm = restart(A, rhs, P, x);
            } while (iter < prm.maxiter && res_norm > eps);

            return boost::make_tuple(iter, res_norm / norm_rhs);
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

    public:
        params prm;

    private:
        size_t n;

        mutable boost::multi_array<coef_type, 2> H;
        mutable std::vector<coef_type> s, cs, sn, y;
        boost::shared_ptr<vector> r, w;
        std::vector< boost::shared_ptr<vector> > v;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return std::abs(sqrt(inner_product(x, x)));
        }

        static void apply_plane_rotation(
                coef_type &dx, coef_type &dy, coef_type cs, coef_type sn
                )
        {
            coef_type tmp = cs * dx + sn * dy;
            dy = -sn * dx + cs * dy;
            dx = tmp;
        }

        static void generate_plane_rotation(
                coef_type dx, coef_type dy, coef_type &cs, coef_type &sn
                )
        {
            if (math::is_zero(dy)) {
                cs = 1;
                sn = 0;
            } else if (std::abs(dy) > std::abs(dx)) {
                coef_type tmp = dx / dy;
                sn = math::inverse(sqrt(math::identity<coef_type>() + tmp * tmp));
                cs = tmp * sn;
            } else {
                coef_type tmp = dy / dx;
                cs = math::inverse(sqrt(math::identity<coef_type>() + tmp * tmp));
                sn = tmp * cs;
            }
        }

        template <class Vec>
        void update(Vec &x, int k) const {
            boost::range::copy(s, y.begin());

            for (int i = k; i >= 0; --i) {
                y[i] /= H[i][i];
                for (int j = i - 1; j >= 0; --j)
                    y[j] -= H[j][i] * y[i];
            }

            // Unroll the loop
            int j = 0;
            for (; j + 1 <= k; j += 2)
                backend::axpbypcz(y[j], *v[j], y[j+1], *v[j+1], math::identity<scalar_type>(), x);
            for (; j <= k; ++j)
                backend::axpby(y[j], *v[j], math::identity<scalar_type>(), x);
        }

        template <class Matrix, class Precond, class Vec1, class Vec2>
        scalar_type restart(const Matrix &A, const Vec1 &rhs,
                const Precond &P, const Vec2 &x) const
        {
            backend::residual(rhs, A, x, *w);
            P.apply(*w, *r);

            boost::fill(s, 0);
            s[0] = norm(*r);

            if (!math::is_zero(s[0]))
                backend::axpby(math::inverse(s[0]), *r, math::zero<scalar_type>(), *v[0]);

            return std::abs(s[0]);
        }

        template <class Matrix, class Precond>
        scalar_type iteration(const Matrix &A, const Precond &P, int i) const
        {
            backend::spmv(math::identity<scalar_type>(), A, *v[i], math::zero<scalar_type>(), *r);
            P.apply(*r, *w);

            for(int k = 0; k <= i; ++k) {
                H[k][i] = inner_product(*w, *v[k]);
                backend::axpby(-H[k][i], *v[k], math::identity<scalar_type>(), *w);
            }

            H[i+1][i] = norm(*w);

            backend::axpby(math::inverse(H[i+1][i]), *w, math::zero<scalar_type>(), *v[i+1]);

            for(int k = 0; k < i; ++k)
                apply_plane_rotation(H[k][i], H[k+1][i], cs[k], sn[k]);

            generate_plane_rotation(H[i][i], H[i+1][i], cs[i], sn[i]);
            apply_plane_rotation(H[i][i], H[i+1][i], cs[i], sn[i]);
            apply_plane_rotation(s[i], s[i+1], cs[i], sn[i]);

            return std::abs(s[i+1]);
        }
};

} // namespace solver
} // namespace amgcl

#endif
