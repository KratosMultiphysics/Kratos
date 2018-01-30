#ifndef AMGCL_SOLVER_GMRES_HPP
#define AMGCL_SOLVER_GMRES_HPP

/*
The MIT License

Copyright (c) 2012-2017 Denis Demidov <dennis.demidov@gmail.com>

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
#include <algorithm>
#include <cmath>

#include <boost/multi_array.hpp>
#include <boost/tuple/tuple.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/solver/detail/givens_rotations.hpp>
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
            unsigned M;

            /// Maximum number of iterations.
            unsigned maxiter;

            /// Target relative residual error.
            scalar_type tol;

            /// Target absolute residual error.
            scalar_type abstol;

            params(unsigned M = 30, unsigned maxiter = 100, scalar_type tol = 1e-8)
                : M(M), maxiter(maxiter), tol(tol),
                  abstol(std::numeric_limits<scalar_type>::min())
            { }

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, M),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol),
                  AMGCL_PARAMS_IMPORT_VALUE(p, abstol)
            {
                AMGCL_PARAMS_CHECK(p, (M)(maxiter)(tol)(abstol));
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, M);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, abstol);
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
              s(prm.M + 1), cs(prm.M + 1), sn(prm.M + 1),
              r( Backend::create_vector(n, backend_prm) ),
              inner_product(inner_product)
        {
            v.reserve(prm.M + 1);
            for(unsigned i = 0; i <= prm.M; ++i)
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

            scalar_type eps = std::max(prm.tol * norm_rhs, prm.abstol);
            scalar_type norm_r = math::zero<scalar_type>();

            while(true) {
                backend::residual(rhs, A, x, *r);

                // -- Check stopping condition
                if ((norm_r = norm(*r)) < eps || iter >= prm.maxiter)
                    break;

                // -- Inner GMRES iteration
                P.apply(*r, *v[0]);

                std::fill(s.begin(), s.end(), 0);
                s[0] = norm(*v[0]);

                precondition(!math::is_zero(s[0]),
                        "Preconditioner returned a zero vector");

                backend::axpby(math::inverse(s[0]), *v[0], math::zero<scalar_type>(), *v[0]);


                unsigned j = 0;
                while(true) {
                    // -- Arnoldi process
                    //
                    // Build an orthonormal basis V and matrix H such that
                    //     A V_{i-1} = V_{i} H
                    vector &v_new = *v[j+1];
                    backend::spmv(math::identity<scalar_type>(), A, *v[j], math::zero<scalar_type>(), *r);
                    P.apply(*r, v_new);

                    for(unsigned k = 0; k <= j; ++k) {
                        H[k][j] = inner_product(v_new, *v[k]);
                        backend::axpby(-H[k][j], *v[k], math::identity<scalar_type>(), v_new);
                    }
                    H[j+1][j] = norm(v_new);

                    backend::axpby(math::inverse(H[j+1][j]), v_new, math::zero<scalar_type>(), v_new);

                    for(unsigned k = 0; k < j; ++k)
                        detail::apply_plane_rotation(H[k][j], H[k+1][j], cs[k], sn[k]);

                    detail::generate_plane_rotation(H[j][j], H[j+1][j], cs[j], sn[j]);
                    detail::apply_plane_rotation(H[j][j], H[j+1][j], cs[j], sn[j]);
                    detail::apply_plane_rotation(s[j], s[j+1], cs[j], sn[j]);

                    scalar_type inner_res = std::abs(s[j+1]);

                    // Check for termination
                    ++j, ++iter;
                    if (iter >= prm.maxiter || j >= prm.M || inner_res <= eps)
                        break;
                }

                // -- GMRES terminated: eval solution
                for (unsigned i = j; i --> 0; ) {
                    s[i] /= H[i][i];
                    for (unsigned k = 0; k < i; ++k)
                        s[k] -= H[k][i] * s[i];
                }

                backend::lin_comb(j, s, v, math::identity<scalar_type>(), x);
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


        friend std::ostream& operator<<(std::ostream &os, const gmres &s) {
            return os << "gmres(" << s.prm.M << "): " << s.n << " unknowns";
        }
    public:
        params prm;

    private:
        size_t n;

        mutable boost::multi_array<coef_type, 2> H;
        mutable std::vector<coef_type> s, cs, sn;
        boost::shared_ptr<vector> r;
        std::vector< boost::shared_ptr<vector> > v;

        InnerProduct inner_product;

        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return std::abs(sqrt(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl

#endif
