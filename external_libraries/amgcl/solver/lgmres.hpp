#ifndef AMGCL_SOLVER_LGMRES_HPP
#define AMGCL_SOLVER_LGMRES_HPP

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
 * \file   lgmres.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  LGMRES method.
 *
 * Ported from scipy lgmres. The original code came with the following license:
 * \verbatim
   Copyright (c) 2001, 2002 Enthought, Inc.
   All rights reserved.

   Copyright (c) 2003-2016 SciPy Developers.
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

     a. Redistributions of source code must retain the above copyright notice,
        this list of conditions and the following disclaimer.
     b. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
     c. Neither the name of Enthought nor the names of the SciPy Developers
        may be used to endorse or promote products derived from this software
        without specific prior written permission.


   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
   BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
   OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   THE POSSIBILITY OF SUCH DAMAGE.
 * \endverbatim
 */

#include <vector>
#include <cmath>

#include <boost/multi_array.hpp>
#include <boost/circular_buffer.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/range/algorithm.hpp>

#include <amgcl/backend/interface.hpp>
#include <amgcl/solver/detail/default_inner_product.hpp>
#include <amgcl/solver/detail/givens_rotations.hpp>
#include <amgcl/util.hpp>

namespace amgcl {

namespace precond {

enum type {
    left,
    right
};

inline std::ostream& operator<<(std::ostream &os, type p) {
    switch (p) {
        case left:
            return os << "left";
        case right:
            return os << "right";
        default:
            return os << "???";
    }
}

inline std::istream& operator>>(std::istream &in, type &p) {
    std::string val;
    in >> val;

    if (val == "left")
        p = left;
    else if (val == "right")
        p = right;
    else
        throw std::invalid_argument("Invalid preconditioning kind");

    return in;
}

} // precond

namespace solver {

/** "Loose" GMRES.
 * \rst
 * The LGMRES algorithm [BaJM05]_  is designed to avoid some problems
 * in the convergence in restarted GMRES, and often converges in fewer
 * iterations.
 * \endrst
 */
template <
    class Backend,
    class InnerProduct = detail::default_inner_product
    >
class lgmres {
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
            /// Preconditioning kind (left/right).
            precond::type pside;

            /// Number of inner GMRES iterations per each outer iteration.
            unsigned M;

            /// Number of vectors to carry between inner GMRES iterations.
            /**
             * According to [BaJM05], good values are in the range of 1...3.
             * However, note that if you want to use the additional vectors to
             * accelerate solving multiple similar problems, larger values may
             * be beneficial.
             */
            unsigned K;

            /// Reset augmented vectors between solves.
            /** If the solver is used to repeatedly solve similar problems,
             *  then keeping the augmented vectors between solves may speed up
             *  subsequent solves.
             *  This flag, when set, resets the augmented vectors at the
             *  beginning of each solve.
             */
            bool always_reset;

            /// Whether LGMRES should store also A*v in addition to vectors `v`.
            bool store_Av;

            /// Maximum number of iterations.
            size_t maxiter;

            /// Target relative residual error.
            scalar_type tol;

            /// Target absolute residual error.
            scalar_type abstol;

            params()
                : pside(precond::right), M(30), K(3), always_reset(true),
                  store_Av(true), maxiter(100), tol(1e-8),
                  abstol(std::numeric_limits<scalar_type>::min())
            { }

            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, pside),
                  AMGCL_PARAMS_IMPORT_VALUE(p, M),
                  AMGCL_PARAMS_IMPORT_VALUE(p, K),
                  AMGCL_PARAMS_IMPORT_VALUE(p, always_reset),
                  AMGCL_PARAMS_IMPORT_VALUE(p, store_Av),
                  AMGCL_PARAMS_IMPORT_VALUE(p, maxiter),
                  AMGCL_PARAMS_IMPORT_VALUE(p, tol),
                  AMGCL_PARAMS_IMPORT_VALUE(p, abstol)
            {
                AMGCL_PARAMS_CHECK(p, (pside)(M)(K)(always_reset)(store_Av)(maxiter)(tol)(abstol));
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, pside);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, M);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, K);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, always_reset);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, store_Av);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, maxiter);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, tol);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, abstol);
            }
        } prm;

        /// Preallocates necessary data structures for the system of size \p n.
        lgmres(
                size_t n,
                const params &prm = params(),
                const backend_params &bprm = backend_params(),
                const InnerProduct &inner_product = InnerProduct()
             )
            : prm(prm), n(n),
              H(boost::extents[prm.M + 1][prm.M]),
              H0(boost::extents[prm.M + 1][prm.M]),
              s(prm.M + 1), cs(prm.M + 1), sn(prm.M + 1),
              r( Backend::create_vector(n, bprm) ),
              ws(prm.M + prm.K), outer_v(prm.K), outer_Av(prm.K),
              inner_product(inner_product)
        {
            if (prm.pside == precond::right)
                t = Backend::create_vector(n, bprm);

            outer_v_data.reserve(prm.K);
            for(unsigned i = 0; i < prm.K; ++i)
                outer_v_data.push_back(Backend::create_vector(n, bprm));

            if (prm.store_Av) {
                y.resize(prm.M + prm.K + 1);
                outer_Av_data.reserve(prm.K);
                for(unsigned i = 0; i < prm.K; ++i)
                    outer_Av_data.push_back(Backend::create_vector(n, bprm));
            }

            vs.reserve(prm.M + prm.K + 1);
            for(unsigned i = 0; i <= prm.M + prm.K; ++i)
                vs.push_back(Backend::create_vector(n, bprm));
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
            if (prm.always_reset) {
                outer_v.clear();
                outer_Av.clear();
            }

            scalar_type norm_rhs = norm(rhs);
            if (norm_rhs < amgcl::detail::eps<scalar_type>(n)) {
                backend::clear(x);
                return boost::make_tuple(0, norm_rhs);
            }

            scalar_type norm_r = math::zero<scalar_type>(), norm_v0;
            scalar_type eps = std::max(prm.tol * norm_rhs, prm.abstol);

            unsigned iter = 0, n_outer = 0;
            while(true) {
                backend::residual(rhs, A, x, *r);

                // -- Check stopping condition
                if ((norm_r = norm(*r)) < eps || iter >= prm.maxiter)
                    break;

                // -- Inner LGMRES iteration
                if (prm.pside == precond::left) {
                    P.apply(*r, *vs[0]);
                    norm_v0 = norm(*vs[0]);

                    precondition(!math::is_zero(norm_v0),
                            "Preconditioner returned a zero vector");

                    backend::axpby(math::inverse(norm_v0), *vs[0],
                            math::zero<scalar_type>(), *vs[0]);
                } else {
                    norm_v0 = norm_r;
                    backend::axpby(math::inverse(norm_r), *r,
                            math::zero<scalar_type>(), *vs[0]);
                }

                boost::fill(s, 0);
                s[0] = norm_v0;

                unsigned j = 0;
                while(true) {
                    // -- Arnoldi process:
                    //
                    // Build an orthonormal basis V and matrices W and H such that
                    //     A W = V H
                    // Columns of W, V, and H are stored in `ws`, `vs` and `hs`.
                    //
                    // The first column of V is always the residual vector,
                    // `vs0`; V has *one more column* than the other of the
                    // three matrices.
                    //
                    // The other columns in V are built by feeding in, one by
                    // one, some vectors `z` and orthonormalizing them against
                    // the basis so far. The trick here is to feed in first
                    // some augmentation vectors, before starting to construct
                    // the Krylov basis on `v0`.
                    //
                    // It was shown in [BaJM05] that a good choice (the LGMRES
                    // choice) for these augmentation vectors are the `dx`
                    // vectors obtained from a couple of the previous restart
                    // cycles.
                    //
                    // Note especially that while `vs0` is always the first
                    // column in V, there is no reason why it should also be
                    // the first column in W. (In fact, below `vs0` comes in W
                    // only after the augmentation vectors.)
                    //
                    // The rest of the algorithm then goes as in GMRES, one
                    // solves a minimization problem in the smaller subspace
                    // spanned by W (range) and V (image).

                    vector &v_new = *vs[j+1];

                    boost::shared_ptr<vector> z;
                    if (j < outer_v.size()) {
                        z = outer_v[j];
                    } else if (j == outer_v.size()) {
                        z = vs[0];
                    } else {
                        z = vs[j];
                    }

                    ws[j] = z;

                    if (j < outer_Av.size()) {
                        backend::copy(*outer_Av[j], v_new);
                    } else {
                        if (prm.pside == precond::left) {
                            backend::spmv(math::identity<scalar_type>(), A, *z,
                                    math::zero<scalar_type>(), *r);
                            P.apply(*r, v_new);
                        } else {
                            P.apply(*z, *r);
                            backend::spmv(math::identity<scalar_type>(), A, *r,
                                    math::zero<scalar_type>(), v_new);
                        }
                    }

                    for(unsigned k = 0; k <= j; ++k) {
                        H0[k][j] = H[k][j] = inner_product(v_new, *vs[k]);
                        backend::axpby(-H[k][j], *vs[k], math::identity<scalar_type>(), v_new);
                    }
                    H0[j+1][j] = H[j+1][j] = norm(v_new);

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

                vector &dx = *r;
                backend::lin_comb(j, s, ws, math::zero<scalar_type>(), dx);

                // -- Apply step
                if (prm.pside == precond::left) {
                    backend::axpby(math::identity<scalar_type>(), dx, math::identity<scalar_type>(), x);
                } else {
                    P.apply(dx, *t);
                    backend::axpby(math::identity<scalar_type>(), *t, math::identity<scalar_type>(), x);
                }

                // -- Store LGMRES augmented vectors
                scalar_type norm_dx = norm(dx);

                if(prm.K > 0 && !math::is_zero(norm_dx)) {
                    unsigned outer_slot = n_outer % prm.K;
                    ++n_outer;

                    norm_dx = math::inverse(norm_dx);
                    backend::axpby(norm_dx, dx, math::zero<scalar_type>(), *outer_v_data[outer_slot]);
                    outer_v.push_back(outer_v_data[outer_slot]);

                    if (prm.store_Av) {
                        // y = H0 * s, AW = Vy
                        for(unsigned k = 0; k <= j; ++k) {
                            coef_type sum = math::zero<coef_type>();
                            for(unsigned i = k ? k-1 : 0; i < j; ++i)
                                sum += H0[k][i] * s[i];
                            y[k] = sum * norm_dx;
                        }

                        backend::lin_comb(j+1, y, vs, math::zero<scalar_type>(),
                                *outer_Av_data[outer_slot]);
                        outer_Av.push_back(outer_Av_data[outer_slot]);
                    }
                }
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


        friend std::ostream& operator<<(std::ostream &os, const lgmres &s) {
            return os << "lgmres(" << s.prm.M << "," << s.prm.K << "): " << s.n << " unknowns";
        }
    private:
        size_t n;

        mutable boost::multi_array<coef_type, 2> H, H0;
        mutable std::vector<coef_type> s, cs, sn, y;
        boost::shared_ptr<vector> r;
        mutable std::vector< boost::shared_ptr<vector> > vs, ws;
        mutable boost::shared_ptr<vector> t;
        mutable std::vector< boost::shared_ptr<vector> > outer_v_data, outer_Av_data;
        mutable boost::circular_buffer< boost::shared_ptr<vector> > outer_v;
        mutable boost::circular_buffer< boost::shared_ptr<vector> > outer_Av;


        InnerProduct inner_product;


        template <class Vec>
        scalar_type norm(const Vec &x) const {
            return std::abs(sqrt(inner_product(x, x)));
        }
};

} // namespace solver
} // namespace amgcl

#endif
