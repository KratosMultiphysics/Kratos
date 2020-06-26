#ifndef AMGCL_RELAXATION_CHEBYSHEV_HPP
#define AMGCL_RELAXATION_CHEBYSHEV_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>
Copyright (c) 2019 Peter Gamnitzer, UIBK (University of Innsbruck)

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
 * \file   amgcl/relaxation/chebyshev.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Chebyshev polynomial smoother.
 *
 * Implements Algorithm 1 from
 * P. Ghysels, P. Kłosiewicz, and W. Vanroose.
 * "Improving the arithmetic intensity of multigrid with the help of polynomial smoothers".
 * Numer. Linear Algebra Appl. 2012;19:253-267. DOI: 10.1002/nla.1808
 */

#include <vector>
#include <cmath>

#include <amgcl/detail/inverse.hpp>
#include <amgcl/util.hpp>

namespace amgcl {
namespace relaxation {

/// Chebyshev polynomial smoother.
/**
 * \param Backend Backend for temporary structures allocation.
 * \ingroup relaxation
 */
template <class Backend>
class chebyshev {
    public:
        typedef typename Backend::value_type value_type;
        typedef typename Backend::vector     vector;

        typedef typename math::scalar_of<value_type>::type scalar_type;

        /// Relaxation parameters.
        struct params {
            /// Chebyshev polynomial degree.
            unsigned degree;

            /// highest eigen value safety upscaling.
            // use boosting factor for a more conservative upper bound estimate
            // See: Adams, Brezina, Hu, Tuminaro,
            //      PARALLEL MULTIGRID SMOOTHING: POLYNOMIAL VERSUS
            //      GAUSS-SEIDEL, J. Comp. Phys. 188 (2003) 593-610.
            //
            float higher;

            /// Lowest-to-highest eigen value ratio.
            float lower;

            // Number of power iterations to apply for the spectral radius
            // estimation. When 0, use Gershgorin disk theorem to estimate
            // spectral radius.
            int power_iters;

            // Scale the system matrix
            bool scale;

            params()
                : degree(5), higher(1.0f), lower(1.0f / 30), power_iters(0),
                  scale(false)
            {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, degree),
                  AMGCL_PARAMS_IMPORT_VALUE(p, higher),
                  AMGCL_PARAMS_IMPORT_VALUE(p, lower),
                  AMGCL_PARAMS_IMPORT_VALUE(p, power_iters),
                  AMGCL_PARAMS_IMPORT_VALUE(p, scale)
            {
                check_params(p, {"degree", "higher", "lower", "power_iters", "scale"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, degree);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, higher);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, lower);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, power_iters);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, scale);
            }
#endif
        } prm;

        /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
        template <class Matrix>
        chebyshev(
                const Matrix &A, const params &prm,
                const typename Backend::params &backend_prm
            ) : prm(prm),
                p( Backend::create_vector(rows(A), backend_prm) ),
                r( Backend::create_vector(rows(A), backend_prm) )
        {
            scalar_type hi, lo;

            if (prm.scale) {
                M  = Backend::copy_vector( diagonal(A, /*invert*/true), backend_prm );
                hi = backend::spectral_radius<true>(A, prm.power_iters);
            } else {
                hi = backend::spectral_radius<false>(A, prm.power_iters);
            }

            lo = hi * prm.lower;
            hi *= prm.higher;

            // Centre of ellipse containing the eigenvalues of A:
            d = 0.5 * (hi + lo);

            // Semi-major axis of ellipse containing the eigenvalues of A:
            c = 0.5 * (hi - lo);
        }

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_pre(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP&
                ) const
        {
            solve(A, rhs, x);
        }

        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_post(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP&
                ) const
        {
            solve(A, rhs, x);
        }

        template <class Matrix, class VectorRHS, class VectorX>
        void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
        {
            backend::clear(x);
            solve(A, rhs, x);
        }

        size_t bytes() const {
            size_t b = backend::bytes(*p) + backend::bytes(*r);
            if (prm.scale) b += backend::bytes(*M);
            return b;
        }

    private:
        std::shared_ptr<typename Backend::matrix_diagonal> M;
        mutable std::shared_ptr<vector> p, r;

        scalar_type c, d;

        template <class Matrix, class VectorB, class VectorX>
        void solve(const Matrix &A, const VectorB &b, VectorX &x) const
        {
            static const scalar_type one  = math::identity<scalar_type>();
            static const scalar_type zero = math::zero<scalar_type>();

            scalar_type alpha = zero, beta = zero;

            for (unsigned k = 0; k < prm.degree; ++k) {
                backend::residual(b, A, x, *r);

                if (prm.scale) backend::vmul(one, *M, *r, zero, *r);

                if (k == 0) {
                    alpha = math::inverse(d);
                    beta  = zero;
                } else if (k == 1) {
                    alpha = 2 * d * math::inverse(2 * d * d - c * c);
                    beta  = alpha * d - one;
                } else {
                    alpha = math::inverse(d - 0.25 * alpha * c * c);
                    beta  = alpha * d - one;
                }

                backend::axpby(alpha, *r, beta, *p);
                backend::axpby(one, *p, one, x);
            }
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
