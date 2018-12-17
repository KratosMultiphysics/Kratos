#ifndef AMGCL_RELAXATION_CHEBYSHEV_HPP
#define AMGCL_RELAXATION_CHEBYSHEV_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

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

            /// Lowest-to-highest eigen value ratio.
            float lower;

            // Number of power iterations to apply for the spectral radius
            // estimation. When 0, use Gershgorin disk theorem to estimate
            // spectral radius.
            int power_iters;

            params() : degree(5), lower(1.0f / 30), power_iters(0) {}

#ifndef AMGCL_NO_BOOST
            params(const boost::property_tree::ptree &p)
                : AMGCL_PARAMS_IMPORT_VALUE(p, degree),
                  AMGCL_PARAMS_IMPORT_VALUE(p, lower),
                  AMGCL_PARAMS_IMPORT_VALUE(p, power_iters)
            {
                check_params(p, {"degree", "lower", "power_iters"});
            }

            void get(boost::property_tree::ptree &p, const std::string &path) const {
                AMGCL_PARAMS_EXPORT_VALUE(p, path, degree);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, lower);
                AMGCL_PARAMS_EXPORT_VALUE(p, path, power_iters);
            }
#endif
        } prm;

        /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
        template <class Matrix>
        chebyshev(
                const Matrix &A, const params &prm,
                const typename Backend::params &backend_prm
            ) : C( prm.degree ),
                p( Backend::create_vector(rows(A), backend_prm) ),
                q( Backend::create_vector(rows(A), backend_prm) )
        {
            scalar_type hi = backend::spectral_radius</*scale=*/false>(A, prm.power_iters);
            scalar_type lo = hi * prm.lower;

            // Chebyshev polynomial roots on the interval [lo, hi].
            std::vector<scalar_type> roots(prm.degree);
            for(unsigned i = 0; i < prm.degree; ++i) {
                scalar_type pi   = static_cast<scalar_type>(3.14159265358979323846);
                scalar_type half = static_cast<scalar_type>(0.5);

                roots[i] = lo + half * (hi - lo) * (1 + cos( pi * ( i + half ) / prm.degree));
            }

            // Construct linear system to determine Chebyshev coefficients.
            multi_array<scalar_type, 2> S(prm.degree, prm.degree);
            std::vector<scalar_type> buf(prm.degree * prm.degree);

            std::vector<scalar_type> rhs(prm.degree);
            for(unsigned i = 0; i < prm.degree; ++i) {
                scalar_type x = roots[i];
                scalar_type x_to_j = 1;
                for(unsigned j = 0; j < prm.degree; ++j) {
                    S(i,j) = x_to_j;
                    x_to_j *= x;
                }
                rhs[i] = -x_to_j;
            }

            // Invert S, compute coefficients.
            amgcl::detail::inverse(prm.degree, S.data(), buf.data());

            scalar_type const_c = 1;
            for(unsigned i = 0; i < prm.degree; ++i) {
                scalar_type c = 0;
                for(unsigned j = 0; j < prm.degree; ++j)
                    c += S(i,j) * rhs[j];
                if (i == 0)
                    const_c = c;
                else
                    C[prm.degree - i] = -c / const_c;
            }
            C[0] = -1 / const_c;
        }

        /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_pre(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
                ) const
        {
            static const scalar_type one  = math::identity<scalar_type>();

            backend::residual(rhs, A, x, tmp);
            solve(A, tmp, *p);
            backend::axpby(one, *p, one, x);
        }

        /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
        template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
        void apply_post(
                const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
                ) const
        {
            static const scalar_type one  = math::identity<scalar_type>();

            backend::residual(rhs, A, x, tmp);
            solve(A, tmp, *p);
            backend::axpby(one, *p, one, x);
        }

        /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
        template <class Matrix, class VectorRHS, class VectorX>
        void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
        {
            solve(A, rhs, x);
        }

        size_t bytes() const {
            return
                backend::bytes(C) +
                backend::bytes(*p) +
                backend::bytes(*q);
        }

    private:
        std::vector<scalar_type> C;
        mutable std::shared_ptr<vector> p, q;

        template <class Matrix, class VectorRHS, class VectorX>
        void solve(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
        {
            static const scalar_type one  = math::identity<scalar_type>();
            static const scalar_type zero = math::zero<scalar_type>();

            backend::axpby(C[0], rhs, zero, x);

            for(auto c = C.begin() + 1; c != C.end(); ++c)
            {
                backend::spmv(one, A, x, zero, *q);
                backend::axpbypcz(*c, rhs, one, *q, zero, x);
            }
        }
};

} // namespace relaxation
} // namespace amgcl

#endif
