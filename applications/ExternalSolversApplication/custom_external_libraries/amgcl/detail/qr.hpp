#ifndef AMGCL_DETAIL_QR_HPP
#define AMGCL_DETAIL_QR_HPP

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
 * \file   amgcl/detail/qr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  QR decomposition.
 */

#include <vector>
#include <cmath>

#include <amgcl/util.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace detail {

/// In-place QR decomposition.
template <typename value_type>
class QR {
    public:
        QR() : m(0), n(0) {}

        void compute(unsigned rows, unsigned cols, value_type *A, bool needQ = true) {
            m = rows;
            n = cols;
            p = std::min(m, n);

            r = A;

            scalar_type absV0, normV, beta;
            value_type alpha;

            if (diagR.size() < p) diagR.resize(p);
            if (betaR.size() < p) betaR.resize(p);

            for(unsigned k = 0; k < p; ++k) {
                // Form k-th Housholder vector.
                normV = 0;
                for(unsigned i = k; i < m; ++i)
                    normV += sqr(math::norm(r[i*p+k]));
                normV = sqrt(normV);

                absV0 = math::norm(r[k*p+k]);
                alpha = -normV / absV0 * r[k*p+k];
                beta = 1 / (normV * (normV + absV0));
                r[k*p+k] -= alpha;


                // Apply transformation to remaining columns.
                for(unsigned j = k + 1; j < n; ++j) {
                    value_type s = math::zero<value_type>();

                    for(unsigned i = k; i < m; ++i)
                        s += math::adjoint(r[i*p+k]) * r[i*p+j];

                    s *= beta;

                    for(unsigned i = k; i < m; ++i)
                        r[i*p+j] -= s * r[i*p+k];
                }

                diagR[k] = alpha;
                betaR[k] = beta;
            }

            if (needQ) {
                if (q.size() < m * p) q.resize(m * p);

                for(unsigned k=p; k --> 0;) {
                    q[k*p+k] = (1 - betaR[k] * sqr(math::norm(r[k*p+k]))) *
                        math::identity<value_type>();

                    for(unsigned i = k+1; i < m; ++i)
                        q[i*p+k] = -betaR[k] * r[i*p+k] * math::adjoint(r[k*p+k]);

                    for(unsigned j = k+1; j < p; ++j) {
                        value_type s = math::zero<value_type>();

                        for(unsigned i = k; i < m; ++i)
                            s += math::adjoint(r[i*p+k]) * q[i*p+j];

                        s *= betaR[k];

                        for(unsigned i = k; i < m; ++i)
                            q[i*p+j] -= s * r[i*p+k];
                    }
                }
            }
        }

        value_type R(unsigned i, unsigned j) const {
            if (j <  i) return math::zero<value_type>();
            if (j == i) return diagR[i];
            return r[i * p + j];
        }

        value_type Q(unsigned i, unsigned j) const {
            return q[i*p+j];
        }

        void solve(value_type *f, value_type *x) const {
            for(unsigned k = 0; k < n; ++k) {
                value_type s = math::zero<value_type>();

                for(unsigned i = k; i < m; ++i)
                    s += math::adjoint(r[i*p+k]) * f[i];

                s *= betaR[k];

                for(unsigned i = k; i < m; ++i)
                    f[i] -= s * r[i*p+k];
            }

            std::copy(f, f + n, x);

            for(unsigned k=n; k --> 0;) {
                x[k] = math::inverse(diagR[k]) * x[k];

                for(unsigned i = 0; i < k; ++i)
                    x[i] -= r[i*p+k] * x[k];
            }
        }
    private:
        typedef typename math::scalar_of<value_type>::type scalar_type;

        static scalar_type sqr(scalar_type x) { return x * x; }

        unsigned m, n, p;

        value_type *r;
        std::vector<value_type>  q;
        std::vector<value_type>  diagR;
        std::vector<scalar_type> betaR;
};

} // namespace detail
} // namespace amgcl

#endif
