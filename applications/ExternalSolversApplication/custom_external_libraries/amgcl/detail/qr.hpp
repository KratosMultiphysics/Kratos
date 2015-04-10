#ifndef AMGCL_DETAIL_QR_HPP
#define AMGCL_DETAIL_QR_HPP

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
 * \file   amgcl/detail/qr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  QR decomposition.
 */

#include <vector>
#include <cmath>
#include <amgcl/util.hpp>

namespace amgcl {
namespace detail {

/// In-place QR decomposition.
template <typename value_type>
class QR {
    public:
        QR() : m(0), n(0) {}

        void compute(
                size_t rows, size_t cols, value_type *A,
                bool needQ = true
                )
        {
            precondition(rows >= cols,
                    "QR decomposion of wide matrices is not supported"
                    );

            m = rows;
            n = cols;
            r = A;

            resize(d, n);

            for(size_t j = 0; j < n; ++j) {
                value_type s = 0;
                for(size_t i = j; i < m; ++i) s += r[i*n+j] * r[i*n+j];
                s = sqrt(s);
                d[j] = r[j*n+j] > 0 ? -s : s;
                value_type fak = sqrt(s * (s + fabs(r[j*n+j])));
                r[j*n+j] -= d[j];
                for(size_t k = j; k < m; ++k) r[k*n+j] /= fak;
                for(size_t i = j + 1; i < n; ++i) {
                    value_type s = 0;
                    for(size_t k = j; k < m; ++k) s += r[k*n+j] * r[k*n+i];
                    for(size_t k = j; k < m; ++k) r[k*n+i] -= r[k*n+j] * s;
                }
            }

            sign = r[0] < 0 ? -1 : 1;

            if (needQ) {
                resize(q, n * m);
                for(size_t i = 0; i < n; ++i) {
                    for(size_t j = 0; j < m; ++j) q[j*n+i] = (j == i);

                    for(size_t j = n; j-- > 0; ) {
                        value_type s = 0;

                        for(size_t k = j; k < m; ++k) s += r[k*n+j] * q[k*n+i];
                        for(size_t k = j; k < m; ++k) q[k*n+i] -= r[k*n+j] * s;
                    }
                }
            }
        }

        value_type R(size_t i, size_t j) const {
            if (i > j)
                return value_type(0);
            if (i == j)
                return sign * d[i];
            else
                return sign * r[i*n+j];
        }

        value_type Q(size_t i, size_t j) const {
            return sign * q[i*n+j];
        }

        void solve(value_type *f, value_type *x) const {
            // x = Q' * x;
            for(size_t j = 0; j < n; ++j) {
                value_type s = 0;

                for(size_t k = j; k < m; ++k) s += r[k*n+j] * f[k];
                for(size_t k = j; k < m; ++k) f[k] -= r[k*n+j] * s;
            }

            // x = R^-1 x
            for (size_t i = n; i --> 0; ) {
                value_type sum = x[i] = f[i];
                for (size_t j = i + 1; j < n; j++) sum -= r[i*n+j] * x[j];
                x[i] = sum / d[i];
            }
        }
    private:
        size_t m, n;
        value_type sign;

        value_type *r;

        std::vector<value_type> d;
        std::vector<value_type> q;

        void resize(std::vector<value_type> &v, size_t size) {
            if (v.size() < size) v.resize(size);
        }
};

} // namespace detail
} // namespace amgcl

#endif
