#ifndef AMGCL_DETAIL_INVERSE_HPP
#define AMGCL_DETAIL_INVERSE_HPP

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
 * \file   amgcl/detail/inverse.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Inverse of a dense matrix.
 */

#include <vector>
#include <algorithm>
#include <amgcl/util.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace detail {

    template <typename value_type>
    static void inverse(int n, value_type *A, value_type *t) {
        // Perform LU-factorization of A in-place
        for(int k = 0; k < n; ++k) {
            value_type d = math::inverse(A[k*n+k]);
            assert(!math::is_zero(d));
            for(int i = k+1; i < n; ++i) {
                A[i*n+k] *= d;
                for(int j = k+1; j < n; ++j)
                    A[i*n+j] -= A[i*n+k] * A[k*n+j];
            }
            A[k*n+k] = d;
        }

        // Invert identity matrix in-place to get the solution.
        for(int k = 0; k < n; ++k) {
            // Lower triangular solve:
            for(int i = 0; i < n; ++i) {
                value_type b = (i == k) ? math::identity<value_type>() : math::zero<value_type>();
                for(int j = 0; j < i; ++j)
                    b -= A[i*n+j] * t[j*n+k];
                t[i*n+k] = b;
            }

            // Upper triangular solve:
            for(int i = n; i --> 0; ) {
                for(int j = i+1; j < n; ++j)
                    t[i*n+k] -= A[i*n+j] * t[j*n+k];
                t[i*n+k] *= A[i*n+i];
            }
        }

        std::copy(t, t + n * n, A);
    }

} // namespace detail
} // namespace amgcl

#endif
