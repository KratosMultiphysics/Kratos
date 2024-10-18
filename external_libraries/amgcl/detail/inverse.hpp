#ifndef AMGCL_DETAIL_INVERSE_HPP
#define AMGCL_DETAIL_INVERSE_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

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

#include <algorithm>
#include <cassert>
#include <numeric>
#include <utility>
#include <amgcl/util.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace detail {

    template <typename value_type>
    static void inverse(int n, value_type *A, value_type *t, int *p) {
        std::iota(p, p + n, 0);

        // Perform LU-factorization of A in-place
        for(int col = 0; col < n; ++col) {
            // Find pivot element
            int pivot_i = col;
            using mag_type = typename math::scalar_of<value_type>::type;
            mag_type pivot_mag = math::zero<mag_type>();
            for (int i = col; i < n; ++i) {
                int row = p[i];
                mag_type mag = math::norm(A[row*n+col]);
                if (mag > pivot_mag) {
                    pivot_mag = mag;
                    pivot_i = i;
                }
            }
            std::swap(p[col], p[pivot_i]);
            int pivot_row = p[col];
            // We have found pivot element, perform Gauss elimination
            value_type d = math::inverse(A[pivot_row*n+col]);
            assert(!math::is_zero(d));
            for (int i = col+1; i < n; ++i) {
                int row = p[i];
                A[row*n+col] *= d;
                for(int j = col+1; j < n; ++j)
                    A[row*n+j] -= A[row*n+col] * A[pivot_row*n+j];
            }
            A[pivot_row*n+col] = d;
        }

        // Invert identity matrix in-place to get the solution.
        for(int k = 0; k < n; ++k) {
            // Lower triangular solve:
            for(int i = 0; i < n; ++i) {
                int row = p[i];
                value_type b = (row == k) ? math::identity<value_type>() : math::zero<value_type>();
                for(int j = 0; j < i; ++j)
                    b -= A[row*n+j] * t[j*n+k];
                t[i*n+k] = b;
            }

            // Upper triangular solve:
            for(int i = n; i --> 0; ) {
                int row = p[i];
                for(int j = i+1; j < n; ++j)
                    t[i*n+k] -= A[row*n+j] * t[j*n+k];
                t[i*n+k] *= A[row*n+i];
            }
        }

        std::copy(t, t + n * n, A);
    }

} // namespace detail
} // namespace amgcl

#endif
