#ifndef AMGCL_CHEBYSHEV_HPP
#define AMGCL_CHEBYSHEV_HPP

/*
The MIT License

Copyright (c) 2012-2014 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   chebyshev.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Find coefficients of a chebyshev polynomial.
 */

#include <vector>
#include <stdexcept>
#include <cmath>
#include <boost/multi_array.hpp>
#include <boost/math/constants/constants.hpp>
#include <amgcl/spmat.hpp>

namespace amgcl {

/// Returns negated coefficients of a Chebyshev polynomial of given \p degree on interval [\p a, \p b].
/**
 * \note The last coefficient (at the 0-th order term) is not included.
 */
template <typename real>
std::vector<real> chebyshev_coefficients(unsigned degree, real a, real b) {
    assert(a > 0 && b > a &&
        "Invalid interval for chebyshev_coefficients");

    std::vector<real> C(degree);

    // Chebyshev polynomial roots on the interval [a, b].
    std::vector<real> roots(degree);
    for(unsigned i = 0; i < degree; ++i) {
        using boost::math::constants::pi;
        using boost::math::constants::half;

        roots[i] = a + half<real>() * ( b - a ) * (
                1 + cos( pi<real>() * ( i + half<real>() ) / degree ) );
    }

    // Construct linear system to determine Chebyshev coefficients.
    boost::multi_array<real, 2> A(boost::extents[degree][degree]);
    std::vector<real> rhs(degree);
    for(unsigned i = 0; i < degree; ++i) {
        real x = roots[i];
        real x_to_j = 1;
        for(unsigned j = 0; j < degree; ++j) {
            A[i][j] = x_to_j;
            x_to_j *= x;
        }
        rhs[i] = -x_to_j;
    }

    // Invert A, compute coefficients.
    sparse::gaussj(degree, A.data());

    real const_c;
    for(unsigned i = 0; i < degree; ++i) {
        real c = 0;
        for(unsigned j = 0; j < degree; ++j)
            c += A[i][j] * rhs[j];

        if (i == 0)
            const_c = c;
        else
            C[degree - i] = -c / const_c;
    }
    C[0] = -1 / const_c;

    return C;
}

/// Estimate the spectral radius of a given matric \p A by application of Gershgorin theorem.
template <typename spmat>
typename sparse::matrix_value<spmat>::type spectral_radius(const spmat &A) {
    typedef typename sparse::matrix_index<spmat>::type index_t;
    typedef typename sparse::matrix_value<spmat>::type value_t;

    const index_t n   = sparse::matrix_rows(A);

    const index_t *Arow = sparse::matrix_outer_index(A);
    const index_t *Acol = sparse::matrix_inner_index(A);
    const value_t *Aval = sparse::matrix_values(A);

    bool first = true;
    value_t emax;

    for(index_t i = 0; i < n; ++i) {
        value_t hi = 0;

        for(index_t j = Arow[i]; j < Arow[i + 1]; ++j)
            hi += std::fabs( Aval[j] );

        if (first) {
            emax = hi;
            first = false;
        } else {
            emax = std::max(emax, hi);
        }
    }

    return emax;
}


} // namespace amgcl

#endif
