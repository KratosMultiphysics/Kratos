#ifndef AMGCL_CG_HPP
#define AMGCL_CG_HPP

/*
The MIT License

Copyright (c) 2012 Denis Demidov <ddemidov@ksu.ru>

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
 * \file   cg.hpp
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  Conjugate Gradient method.
 *
 * Implementation is based on \ref Templates_1994 "Barrett (1994)"
 */

#include <utility>

#include <amgcl/common.hpp>

namespace amgcl {

/**
 * \defgroup iterative Iterative solvers
 * \brief List of provided iterative solvers.
 */

/// Conjugate Gradient method.
/**
 * Implementation is based on \ref Templates_1994 "Barrett (1994)"
 *
 * \param A   The system matrix.
 * \param rhs The right-hand side.
 * \param P   The preconditioner. Should provide apply(rhs, x) method.
 * \param x   The solution. Contains an initial approximation on input, and
 *            the approximated solution on output.
 * \param prm The control parameters.
 *
 * \returns a pair containing number of iterations made and precision
 * achieved.
 *
 * \ingroup iterative
 */
template <class matrix, class vector, class precond>
std::pair< int, typename value_type<vector>::type >
solve(const matrix &A, const vector &rhs, const precond &P, vector &x, cg_tag prm = cg_tag())
{
    typedef typename value_type<vector>::type value_t;

    const size_t n = x.size();

    vector r(n), s(n), p(n), q(n);
    residual(A, x, rhs, r);

    value_t rho1 = 0, rho2 = 0;
    value_t norm_of_rhs = norm(rhs);

    if (norm_of_rhs == 0) {
        clear(x);
        return std::make_pair(0, norm_of_rhs);
    }

    int     iter = 0;
    value_t res;

    for(; (res = norm(r) / norm_of_rhs) > prm.tol && iter < prm.maxiter; ++iter)
    {
        clear(s);
        P.apply(r, s);

        rho2 = rho1;
        rho1 = inner_prod(r, s);

        if (iter)
            p = s + (rho1 / rho2) * p;
        else
            p = s;

        axpy(A, p, q);

        value_t alpha = rho1 / inner_prod(q, p);

        x += alpha * p;
        r -= alpha * q;
    }

    return std::make_pair(iter, res);
}

} // namespace amgcl

#endif
