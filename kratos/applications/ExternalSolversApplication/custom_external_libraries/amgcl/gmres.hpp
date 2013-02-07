#ifndef AMGCL_GMRES_HPP
#define AMGCL_GMRES_HPP

/*
The MIT License

Copyright (c) 2012-2013 Denis Demidov <ddemidov@ksu.ru>

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
 * \author Denis Demidov <ddemidov@ksu.ru>
 * \brief  GMRES method.
 *
 * Implementation is based on \ref Templates_1994 "Barrett (1994)"
 */

#include <vector>
#include <utility>
#include <algorithm>

#include <amgcl/common.hpp>

namespace amgcl {

namespace gmres_ops {

template <typename value_t>
void apply_plane_rotation(value_t &dx, value_t &dy, value_t cs, value_t sn) {
    value_t tmp = cs * dx + sn * dy;
    dy = -sn * dx + cs * dy;
    dx = tmp;
}

template <typename value_t>
void generate_plane_rotation(value_t dx, value_t dy, value_t &cs, value_t &sn) {
    if (dy == 0) {
        cs = 1;
        sn = 0;
    } else if (fabs(dy) > fabs(dx)) {
        value_t tmp = dx / dy;
        sn = 1 / sqrt(1 + tmp * tmp);
        cs = tmp * sn;
    } else {
        value_t tmp = dy / dx;
        cs = 1 / sqrt(1 + tmp * tmp);
        sn = tmp * cs;
    }
}

template <class GMRES, class Vector>
void update(GMRES &gmres, Vector &x, int k) {
    std::copy(gmres.s.begin(), gmres.s.end(), gmres.y.begin());

    for (int i = k; i >= 0; --i) {
        gmres.y[i] /= gmres.H[i * gmres.M + i];
        for (int j = i - 1; j >= 0; --j)
            gmres.y[j] -= gmres.H[j * gmres.M + i] * gmres.y[i];
    }

    for (int j = 0; j <= k; j++)
        x += gmres.y[j] * gmres.v[j];
}

template <class GMRES, class matrix, class Vector, class precond>
typename value_type<Vector>::type restart(GMRES &gmres,
        const matrix &A, const Vector &rhs, const precond &P, const Vector &x
        )
{
    typedef typename value_type<Vector>::type value_t;

    residual(A, x, rhs, gmres.w);
    clear(gmres.r);
    P.apply(gmres.w, gmres.r);

    gmres.s[0] = norm(gmres.r);
    gmres.v[0] = gmres.r / gmres.s[0];

    std::fill(gmres.s.begin() + 1, gmres.s.end(), value_t());

    return gmres.s[0];
}

template <class GMRES, class matrix, class precond>
typename GMRES::value_t iteration(
        GMRES &gmres, const matrix &A, const precond &P, int i)
{
    axpy(A, gmres.v[i], gmres.r);
    clear(gmres.w);
    P.apply(gmres.r, gmres.w);

    for(int k = 0; k <= i; ++k) {
        gmres.H[k * gmres.M + i] = inner_prod(gmres.w, gmres.v[k]);
        gmres.w -= gmres.H[k * gmres.M + i] * gmres.v[k];
    }

    gmres.H[(i+1) * gmres.M + i] = norm(gmres.w);

    gmres.v[i+1] = gmres.w / gmres.H[(i+1) * gmres.M + i];

    for(int k = 0; k < i; ++k)
        apply_plane_rotation(gmres.H[k * gmres.M + i], gmres.H[(k+1) * gmres.M + i], gmres.cs[k], gmres.sn[k]);

    generate_plane_rotation(gmres.H[i * gmres.M + i], gmres.H[(i+1) * gmres.M + i], gmres.cs[i], gmres.sn[i]);
    apply_plane_rotation(gmres.H[i * gmres.M + i], gmres.H[(i+1) * gmres.M + i], gmres.cs[i], gmres.sn[i]);
    apply_plane_rotation(gmres.s[i], gmres.s[i+1], gmres.cs[i], gmres.sn[i]);

    return fabs(gmres.s[i + 1]);
}

} // namespace gmres_ops

template <class Vector>
struct gmres_data {
    typedef typename value_type<Vector>::type value_t;

    int M;
    std::vector<value_t> H, s, cs, sn, y;
    Vector r, w;
    std::vector<Vector> v;

    gmres_data(int M, size_t n)
        : M(M), H(M * (M + 1)), s(M + 1), cs(M + 1), sn(M + 1), y(M + 1),
          r(n), w(n), v(M + 1)
    {
        for(BOOST_AUTO(vp, v.begin()); vp != v.end(); ++vp) vp->resize(n);
    }

    void update(Vector &x, int k) {
        gmres_ops::update(*this, x, k);
    }

    template <class matrix, class precond>
    value_t restart(
            const matrix &A, const Vector &rhs, const precond &P, const Vector &x
            )
    {
        return gmres_ops::restart(*this, A, rhs, P, x);
    }

    template <class matrix, class precond>
    value_t iteration(const matrix &A, const precond &P, int i) {
        return gmres_ops::iteration(*this, A, P, i);
    }
};

/// GMRES method.
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
solve(const matrix &A, const vector &rhs, const precond &P, vector &x, gmres_tag prm = gmres_tag())
{
    typedef typename value_type<vector>::type value_t;
    const size_t n = x.size();

    gmres_data<vector> gmres(prm.M, n);

    int     iter = 0;
    value_t res;

    do {
        res = gmres.restart(A, rhs, P, x);

        for(int i = 0; i < prm.M && iter < prm.maxiter; ++i, ++iter) {
            res = gmres.iteration(A, P, i);

	    if (res < prm.tol) {
                gmres.update(x, i);
		return std::make_pair(iter + 1, res);
	    };
	}

        gmres.update(x, prm.M - 1);
    } while (iter < prm.maxiter && res > prm.tol);

    return std::make_pair(iter, res);
}

} // namespace amgcl

#endif
