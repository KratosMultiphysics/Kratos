/*
    enoki/special.h -- Special functions: Bessel functions, Elliptic
    and exponential integrals, etc. (still incomplete)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <enoki/array.h>

#pragma once

NAMESPACE_BEGIN(enoki)

/// Evaluates a series of Chebyshev polynomials at argument x/2.
template <typename T, typename T2, size_t Size,
          typename Expr = expr_t<T>> Expr chbevl(const T &x, T2 (&coeffs)[Size]) {
    using Scalar = scalar_t<Expr>;

    Expr b0 = Scalar(coeffs[0]);
    Expr b1 = Scalar(0);
    Expr b2;

    ENOKI_UNROLL for (size_t i = 0; i < Size; ++i) {
        b2 = b1;
        b1 = b0;
        b0 = fmsub(x, b1, b2 - Scalar(coeffs[i]));
    }

    return (b0 - b2) * Scalar(0.5f);
}

template <typename T, enable_if_not_array_t<T> = 0> T erf(const T &x) {
    return std::erf(x);
}


template <typename T, enable_if_not_array_t<T> = 0> T erfc(const T &x) {
    return std::erfc(x);
}

template <typename T, bool Recurse = true, typename Expr = expr_t<T>,
          enable_if_array_t<T> = 0>
Expr erfc(const T &x);

template <typename T, bool Recurse = true, typename Expr = expr_t<T>,
          enable_if_array_t<T> = 0>
Expr erf(const T &x);

template <typename T, bool Recurse, typename Expr, enable_if_array_t<T>>
Expr erfc(const T &x) {
    constexpr bool Single = std::is_same_v<scalar_t<T>, float>;
    using Scalar = scalar_t<T>;

    Expr r;
    Expr xa = abs(x),
         z  = exp(-x*x);

    auto erf_mask   = xa < Scalar(1),
         large_mask = xa > Scalar(Single ? 2 : 8);

    ENOKI_MARK_USED(erf_mask);

    if constexpr (Single) {
        Expr q  = rcp(xa),
             y  = q*q, p_small, p_large;

        if (is_cuda_array_v<Expr> || !all_nested(large_mask))
            p_small = poly8(y, 5.638259427386472e-1, -2.741127028184656e-1,
                               3.404879937665872e-1, -4.944515323274145e-1,
                               6.210004621745983e-1, -5.824733027278666e-1,
                               3.687424674597105e-1, -1.387039388740657e-1,
                               2.326819970068386e-2);

        if (is_cuda_array_v<Expr> || any_nested(large_mask))
            p_large = poly7(y, 5.641895067754075e-1, -2.820767439740514e-1,
                               4.218463358204948e-1, -1.015265279202700e+0,
                               2.921019019210786e+0, -7.495518717768503e+0,
                               1.297719955372516e+1, -1.047766399936249e+1);
        r = z * q * select(large_mask, p_large, p_small);
    } else {
        Expr p_small, p_large, q_small, q_large;

        if (is_cuda_array_v<Expr> || !all_nested(large_mask)) {
            p_small = poly8(xa, 5.57535335369399327526e2, 1.02755188689515710272e3,
                                9.34528527171957607540e2, 5.26445194995477358631e2,
                                1.96520832956077098242e2, 4.86371970985681366614e1,
                                7.46321056442269912687e0, 5.64189564831068821977e-1,
                                2.46196981473530512524e-10);

            q_small = poly8(xa, 5.57535340817727675546e2, 1.65666309194161350182e3,
                                2.24633760818710981792e3, 1.82390916687909736289e3,
                                9.75708501743205489753e2, 3.54937778887819891062e2,
                                8.67072140885989742329e1, 1.32281951154744992508e1,
                                1.00000000000000000000e0);
        }


        if (is_cuda_array_v<Expr> || any_nested(large_mask)) {
            p_large = poly5(xa, 2.97886665372100240670e0, 7.40974269950448939160e0,
                                6.16021097993053585195e0, 5.01905042251180477414e0,
                                1.27536670759978104416e0, 5.64189583547755073984e-1);

            q_large = poly6(xa, 3.36907645100081516050e0, 9.60896809063285878198e0,
                                1.70814450747565897222e1, 1.20489539808096656605e1,
                                9.39603524938001434673e0, 2.26052863220117276590e0,
                                1.00000000000000000000e0);
        }

        r = (z * select(large_mask, p_large, p_small)) /
                 select(large_mask, q_large, q_small);

        r &= neq(z, zero<Expr>());
    }

    r[x < Scalar(0)] = Scalar(2) - r;

    if constexpr (Recurse) {
        if (ENOKI_UNLIKELY(is_cuda_array_v<Expr> || any_nested(erf_mask)))
            r[erf_mask] = Scalar(1) - erf<T, false>(x);
    }
    return r;
}

template <typename T, bool Recurse, typename Expr, enable_if_array_t<T>>
Expr erf(const T &x) {
    using Scalar = scalar_t<T>;

    Expr r;
    auto erfc_mask = abs(x) > Scalar(1);
    ENOKI_MARK_USED(erfc_mask);

    Expr z = x * x;

    constexpr bool Single = std::is_same_v<scalar_t<T>, float>;
    if constexpr (Single) {
        r = poly6(z, 1.128379165726710e+0, -3.761262582423300e-1,
                     1.128358514861418e-1, -2.685381193529856e-2,
                     5.188327685732524e-3, -8.010193625184903e-4,
                     7.853861353153693e-5);
    } else {
        r = poly4(z, 5.55923013010394962768e4, 7.00332514112805075473e3,
                     2.23200534594684319226e3, 9.00260197203842689217e1,
                     9.60497373987051638749e0) /
            poly5(z, 4.92673942608635921086e4, 2.26290000613890934246e4,
                     4.59432382970980127987e3, 5.21357949780152679795e2,
                     3.35617141647503099647e1, 1.00000000000000000000e0);
    }

    r *= x;

    if constexpr (Recurse) {
        if (ENOKI_UNLIKELY(is_cuda_array_v<Expr> || any_nested(erfc_mask)))
            r[erfc_mask] = Scalar(1) - erfc<T, false>(x);
    }

    return r;
}


/// Modified Bessel function of the first kind, order zero (exponentially scaled)
template <typename T, typename Expr = expr_t<T>> Expr i0e(const T &x_) {
    using Scalar = scalar_t<T>;

    /* Chebyshev coefficients for exp(-x) I0(x)
     * in the interval [0,8].
     *
     * lim(x->0) { exp(-x) I0(x) } = 1.
     */

    static Scalar A[] = {
        Scalar(-1.30002500998624804212E-8), Scalar(6.04699502254191894932E-8),
        Scalar(-2.67079385394061173391E-7), Scalar(1.11738753912010371815E-6),
        Scalar(-4.41673835845875056359E-6), Scalar(1.64484480707288970893E-5),
        Scalar(-5.75419501008210370398E-5), Scalar(1.88502885095841655729E-4),
        Scalar(-5.76375574538582365885E-4), Scalar(1.63947561694133579842E-3),
        Scalar(-4.32430999505057594430E-3), Scalar(1.05464603945949983183E-2),
        Scalar(-2.37374148058994688156E-2), Scalar(4.93052842396707084878E-2),
        Scalar(-9.49010970480476444210E-2), Scalar(1.71620901522208775349E-1),
        Scalar(-3.04682672343198398683E-1), Scalar(6.76795274409476084995E-1)
    };


    /* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
     * in the inverted interval [8,infinity].
     *
     * lim(x->inf) { exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
     */

    static Scalar B[] = {
        Scalar(3.39623202570838634515E-9), Scalar(2.26666899049817806459E-8),
        Scalar(2.04891858946906374183E-7), Scalar(2.89137052083475648297E-6),
        Scalar(6.88975834691682398426E-5), Scalar(3.36911647825569408990E-3),
        Scalar(8.04490411014108831608E-1)
    };


    Expr x = abs(x_);

    auto mask_big = x > Scalar(8);

    Expr r_big, r_small;

    if (!all_nested(mask_big))
        r_small = chbevl(fmsub(x, Expr(Scalar(0.5)), Expr(Scalar(2))), A);

    if (any_nested(mask_big))
        r_big = chbevl(fmsub(Expr(Scalar(32)), rcp(x), Expr(Scalar(2))), B) *
                rsqrt(x);

    return select(mask_big, r_big, r_small);
}

// Inverse real error function approximation based on on "Approximating the
// erfinv function" by Mark Giles
template <typename T, typename Expr = expr_t<T>> Expr erfinv(const T &x_) {
    using Scalar = scalar_t<T>;

    Expr x(x_);
    Expr w = -log((Expr(Scalar(1)) - x) * (Expr(Scalar(1)) + x));

    Expr w1 = w - Scalar(2.5);
    Expr w2 = sqrt(w) - Scalar(3);

    Expr p1 = poly8(w1,
         1.50140941,     0.246640727,
        -0.00417768164, -0.00125372503,
         0.00021858087, -4.39150654e-06,
        -3.5233877e-06,  3.43273939e-07,
         2.81022636e-08);

    Expr p2 = poly8(w2,
         2.83297682,     1.00167406,
         0.00943887047, -0.0076224613,
         0.00573950773, -0.00367342844,
         0.00134934322,  0.000100950558,
        -0.000200214257);

    return select(w < Scalar(5), p1, p2) * x;
}

/// Evaluates Dawson's integral (e^(-x^2) \int_0^x e^(y^2) dy)
template <typename T, typename Expr = expr_t<T>> Expr dawson(const T &x) {
    // Rational minimax approximation to Dawson's integral with relative
    // error < 1e-6 on the real number line. July 2017, Wenzel Jakob

    Expr x2 = x*x;
    Expr num = poly6(x2, 1.00000080272429,9.18170212243285e-2,
                         4.25835373536124e-2, 6.0536496345054e-3,
                         9.88555033724111e-4, 3.64943550840577e-5,
                         1.55942290996993e-5);

    Expr denom = poly7(x2, 1.0, 7.58517175815194e-1,
                           2.81364355593059e-1, 6.81783097841267e-2,
                           1.13586116798019e-2, 1.92020805811771e-3,
                           5.74217664074868e-5, 3.11884331363595e-5);

    return num / denom * x;
}

/// Imaginary component of the error function
template <typename T, typename Expr = expr_t<T>> Expr erfi(const T &x) {
    using Scalar = scalar_t<T>;

    return Scalar(M_2_SQRTPI) * dawson(x) * exp(x * x);
}

/// Natural logarithm of the Gamma function
template <typename Value> Value lgamma(Value x_) {
    using Mask = mask_t<Value>;
    using Scalar = scalar_t<Value>;

    // 'g' and 'n' parameters of the Lanczos approximation
    // See mrob.com/pub/ries/lanczos-gamma.html
    const int n = 6;
    const Scalar g = 5.0f;
    const Scalar log_sqrt2pi = Scalar(0.91893853320467274178);
    const Scalar coeff[n + 1] = { (Scalar)  1.000000000190015, (Scalar) 76.18009172947146,
                                  (Scalar) -86.50532032941677, (Scalar) 24.01409824083091,
                                  (Scalar) -1.231739572450155, (Scalar) 0.1208650973866179e-2,
                                  (Scalar) -0.5395239384953e-5 };

    // potentially reflect using gamma(x) = pi / (sin(pi*x) * gamma(1-x))
    Mask reflect = x_ < .5f;

    Value x = select(reflect, -x_, x_ - 1.f),
          b = x + g + .5f; // base

    Value sum = 0;
    for (int i = n; i >= 1; --i)
        sum += coeff[i] / (x + Scalar(i));
    sum += coeff[0];

    // gamma(x) = sqrt(2*pi) * sum * b^(x + .5) / exp(b)
    Value result = ((log_sqrt2pi + log(sum)) - b) + log(b) * (x + .5f);

    if (is_cuda_array_v<Value> || any_nested(reflect)) {
        masked(result, reflect) = log(abs(Scalar(M_PI) / sin(Scalar(M_PI) * x_))) - result;
        masked(result, reflect && eq(x_, round(x_))) = std::numeric_limits<Scalar>::infinity();
    }

    return result;
}

/// Gamma function
template <typename Value> Value tgamma(Value x) { return exp(lgamma(x)); }

/**
 * Computes a Carlson integral of the form
 *
 * R_F(X, Y, Z) = 1/2 * \int_{0}^\infty ((t + x) (t + y) (t + z))^(-1/2) dt
 *
 * Based on
 *
 *   Computing elliptic integrals by duplication
 *   B. C. Carlson
 *   Numerische Mathematik, March 1979, Volume 33, Issue 1
 */
template <typename Vector3,
          typename Value = value_t<Vector3>,
          typename Scalar = scalar_t<Vector3>>
Value carlson_rf(Vector3 xyz) {
    static_assert(
        Vector3::Size == 3,
        "carlson_rf(): Expected a three-dimensional input vector (x, y, z)");
    assert(all_nested(xyz.x() >= Scalar(0) && xyz.y() > Scalar(0) && xyz.z() > Scalar(0)));

    Vector3 XYZ;
    Value mu_inv;
    mask_t<Value> active = true;
    int iterations = 0;

    while (true) {
        Vector3 sqrt_xyz = sqrt(xyz);
        Value lambda = dot(shuffle<1, 2, 0>(sqrt_xyz), sqrt_xyz);
        Value mu = hsum(xyz) * Scalar(1.0 / 3.0);
        mu_inv = rcp(mu);
        XYZ = fnmadd(xyz, mu_inv, Scalar(1));
        Value eps = hmax(abs(XYZ));
        active &= eps > Scalar(std::is_same_v<Scalar, double>
                                   ? 0.0024608
                                   : 0.070154); // eps ^ (1/6)

        if (none(active) || ++iterations == 10)
            break;

        xyz[mask_t<Vector3>(active)] = (xyz + lambda) * Scalar(0.25);
    }

    /* Use recurrences for cheaper polynomial evaluation. Based
       on Numerical Recipes (3rd ed) by Press, Teukolsky,
       Vetterling, and Flannery */

    Value e2 = XYZ.x() * XYZ.y() - XYZ.z() * XYZ.z(),
          e3 = hprod(XYZ),
          er = (Scalar(1.0 / 24.0) * e2 - Scalar(1.0 / 10.0) -
                Scalar(3.0 / 44.0) * e3) * e2 + Scalar(1.0 / 14.0) * e3;

    return sqrt(mu_inv) * (Scalar(1) + er);
}

/**
 * Computes a Carlson integral of the form
 *
 * R_D(x, y, z) = 3/2 * \int_{0}^\infty (t + x)^(-1/2) (t + y)^(-1/2) (t + z)^(-3/2) dt
 *
 * Based on
 *
 *   Computing elliptic integrals by duplication
 *   B. C. Carlson
 *   Numerische Mathematik, March 1979, Volume 33, Issue 1
 */
template <typename Vector3,
          typename Value = value_t<Vector3>,
          typename Scalar = scalar_t<Vector3>>
Value carlson_rd(Vector3 xyz) {
    static_assert(
        Vector3::Size == 3,
        "carlson_rd(): Expected a three-dimensional input vector (x, y, z)");
    assert(all_nested(xyz.x() >= Scalar(0) && xyz.y() > Scalar(0) && xyz.z() > Scalar(0)));

    Vector3 XYZ;
    Value mu_inv;
    mask_t<Value> active = true;
    int iterations = 0;
    Value sum = 0;
    Value num = 1;
    const Vector3 W(Scalar(1.0 / 5.0), Scalar(1.0 / 5.0), Scalar(3.0 / 5.0));

    while (true) {
        Vector3 sqrt_xyz = sqrt(xyz);
        Value lambda = dot(shuffle<1, 2, 0>(sqrt_xyz), sqrt_xyz);
        Value mu = hsum(xyz * W);
        mu_inv = rcp(mu);
        XYZ = fnmadd(xyz, mu_inv, Scalar(1));
        Value eps = hmax(abs(XYZ));
        active &= eps > Scalar(std::is_same_v<Scalar, double>
                                   ? (0.0024608 * 0.6)
                                   : (0.070154 * 0.6)); // eps ^ (1/6) * 0.6

        if (none(active) || ++iterations == 10)
            break;

        masked(sum, active) += num / (sqrt(xyz.z()) * (xyz.z() + lambda));
        masked(num, active) *= Scalar(0.25f);
        masked(xyz, mask_t<Vector3>(active)) = (xyz + lambda) * Scalar(0.25f);
    }

    /* Use recurrences for cheaper polynomial evaluation. Based
       on Numerical Recipes (3rd ed) by Press, Teukolsky,
       Vetterling, and Flannery */

    Value z  = XYZ.z(),
          ea = XYZ.x() * XYZ.y(),
          eb = z * z,
          ec = ea - eb,
          ed = fnmadd(Scalar(6), eb, ea),
          ee = fmadd(ec, Scalar(2), ed);

    Value p = ed * (-Scalar(3.0 / 14.0) + Scalar(9.0 / 88.0) * ed -
                    Scalar(1.0 / 4.0) * z * ee) +
              z * (Scalar(1.0 / 6.0) * ee + z *
                    (-Scalar(9.0 / 22.0) * ec + z * Scalar(3.0 / 26.0) * ea));

    return Scalar(3) * sum + num * mu_inv * sqrt(mu_inv) * (Scalar(1.0) + p);
}

/**
 * Computes a Carlson integral of the form
 *
 * R_C(x, y) = 1/2 * \int_{0}^\infty (t + x)^(-1/2) (t + y)^-1 dt
 *
 * Based on
 *
 *   Computing elliptic integrals by duplication
 *   B. C. Carlson
 *   Numerische Mathematik, March 1979, Volume 33, Issue 1
 */
template <typename Vector2,
          typename Value = value_t<Vector2>,
          typename Scalar = scalar_t<Vector2>>
Value carlson_rc(Vector2 xy) {
    static_assert(
        Vector2::Size == 2,
        "carlson_rc(): Expected a two-dimensional input vector (x, y)");
    assert(all(xy.x() >= Scalar(0) && xy.y() > Scalar(0)));

    mask_t<Value> active = true;
    Value inv_mu, s;
    int iterations = 0;

    while (true) {
        Value lambda = hprod(sqrt(xy));
        lambda += lambda + xy.y();
        Value mu = fmadd(xy.x(), Scalar(1.0 / 3.0), xy.y() * Scalar(2.0 / 3.0));
        inv_mu = rcp(mu);
        s = (xy.y() - mu) * inv_mu;

        active &= abs(s) > Scalar(std::is_same_v<Scalar, double>
                                   ? (0.0024608 * 0.48)
                                   : (0.070154 * 0.48)); // eps ^ (1/6) * 0.48

        if (none(active) || ++iterations == 10)
            break;

        masked(xy, mask_t<Vector2>(active)) = (xy + lambda) * Scalar(0.25f);
    }

    /* Use recurrences for cheaper polynomial evaluation. Based
       on Numerical Recipes (3rd ed) by Press, Teukolsky,
       Vetterling, and Flannery */

    return sqrt(inv_mu) * (Scalar(1) + s * s *
              (Scalar(0.3) + s * (Scalar(1.0 / 7.0) +
               s * (Scalar(0.375) + s * Scalar(9.0 / 22.0)))));
}

/**
 * Computes a Carlson integral of the form
 *
 * R_J(x, y, z, rho) = 3/2 * \int_{0}^\infty ((t + x) (t + y) (t + z))^(-1/2) (t+rho)^(-1) dt
 *
 * Based on
 *
 *   Computing elliptic integrals by duplication
 *   B. C. Carlson
 *   Numerische Mathematik, March 1979, Volume 33, Issue 1
 */
template <typename Vector4,
          typename Value = value_t<Vector4>,
          typename Vector2 = Array<Value, 2>,
          typename Scalar = scalar_t<Vector4>>
Value carlson_rj(Vector4 xyzr) {
    static_assert(
        Vector4::Size == 4,
        "carlson_rj(): Expected a four-dimensional input vector (x, y, z, rho)");
    assert(all(xyzr.x() >= Scalar(0) && xyzr.y() > Scalar(0) && xyzr.z() > Scalar(0) && xyzr.w() > Scalar(0)));

    Vector4 XYZR;
    Value mu_inv;
    mask_t<Value> active = true;
    int iterations = 0;
    Value sum = 0;
    Value num = 1;

    while (true) {
        auto xyz = head<3>(xyzr);
        auto rho = xyzr.w();
        auto sqrt_xyz = sqrt(xyz);
        Value lambda = dot(shuffle<1, 2, 0>(sqrt_xyz), sqrt_xyz);

        Value mu = (hsum(xyzr) + rho) * Scalar(1.0 / 5.0);
        mu_inv = rcp(mu);
        XYZR = fnmadd(xyzr, mu_inv, Scalar(1));
        Value eps = hmax(abs(XYZR));
        active &= eps > Scalar(std::is_same_v<Scalar, double>
                                   ? (0.0024608 * 0.6)
                                   : (0.070154 * 0.6)); // eps ^ (1/6) * 0.6

        Value alpha = rho * hsum(sqrt(xyz)) + sqrt(hprod(xyz));
        alpha *= alpha;
        Value beta = rho * (rho + lambda) * (rho + lambda);

        if (none(active) || ++iterations == 10)
            break;

        masked(sum, active) += num * carlson_rc(Vector2(alpha, beta));
        masked(num, active) *= Scalar(0.25f);
        masked(xyzr, mask_t<Vector4>(active)) = (xyzr + lambda) * Scalar(0.25f);
    }

    /* Use recurrences for cheaper polynomial evaluation. Based
       on Numerical Recipes (3rd ed) by Press, Teukolsky,
       Vetterling, and Flannery */

    Value ea = XYZR.x() * (XYZR.y() + XYZR.z()) + XYZR.y() * XYZR.z(),
          eb = XYZR.x() * XYZR.y() * XYZR.z(),
          R  = XYZR.w(),
          ec = R * R,
          ed = ea - Scalar(3) * ec,
          ee = eb + Scalar(2) * R * (ea - ec);

    return Scalar(3) * sum +
           num * mu_inv * sqrt(mu_inv) *
               (Scalar(1) +
                ed * (-Scalar(3.0 / 14.0) + Scalar(9.0 / 88.0) * ed -
                      Scalar(9.0 / 52.0) * ee) +
                eb * (Scalar(1.0 / 6.0) +
                      R * (-Scalar(3.0 / 11.0) + R * Scalar(3.0 / 26.0))) +
                R * ea * (Scalar(1.0 / 3.0) - R * Scalar(3.0 / 22.0)) -
                Scalar(1.0 / 3.0) * R * ec);
}

// -----------------------------------------------------------------------
//! @{ \name Complete and incomplete elliptic integrals
//! Caution: the 'k' factor is squared in the elliptic integral, which
//! differs from the convention of Mathematica's EllipticK etc.
// -----------------------------------------------------------------------

/// Complete elliptic integral of the first kind
template <typename K, typename Value = expr_t<K>,
          typename Scalar = scalar_t<Value>,
          typename Vector3 = Array<Value, 3>>
Value comp_ellint_1(K k) {
    return carlson_rf(Vector3(Scalar(0), Scalar(1) - k * k, Scalar(1)));
}


/// Incomplete elliptic integral of the first kind
template <typename Phi, typename K,
          typename Value = expr_t<Phi, K>,
          typename Scalar = scalar_t<Value>,
          typename Vector3 = Array<Value, 3>>
Value ellint_1(Phi phi_, K k) {
    Value phi = phi_,
          n = floor(fmadd(phi, Scalar(1.0 / M_PI), Scalar(.5f))),
          result = 0;

    if (ENOKI_UNLIKELY(any(neq(n, Scalar(0))))) {
        result = comp_ellint_1(k) * n * Scalar(2);
        phi = fnmadd(n, Scalar(M_PI), phi);
    }

    auto [sin_phi, cos_phi] = sincos(phi);
    Vector3 xyz(cos_phi * cos_phi, Scalar(1) - k * k * sin_phi * sin_phi,
                Scalar(1));
    result += sin_phi * carlson_rf(xyz);

    return result;
}

/// Complete elliptic integral of the second kind
template <typename K, typename Value = expr_t<K>,
          typename Scalar = scalar_t<Value>,
          typename Vector3 = Array<Value, 3>>
Value comp_ellint_2(K k) {
    auto k2 = k*k;
    Vector3 xyz(Scalar(0), Scalar(1) - k2, Scalar(1));
    return carlson_rf(xyz) - Scalar(1.0 / 3.0) * k2 * carlson_rd(xyz);
}

/// Incomplete elliptic integral of the second kind
template <typename Phi, typename K,
          typename Value = expr_t<Phi, K>,
          typename Scalar = scalar_t<Value>,
          typename Vector3 = Array<Value, 3>>
Value ellint_2(Phi phi_, K k) {
    Value phi = phi_,
          k2 = k*k,
          n = floor(fmadd(phi, Scalar(1.0 / M_PI), Scalar(.5f))),
          result = 0;

    if (ENOKI_UNLIKELY(any(neq(n, Scalar(0))))) {
        result = comp_ellint_2(k) * n * Scalar(2);
        phi = fnmadd(n, Scalar(M_PI), phi);
    }

    auto [sin_phi, cos_phi] = sincos(phi);
    auto sin_phi_k_2 = sin_phi * sin_phi * k2;
    Vector3 xyz(cos_phi * cos_phi, Scalar(1) - sin_phi_k_2, Scalar(1));
    result += sin_phi * (carlson_rf(xyz) -
                         Scalar(1.0 / 3.0) * sin_phi_k_2 * carlson_rd(xyz));

    return result;
}

/// Complete elliptic integral of the third kind
template <typename K, typename Nu,
          typename Value = expr_t<K, Nu>,
          typename Scalar = scalar_t<Value>,
          typename Vector4 = Array<Value, 4>>
Value comp_ellint_3(K k, Nu nu) {
    auto k2 = k*k;
    Vector4 xyzr(Scalar(0), Scalar(1) - k2, Scalar(1), Scalar(1) + nu);
    return carlson_rf(head<3>(xyzr)) -
           Scalar(1.0 / 3.0) * nu * carlson_rj(xyzr);
}

/// Incomplete elliptic integral of the third kind
template <typename Phi, typename K, typename Nu,
          typename Value = expr_t<Phi, K, Nu>,
          typename Scalar = scalar_t<Value>,
          typename Vector4 = Array<Value, 4>>
Value ellint_3(Phi phi_, K k, Nu nu) {
    Value phi = phi_,
          k2 = k*k,
          n = floor(fmadd(phi, Scalar(1.0 / M_PI), Scalar(.5f))),
          result = 0;

    if (ENOKI_UNLIKELY(any(neq(n, Scalar(0))))) {
        result = comp_ellint_3(k, nu) * n * Scalar(2);
        phi = fnmadd(n, Scalar(M_PI), phi);
    }


    auto [sin_phi, cos_phi] = sincos(phi);
    auto sin_phi_2 = sin_phi * sin_phi;
    Vector4 xyzr(cos_phi * cos_phi, Scalar(1) - k2 * sin_phi_2, Scalar(1),
                 Scalar(1) + nu * sin_phi_2);
    result += sin_phi * (carlson_rf(head<3>(xyzr)) -
                         Scalar(1.0 / 3.0) * nu * sin_phi_2 * carlson_rj(xyzr));

    return result;
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
