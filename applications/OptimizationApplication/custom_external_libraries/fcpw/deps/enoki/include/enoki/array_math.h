/*
    enoki/array_math.h -- Mathematical support library

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <enoki/array_generic.h>

#pragma once

NAMESPACE_BEGIN(enoki)

// -----------------------------------------------------------------------
//! @{ \name Polynomial evaluation with short dependency chains and
//           fused multply-adds based on Estrin's scheme
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T = expr_t<T1>, typename S = scalar_t<T1>>
ENOKI_INLINE T poly2(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2) {
    T x2 = x * x;
    return fmadd(x2, S(c2), fmadd(x, S(c1), S(c0)));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly3(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3) {
    T x2 = x * x;
    return fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0)));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly4(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4) {
    T x2 = x * x, x4 = x2 * x2;
    return fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0)) + S(c4) * x4);
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly5(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4, const T2 &c5) {
    T x2 = x * x, x4 = x2 * x2;
    return fmadd(x2, fmadd(x, S(c3), S(c2)),
                     fmadd(x4, fmadd(x, S(c5), S(c4)), fmadd(x, S(c1), S(c0))));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly6(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4, const T2 &c5, const T2 &c6) {
    T x2 = x * x, x4 = x2 * x2;
    return fmadd(x4, fmadd(x2, S(c6), fmadd(x, S(c5), S(c4))),
                     fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0))));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly7(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4, const T2 &c5, const T2 &c6,
                     const T2 &c7) {
    T x2 = x * x, x4 = x2 * x2;
    return fmadd(x4, fmadd(x2, fmadd(x, S(c7), S(c6)), fmadd(x, S(c5), S(c4))),
                     fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0))));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly8(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4, const T2 &c5, const T2 &c6,
                     const T2 &c7, const T2 &c8) {
    T x2 = x * x, x4 = x2 * x2, x8 = x4 * x4;
    return fmadd(x4, fmadd(x2, fmadd(x, S(c7), S(c6)), fmadd(x, S(c5), S(c4))),
                     fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0)) + S(c8) * x8));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly9(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                     const T2 &c3, const T2 &c4, const T2 &c5, const T2 &c6,
                     const T2 &c7, const T2 &c8, const T2 &c9) {
    T x2 = x * x, x4 = x2 * x2, x8 = x4 * x4;
    return fmadd(x8, fmadd(x, S(c9), S(c8)),
                     fmadd(x4, fmadd(x2, fmadd(x, S(c7), S(c6)), fmadd(x, S(c5), S(c4))),
                               fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0)))));
}

template <typename T1, typename T2, typename T = expr_t<T1>,
          typename S = scalar_t<T1>>
ENOKI_INLINE T poly10(const T1 &x, const T2 &c0, const T2 &c1, const T2 &c2,
                      const T2 &c3, const T2 &c4, const T2 &c5, const T2 &c6,
                      const T2 &c7, const T2 &c8, const T2 &c9, const T2 &c10) {
    T x2 = x * x, x4 = x2 * x2, x8 = x4 * x4;
    return fmadd(x8, fmadd(x2, S(c10), fmadd(x, S(c9), S(c8))),
                     fmadd(x4, fmadd(x2, fmadd(x, S(c7), S(c6)), fmadd(x, S(c5), S(c4))),
                               fmadd(x2, fmadd(x, S(c3), S(c2)), fmadd(x, S(c1), S(c0)))));
}

//! @}
// -----------------------------------------------------------------------

#define ENOKI_UNARY_OPERATION(name, expr)                                      \
    namespace detail {                                                         \
        template <typename T>                                                  \
        using has_##name = decltype(std::declval<T>().name##_());              \
        template <typename T>                                                  \
        constexpr bool has_##name##_v = is_detected_v<has_##name, T>;          \
        template <typename Value, typename Scalar = scalar_t<Value>,           \
                  typename Mask = mask_t<Value>,                               \
                  bool Single = std::is_same_v<float, Scalar>>                 \
        Value name##_impl(const Value &);                                      \
    }                                                                          \
    template <typename T> auto name(const T &x) {                              \
        using E = expr_t<T>;                                                   \
        using Value = value_t<E>;                                              \
        if constexpr (detail::has_##name##_v<E>) {                             \
            return ((const E &) x).name##_();                                  \
        } else if constexpr (is_recursive_array_v<E>) {                        \
            return E(name(low(x)), name(high(x)));                             \
        } else if constexpr (is_dynamic_array_v<E> &&                          \
                            !is_diff_array_v<E> &&                             \
                            !is_cuda_array_v<E>) {                             \
            E r = empty<E>(x.size());                                          \
            auto pr = r.packet_ptr();                                          \
            auto px = x.packet_ptr();                                          \
            for (size_t i = 0, n = r.packets(); i < n; ++i, ++pr, ++px)        \
                *pr = name(*px);                                               \
            return r;                                                          \
        } else if constexpr (array_depth_v<E> > 1) {                           \
            E r;                                                               \
            ENOKI_CHKSCALAR(#name);                                            \
            for (size_t i = 0; i < x.size(); ++i)                              \
                r.coeff(i) = name(x.coeff(i));                                 \
            return r;                                                          \
        } else if constexpr (is_array_v<E>) {                                  \
            return detail::name##_impl((const E &) x);                         \
        } else {                                                               \
            return expr;                                                       \
        }                                                                      \
    }                                                                          \
    template <typename Value, typename Scalar, typename Mask, bool Single>     \
    ENOKI_INLINE Value enoki::detail::name##_impl(const Value &x)

#define ENOKI_UNARY_OPERATION_PAIR(name, expr)                                 \
    namespace detail {                                                         \
        template <typename T>                                                  \
        using has_##name = decltype(std::declval<T>().name##_());              \
        template <typename T>                                                  \
        constexpr bool has_##name##_v = is_detected_v<has_##name, T>;          \
        template <typename Value, typename Scalar = scalar_t<Value>,           \
                  typename Mask = mask_t<Value>,                               \
                  bool Single = std::is_same_v<float, Scalar>>                 \
        std::pair<Value, Value> name##_impl(const Value &);                    \
    }                                                                          \
    template <typename T> auto name(const T &x) {                              \
        using E = expr_t<T>;                                                   \
        using Value = value_t<E>;                                              \
        if constexpr (detail::has_##name##_v<E>) {                             \
            return ((const E &) x).name##_();                                  \
        } else if constexpr (is_recursive_array_v<E>) {                        \
            auto l = name(low(x));                                             \
            auto h = name(high(x));                                            \
            return std::pair<E, E>(E(l.first, h.first),                        \
                                   E(l.second, h.second));                     \
        } else if constexpr (is_dynamic_array_v<E> &&                          \
                            !is_cuda_array_v<E> &&                             \
                            !is_diff_array_v<E>) {                             \
            std::pair<E, E> r(empty<E>(x.size()), empty<E>(x.size()));         \
            auto pr0 = r.first.packet_ptr(),                                   \
                 pr1 = r.second.packet_ptr();                                  \
            auto px = x.packet_ptr();                                          \
            for (size_t i = 0, n = x.packets();                                \
                 i < n; ++i, ++pr0, ++pr1, ++px)                               \
                std::tie(*pr0, *pr1) = name(*px);                              \
            return r;                                                          \
        } else if constexpr (array_depth_v<E> > 1) {                           \
            std::pair<E, E> r;                                                 \
            ENOKI_CHKSCALAR(#name);                                            \
            for (size_t i = 0; i < x.size(); ++i)                              \
                std::tie(r.first.coeff(i),                                     \
                         r.second.coeff(i)) = name(x.coeff(i));                \
            return r;                                                          \
        } else if constexpr (is_array_v<E>) {                                  \
            return detail::name##_impl((const E &) x);                         \
        } else {                                                               \
            return expr;                                                       \
        }                                                                      \
                                                                               \
    }                                                                          \
    template <typename Value, typename Scalar, typename Mask, bool Single>     \
    ENOKI_INLINE std::pair<Value, Value> enoki::detail::name##_impl(const Value &x)


#define ENOKI_BINARY_OPERATION(name, expr)                                     \
    namespace detail {                                                         \
        template <typename T>                                                  \
        using has_##name = decltype(std::declval<T>()                          \
                                        .name##_(std::declval<T>()));          \
        template <typename T>                                                  \
        constexpr bool has_##name##_v = is_detected_v<has_##name, T>;          \
        template <typename Value, typename Scalar = scalar_t<Value>,           \
                  typename Mask = mask_t<Value>,                               \
                  bool Single = std::is_same_v<float, Scalar>>                 \
        Value name##_impl(const Value &, const Value &);                       \
    }                                                                          \
    template <typename T1, typename T2> auto name(const T1 &x, const T2 &y) {  \
        using E = expr_t<T1, T2>;                                              \
        using Value = value_t<E>;                                              \
        if constexpr (detail::has_##name##_v<E>) {                             \
            return ((const E &) x).name##_((const E &) y);                     \
        } else if constexpr (is_recursive_array_v<E>) {                        \
            return E(name(low(x), low(y)), name(high(x), high(y)));            \
        } else if constexpr (!std::is_same_v<T1, E> ||                         \
                             !std::is_same_v<T2, E>) {                         \
            return name((const E& ) x, (const E &) y);                         \
        } else if constexpr (is_dynamic_array_v<E> &&                          \
                            !is_cuda_array_v<E> &&                             \
                            !is_diff_array_v<E>) {                             \
            E r;                                                               \
            r.resize_like(x, y);                                               \
            size_t xs = x.size() == 1 ? 0 : 1,                                 \
                   ys = y.size() == 1 ? 0 : 1;                                 \
            auto pr = r.packet_ptr();                                          \
            auto px = x.packet_ptr();                                          \
            auto py = y.packet_ptr();                                          \
            for (size_t i = 0, n = r.packets(); i < n;                         \
                 ++i, pr += 1, px += xs, py += ys)                             \
                *pr = name(*px, *py);                                          \
            return r;                                                          \
        } else if constexpr (array_depth_v<E> > 1) {                           \
            assert(x.size() == y.size());                                      \
            E r;                                                               \
            ENOKI_CHKSCALAR(#name);                                            \
            for (size_t i = 0; i < x.size(); ++i)                              \
                r.coeff(i) = name(x.coeff(i), y.coeff(i));                     \
            return r;                                                          \
        } else if constexpr (is_array_v<E>) {                                  \
            return detail::name##_impl((const E &) x, (const E &) y);          \
        } else {                                                               \
            return expr;                                                       \
        }                                                                      \
                                                                               \
    }                                                                          \
    template <typename Value, typename Scalar, typename Mask, bool Single>     \
    ENOKI_INLINE Value enoki::detail::name##_impl(const Value &x, const Value &y)


// -----------------------------------------------------------------------
//! @{ \name Trigonometric functions and their inverses
// -----------------------------------------------------------------------

namespace detail {
    template <bool Sin, bool Cos, typename Value>
    ENOKI_INLINE void sincos_approx(const Value &x, Value &s_out, Value &c_out) {
        using Scalar = scalar_t<Value>;
        constexpr bool Single = std::is_same_v<Scalar, float>;
        using IntArray = int_array_t<Value>;
        using Int = scalar_t<IntArray>;
        using Mask = mask_t<Value>;
        ENOKI_MARK_USED(s_out);
        ENOKI_MARK_USED(c_out);

        /* Joint sine & cosine function approximation based on CEPHES.
           Excellent accuracy in the domain |x| < 8192

           Redistributed under a BSD license with permission of the author, see
           https://github.com/deepmind/torch-cephes/blob/master/LICENSE.txt

         - sin (in [-8192, 8192]):
           * avg abs. err = 6.61896e-09
           * avg rel. err = 1.37888e-08
              -> in ULPs  = 0.166492
           * max abs. err = 5.96046e-08
             (at x=-8191.31)
           * max rel. err = 1.76826e-06
             -> in ULPs   = 19
             (at x=-6374.29)

         - cos (in [-8192, 8192]):
           * avg abs. err = 6.59965e-09
           * avg rel. err = 1.37432e-08
              -> in ULPs  = 0.166141
           * max abs. err = 5.96046e-08
             (at x=-8191.05)
           * max rel. err = 3.13993e-06
             -> in ULPs   = 47
             (at x=-6199.93)
        */

        Value xa = abs(x);

        /* Scale by 4/Pi and get the integer part */
        IntArray j(xa * Scalar(1.2732395447351626862));

        /* Map zeros to origin; if (j & 1) j += 1 */
        j = (j + Int(1)) & Int(~1u);

        /* Cast back to a floating point value */
        Value y(j);

        /* Determine sign of result */
        Value sign_sin, sign_cos;
        constexpr size_t Shift = sizeof(Scalar) * 8 - 3;

        if constexpr (Sin)
            sign_sin = reinterpret_array<Value>(sl<Shift>(j)) ^ x;

        if constexpr (Cos)
            sign_cos = reinterpret_array<Value>(sl<Shift>(~(j - Int(2))));

        /* Extended precision modular arithmetic */
        if constexpr (Single) {
            y = xa - y * Scalar(0.78515625)
                   - y * Scalar(2.4187564849853515625e-4)
                   - y * Scalar(3.77489497744594108e-8);
        } else {
            y = xa - y * Scalar(7.85398125648498535156e-1)
                   - y * Scalar(3.77489470793079817668e-8)
                   - y * Scalar(2.69515142907905952645e-15);
        }

        Value z = y * y, s, c;
        z |= eq(xa, std::numeric_limits<Scalar>::infinity());

        if constexpr (Single) {
            s = poly2(z, -1.6666654611e-1,
                          8.3321608736e-3,
                         -1.9515295891e-4) * z;

            c = poly2(z,  4.166664568298827e-2,
                         -1.388731625493765e-3,
                          2.443315711809948e-5) * z;
        } else {
            s = poly5(z, -1.66666666666666307295e-1,
                          8.33333333332211858878e-3,
                         -1.98412698295895385996e-4,
                          2.75573136213857245213e-6,
                         -2.50507477628578072866e-8,
                          1.58962301576546568060e-10) * z;

            c = poly5(z,  4.16666666666665929218e-2,
                         -1.38888888888730564116e-3,
                          2.48015872888517045348e-5,
                         -2.75573141792967388112e-7,
                          2.08757008419747316778e-9,
                         -1.13585365213876817300e-11) * z;
        }

        s = fmadd(s, y, y);
        c = fmadd(c, z, fmadd(z, Scalar(-0.5), Scalar(1)));

        Mask polymask(eq(j & Int(2), zero<IntArray>()));

        if constexpr (Sin)
            s_out = mulsign(select(polymask, s, c), sign_sin);

        if constexpr (Cos)
            c_out = mulsign(select(polymask, c, s), sign_cos);
    }

    template <bool Tan, typename Value>
    ENOKI_INLINE auto tancot_approx(const Value &x) {
        using Scalar = scalar_t<Value>;
        constexpr bool Single = std::is_same_v<Scalar, float>;
        using IntArray = int_array_t<Value>;
        using Int = scalar_t<IntArray>;

        /*
         - tan (in [-8192, 8192]):
           * avg abs. err = 4.63693e-06
           * avg rel. err = 3.60191e-08
              -> in ULPs  = 0.435442
           * max abs. err = 0.8125
             (at x=-6199.93)
           * max rel. err = 3.12284e-06
             -> in ULPs   = 30
             (at x=-7406.3)
        */

        Value xa = abs(x);

        /* Scale by 4/Pi and get the integer part */
        IntArray j(xa * Scalar(1.2732395447351626862));

        /* Map zeros to origin; if (j & 1) j += 1 */
        j = (j + Int(1)) & Int(~1u);

        /* Cast back to a floating point value */
        Value y(j);

        /* Extended precision modular arithmetic */
        if constexpr (Single) {
            y = xa - y * Scalar(0.78515625)
                   - y * Scalar(2.4187564849853515625e-4)
                   - y * Scalar(3.77489497744594108e-8);
        } else {
            y = xa - y * Scalar(7.85398125648498535156e-1)
                   - y * Scalar(3.77489470793079817668e-8)
                   - y * Scalar(2.69515142907905952645e-15);
        }

        Value z = y * y;
        z |= eq(xa, std::numeric_limits<Scalar>::infinity());

        Value r;
        if constexpr (Single) {
            r = poly5(z, 3.33331568548e-1,
                         1.33387994085e-1,
                         5.34112807005e-2,
                         2.44301354525e-2,
                         3.11992232697e-3,
                         9.38540185543e-3);
        } else {
            r = poly2(z, -1.79565251976484877988e7,
                          1.15351664838587416140e6,
                         -1.30936939181383777646e4) /
                poly4(z, -5.38695755929454629881e7,
                          2.50083801823357915839e7,
                         -1.32089234440210967447e6,
                          1.36812963470692954678e4,
                          1.00000000000000000000e0);
        }

        r = fmadd(r, z * y, y);

        auto recip_mask = Tan ? neq(j & Int(2), Int(0)) :
                                 eq(j & Int(2), Int(0));
        r[xa < Scalar(1e-4)] = y;
        r[recip_mask] = rcp(r);

        Value sign = reinterpret_array<Value>(sl<sizeof(Scalar) * 8 - 2>(j)) ^ x;

        return mulsign(r, sign);
    }
}

ENOKI_UNARY_OPERATION(sin, std::sin(x)) {
    Value r;
    detail::sincos_approx<true, false>(x, r, r);
    return r;
}

ENOKI_UNARY_OPERATION(cos, std::cos(x)) {
    Value r;
    detail::sincos_approx<false, true>(x, r, r);
    return r;
}

ENOKI_UNARY_OPERATION_PAIR(sincos, std::make_pair(std::sin(x), std::cos(x))) {
    Value s, c;
    detail::sincos_approx<true, true>(x, s, c);
    return std::make_pair(s, c);
}

template <typename T> auto csc(const T &a) { return rcp(sin(a)); }
template <typename T> auto sec(const T &a) { return rcp(cos(a)); }

ENOKI_UNARY_OPERATION(tan, std::tan(x)) {
    return detail::tancot_approx<true>(x);
}

ENOKI_UNARY_OPERATION(cot, 1 / std::tan(x)) {
    return detail::tancot_approx<false>(x);
}

ENOKI_UNARY_OPERATION(asin, std::asin(x)) {
    /*
       Arc sine function approximation based on CEPHES.

     - asin (in [-1, 1]):
       * avg abs. err = 2.25422e-08
       * avg rel. err = 2.85777e-08
          -> in ULPs  = 0.331032
       * max abs. err = 1.19209e-07
         (at x=-0.999998)
       * max rel. err = 2.27663e-07
         -> in ULPs   = 2
         (at x=-0.841416)
    */

    Value xa          = abs(x),
          x2          = sqr(x),
          r;

    if constexpr (Single) {
        Mask mask_big = xa > Scalar(0.5);

        Value x1 = Scalar(0.5) * (Scalar(1) - xa);
        Value x3 = select(mask_big, x1, x2);
        Value x4 = select(mask_big, sqrt(x1), xa);

        Value z1 = poly4(x3, 1.6666752422e-1f,
                            7.4953002686e-2f,
                            4.5470025998e-2f,
                            2.4181311049e-2f,
                            4.2163199048e-2f);

        z1 = fmadd(z1, x3*x4, x4);

        r = select(mask_big, Scalar(M_PI_2) - (z1 + z1), z1);
    } else {
        constexpr bool IsCuda = is_cuda_array_v<Value>;
        Mask mask_big = xa > Scalar(0.625);

        if (IsCuda || any_nested(mask_big)) {
            const Scalar pio4 = Scalar(0.78539816339744830962);
            const Scalar more_bits = Scalar(6.123233995736765886130e-17);

            /* arcsin(1-x) = pi/2 - sqrt(2x)(1+R(x))  */
            Value zz = Scalar(1) - xa;
            Value p = poly4(zz, 2.853665548261061424989e1,
                               -2.556901049652824852289e1,
                                6.968710824104713396794e0,
                               -5.634242780008963776856e-1,
                                2.967721961301243206100e-3) /
                      poly4(zz, 3.424398657913078477438e2,
                               -3.838770957603691357202e2,
                                1.470656354026814941758e2,
                               -2.194779531642920639778e1,
                                1.000000000000000000000e0) * zz;
            zz = sqrt(zz + zz);
            Value z = pio4 - zz;
            r[mask_big] = z - fmsub(zz, p, more_bits) + pio4;
        }

        if (IsCuda || !all_nested(mask_big)) {
            Value z = poly5(x2, -8.198089802484824371615e0,
                                 1.956261983317594739197e1,
                                -1.626247967210700244449e1,
                                 5.444622390564711410273e0,
                                -6.019598008014123785661e-1,
                                 4.253011369004428248960e-3) /
                      poly5(x2, -4.918853881490881290097e1,
                                 1.395105614657485689735e2,
                                -1.471791292232726029859e2,
                                 7.049610280856842141659e1,
                                -1.474091372988853791896e1,
                                 1.000000000000000000000e0) * x2;
            z = fmadd(xa, z, xa);
            z = select(xa < Scalar(1e-8), xa, z);
            r[~mask_big] = z;
        }
    }
    return copysign(r, x);
}

ENOKI_UNARY_OPERATION(acos, std::acos(x)) {
    /*
       Arc cosine function approximation based on CEPHES.

     - acos (in [-1, 1]):
       * avg abs. err = 4.72002e-08
       * avg rel. err = 2.85612e-08
          -> in ULPs  = 0.33034
       * max abs. err = 2.38419e-07
         (at x=-0.99999)
       * max rel. err = 1.19209e-07
         -> in ULPs   = 1
         (at x=-0.99999)
    */

    if constexpr (Single) {
        Value xa = abs(x), x2 = sqr(x);

        Mask mask_big = xa > Scalar(0.5);

        Value x1 = Scalar(0.5) * (Scalar(1) - xa);
        Value x3 = select(mask_big, x1, x2);
        Value x4 = select(mask_big, sqrt(x1), xa);

        Value z1 = poly4(x3, 1.666675242e-1f,
                             7.4953002686e-2f,
                             4.5470025998e-2f,
                             2.4181311049e-2f,
                             4.2163199048e-2f);

        z1 = fmadd(z1, x3 * x4, x4);
        Value z2 = z1 + z1;
        z2 = select(x < Scalar(0), Scalar(M_PI) - z2, z2);

        Value z3 = Scalar(M_PI_2) - copysign(z1, x);
        return select(mask_big, z2, z3);
    } else {
        const Scalar pio4 = Scalar(0.78539816339744830962);
        const Scalar more_bits = Scalar(6.123233995736765886130e-17);
        const Scalar h = Scalar(0.5);

        Mask mask = x > h;

        Value y = asin(select(mask, sqrt(fnmadd(h, x, h)), x));
        return select(mask, y + y, pio4 - y + more_bits + pio4);
    }
}

ENOKI_BINARY_OPERATION(atan2, std::atan2(x, y)) {
    /*
       MiniMax fit by Wenzel Jakob, May 2016

     - atan2() tested via atan() (in [-1, 1]):
       * avg abs. err = 1.81543e-07
       * avg rel. err = 4.15224e-07
          -> in ULPs  = 4.9197
       * max abs. err = 5.96046e-07
         (at x=-0.976062)
       * max rel. err = 7.73931e-07
         -> in ULPs   = 12
         (at x=-0.015445)
    */
    Value x_ = y,
          y_ = x,
          abs_x      = abs(x_),
          abs_y      = abs(y_),
          min_val    = min(abs_y, abs_x),
          max_val    = max(abs_x, abs_y),
          scale      = Scalar(1) / max_val,
          scaled_min = min_val * scale,
          z          = scaled_min * scaled_min;

    // How to find these:
    // f[x_] = MiniMaxApproximation[ArcTan[Sqrt[x]]/Sqrt[x],
    //         {x, {1/10000, 1}, 6, 0}, WorkingPrecision->20][[2, 1]]

    Value t;
    if constexpr (Single) {
        t = poly6(z, 0.99999934166683966009,
                    -0.33326497518773606976,
                    +0.19881342388439013552,
                    -0.13486708938456973185,
                    +0.083863120428809689910,
                    -0.037006525670417265220,
                     0.0078613793713198150252);
    } else {
        t = poly6(z, 9.9999999999999999419e-1,
                     2.50554429737833465113e0,
                     2.28289058385464073556e0,
                     9.20960512187107069075e-1,
                     1.59189681028889623410e-1,
                     9.35911604785115940726e-3,
                     8.07005540507283419124e-5) /
            poly6(z, 1.00000000000000000000e0,
                     2.83887763071166519407e0,
                     3.02918312742541450749e0,
                     1.50576983803701596773e0,
                     3.49719171130492192607e-1,
                     3.29968942624402204199e-2,
                     8.26619391703564168942e-4);
    }

    t = t * scaled_min;

    t = select(abs_y > abs_x, Scalar(M_PI_2) - t, t);
    t = select(x_ < zero<Value>(), Scalar(M_PI) - t, t);
    Value r = select(y_ < zero<Value>(), -t, t);
    r &= neq(max_val, Scalar(0));
    return r;
}

ENOKI_UNARY_OPERATION(atan, std::atan(x)) {
    return atan2(x, Scalar(1));
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Exponential function, logarithm, power
// -----------------------------------------------------------------------

ENOKI_BINARY_OPERATION(ldexp, detail::ldexp_scalar(x, y)) {
    return x * reinterpret_array<Value>(
        sl<Single ? 23 : 52>(int_array_t<Value>(y) + (Single ? 0x7f : 0x3ff)));
}

ENOKI_UNARY_OPERATION_PAIR(frexp, detail::frexp_scalar(x)) {
    using IntArray = int_array_t<Value>;
    using Int = scalar_t<IntArray>;
    using IntMask = mask_t<IntArray>;

    const IntArray
        exponent_mask(Int(Single ? 0x7f800000ull : 0x7ff0000000000000ull)),
        mantissa_sign_mask(Int(Single ? ~0x7f800000ull : ~0x7ff0000000000000ull)),
        bias(Int(Single ? 0x7f : 0x3ff));

    IntArray xi = reinterpret_array<IntArray>(x);
    IntArray exponent_bits = xi & exponent_mask;

    /* Detect zero/inf/NaN */
    IntMask is_normal =
        IntMask(neq(x, zero<Value>())) &
        neq(exponent_bits, exponent_mask);

    IntArray exponent_i = (sr<Single ? 23 : 52>(exponent_bits)) - bias;

    IntArray mantissa = (xi & mantissa_sign_mask) |
                         IntArray(memcpy_cast<Int>(Scalar(.5f)));

    return std::make_pair(
        reinterpret_array<Value>(select(is_normal, mantissa, xi)),
        Value(exponent_i & is_normal)
    );
}

ENOKI_UNARY_OPERATION(exp, std::exp(x)) {
    /* Exponential function approximation based on CEPHES

       Redistributed under a BSD license with permission of the author, see
       https://github.com/deepmind/torch-cephes/blob/master/LICENSE.txt

     - exp (in [-20, 30]):
       * avg abs. err = 7155.01
       * avg rel. err = 2.35929e-08
          -> in ULPs  = 0.273524
       * max abs. err = 1.04858e+06
         (at x=29.8057)
       * max rel. err = 1.192e-07
         -> in ULPs   = 1
         (at x=-19.9999)
    */

    const Scalar inf = std::numeric_limits<Scalar>::infinity();
    const Scalar max_range = Scalar(Single ? +88.3762588501 : +7.0943613930310391424428e2);
    const Scalar min_range = Scalar(Single ? -88.3762588501 : -7.0943613930310391424428e2);

    Mask mask_overflow  = x > max_range,
         mask_underflow = x < min_range;

    /* Valueess e^x = e^g 2^n
         = e^g e^(n loge(2))
         = e^(g + n loge(2))
    */
    Value n = floor(fmadd(Scalar(1.4426950408889634073599), x, Scalar(0.5)));
    Value xr = x;
    if constexpr (Single) {
        xr = fnmadd(n, Scalar(0.693359375), xr);
        xr = fnmadd(n, Scalar(-2.12194440e-4), xr);
    } else {
        xr = fnmadd(n, Scalar(6.93145751953125e-1), xr);
        xr = fnmadd(n, Scalar(1.42860682030941723212e-6), xr);
    }

    Value z = sqr(xr);

    if constexpr (Single) {
        z = poly5(xr, 5.0000001201e-1, 1.6666665459e-1,
                      4.1665795894e-2, 8.3334519073e-3,
                      1.3981999507e-3, 1.9875691500e-4);
        z = fmadd(z, xr * xr, xr + Scalar(1));
    } else {
        /* Rational approximation for exponential
           of the fractional part:
              e^x = 1 + 2x P(x^2) / (Q(x^2) - P(x^2))
         */
        Value p = poly2(z, 9.99999999999999999910e-1,
                           3.02994407707441961300e-2,
                           1.26177193074810590878e-4) * xr;

        Value q = poly3(z, 2.00000000000000000009e0,
                           2.27265548208155028766e-1,
                           2.52448340349684104192e-3,
                           3.00198505138664455042e-6);

        Value pq = p / (q-p);
        z = pq + pq + Scalar(1);
    }

    return select(mask_overflow, inf,
                  select(mask_underflow, zero<Value>(), ldexp(z, n)));
}

ENOKI_UNARY_OPERATION(log, std::log(x)) {
    /* Logarithm function approximation based on CEPHES

       Redistributed under a BSD license with permission of the author, see
       https://github.com/deepmind/torch-cephes/blob/master/LICENSE.txt

     - log (in [1e-20, 1000]):
       * avg abs. err = 8.8672e-09
       * avg rel. err = 1.57541e-09
          -> in ULPs  = 0.020038
       * max abs. err = 4.76837e-07
         (at x=54.7661)
       * max rel. err = 1.19194e-07
         -> in ULPs   = 1
         (at x=0.021)
    */

    using UInt = scalar_t<int_array_t<Value>>;

    /* Catch negative and NaN values */
    Mask valid_mask = x >= Scalar(0);
    Value input = x, xm;

    /* The frexp in array_base.h does not handle denormalized numbers,
       cut them off. The AVX512 backend does support them, however. */
    if constexpr (!has_avx512f) {
        Scalar limit = memcpy_cast<Scalar>(
            UInt(Single ? 0x00800000u : 0x0010000000000000ull));
        xm = max(x, limit);
    } else {
        xm = x;
    }

    Value e;
    std::tie(xm, e) = frexp(x);

    const Scalar sqrt_half = Scalar(0.70710678118654752440);
    Mask mask_e_big = abs(e) > Scalar(2);
    Mask mask_ge_inv_sqrt2 = xm >= sqrt_half;
    ENOKI_MARK_USED(mask_e_big);

    e[mask_ge_inv_sqrt2] += Scalar(1);

    Value r;
    if constexpr (Single) {
        xm += (xm & ~mask_ge_inv_sqrt2) - Scalar(1);

        Value z = xm * xm;
        Value y = poly8(xm, 3.3333331174e-1, -2.4999993993e-1,
                            2.0000714765e-1, -1.6668057665e-1,
                            1.4249322787e-1, -1.2420140846e-1,
                            1.1676998740e-1, -1.1514610310e-1,
                            7.0376836292e-2);

        y *= xm * z;

        y = fmadd(e, Scalar(-2.12194440e-4), y);
        z = fmadd(z, Scalar(-0.5), xm + y);
        r = fmadd(e, Scalar(0.693359375), z);
    } else {
        constexpr bool IsCuda = is_cuda_array_v<Value>;
        const Scalar half = Scalar(0.5);
        Value r_big, r_small;

        if (IsCuda || any_nested(mask_e_big)) {
            /* logarithm using log(x) = z + z**3 P(z)/Q(z), where z = 2(x-1)/x+1) */
            Value z = xm - half;

            z[mask_ge_inv_sqrt2] -= half;

            Value y = half * select(mask_ge_inv_sqrt2, xm, z) + half;
            Value x2 = z / y;

            z = x2 * x2;
            z = x2 * (z * poly2(z, -6.41409952958715622951e1,
                                    1.63866645699558079767e1,
                                   -7.89580278884799154124e-1) /
                          poly3(z, -7.69691943550460008604e2,
                                    3.12093766372244180303e2,
                                   -3.56722798256324312549e1,
                                    1.00000000000000000000e0));

            r_big = fnmadd(e, Scalar(2.121944400546905827679e-4), z) + x2;
        }

        if (IsCuda || !all_nested(mask_e_big)) {
            /* logarithm using log(1+x) = x - .5x**2 + x**3 P(x)/Q(x) */
            Value x2 = select(mask_ge_inv_sqrt2, xm, xm + xm) - Scalar(1);

            Value z = x2*x2;
            Value y = x2 * (z * poly5(x2, 7.70838733755885391666e0,
                                          1.79368678507819816313e1,
                                          1.44989225341610930846e1,
                                          4.70579119878881725854e0,
                                          4.97494994976747001425e-1,
                                          1.01875663804580931796e-4) /
                                poly5(x2, 2.31251620126765340583e1,
                                          7.11544750618563894466e1,
                                          8.29875266912776603211e1,
                                          4.52279145837532221105e1,
                                          1.12873587189167450590e1,
                                          1.00000000000000000000e0));

            y = fnmadd(e, Scalar(2.121944400546905827679e-4), y);

            r_small = x2 + fnmadd(half, z, y);
        }

        r = select(mask_e_big, r_big, r_small);
        r = fmadd(e, Scalar(0.693359375), r);
    }

    /* Handle a few special cases */
    const Scalar n_inf(-std::numeric_limits<Scalar>::infinity());
    const Scalar p_inf(std::numeric_limits<Scalar>::infinity());

    r[eq(input, p_inf)] = p_inf;
    r[eq(input, Scalar(0))] = n_inf;

    return r | ~valid_mask;
}

ENOKI_UNARY_OPERATION(cbrt, std::cbrt(x)) {
    /* Cubic root approximation based on CEPHES

       Redistributed under a BSD license with permission of the author, see
       https://github.com/deepmind/torch-cephes/blob/master/LICENSE.txt

     - cbrt (in [-10, 10]):
       * avg abs. err = 2.91027e-17
       * avg rel. err = 1.79292e-17
          -> in ULPs  = 0.118351
       * max abs. err = 4.44089e-16
         (at x=-9.99994)
       * max rel. err = 2.22044e-16
         -> in ULPs   = 1
         (at x=-9.99994)
    */

    const Scalar CBRT2 = Scalar(1.25992104989487316477),
                 CBRT4 = Scalar(1.58740105196819947475),
                 THIRD = Scalar(1.0 / 3.0);

    Value xa = abs(x);

    auto [xm, xe] = frexp(xa);
    xe += Scalar(1);

    Value xea = abs(xe),
          xea1 = floor(xea * THIRD),
          rem = fnmadd(xea1, Scalar(3), xea);

    /* Approximate cube root of number between .5 and 1,
       peak relative error = 9.2e-6 */
    xm = poly4(xm, 0.40238979564544752126924,
                   1.1399983354717293273738,
                  -0.95438224771509446525043,
                   0.54664601366395524503440,
                  -0.13466110473359520655053);

    Value f1 = select(xe >= Scalar(0), Value(CBRT2), Value(Scalar(1) / CBRT2)),
          f2 = select(xe >= Scalar(0), Value(CBRT4), Value(Scalar(1) / CBRT4)),
          f  = select(eq(rem, 1.f), f1, f2);

    xm[neq(rem, 0.f)] *= f;

    Value r = ldexp(xm, mulsign(xea1, xe));
    r = mulsign(r, x);

    // Newton iteration
    r -= (r - (x / sqr(r))) * THIRD;

    if constexpr (!Single)
        r -= (r - (x / sqr(r))) * THIRD;

    return select(isfinite(x), r, x);
}

ENOKI_BINARY_OPERATION(pow, std::pow(x, y)) {
    return exp(log(x) * y);
}

template <typename T, typename E = expr_t<T>>
ENOKI_INLINE E pow(const T &x_, const int &y) {
    int n = std::abs(y);
    E result(1.f), x(x_);

    while (n > 0) {
        if (n & 1)
            result *= x;
        x *= x;
        n /= 2;
    }

    return (y >= 0) ? result : rcp(result);
}

template <typename T, typename E = expr_t<T, float>,
          enable_if_t<is_array_v<T>> = 0>
ENOKI_INLINE E pow(const T &x, const float &y) {
    if (enoki::round(y) == y)
        return enoki::pow(E(x), (int) y);
    else
        return enoki::pow(E(x), E(y));
}

template <typename T, typename E = expr_t<T, double>,
          enable_if_t<is_array_v<T>> = 0>
ENOKI_INLINE E pow(const T &x, const double &y) {
    if (enoki::round(y) == y)
        return enoki::pow(E(x), (int) y);
    else
        return enoki::pow(E(x), E(y));
}

// -----------------------------------------------------------------------
//! @{ \name Hyperbolic and inverse hyperbolic functions
// -----------------------------------------------------------------------

ENOKI_UNARY_OPERATION(sinh, std::sinh(x)) {
    /*
     - sinh (in [-10, 10]):
       * avg abs. err = 2.92524e-05
       * avg rel. err = 2.80831e-08
          -> in ULPs  = 0.336485
       * max abs. err = 0.00195312
         (at x=-9.99894)
       * max rel. err = 2.36862e-07
         -> in ULPs   = 3
         (at x=-9.69866)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    Value xa = abs(x),
          r_small, r_big;

    Mask mask_big = xa > Scalar(1);

    if (IsCuda || any_nested(mask_big)) {
        Value exp0 = exp(x),
              exp1 = rcp(exp0);

        r_big = (exp0 - exp1) * Scalar(0.5);
    }

    if (IsCuda || !all_nested(mask_big)) {
        Value x2 = x * x;

        if constexpr (Single) {
            r_small = fmadd(poly2(x2, 1.66667160211e-1,
                                      8.33028376239e-3,
                                      2.03721912945e-4),
                            x2 * x, x);
        } else {
            r_small = fmadd(poly3(x2, -3.51754964808151394800e5,
                                      -1.15614435765005216044e4,
                                      -1.63725857525983828727e2,
                                      -7.89474443963537015605e-1) /
                            poly3(x2, -2.11052978884890840399e6,
                                       3.61578279834431989373e4,
                                      -2.77711081420602794433e2,
                                       1.00000000000000000000e0),
                            x2 * x, x);
        }
    }

    return select(mask_big, r_big, r_small);
}

ENOKI_UNARY_OPERATION(cosh, std::cosh(x)) {
    /*
     - cosh (in [-10, 10]):
       * avg abs. err = 4.17738e-05
       * avg rel. err = 3.15608e-08
          -> in ULPs  = 0.376252
       * max abs. err = 0.00195312
         (at x=-9.99894)
       * max rel. err = 2.38001e-07
         -> in ULPs   = 3
         (at x=-9.70164)
    */

    Value exp0 = exp(x),
          exp1 = rcp(exp0);

    return (exp0 + exp1) * Scalar(.5f);
}

ENOKI_UNARY_OPERATION_PAIR(sincosh, std::make_pair(std::sinh(x), std::cosh(x))) {
    /*
     - sinh (in [-10, 10]):
       * avg abs. err = 2.92524e-05
       * avg rel. err = 2.80831e-08
          -> in ULPs  = 0.336485
       * max abs. err = 0.00195312
         (at x=-9.99894)
       * max rel. err = 2.36862e-07
         -> in ULPs   = 3
         (at x=-9.69866)

     - cosh (in [-10, 10]):
       * avg abs. err = 4.17738e-05
       * avg rel. err = 3.15608e-08
          -> in ULPs  = 0.376252
       * max abs. err = 0.00195312
         (at x=-9.99894)
       * max rel. err = 2.38001e-07
         -> in ULPs   = 3
         (at x=-9.70164)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    const Scalar half = Scalar(0.5);

    Value xa    = abs(x),
          exp0  = exp(x),
          exp1  = rcp(exp0),
          r_big = (exp0 - exp1) * half,
          r_small;

    Mask mask_big = xa > Scalar(1);

    if (IsCuda || !all_nested(mask_big)) {
        Value x2 = x * x;

        if constexpr (Single) {
            r_small = fmadd(poly2(x2, 1.66667160211e-1,
                                      8.33028376239e-3,
                                      2.03721912945e-4),
                            x2 * x, x);
        } else {
            r_small = fmadd(poly3(x2, -3.51754964808151394800e5,
                                      -1.15614435765005216044e4,
                                      -1.63725857525983828727e2,
                                      -7.89474443963537015605e-1) /
                            poly3(x2, -2.11052978884890840399e6,
                                       3.61578279834431989373e4,
                                      -2.77711081420602794433e2,
                                       1.00000000000000000000e0),
                            x2 * x, x);
        }
    }

    return std::make_pair(
        select(mask_big, r_big, r_small),
        half * (exp0 + exp1)
    );
}

ENOKI_UNARY_OPERATION(tanh, std::tanh(x)) {
    /*
       Hyperbolic tangent function approximation based on CEPHES.

     - tanh (in [-10, 10]):
       * avg abs. err = 4.44655e-08
       * avg rel. err = 4.58074e-08
          -> in ULPs  = 0.698044
       * max abs. err = 3.57628e-07
         (at x=-2.12867)
       * max rel. err = 4.1006e-07
         -> in ULPs   = 6
         (at x=-2.12867)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    Value r_big, r_small;

    Mask mask_big = abs(x) >= Scalar(0.625);

    if (IsCuda || !all_nested(mask_big)) {
        Value x2 = x*x;

        if constexpr (Single) {
            r_small = poly4(x2, -3.33332819422e-1,
                                 1.33314422036e-1,
                                -5.37397155531e-2,
                                 2.06390887954e-2,
                                -5.70498872745e-3);
        } else {
            r_small = poly2(x2, -1.61468768441708447952e3,
                                -9.92877231001918586564e1,
                                -9.64399179425052238628e-1) /
                      poly3(x2,  4.84406305325125486048e3,
                                 2.23548839060100448583e3,
                                 1.12811678491632931402e2,
                                 1.00000000000000000000e0);
        }

        r_small = fmadd(r_small, x2 * x, x);
    }

    if (IsCuda || any_nested(mask_big)) {
        Value e  = exp(x + x),
              e2 = rcp(e + Scalar(1));
        r_big = Scalar(1) - (e2 + e2);
    }

    return select(mask_big, r_big, r_small);
}

template <typename T> auto csch(const T &a) { return rcp(sinh(a)); }
template <typename T> auto sech(const T &a) { return rcp(cosh(a)); }
template <typename T> auto coth(const T &a) { return rcp(tanh(a)); }

ENOKI_UNARY_OPERATION(asinh, std::asinh(x)) {
    /*
       Hyperbolic arc sine function approximation based on CEPHES.

     - asinh (in [-10, 10]):
       * avg abs. err = 2.75626e-08
       * avg rel. err = 1.51762e-08
          -> in ULPs  = 0.178341
       * max abs. err = 2.38419e-07
         (at x=-10)
       * max rel. err = 1.71857e-07
         -> in ULPs   = 2
         (at x=-1.17457)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    Value x2 = x*x,
          xa = abs(x),
          r_big, r_small;

    Mask mask_big  = xa >= Scalar(Single ? 0.51 : 0.533),
         mask_huge = xa >= Scalar(Single ? 1e10 : 1e20);

    if (IsCuda || !all_nested(mask_big)) {
        if constexpr (Single) {
            r_small = poly3(x2, -1.6666288134e-1,
                                 7.4847586088e-2,
                                -4.2699340972e-2,
                                 2.0122003309e-2);
        } else {
            r_small = poly4(x2, -5.56682227230859640450e0,
                                -9.09030533308377316566e0,
                                -4.37390226194356683570e0,
                                -5.91750212056387121207e-1,
                                -4.33231683752342103572e-3) /
                      poly4(x2, 3.34009336338516356383e1,
                                6.95722521337257608734e1,
                                4.86042483805291788324e1,
                                1.28757002067426453537e1,
                                1.00000000000000000000e0);
        }
        r_small = fmadd(r_small, x2 * x, x);
    }

    if (IsCuda || any_nested(mask_big)) {
        r_big = log(xa + (sqrt(x2 + Scalar(1)) & ~mask_huge));
        r_big[mask_huge] += Scalar(M_LN2);
        r_big = copysign(r_big, x);
    }

    return select(mask_big, r_big, r_small);
}

ENOKI_UNARY_OPERATION(acosh, std::acosh(x)) {
    /*
       Hyperbolic arc cosine function approximation based on CEPHES.

     - acosh (in [-10, 10]):
       * avg abs. err = 2.8897e-08
       * avg rel. err = 1.49658e-08
          -> in ULPs  = 0.175817
       * max abs. err = 2.38419e-07
         (at x=3.76221)
       * max rel. err = 2.35024e-07
         -> in ULPs   = 3
         (at x=1.02974)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    Value x1 = x - Scalar(1),
         r_big, r_small;

    Mask mask_big  = x1 >= Scalar(0.49),
         mask_huge = x1 >= Scalar(1e10);

    if (IsCuda || !all_nested(mask_big)) {
        if constexpr (Single) {
            r_small = poly4(x1,  1.4142135263e+0,
                                -1.1784741703e-1,
                                 2.6454905019e-2,
                                -7.5272886713e-3,
                                 1.7596881071e-3);
        } else {
            r_small = poly4(x1, 1.10855947270161294369E5,
                                1.08102874834699867335E5,
                                3.43989375926195455866E4,
                                3.94726656571334401102E3,
                                1.18801130533544501356E2) /
                      poly5(x1, 7.83869920495893927727E4,
                                8.29725251988426222434E4,
                                2.97683430363289370382E4,
                                4.15352677227719831579E3,
                                1.86145380837903397292E2,
                                1.00000000000000000000E0);
        }

        r_small *= sqrt(x1);
        r_small |= x1 < zero<Value>();
    }

    if (IsCuda || any_nested(mask_big)) {
        r_big = log(x + (sqrt(fmsub(x, x, Scalar(1))) & ~mask_huge));
        r_big[mask_huge] += Scalar(M_LN2);
    }

    return select(mask_big, r_big, r_small);
}

ENOKI_UNARY_OPERATION(atanh, std::atanh(x)) {
    /*
       Hyperbolic arc tangent function approximation based on CEPHES.


     - acosh (in [-10, 10]):
       * avg abs. err = 9.87529e-09
       * avg rel. err = 1.52741e-08
          -> in ULPs  = 0.183879
       * max abs. err = 2.38419e-07
         (at x=-0.998962)
       * max rel. err = 1.19209e-07
         -> in ULPs   = 1
         (at x=-0.998962)
    */

    constexpr bool IsCuda = is_cuda_array_v<Value>;

    Value xa = abs(x),
          r_big, r_small;

    Mask mask_big  = xa >= Scalar(0.5);

    if (IsCuda || !all_nested(mask_big)) {
        Value x2 = x*x;
        if constexpr (Single) {
            r_small = poly4(x2, 3.33337300303e-1,
                                1.99782164500e-1,
                                1.46691431730e-1,
                                8.24370301058e-2,
                                1.81740078349e-1);
        } else {
            r_small = poly4(x2, -3.09092539379866942570e1,
                                 6.54566728676544377376e1,
                                -4.61252884198732692637e1,
                                 1.20426861384072379242e1,
                                -8.54074331929669305196e-1) /
                      poly5(x2, -9.27277618139601130017e1,
                                 2.52006675691344555838e2,
                                -2.49839401325893582852e2,
                                 1.08938092147140262656e2,
                                -1.95638849376911654834e1,
                                 1.00000000000000000000e0);
        }
        r_small = fmadd(r_small, x2*x, x);
    }

    if (IsCuda || any_nested(mask_big)) {
        r_big = log((Scalar(1) + xa) / (Scalar(1) - xa)) * Scalar(0.5);
        r_big = copysign(r_big, x);
    }

    return select(mask_big, r_big, r_small);
}

/// Linearly interpolate between 'a' and 'b', using 't'
template <typename Value1, typename Value2, typename Value3>
auto lerp(const Value1 &a, const Value2 &b, const Value3 &t) {
    return fmadd(b, t, fnmadd(a, t, a));
}

/// Clamp the value 'value' to the range [min, max]
template <typename Value1, typename Value2, typename Value3>
auto clamp(const Value1 &value, const Value2 &min, const Value3 &max) {
    return enoki::max(enoki::min(value, max), min);
}

/// Compute the hypotenuse of 'a' and 'b', while avoiding under/overflow
template <typename T1, typename T2>
ENOKI_INLINE auto hypot(const T1 &a, const T2 &b) {
    auto abs_a  = abs(a);
    auto abs_b  = abs(b);
    auto maxval = max(abs_a, abs_b),
         minval = min(abs_a, abs_b),
         ratio  = minval / maxval;

    using Scalar = scalar_t<decltype(ratio)>;
    const Scalar inf = std::numeric_limits<Scalar>::infinity();

    return select(
        (abs_a < inf) && (abs_b < inf) && (ratio < inf),
        maxval * sqrt(Scalar(1) + sqr(ratio)),
        abs_a + abs_b
    );
}

ENOKI_BINARY_OPERATION(fmod, std::fmod(x, y)) {
    return fnmadd(trunc(x / y), y, x);
}

// -----------------------------------------------------------------------
//! @{ \name "Safe" functions that avoid domain errors due to rounding
// -----------------------------------------------------------------------

template <typename T> ENOKI_INLINE auto safe_sqrt(const T &a) {
    return sqrt(max(a, zero<T>()));
}

template <typename T> ENOKI_INLINE auto safe_rsqrt(const T &a) {
    return rsqrt(max(a, zero<T>()));
}

template <typename T> ENOKI_INLINE auto safe_asin(const T &a) {
    return asin(min(T(1), max(T(-1), a)));
}

template <typename T> ENOKI_INLINE auto safe_acos(const T &a) {
    return acos(min(T(1), max(T(-1), a)));
}

/**
 * \brief Numerically well-behaved routine for computing the angle
 * between two unit direction vectors
 *
 * This should be used wherever one is tempted to compute the
 * arc cosine of a dot product.
 *
 * By Don Hatch at http://www.plunk.org/~hatch/rightway.php
 */
template <typename T, typename Expr = expr_t<value_t<T>>>
Expr unit_angle(const T &a, const T &b) {
    Expr dot_uv = dot(a, b),
         temp = 2.f * asin(.5f * norm(b - mulsign(a, dot_uv)));
    return select(dot_uv >= 0, temp, scalar_t<Expr>(M_PI) - temp);
}

/**
 * \brief Numerically well-behaved routine for computing the angle
 * between the unit direction vector 'v' and the z-axis
 *
 * This should be used wherever one is tempted to compute
 * std::acos(v.z())
 *
 * By Don Hatch at http://www.plunk.org/~hatch/rightway.php
 */
template <typename T, typename Expr = expr_t<value_t<T>>>
Expr unit_angle_z(const T &v) {
    static_assert(T::Size == 3, "unit_angle_z(): input is not a 3D vector");
    Expr temp = 2.f * asin(.5f * sqrt(sqr(v.x()) + sqr(v.y()) +
                                      sqr(v.z() - copysign(Expr(1.f), v.z()))));
    return select(v.z() >= 0, temp, scalar_t<Expr>(M_PI) - temp);
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Floating point manipulation routines
// -----------------------------------------------------------------------

template <typename Value, typename Expr = expr_t<Value>>
ENOKI_INLINE Expr prev_float(const Value &value) {
    using Int = int_array_t<Expr>;
    using IntScalar = scalar_t<Int>;

    const Int exponent_mask = sizeof(IntScalar) == 4
                                  ? IntScalar(0x7f800000)
                                  : IntScalar(0x7ff0000000000000ll);

    const Int pos_denorm = sizeof(IntScalar) == 4
                              ? IntScalar(0x80000001)
                              : IntScalar(0x8000000000000001ll);

    Int i = reinterpret_array<Int>(value);

    auto is_nan_inf = eq(i & exponent_mask, exponent_mask);
    auto is_pos_0   = eq(i, 0);
    auto is_gt_0    = i >= 0;
    auto is_special = is_nan_inf | is_pos_0;

    Int j1 = i + select(is_gt_0, Int(-1), Int(1));
    Int j2 = select(is_pos_0, pos_denorm, i);

    return reinterpret_array<Expr>(select(is_special, j2, j1));
}

template <typename Value, typename Expr = expr_t<Value>>
ENOKI_INLINE Expr next_float(const Value &value) {
    using Int = int_array_t<Expr>;
    using IntScalar = scalar_t<Int>;

    const Int exponent_mask = sizeof(IntScalar) == 4
                                  ? IntScalar(0x7f800000)
                                  : IntScalar(0x7ff0000000000000ll);

    const Int sign_mask = sizeof(IntScalar) == 4
                              ? IntScalar(0x80000000)
                              : IntScalar(0x8000000000000000ll);

    Int i = reinterpret_array<Int>(value);

    auto is_nan_inf = eq(i & exponent_mask, exponent_mask);
    auto is_neg_0   = eq(i, sign_mask);
    auto is_gt_0    = i >= 0;
    auto is_special = is_nan_inf | is_neg_0;

    Int j1 = i + select(is_gt_0, Int(1), Int(-1));
    Int j2 = select(is_neg_0, Int(1), i);

    return reinterpret_array<Expr>(select(is_special, j2, j1));
}

template <typename Arg> auto isdenormal(const Arg &a) {
    return abs(a) < std::numeric_limits<scalar_t<Arg>>::min() &&
           neq(a, zero<Arg>());
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
