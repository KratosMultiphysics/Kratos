/*
    enoki/array_fallbacks.h -- Scalar fallback implementations of various
    operations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_intrin.h>

NAMESPACE_BEGIN(enoki)
NAMESPACE_BEGIN(detail)

/// Reciprocal (scalar fallback)
template <typename T> ENOKI_INLINE T rcp_scalar(const T &a) {
#if defined(ENOKI_X86_AVX512ER)
    if (std::is_same_v<T, float>) {
        __m128 v = _mm_set_ss((float) a);
        return T(_mm_cvtss_f32(_mm_rcp28_ss(v, v))); /* rel error < 2^-28 */
    }
#endif

    if constexpr (std::is_same_v<T, float>) {
#if defined(ENOKI_X86_SSE42)
        __m128 v = _mm_set_ss((float) a), r;

        #if defined(ENOKI_X86_AVX512F)
            r = _mm_rcp14_ss(v, v); /* rel error < 2^-14 */
        #else
            r = _mm_rcp_ss(v);      /* rel error < 1.5*2^-12 */
        #endif

        /* Refine using one Newton-Raphson iteration */
        __m128 ro = r;

        __m128 t0 = _mm_add_ss(r, r);
        __m128 t1 = _mm_mul_ss(r, v);

        #if defined(ENOKI_X86_FMA)
            r = _mm_fnmadd_ss(r, t1, t0);
        #else
            r = _mm_sub_ss(t0, _mm_mul_ss(r, t1));
        #endif

        #if defined(ENOKI_X86_AVX512F)
            (void) ro;
            r = _mm_fixupimm_ss(r, v, _mm_set1_epi32(0x0087A622), 0);
        #else
            r = _mm_blendv_ps(r, ro, t1); /* mask bit is '1' iff t1 == nan */
        #endif

        return T(_mm_cvtss_f32(r));
#elif defined(ENOKI_ARM_NEON) && defined(ENOKI_ARM_64)
        float v = (float) a;
        float r = vrecpes_f32(v);
        r *= vrecpss_f32(r, v);
        r *= vrecpss_f32(r, v);
        return T(r);
#endif
    }

#if defined(ENOKI_X86_AVX512F) || defined(ENOKI_X86_AVX512ER)
    if constexpr (std::is_same_v<T, double>) {
        __m128d v = _mm_set_sd((double) a), r;

        #if defined(ENOKI_X86_AVX512ER)
            r = _mm_rcp28_sd(v, v);  /* rel error < 2^-28 */
        #elif defined(ENOKI_X86_AVX512F)
            r = _mm_rcp14_sd(v, v);  /* rel error < 2^-14 */
        #endif

        __m128d ro = r, t0, t1;

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            t0 = _mm_add_sd(r, r);
            t1 = _mm_mul_sd(r, v);

            #if defined(ENOKI_X86_FMA)
                r = _mm_fnmadd_sd(t1, r, t0);
            #else
                r = _mm_sub_sd(t0, _mm_mul_sd(r, t1));
            #endif
        }

        r = _mm_blendv_pd(r, ro, t1); /* mask bit is '1' iff t1 == nan */

        return T(_mm_cvtsd_f64(r));
    }
#endif

    return T(1) / a;
}

/// Reciprocal square root (scalar fallback)
template <typename T> ENOKI_INLINE T rsqrt_scalar(const T &a) {
#if defined(ENOKI_X86_AVX512ER)
    if (std::is_same_v<T, float>) {
        __m128 v = _mm_set_ss((float) a);
        return T(_mm_cvtss_f32(_mm_rsqrt28_ss(v, v))); /* rel error < 2^-28 */
    }
#endif

    if constexpr (std::is_same_v<T, float>) {
#if defined(ENOKI_X86_SSE42)
        __m128 v = _mm_set_ss((float) a), r;
        #if defined(ENOKI_X86_AVX512F)
            r = _mm_rsqrt14_ss(v, v);  /* rel error < 2^-14 */
        #else
            r = _mm_rsqrt_ss(v);       /* rel error < 1.5*2^-12 */
        #endif

        /* Refine using one Newton-Raphson iteration */
        const __m128 c0 = _mm_set_ss(0.5f),
                     c1 = _mm_set_ss(3.0f);

        __m128 t0 = _mm_mul_ss(r, c0),
               t1 = _mm_mul_ss(r, v),
               ro = r;

        #if defined(ENOKI_X86_FMA)
            r = _mm_mul_ss(_mm_fnmadd_ss(t1, r, c1), t0);
        #else
            r = _mm_mul_ss(_mm_sub_ss(c1, _mm_mul_ss(t1, r)), t0);
        #endif

        #if defined(ENOKI_X86_AVX512F)
            (void) ro;
            r = _mm_fixupimm_ss(r, v, _mm_set1_epi32(0x0383A622), 0);
        #else
            r = _mm_blendv_ps(r, ro, t1); /* mask bit is '1' iff t1 == nan */
        #endif

        return T(_mm_cvtss_f32(r));
#elif defined(ENOKI_ARM_NEON) && defined(ENOKI_ARM_64)
        float v = (float) a;
        float r = vrsqrtes_f32(v);
        r *= vrsqrtss_f32(r*r, v);
        r *= vrsqrtss_f32(r*r, v);
        return r;
#endif
    }

#if defined(ENOKI_X86_AVX512F) || defined(ENOKI_X86_AVX512ER)
    if constexpr (std::is_same_v<T, double>) {
        __m128d v = _mm_set_sd((double) a), r;

        #if defined(ENOKI_X86_AVX512ER)
            r = _mm_rsqrt28_sd(v, v);  /* rel error < 2^-28 */
        #elif defined(ENOKI_X86_AVX512F)
            r = _mm_rsqrt14_sd(v, v);  /* rel error < 2^-14 */
        #endif

        const __m128d c0 = _mm_set_sd(0.5),
                      c1 = _mm_set_sd(3.0);

        __m128d ro = r, t0, t1;

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            t0 = _mm_mul_sd(r, c0);
            t1 = _mm_mul_sd(r, v);

            #if defined(ENOKI_X86_FMA)
                r = _mm_mul_sd(_mm_fnmadd_sd(t1, r, c1), t0);
            #else
                r = _mm_mul_sd(_mm_sub_sd(c1, _mm_mul_sd(t1, r)), t0);
            #endif
        }

        r = _mm_blendv_pd(r, ro, t1); /* mask bit is '1' iff t1 == nan */

        return T(_mm_cvtsd_f64(r));
    }
#endif

    return T(1) / std::sqrt(a);
}

template <typename T> ENOKI_INLINE T popcnt_scalar(T v) {
    static_assert(std::is_integral_v<T>, "popcnt(): requires an integer argument!");
#if defined(ENOKI_X86_SSE42)
    if constexpr (sizeof(T) <= 4) {
        return (T) _mm_popcnt_u32((unsigned int) v);
    } else {
        #if defined(ENOKI_X86_64)
            return (T) _mm_popcnt_u64((unsigned long long) v);
        #else
            unsigned long long v_ = (unsigned long long) v;
            unsigned int lo = (unsigned int) v_;
            unsigned int hi = (unsigned int) (v_ >> 32);
            return (T) (_mm_popcnt_u32(lo) + _mm_popcnt_u32(hi));
        #endif
    }
#elif defined(_MSC_VER)
    if constexpr (sizeof(T) <= 4) {
        uint32_t w = (uint32_t) v;
        w -= (w >> 1) & 0x55555555;
        w = (w & 0x33333333) + ((w >> 2) & 0x33333333);
        w = (w + (w >> 4)) & 0x0F0F0F0F;
        w = (w * 0x01010101) >> 24;
        return (T) w;
    } else {
        uint64_t w = (uint64_t) v;
        w -= (w >> 1) & 0x5555555555555555ull;
        w = (w & 0x3333333333333333ull) + ((w >> 2) & 0x3333333333333333ull);
        w = (w + (w >> 4)) & 0x0F0F0F0F0F0F0F0Full;
        w = (w * 0x0101010101010101ull) >> 56;
        return (T) w;
    }
#else
    if constexpr (sizeof(T) <= 4)
        return (T) __builtin_popcount((unsigned int) v);
    else
        return (T) __builtin_popcountll((unsigned long long) v);
#endif
}

template <typename T> ENOKI_INLINE T lzcnt_scalar(T v) {
    static_assert(std::is_integral_v<T>, "lzcnt(): requires an integer argument!");
#if defined(ENOKI_X86_AVX2)
    if constexpr (sizeof(T) <= 4) {
        return (T) _lzcnt_u32((unsigned int) v);
    } else {
        #if defined(ENOKI_X86_64)
            return (T) _lzcnt_u64((unsigned long long) v);
        #else
            unsigned long long v_ = (unsigned long long) v;
            unsigned int lo = (unsigned int) v_;
            unsigned int hi = (unsigned int) (v_ >> 32);
            return (T) (hi != 0 ? _lzcnt_u32(hi) : (_lzcnt_u32(lo) + 32));
        #endif
    }
#elif defined(_MSC_VER)
    unsigned long result;
    if constexpr (sizeof(T) <= 4) {
        _BitScanReverse(&result, (unsigned long) v);
        return (v != 0) ? (31 - result) : 32;
    } else {
        _BitScanReverse64(&result, (unsigned long long) v);
        return (v != 0) ? (63 - result) : 64;
    }
#else
    if constexpr (sizeof(T) <= 4)
        return (T) (v != 0 ? __builtin_clz((unsigned int) v) : 32);
    else
        return (T) (v != 0 ? __builtin_clzll((unsigned long long) v) : 64);
#endif
}

template <typename T> ENOKI_INLINE T tzcnt_scalar(T v) {
    static_assert(std::is_integral_v<T>, "tzcnt(): requires an integer argument!");
#if defined(ENOKI_X86_AVX2)
    if (sizeof(T) <= 4)
        return (T) _tzcnt_u32((unsigned int) v);
    #if defined(ENOKI_X86_64)
        return (T) _tzcnt_u64((unsigned long long) v);
    #else
        unsigned long long v_ = (unsigned long long) v;
        unsigned int lo = (unsigned int) v_;
        unsigned int hi = (unsigned int) (v_ >> 32);
        return (T) (lo != 0 ? _tzcnt_u32(lo) : (_tzcnt_u32(hi) + 32));
    #endif
#elif defined(_MSC_VER)
    unsigned long result;
    if (sizeof(T) <= 4) {
        _BitScanForward(&result, (unsigned long) v);
        return (v != 0) ? result : 32;
    } else {
        _BitScanForward64(&result, (unsigned long long) v);
        return (v != 0) ? result: 64;
    }
#else
    if (sizeof(T) <= 4)
        return (T) (v != 0 ? __builtin_ctz((unsigned int) v) : 32);
    else
        return (T) (v != 0 ? __builtin_ctzll((unsigned long long) v) : 64);
#endif
}

template <typename T1, typename T2>
ENOKI_INLINE T1 ldexp_scalar(const T1 &a1, const T2 &a2) {
#if defined(ENOKI_X86_AVX512F)
    if constexpr (std::is_same_v<T1, float>) {
        __m128 v1 = _mm_set_ss((float) a1),
               v2 = _mm_set_ss((float) a2);
        return T1(_mm_cvtss_f32(_mm_scalef_ss(v1, v2)));
    } else if constexpr (std::is_same_v<T1, double>) {
        __m128d v1 = _mm_set_sd((double) a1),
                v2 = _mm_set_sd((double) a2);
        return T1(_mm_cvtsd_f64(_mm_scalef_sd(v1, v2)));
    } else {
        return std::ldexp(a1, int(a2));
    }
#else
    return std::ldexp(a1, int(a2));
#endif
}

/// Break floating-point number into normalized fraction and power of 2 (scalar fallback)
template <typename T>
ENOKI_INLINE std::pair<T, T> frexp_scalar(const T &a) {
#if defined(ENOKI_X86_AVX512F)
    if constexpr (std::is_same_v<T, float>) {
        __m128 v = _mm_set_ss((float) a);
        return std::make_pair(
            T(_mm_cvtss_f32(_mm_getmant_ss(v, v, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src))),
            T(_mm_cvtss_f32(_mm_getexp_ss(v, v))));
    } else if constexpr (std::is_same_v<T, double>) {
        __m128d v = _mm_set_sd((double) a);
        return std::make_pair(
            T(_mm_cvtsd_f64(_mm_getmant_sd(v, v, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src))),
            T(_mm_cvtsd_f64(_mm_getexp_sd(v, v))));
    } else {
        int tmp;
        T result = std::frexp(a, &tmp);
        return std::make_pair(result, T(tmp) - T(1));
    }
#else
    int tmp;
    T result = std::frexp(a, &tmp);
    return std::make_pair(result, T(tmp) - T(1));
#endif
}

ENOKI_INLINE int32_t mulhi_scalar(int32_t x, int32_t y) {
    int64_t rl = (int64_t) x * (int64_t) y;
    return (int32_t) (rl >> 32);
}

ENOKI_INLINE uint32_t mulhi_scalar(uint32_t x, uint32_t y) {
    uint64_t rl = (uint64_t) x * (uint64_t) y;
    return (uint32_t) (rl >> 32);
}

ENOKI_INLINE uint64_t mulhi_scalar(uint64_t x, uint64_t y) {
#if defined(_MSC_VER) && defined(ENOKI_X86_64)
    return __umulh(x, y);
#elif defined(__SIZEOF_INT128__)
    __uint128_t rl = (__uint128_t) x * (__uint128_t) y;
    return (uint64_t)(rl >> 64);
#else
    // full 128 bits are x0 * y0 + (x0 * y1 << 32) + (x1 * y0 << 32) + (x1 * y1 << 64)
    const uint32_t mask = 0xFFFFFFFF;
    const uint32_t x0 = (uint32_t) (x & mask), x1 = (uint32_t) (x >> 32);
    const uint32_t y0 = (uint32_t) (y & mask), y1 = (uint32_t) (y >> 32);
    const uint32_t x0y0_hi = mulhi_scalar(x0, y0);
    const uint64_t x0y1 = x0 * (uint64_t) y1;
    const uint64_t x1y0 = x1 * (uint64_t) y0;
    const uint64_t x1y1 = x1 * (uint64_t) y1;
    const uint64_t temp = x1y0 + x0y0_hi;
    const uint64_t temp_lo = temp & mask, temp_hi = temp >> 32;

    return x1y1 + temp_hi + ((temp_lo + x0y1) >> 32);
#endif
}

ENOKI_INLINE int64_t mulhi_scalar(int64_t x, int64_t y) {
#if defined(_MSC_VER) && defined(_M_X64)
    return __mulh(x, y);
#elif defined(__SIZEOF_INT128__)
    __int128_t rl = (__int128_t) x * (__int128_t) y;
    return (int64_t)(rl >> 64);
#else
    // full 128 bits are x0 * y0 + (x0 * y1 << 32) + (x1 * y0 << 32) + (x1 * y1 << 64)
    const uint32_t mask = 0xFFFFFFFF;
    const uint32_t x0 = (uint32_t) (x & mask), y0 = (uint32_t) (y & mask);
    const int32_t x1 = (int32_t) (x >> 32), y1 = (int32_t) (y >> 32);
    const uint32_t x0y0_hi = mulhi_scalar(x0, y0);
    const int64_t t = x1 * (int64_t) y0 + x0y0_hi;
    const int64_t w1 = x0 * (int64_t) y1 + (t & mask);

    return x1 * (int64_t) y1 + (t >> 32) + (w1 >> 32);
#endif
}

template <typename T> ENOKI_INLINE T abs_scalar(const T &a) {
    if constexpr (std::is_signed_v<T>)
        return std::abs(a);
    else
        return a;
}

template <typename T1, typename T2, typename T3,
          typename E = expr_t<T1, T2, T3>> ENOKI_INLINE E fmadd_scalar(T1 a1, T2 a2, T3 a3) {
#if defined(ENOKI_X86_FMA) || defined(ENOKI_ARM_FMA)
    if constexpr (std::is_floating_point_v<E>)
        return (E) std::fma((E) a1, (E) a2, (E) a3);
#endif
    return (E) a1 * (E) a2 + (E) a3;
}

template <typename T, typename Arg>
T ceil2int_scalar(Arg x) {
#if defined(ENOKI_X86_AVX512F)
    if constexpr (std::is_same_v<Arg, float>) {
        __m128 y = _mm_set_ss(x);
        if constexpr (sizeof(T) == 4) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundss_i32(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundss_u32(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
        } else if constexpr (sizeof(T) == 8) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundss_i64(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundss_u64(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
        }
    } else if constexpr (std::is_same_v<Arg, double>) {
        __m128d y = _mm_set_sd(x);
        if constexpr (sizeof(T) == 4) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundsd_i32(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundsd_u32(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
        } else if constexpr (sizeof(T) == 8) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundsd_i64(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundsd_u64(y, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC);
        }
    }
#endif
    return T(std::ceil(x));
}

template <typename T, typename Arg>
T floor2int_scalar(Arg x) {
#if defined(ENOKI_X86_AVX512F)
    if constexpr (std::is_same_v<Arg, float>) {
        __m128 y = _mm_set_ss(x);
        if constexpr (sizeof(T) == 4) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundss_i32(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundss_u32(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
        } else if constexpr (sizeof(T) == 8) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundss_i64(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundss_u64(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
        }
    } else if constexpr (std::is_same_v<Arg, double>) {
        __m128d y = _mm_set_sd(x);
        if constexpr (sizeof(T) == 4) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundsd_i32(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundsd_u32(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
        } else if constexpr (sizeof(T) == 8) {
            if constexpr (std::is_signed_v<T>)
                return _mm_cvt_roundsd_i64(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
            else
                return _mm_cvt_roundsd_u64(y, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);
        }
    }
#endif
    return T(std::floor(x));
}

template <typename T> auto or_(const T &a1, const T &a2) {
    using Int = int_array_t<T, false>;

    if constexpr (is_array_v<T> || std::is_integral_v<T>)
        return a1 | a2;
    else
        return memcpy_cast<T>(memcpy_cast<Int>(a1) | memcpy_cast<Int>(a2));
}

template <typename T> auto and_(const T &a1, const T &a2) {
    using Int = int_array_t<T, false>;

    if constexpr (is_array_v<T> || std::is_integral_v<T>)
        return a1 & a2;
    else
        return memcpy_cast<T>(memcpy_cast<Int>(a1) & memcpy_cast<Int>(a2));
}

template <typename T> auto andnot_(const T &a1, const T &a2) {
    using Int = int_array_t<T, false>;

    if constexpr (is_array_v<T>)
        return andnot(a1, a2);
    else if constexpr (std::is_same_v<T, bool>)
        return a1 && !a2;
    else if constexpr (std::is_integral_v<T>)
        return a1 & ~a2;
    else
        return memcpy_cast<T>(memcpy_cast<Int>(a1) & ~memcpy_cast<Int>(a2));
}

template <typename T> auto xor_(const T &a1, const T &a2) {
    using Int = int_array_t<T, false>;

    if constexpr (is_array_v<T> || std::is_integral_v<T>)
        return a1 ^ a2;
    else
        return memcpy_cast<T>(memcpy_cast<Int>(a1) ^ memcpy_cast<Int>(a2));
}

template <typename T, enable_if_t<!std::is_same_v<T, bool>> = 0> auto or_(const T &a, const bool &b) {
    using Scalar = scalar_t<T>;
    using Int    = int_array_t<Scalar>;
    return or_(a, b ? memcpy_cast<Scalar>(Int(-1)) : memcpy_cast<Scalar>(Int(0)));
}

template <typename T, enable_if_t<!std::is_same_v<T, bool>> = 0> auto and_(const T &a, const bool &b) {
    using Scalar = scalar_t<T>;
    using Int    = int_array_t<Scalar>;
    return and_(a, b ? memcpy_cast<Scalar>(Int(-1)) : memcpy_cast<Scalar>(Int(0)));
}

template <typename T, enable_if_t<!std::is_same_v<T, bool>> = 0> auto andnot_(const T &a, const bool &b) {
    using Scalar = scalar_t<T>;
    using Int    = int_array_t<Scalar>;
    return andnot_(a, b ? memcpy_cast<Scalar>(Int(-1)) : memcpy_cast<Scalar>(Int(0)));
}

template <typename T, enable_if_t<!std::is_same_v<T, bool>> = 0> auto xor_(const T &a, const bool &b) {
    using Scalar = scalar_t<T>;
    using Int    = int_array_t<Scalar>;
    return xor_(a, b ? memcpy_cast<Scalar>(Int(-1)) : memcpy_cast<Scalar>(Int(0)));
}

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
auto or_(const T1 &a1, const T2 &a2) { return a1 | a2; }

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
auto and_(const T1 &a1, const T2 &a2) { return a1 & a2; }

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
auto andnot_(const T1 &a1, const T2 &a2) { return andnot(a1, a2); }

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
auto xor_(const T1 &a1, const T2 &a2) { return a1 ^ a2; }

NAMESPACE_END(detail)
NAMESPACE_END(enoki)
