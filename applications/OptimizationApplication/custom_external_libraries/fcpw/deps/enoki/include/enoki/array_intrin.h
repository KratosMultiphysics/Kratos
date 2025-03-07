/*
    enoki/array_kmask.h -- Hardware-specific intrinsics and compatibility
    wrappers

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once
#include <enoki/fwd.h>

#if defined(ENOKI_X86_64) || defined(ENOKI_X86_32)
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic push
#    pragma GCC diagnostic ignored "-Wconversion"
#    pragma GCC diagnostic ignored "-Wuninitialized"
#    pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#  endif
#  include <immintrin.h>
#  if defined(__GNUC__) && !defined(__clang__)
#    pragma GCC diagnostic pop
#  endif
#endif

#if defined(ENOKI_ARM_NEON)
#  include <arm_neon.h>
#endif

#if defined(_MSC_VER)
#  include <intrin.h>
#endif


NAMESPACE_BEGIN(enoki)

// -----------------------------------------------------------------------
//! @{ \name Available instruction sets
// -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512F)
    static constexpr bool has_avx512f = true;
#else
    static constexpr bool has_avx512f = false;
#endif

#if defined(ENOKI_X86_AVX512CD)
    static constexpr bool has_avx512cd = true;
#else
    static constexpr bool has_avx512cd = false;
#endif

#if defined(ENOKI_X86_AVX512DQ)
    static constexpr bool has_avx512dq = true;
#else
    static constexpr bool has_avx512dq = false;
#endif

#if defined(ENOKI_X86_AVX512VL)
    static constexpr bool has_avx512vl = true;
#else
    static constexpr bool has_avx512vl = false;
#endif

#if defined(ENOKI_X86_AVX512BW)
    static constexpr bool has_avx512bw = true;
#else
    static constexpr bool has_avx512bw = false;
#endif

#if defined(ENOKI_X86_AVX512PF)
    static constexpr bool has_avx512pf = true;
#else
    static constexpr bool has_avx512pf = false;
#endif

#if defined(ENOKI_X86_AVX512ER)
    static constexpr bool has_avx512er = true;
#else
    static constexpr bool has_avx512er = false;
#endif

#if defined(__AVX512VBMI__)
    static constexpr bool has_avx512vbmi = true;
#else
    static constexpr bool has_avx512vbmi = false;
#endif

#if defined(ENOKI_X86_AVX512VPOPCNTDQ)
    static constexpr bool has_avx512vpopcntdq = true;
#else
    static constexpr bool has_avx512vpopcntdq = false;
#endif

#if defined(ENOKI_X86_AVX2)
    static constexpr bool has_avx2 = true;
#else
    static constexpr bool has_avx2 = false;
#endif

#if defined(ENOKI_X86_FMA) || defined(ENOKI_ARM_FMA)
    static constexpr bool has_fma = true;
#else
    static constexpr bool has_fma = false;
#endif

#if defined(ENOKI_X86_F16C)
    static constexpr bool has_f16c = true;
#else
    static constexpr bool has_f16c = false;
#endif

#if defined(ENOKI_X86_AVX)
    static constexpr bool has_avx = true;
#else
    static constexpr bool has_avx = false;
#endif

#if defined(ENOKI_X86_SSE42)
    static constexpr bool has_sse42 = true;
#else
    static constexpr bool has_sse42 = false;
#endif

#if defined(ENOKI_X86_32)
    static constexpr bool has_x86_32 = true;
#else
    static constexpr bool has_x86_32 = false;
#endif

#if defined(ENOKI_X86_64)
    static constexpr bool has_x86_64 = true;
#else
    static constexpr bool has_x86_64 = false;
#endif

#if defined(ENOKI_ARM_NEON)
    static constexpr bool has_neon = true;
#else
    static constexpr bool has_neon = false;
#endif

#if defined(ENOKI_ARM_32)
    static constexpr bool has_arm_32 = true;
#else
    static constexpr bool has_arm_32 = false;
#endif

#if defined(ENOKI_ARM_64)
    static constexpr bool has_arm_64 = true;
#else
    static constexpr bool has_arm_64 = false;
#endif

static constexpr bool has_vectorization = has_sse42 || has_neon;

//! @}
// -----------------------------------------------------------------------

#if defined(ENOKI_X86_SSE42)
/// Flush denormalized numbers to zero
inline void set_flush_denormals(bool value) {
    _MM_SET_FLUSH_ZERO_MODE(value ? _MM_FLUSH_ZERO_ON : _MM_FLUSH_ZERO_OFF);
    _MM_SET_DENORMALS_ZERO_MODE(value ? _MM_DENORMALS_ZERO_ON : _MM_DENORMALS_ZERO_OFF);
}

inline bool flush_denormals() {
    return _MM_GET_FLUSH_ZERO_MODE() == _MM_FLUSH_ZERO_ON;
}

#else
inline void set_flush_denormals(bool) { }
inline bool flush_denormals() { return false; }
#endif

struct scoped_flush_denormals {
public:
    scoped_flush_denormals(bool value) {
        m_old_value = flush_denormals();
        set_flush_denormals(value);

    }

    ~scoped_flush_denormals() {
        set_flush_denormals(m_old_value);
    }
private:
    bool m_old_value;
};

NAMESPACE_BEGIN(detail)

// -----------------------------------------------------------------------
//! @{ \name Helper routines to merge smaller arrays into larger ones
// -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX)
ENOKI_INLINE __m256 concat(__m128 l, __m128 h) {
    return _mm256_insertf128_ps(_mm256_castps128_ps256(l), h, 1);
}

ENOKI_INLINE __m256d concat(__m128d l, __m128d h) {
    return _mm256_insertf128_pd(_mm256_castpd128_pd256(l), h, 1);
}

ENOKI_INLINE __m256i concat(__m128i l, __m128i h) {
    return _mm256_insertf128_si256(_mm256_castsi128_si256(l), h, 1);
}
#endif

#if defined(ENOKI_X86_AVX512F)
ENOKI_INLINE __m512 concat(__m256 l, __m256 h) {
    #if defined(ENOKI_X86_AVX512DQ)
        return _mm512_insertf32x8(_mm512_castps256_ps512(l), h, 1);
    #else
        return _mm512_castpd_ps(
            _mm512_insertf64x4(_mm512_castps_pd(_mm512_castps256_ps512(l)),
                               _mm256_castps_pd(h), 1));
    #endif
}

ENOKI_INLINE __m512d concat(__m256d l, __m256d h) {
    return _mm512_insertf64x4(_mm512_castpd256_pd512(l), h, 1);
}

ENOKI_INLINE __m512i concat(__m256i l, __m256i h) {
    return _mm512_inserti64x4(_mm512_castsi256_si512(l), h, 1);
}
#endif

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Mask conversion routines for various platforms
// -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX)
ENOKI_INLINE __m256i mm256_cvtepi32_epi64(__m128i x) {
#if defined(ENOKI_X86_AVX2)
    return _mm256_cvtepi32_epi64(x);
#else
    /* This version is only suitable for mask conversions */
    __m128i xl = _mm_shuffle_epi32(x, _MM_SHUFFLE(1, 1, 0, 0));
    __m128i xh = _mm_shuffle_epi32(x, _MM_SHUFFLE(3, 3, 2, 2));
    return detail::concat(xl, xh);
#endif
}

ENOKI_INLINE __m128i mm256_cvtepi64_epi32(__m256i x) {
#if defined(ENOKI_X86_AVX512VL)
    return _mm256_cvtepi64_epi32(x);
#else
    __m128i x0 = _mm256_castsi256_si128(x);
    __m128i x1 = _mm256_extractf128_si256(x, 1);
    return _mm_castps_si128(_mm_shuffle_ps(
        _mm_castsi128_ps(x0), _mm_castsi128_ps(x1), _MM_SHUFFLE(2, 0, 2, 0)));
#endif
}

ENOKI_INLINE __m256i mm512_cvtepi64_epi32(__m128i x0, __m128i x1, __m128i x2, __m128i x3) {
    __m128i y0 = _mm_castps_si128(_mm_shuffle_ps(
        _mm_castsi128_ps(x0), _mm_castsi128_ps(x1), _MM_SHUFFLE(2, 0, 2, 0)));
    __m128i y1 = _mm_castps_si128(_mm_shuffle_ps(
        _mm_castsi128_ps(x2), _mm_castsi128_ps(x3), _MM_SHUFFLE(2, 0, 2, 0)));
    return detail::concat(y0, y1);
}

ENOKI_INLINE __m256i mm512_cvtepi64_epi32(__m256i x0, __m256i x1) {
    __m128i y0 = _mm256_castsi256_si128(x0);
    __m128i y1 = _mm256_extractf128_si256(x0, 1);
    __m128i y2 = _mm256_castsi256_si128(x1);
    __m128i y3 = _mm256_extractf128_si256(x1, 1);
    return mm512_cvtepi64_epi32(y0, y1, y2, y3);
}
#endif

#if defined(ENOKI_X86_SSE42)

ENOKI_INLINE __m128i mm256_cvtepi64_epi32(__m128i x0, __m128i x1) {
    return _mm_castps_si128(_mm_shuffle_ps(
        _mm_castsi128_ps(x0), _mm_castsi128_ps(x1), _MM_SHUFFLE(2, 0, 2, 0)));
}

ENOKI_INLINE __m128i mm_cvtsi64_si128(long long a)  {
    #if defined(ENOKI_X86_64)
        return _mm_cvtsi64_si128(a);
    #else
        alignas(16) long long x[2] = { a, 0ll };
        return _mm_load_si128((__m128i *) x);
    #endif
}

ENOKI_INLINE long long mm_cvtsi128_si64(__m128i m)  {
    #if defined(ENOKI_X86_64)
        return _mm_cvtsi128_si64(m);
    #else
        alignas(16) long long x[2];
        _mm_store_si128((__m128i *) x, m);
        return x[0];
    #endif
}

template <int Imm8>
ENOKI_INLINE long long mm_extract_epi64(__m128i m)  {
    #if defined(ENOKI_X86_64)
        return _mm_extract_epi64(m, Imm8);
    #else
        alignas(16) long long x[2];
        _mm_store_si128((__m128i *) x, m);
        return x[Imm8];
    #endif
}

#endif

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(detail)
NAMESPACE_END(enoki)
