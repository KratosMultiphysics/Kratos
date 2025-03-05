/*
    enoki/array_sse42.h -- Packed SIMD array (SSE4.2 specialization)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

NAMESPACE_BEGIN(enoki)
NAMESPACE_BEGIN(detail)

/// Compressed look-up table for the store_compress() operation [256 bytes]
alignas(16) const uint8_t compress_lut_128[16*16] = {
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x00, 0x01, 0x02, 0x03, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x04, 0x05, 0x06, 0x07,
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x08, 0x09, 0x0a, 0x0b, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00, 0x01, 0x02, 0x03,
    0x08, 0x09, 0x0a, 0x0b, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
    0x08, 0x09, 0x0a, 0x0b, 0x80, 0x80, 0x80, 0x80, 0x0c, 0x0d, 0x0e, 0x0f,
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
    0x00, 0x01, 0x02, 0x03, 0x0c, 0x0d, 0x0e, 0x0f, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f,
    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00, 0x01, 0x02, 0x03,
    0x04, 0x05, 0x06, 0x07, 0x0c, 0x0d, 0x0e, 0x0f, 0x80, 0x80, 0x80, 0x80,
    0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x80, 0x80, 0x80, 0x80,
    0x80, 0x80, 0x80, 0x80, 0x00, 0x01, 0x02, 0x03, 0x08, 0x09, 0x0a, 0x0b,
    0x0c, 0x0d, 0x0e, 0x0f, 0x80, 0x80, 0x80, 0x80, 0x04, 0x05, 0x06, 0x07,
    0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x80, 0x80, 0x80, 0x80,
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b,
    0x0c, 0x0d, 0x0e, 0x0f
};

template <> struct is_native<float, 4> : std::true_type { } ;
template <> struct is_native<float, 3> : std::true_type { };
template <> struct is_native<double, 2> : std::true_type { };
template <typename Value>    struct is_native<Value, 4, enable_if_int32_t<Value>> : std::true_type { };
template <typename Value>    struct is_native<Value, 3, enable_if_int32_t<Value>> : std::true_type { };
template <typename Value>    struct is_native<Value, 2, enable_if_int64_t<Value>> : std::true_type { };

NAMESPACE_END(detail)

/// Partial overload of StaticArrayImpl using SSE4.2 intrinsics (single precision)
template <bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<float, 4,  IsMask_, Derived_>
  : StaticArrayBase<float, 4, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(float, 4, __m128)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm_set1_ps(value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3)
        : m(_mm_setr_ps(v0, v1, v2, v3)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_F16C)
    ENOKI_CONVERT(half) {
        m = _mm_cvtph_ps(_mm_loadl_epi64((const __m128i *) a.derived().data()));
    }
#endif

    ENOKI_CONVERT(float) : m(a.derived().m) { }
    ENOKI_CONVERT(int32_t) : m(_mm_cvtepi32_ps(a.derived().m)) { }

    ENOKI_CONVERT(uint32_t) {
        #if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            m = _mm_cvtepu32_ps(a.derived().m);
        #else
            int32_array_t<Derived> ai(a);
            Derived result =
                Derived(ai & 0x7fffffff) +
                (Derived(float(1u << 31)) & mask_t<Derived>(sr<31>(ai)));
            m = result.m;
        #endif
    }

#if defined(ENOKI_X86_AVX)
    ENOKI_CONVERT(double) : m(_mm256_cvtpd_ps(a.derived().m)) { }
#else
    ENOKI_CONVERT(double)
        : m(_mm_shuffle_ps(_mm_cvtpd_ps(low(a).m), _mm_cvtpd_ps(high(a).m),
                           _MM_SHUFFLE(1, 0, 1, 0))) { }
#endif

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_CONVERT(int64_t) : m(_mm256_cvtepi64_ps(a.derived().m)) { }
    ENOKI_CONVERT(uint64_t) : m(_mm256_cvtepu64_ps(a.derived().m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        int ival;
        memcpy(&ival, a.derived().data(), 4);
        m = _mm_castsi128_ps(_mm_cvtepi8_epi32(
            _mm_cmpgt_epi8(_mm_cvtsi32_si128(ival), _mm_setzero_si128())));
    }

    ENOKI_REINTERPRET(float) : m(a.derived().m) { }
    ENOKI_REINTERPRET(int32_t) : m(_mm_castsi128_ps(a.derived().m)) { }
    ENOKI_REINTERPRET(uint32_t) : m(_mm_castsi128_ps(a.derived().m)) { }

#if defined(ENOKI_X86_AVX)
    ENOKI_REINTERPRET(double)
        : m(_mm_castsi128_ps(
              detail::mm256_cvtepi64_epi32(_mm256_castpd_si256(a.derived().m)))) { }
#else
    ENOKI_REINTERPRET(double)
        : m(_mm_castsi128_ps(detail::mm256_cvtepi64_epi32(
              _mm_castpd_si128(low(a).m), _mm_castpd_si128(high(a).m)))) { }
#endif

#if defined(ENOKI_X86_AVX2)
    ENOKI_REINTERPRET(uint64_t)
        : m(_mm_castsi128_ps(
              detail::mm256_cvtepi64_epi32(a.derived().m))) { }
    ENOKI_REINTERPRET(int64_t)
        : m(_mm_castsi128_ps(
              detail::mm256_cvtepi64_epi32(a.derived().m))) { }
#else
    ENOKI_REINTERPRET(uint64_t)
        : m(_mm_castsi128_ps(
              detail::mm256_cvtepi64_epi32(low(a).m, high(a).m))) { }
    ENOKI_REINTERPRET(int64_t)
        : m(_mm_castsi128_ps(
              detail::mm256_cvtepi64_epi32(low(a).m, high(a).m))) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(_mm_setr_ps(a1.coeff(0), a1.coeff(1), a2.coeff(0), a2.coeff(1))) { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0), coeff(1)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(2), coeff(3)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const    { return _mm_add_ps(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const    { return _mm_sub_ps(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const    { return _mm_mul_ps(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const    { return _mm_div_ps(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_ps(m, a.k, _mm_set1_ps(memcpy_cast<Value>(int32_t(-1))));
            else
        #endif
        return _mm_or_ps(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_maskz_mov_ps(a.k, m);
            else
        #endif
        return _mm_and_ps(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_xor_ps(m, a.k, m, _mm_set1_ps(memcpy_cast<Value>(int32_t(-1))));
            else
        #endif
        return _mm_xor_ps(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_ps(m, a.k, _mm_setzero_ps());
            else
        #endif
        return _mm_andnot_ps(a.m, m);
    }

    #if defined(ENOKI_X86_AVX512VL)
        #define ENOKI_COMP(name, NAME) mask_t<Derived>::from_k(_mm_cmp_ps_mask(m, a.m, _CMP_##NAME))
    #elif defined(ENOKI_X86_AVX)
        #define ENOKI_COMP(name, NAME) mask_t<Derived>(_mm_cmp_ps(m, a.m, _CMP_##NAME))
    #else
        #define ENOKI_COMP(name, NAME) mask_t<Derived>(_mm_cmp##name##_ps(m, a.m))
    #endif

    ENOKI_INLINE auto lt_ (Ref a) const { return ENOKI_COMP(lt,  LT_OQ);  }
    ENOKI_INLINE auto gt_ (Ref a) const { return ENOKI_COMP(gt,  GT_OQ);  }
    ENOKI_INLINE auto le_ (Ref a) const { return ENOKI_COMP(le,  LE_OQ);  }
    ENOKI_INLINE auto ge_ (Ref a) const { return ENOKI_COMP(ge,  GE_OQ);  }
    ENOKI_INLINE auto eq_ (Ref a) const {
        using Int = int_array_t<Derived>;
        if constexpr (IsMask_)
            return mask_t<Derived>(eq(Int(derived()), Int(a)));
        else
            return ENOKI_COMP(eq, EQ_OQ);
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        using Int = int_array_t<Derived>;
        if constexpr (IsMask_)
            return mask_t<Derived>(neq(Int(derived()), Int(a)));
        else
            return ENOKI_COMP(neq, NEQ_UQ);
    }

    #undef ENOKI_COMP

    ENOKI_INLINE Derived abs_()      const { return _mm_andnot_ps(_mm_set1_ps(-0.f), m); }
    ENOKI_INLINE Derived min_(Ref b) const { return _mm_min_ps(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return _mm_max_ps(b.m, m); }
    ENOKI_INLINE Derived ceil_()     const { return _mm_ceil_ps(m);     }
    ENOKI_INLINE Derived floor_()    const { return _mm_floor_ps(m);    }
    ENOKI_INLINE Derived sqrt_()     const { return _mm_sqrt_ps(m);     }

    ENOKI_INLINE Derived round_() const {
        return _mm_round_ps(m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    ENOKI_INLINE Derived trunc_() const {
        return _mm_round_ps(m, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        #if !defined(ENOKI_X86_AVX512VL)
            return _mm_blendv_ps(f.m, t.m, m.m);
        #else
            return _mm_mask_blend_ps(m.k, f.m, t.m);
        #endif
    }

#if defined(ENOKI_X86_FMA)
    ENOKI_INLINE Derived fmadd_   (Ref b, Ref c) const { return _mm_fmadd_ps   (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsub_   (Ref b, Ref c) const { return _mm_fmsub_ps   (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmadd_  (Ref b, Ref c) const { return _mm_fnmadd_ps  (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmsub_  (Ref b, Ref c) const { return _mm_fnmsub_ps  (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsubadd_(Ref b, Ref c) const { return _mm_fmsubadd_ps(m, b.m, c.m); }
    ENOKI_INLINE Derived fmaddsub_(Ref b, Ref c) const { return _mm_fmaddsub_ps(m, b.m, c.m); }
#endif

    template <int I0, int I1, int I2, int I3>
    ENOKI_INLINE Derived shuffle_() const {
        #if defined(ENOKI_X86_AVX)
            return _mm_permute_ps(m, _MM_SHUFFLE(I3, I2, I1, I0));
        #else
            return _mm_shuffle_ps(m, m, _MM_SHUFFLE(I3, I2, I1, I0));
        #endif
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        #if defined(ENOKI_X86_AVX)
            return _mm_permutevar_ps(m, index.m);
        #else
            return Base::shuffle_(index);
        #endif
    }

#if defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived ldexp_(Ref arg) const { return _mm_scalef_ps(m, arg.m); }

    ENOKI_INLINE std::pair<Derived, Derived> frexp_() const {
        return std::make_pair<Derived, Derived>(
            _mm_getmant_ps(m, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src),
            _mm_getexp_ps(m));
    }
#endif

    ENOKI_INLINE Derived rcp_() const {
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28, use as is (even in non-approximate mode) */
            return _mm512_castps512_ps128(
                _mm512_rcp28_ps(_mm512_castps128_ps512(m)));
        #else
            /* Use best reciprocal approximation available on the current
               hardware and refine */
            __m128 r;
            #if defined(ENOKI_X86_AVX512VL)
                r = _mm_rcp14_ps(m); /* rel error < 2^-14 */
            #else
                r = _mm_rcp_ps(m);   /* rel error < 1.5*2^-12 */
            #endif

            /* Refine using one Newton-Raphson iteration */
            __m128 t0 = _mm_add_ps(r, r),
                   t1 = _mm_mul_ps(r, m),
                   ro = r;
            (void) ro;

            #if defined(ENOKI_X86_FMA)
                r = _mm_fnmadd_ps(t1, r, t0);
            #else
                r = _mm_sub_ps(t0, _mm_mul_ps(r, t1));
            #endif

            #if defined(ENOKI_X86_AVX512VL)
                return _mm_fixupimm_ps(r, m, _mm_set1_epi32(0x0087A622), 0);
            #else
                return _mm_blendv_ps(r, ro, t1); /* mask bit is '1' iff t1 == nan */
            #endif
        #endif
    }

    ENOKI_INLINE Derived rsqrt_() const {
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28, use as is (even in non-approximate mode) */
            return _mm512_castps512_ps128(
                _mm512_rsqrt28_ps(_mm512_castps128_ps512(m)));
        #else
            /* Use best reciprocal square root approximation available
               on the current hardware and refine */
            __m128 r;
            #if defined(ENOKI_X86_AVX512VL)
                r = _mm_rsqrt14_ps(m); /* rel error < 2^-14 */
            #else
                r = _mm_rsqrt_ps(m);   /* rel error < 1.5*2^-12 */
            #endif

            /* Refine using one Newton-Raphson iteration */
            const __m128 c0 = _mm_set1_ps(.5f),
                         c1 = _mm_set1_ps(3.f);

            __m128 t0 = _mm_mul_ps(r, c0),
                   t1 = _mm_mul_ps(r, m),
                   ro = r;
            (void) ro;

            #if defined(ENOKI_X86_FMA)
                r = _mm_mul_ps(_mm_fnmadd_ps(t1, r, c1), t0);
            #else
                r = _mm_mul_ps(_mm_sub_ps(c1, _mm_mul_ps(t1, r)), t0);
            #endif

            #if defined(ENOKI_X86_AVX512VL)
                return _mm_fixupimm_ps(r, m, _mm_set1_epi32(0x0383A622), 0);
            #else
                return _mm_blendv_ps(r, ro, t1); /* mask bit is '1' iff t1 == nan */
            #endif
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op)                                    \
        ENOKI_INLINE Value name##_() const {                                 \
            __m128 t1 = _mm_movehdup_ps(m);                                  \
            __m128 t2 = _mm_##op##_ps(m, t1);                                \
            t1 = _mm_movehl_ps(t1, t2);                                      \
            t2 = _mm_##op##_ss(t2, t1);                                      \
            return _mm_cvtss_f32(t2);                                        \
        }

    ENOKI_HORIZONTAL_OP(hsum, add)
    ENOKI_HORIZONTAL_OP(hprod, mul)
    ENOKI_HORIZONTAL_OP(hmin, min)
    ENOKI_HORIZONTAL_OP(hmax, max)

    #undef ENOKI_HORIZONTAL_OP

    ENOKI_INLINE bool all_()  const { return _mm_movemask_ps(m) == 0xF;}
    ENOKI_INLINE bool any_()  const { return _mm_movemask_ps(m) != 0x0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_ps(m); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    ENOKI_INLINE Value dot_(Ref a) const {
        return _mm_cvtss_f32(_mm_dp_ps(m, a.m, 0b11110001));
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm_mask_mov_ps(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm_mask_add_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm_mask_sub_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm_mask_mul_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mdiv_   (const Derived &a, const Mask &mask) { m = _mm_mask_div_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm_mask_or_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm_mask_and_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm_mask_xor_ps(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
      assert((uintptr_t) ptr % 16 == 0);
        _mm_store_ps((Value *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_store_ps((Value *) ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX)
            _mm_maskstore_ps((Value *) ptr, _mm_castps_si128(mask.m), m);
        #else
            Base::store_(ptr, mask);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm_storeu_ps((Value *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_storeu_ps((Value *) ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX)
            _mm_maskstore_ps((Value *) ptr, _mm_castps_si128(mask.m), m);
        #else
            Base::store_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
      assert((uintptr_t) ptr % 16 == 0);
        return _mm_load_ps((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_load_ps(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX)
            return _mm_maskload_ps((const Value *) ptr, _mm_castps_si128(mask.m));
        #else
            return Base::load_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm_loadu_ps((const Value *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_loadu_ps(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX)
            return _mm_maskload_ps((const Value *) ptr, _mm_castps_si128(mask.m));
        #else
            return Base::load_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm_setzero_ps(); }

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm_mmask_i32gather_ps(_mm_setzero_ps(), mask.k, index.m, (const float *) ptr, Stride);
            else
                return _mm256_mmask_i64gather_ps(_mm_setzero_ps(), mask.k, index.m, (const float *) ptr, Stride);
        #else
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm_mask_i32gather_ps(_mm_setzero_ps(), (const float *) ptr, index.m, mask.m, Stride);
            else
                return _mm256_mask_i64gather_ps(_mm_setzero_ps(), (const float *) ptr, index.m, mask.m, Stride);
        #endif
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm_mask_i32scatter_ps(ptr, mask.k, index.m, m, Stride);
        else
            _mm256_mask_i64scatter_ps(ptr, mask.k, index.m, m, Stride);
    }
#endif

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        #if !defined(ENOKI_X86_AVX512VL)
            unsigned int k = (unsigned int) _mm_movemask_ps(mask.m);
            return coeff((size_t) (detail::tzcnt_scalar(k) & 3));
        #else
            return _mm_cvtss_f32(_mm_mask_compress_ps(_mm_setzero_ps(), mask.k, m));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(float *&ptr, const Mask &mask) const {
        #if !defined(ENOKI_X86_AVX512VL)
            unsigned int k = (unsigned int) _mm_movemask_ps(mask.m);

            /** Fancy LUT-based partitioning algorithm, see
                https://deplinenoise.files.wordpress.com/2015/03/gdc2015_afredriksson_simd.pdf */

            __m128i shuf = _mm_load_si128(((const __m128i *) detail::compress_lut_128) + k),
                    perm = _mm_shuffle_epi8(_mm_castps_si128(m), shuf);

            _mm_storeu_si128((__m128i *) ptr, perm);
        #else
             unsigned int k = (unsigned int) mask.k;
            _mm_storeu_ps(ptr, _mm_mask_compress_ps(_mm_setzero_ps(), mask.k, m));
        #endif

        size_t kn = (size_t) _mm_popcnt_u32(k);
        ptr += kn;
        return kn;
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using SSE4.2 intrinsics (double precision)
template <bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<double, 2, IsMask_, Derived_>
  : StaticArrayBase<double, 2, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(double, 2, __m128d)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm_set1_pd(value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1)
        : m(_mm_setr_pd(v0, v1)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    /* No vectorized conversions from float/[u]int32_t (too small) */

    ENOKI_CONVERT(double) : m(a.derived().m) { }

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_CONVERT(int64_t) : m(_mm_cvtepi64_pd(a.derived().m)) { }
    ENOKI_CONVERT(uint64_t) : m(_mm_cvtepu64_pd(a.derived().m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        int16_t ival;
        memcpy(&ival, a.derived().data(), 2);
        m = _mm_castsi128_pd(_mm_cvtepi8_epi64(_mm_cmpgt_epi8(
            _mm_cvtsi32_si128((int) ival), _mm_setzero_si128())));
    }

    ENOKI_REINTERPRET(float) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, float32[2] -> double[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_castps_pd(_mm_setr_ps(v0, v0, v1, v1));
    }

    ENOKI_REINTERPRET(int32_t) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, int32[2] -> double[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_castsi128_pd(_mm_setr_epi32(v0, v0, v1, v1));
    }

    ENOKI_REINTERPRET(uint32_t) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, uint32[2] -> double[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_castsi128_pd(_mm_setr_epi32((int32_t) v0, (int32_t) v0,
                                            (int32_t) v1, (int32_t) v1));
    }

    ENOKI_REINTERPRET(double) : m(a.derived().m) { }
    ENOKI_REINTERPRET(int64_t) : m(_mm_castsi128_pd(a.derived().m)) { }
    ENOKI_REINTERPRET(uint64_t) : m(_mm_castsi128_pd(a.derived().m)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(_mm_setr_pd(a1.coeff(0), a2.coeff(0))) { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(1)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm_add_pd(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm_sub_pd(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm_mul_pd(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const { return _mm_div_pd(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_pd(m, a.k, _mm_set1_pd(memcpy_cast<Value>(int64_t(-1))));
            else
        #endif
        return _mm_or_pd(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_maskz_mov_pd(a.k, m);
            else
        #endif
        return _mm_and_pd(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_xor_pd(m, a.k, m, _mm_set1_pd(memcpy_cast<Value>(int64_t(-1))));
            else
        #endif
        return _mm_xor_pd(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_pd(m, a.k, _mm_setzero_pd());
            else
        #endif
        return _mm_andnot_pd(a.m, m);
    }

    #if defined(ENOKI_X86_AVX512VL)
        #define ENOKI_COMP(name, NAME) mask_t<Derived>::from_k(_mm_cmp_pd_mask(m, a.m, _CMP_##NAME))
    #elif defined(ENOKI_X86_AVX)
        #define ENOKI_COMP(name, NAME) mask_t<Derived>(_mm_cmp_pd(m, a.m, _CMP_##NAME))
    #else
        #define ENOKI_COMP(name, NAME) mask_t<Derived>(_mm_cmp##name##_pd(m, a.m))
    #endif

    ENOKI_INLINE auto lt_ (Ref a) const { return ENOKI_COMP(lt,  LT_OQ);  }
    ENOKI_INLINE auto gt_ (Ref a) const { return ENOKI_COMP(gt,  GT_OQ);  }
    ENOKI_INLINE auto le_ (Ref a) const { return ENOKI_COMP(le,  LE_OQ);  }
    ENOKI_INLINE auto ge_ (Ref a) const { return ENOKI_COMP(ge,  GE_OQ);  }

    ENOKI_INLINE auto eq_ (Ref a) const {
        using Int = int_array_t<Derived>;
        if constexpr (IsMask_)
            return mask_t<Derived>(eq(Int(derived()), Int(a)));
        else
            return ENOKI_COMP(eq, EQ_OQ);
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        using Int = int_array_t<Derived>;
        if constexpr (IsMask_)
            return mask_t<Derived>(neq(Int(derived()), Int(a)));
        else
            return ENOKI_COMP(neq, NEQ_UQ);
    }

    #undef ENOKI_COMP

    ENOKI_INLINE Derived abs_()      const { return _mm_andnot_pd(_mm_set1_pd(-0.), m); }
    ENOKI_INLINE Derived min_(Ref b) const { return _mm_min_pd(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return _mm_max_pd(b.m, m); }
    ENOKI_INLINE Derived ceil_()     const { return _mm_ceil_pd(m);     }
    ENOKI_INLINE Derived floor_()    const { return _mm_floor_pd(m);    }
    ENOKI_INLINE Derived sqrt_()     const { return _mm_sqrt_pd(m);     }

    ENOKI_INLINE Derived round_() const {
        return _mm_round_pd(m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    ENOKI_INLINE Derived trunc_() const {
        return _mm_round_pd(m, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        #if !defined(ENOKI_X86_AVX512VL)
            return _mm_blendv_pd(f.m, t.m, m.m);
        #else
            return _mm_mask_blend_pd(m.k, f.m, t.m);
        #endif
    }

#if defined(ENOKI_X86_FMA)
    ENOKI_INLINE Derived fmadd_   (Ref b, Ref c) const { return _mm_fmadd_pd   (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsub_   (Ref b, Ref c) const { return _mm_fmsub_pd   (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmadd_  (Ref b, Ref c) const { return _mm_fnmadd_pd  (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmsub_  (Ref b, Ref c) const { return _mm_fnmsub_pd  (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsubadd_(Ref b, Ref c) const { return _mm_fmsubadd_pd(m, b.m, c.m); }
    ENOKI_INLINE Derived fmaddsub_(Ref b, Ref c) const { return _mm_fmaddsub_pd(m, b.m, c.m); }
#endif

    #if defined(ENOKI_X86_AVX)
        #define ENOKI_SHUFFLE_PD(m, flags) _mm_permute_pd(m, flags)
    #else
        #define ENOKI_SHUFFLE_PD(m, flags) _mm_shuffle_pd(m, m, flags)
    #endif

    template <int I0, int I1>
    ENOKI_INLINE Derived shuffle_() const {
        return ENOKI_SHUFFLE_PD(m, (I1 << 1) | I0);
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        #if defined(ENOKI_X86_AVX)
            return _mm_permutevar_pd(m, _mm_slli_epi64(index.m, 1));
        #else
            return Base::shuffle_(index);
        #endif
    }

#if defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived ldexp_(Ref arg) const { return _mm_scalef_pd(m, arg.m); }

    ENOKI_INLINE std::pair<Derived, Derived> frexp_() const {
        return std::make_pair<Derived, Derived>(
            _mm_getmant_pd(m, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src),
            _mm_getexp_pd(m));
    }
#endif

#if defined(ENOKI_X86_AVX512VL) || defined(ENOKI_X86_AVX512ER)
    ENOKI_INLINE Derived rcp_() const {
        /* Use best reciprocal approximation available on the current
           hardware and refine */
        __m128d r;
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28 */
            r = _mm512_castpd512_pd128(
                _mm512_rcp28_pd(_mm512_castpd128_pd512(m)));
        #elif defined(ENOKI_X86_AVX512VL)
            r = _mm_rcp14_pd(m); /* rel error < 2^-14 */
        #endif

        __m128d ro = r, t0, t1;
        (void) ro;

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            t0 = _mm_add_pd(r, r);
            t1 = _mm_mul_pd(r, m);
            r = _mm_fnmadd_pd(t1, r, t0);
        }

        #if defined(ENOKI_X86_AVX512VL)
            return _mm_fixupimm_pd(r, m, _mm_set1_epi32(0x0087A622), 0);
        #else
            return _mm_blendv_pd(r, ro, t1); /* mask bit is '1' iff t1 == nan */
        #endif
    }

    ENOKI_INLINE Derived rsqrt_() const {
        /* Use best reciprocal square root approximation available
           on the current hardware and refine */
        __m128d r;
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28 */
            r = _mm512_castpd512_pd128(
                _mm512_rsqrt28_pd(_mm512_castpd128_pd512(m)));
        #elif defined(ENOKI_X86_AVX512VL)
            r = _mm_rsqrt14_pd(m); /* rel error < 2^-14 */
        #endif

        const __m128d c0 = _mm_set1_pd(0.5),
                      c1 = _mm_set1_pd(3.0);

        __m128d ro = r, t0, t1;
        (void) ro;

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            t0 = _mm_mul_pd(r, c0);
            t1 = _mm_mul_pd(r, m);
            r = _mm_mul_pd(_mm_fnmadd_pd(t1, r, c1), t0);
        }

        #if defined(ENOKI_X86_AVX512VL)
            return _mm_fixupimm_pd(r, m, _mm_set1_epi32(0x0383A622), 0);
        #else
            return _mm_blendv_pd(r, ro, t1); /* mask bit is '1' iff t1 == nan */
        #endif
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op) \
        ENOKI_INLINE Value name##_() const { \
            __m128d t0 = ENOKI_SHUFFLE_PD(m, 1); \
            __m128d t1 = _mm_##op##_sd(t0, m); \
            return  _mm_cvtsd_f64(t1); \
        }

    ENOKI_HORIZONTAL_OP(hsum, add)
    ENOKI_HORIZONTAL_OP(hprod, mul)
    ENOKI_HORIZONTAL_OP(hmin, min)
    ENOKI_HORIZONTAL_OP(hmax, max)

    #undef ENOKI_HORIZONTAL_OP
    #undef ENOKI_SHUFFLE_PD

    ENOKI_INLINE bool all_()  const { return _mm_movemask_pd(m) == 0x3;}
    ENOKI_INLINE bool any_()  const { return _mm_movemask_pd(m) != 0x0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_pd(m); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    ENOKI_INLINE Value dot_(Ref a) const {
        return _mm_cvtsd_f64(_mm_dp_pd(m, a.m, 0b00110001));
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm_mask_mov_pd(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm_mask_add_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm_mask_sub_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm_mask_mul_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mdiv_   (const Derived &a, const Mask &mask) { m = _mm_mask_div_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm_mask_or_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm_mask_and_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm_mask_xor_pd(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        _mm_store_pd((Value *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }
    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_store_pd((Value *) ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX)
            _mm_maskstore_pd((Value *) ptr, _mm_castpd_si128(mask.m), m);
        #else
            Base::store_(ptr, mask);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm_storeu_pd((Value *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_storeu_pd((Value *) ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX)
            _mm_maskstore_pd((Value *) ptr, _mm_castpd_si128(mask.m), m);
        #else
            Base::store_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 16 == 0);
        return _mm_load_pd((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_load_pd(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX)
            return _mm_maskload_pd((const Value *) ptr, _mm_castpd_si128(mask.m));
        #else
            return Base::load_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm_loadu_pd((const Value *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_loadu_pd(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX)
            return _mm_maskload_pd((const Value *) ptr, _mm_castpd_si128(mask.m));
        #else
            return Base::load_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm_setzero_pd(); }

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            return Base::template gather_<Stride>(ptr, index, mask);
        } else {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_mmask_i64gather_pd(_mm_setzero_pd(), mask.k, index.m, (const double *) ptr, Stride);
            #else
                return _mm_mask_i64gather_pd(_mm_setzero_pd(), (const double *) ptr, index.m, mask.m, Stride);
            #endif
        }
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            Base::template scatter_<Stride>(ptr, index, mask);
        else
            _mm_mask_i64scatter_pd(ptr, mask.k, index.m, m, Stride);
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return _mm_cvtsd_f64(_mm_mask_compress_pd(_mm_setzero_pd(), mask.k, m));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(double *&ptr, const Mask &mask) const {
        _mm_storeu_pd(ptr, _mm_mask_compress_pd(_mm_setzero_pd(), mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }
#endif

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using SSE4.2 intrinsics (32 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<Value_, 4, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayBase<Value_, 4, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 4, __m128i)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm_set1_epi32((int32_t) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3)
        : m(_mm_setr_epi32((int32_t) v0, (int32_t) v1, (int32_t) v2, (int32_t) v3)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(float) {
        if constexpr (std::is_signed_v<Value>) {
            m = _mm_cvttps_epi32(a.derived().m);
        } else {
#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            m = _mm_cvttps_epu32(a.derived().m);
#else
            constexpr uint32_t limit = 1u << 31;
            const __m128  limit_f = _mm_set1_ps((float) limit);
            const __m128i limit_i = _mm_set1_epi32((int) limit);

            __m128 v = a.derived().m;

            __m128i mask =
                _mm_castps_si128(_mm_cmpge_ps(v, limit_f));

            __m128i b2 = _mm_add_epi32(
                _mm_cvttps_epi32(_mm_sub_ps(v, limit_f)), limit_i);

            __m128i b1 = _mm_cvttps_epi32(v);

            m = _mm_blendv_epi8(b1, b2, mask);
#endif
        }
    }

    ENOKI_CONVERT(int32_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint32_t) : m(a.derived().m) { }

#if defined(ENOKI_X86_AVX)
    ENOKI_CONVERT(double) {
        if constexpr (std::is_signed_v<Value>) {
            m = _mm256_cvttpd_epi32(a.derived().m);
        } else {
#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            m = _mm256_cvttpd_epu32(a.derived().m);
#else
            ENOKI_TRACK_SCALAR("Constructor (converting, double[4] -> uint32[4])");
            for (size_t i = 0; i < Size; ++i)
                coeff(i) = Value(a.derived().coeff(i));
#endif
        }
    }
#endif

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_CONVERT(int64_t) { m = _mm256_cvtepi64_epi32(a.derived().m); }
    ENOKI_CONVERT(uint64_t) { m = _mm256_cvtepi64_epi32(a.derived().m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        int ival;
        memcpy(&ival, a.derived().data(), 4);
        m = _mm_cvtepi8_epi32(
            _mm_cmpgt_epi8(_mm_cvtsi32_si128(ival), _mm_setzero_si128()));
    }

    ENOKI_REINTERPRET(float) : m(_mm_castps_si128(a.derived().m)) { }
    ENOKI_REINTERPRET(int32_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint32_t) : m(a.derived().m) { }

#if defined(ENOKI_X86_AVX)
    ENOKI_REINTERPRET(double)
        : m(detail::mm256_cvtepi64_epi32(_mm256_castpd_si256(a.derived().m))) { }
#else
    ENOKI_REINTERPRET(double)
        : m(detail::mm256_cvtepi64_epi32(_mm_castpd_si128(low(a).m),
                                         _mm_castpd_si128(high(a).m))) { }
#endif

#if defined(ENOKI_X86_AVX2)
    ENOKI_REINTERPRET(uint64_t)
        : m(detail::mm256_cvtepi64_epi32(a.derived().m)) { }
    ENOKI_REINTERPRET(int64_t)
        : m(detail::mm256_cvtepi64_epi32(a.derived().m)) {}
#else
    ENOKI_REINTERPRET(uint64_t)
        : m(detail::mm256_cvtepi64_epi32(low(a).m, high(a).m)) { }
    ENOKI_REINTERPRET(int64_t)
        : m(detail::mm256_cvtepi64_epi32(low(a).m, high(a).m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(_mm_setr_epi32((int32_t) a1.coeff(0), (int32_t) a1.coeff(1),
                           (int32_t) a2.coeff(0), (int32_t) a2.coeff(1))) { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0), coeff(1)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(2), coeff(3)); }

    //! @}
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm_add_epi32(m, a.m);   }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm_sub_epi32(m, a.m);   }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm_mullo_epi32(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_epi32(m, a.k, _mm_set1_epi32(-1));
            else
        #endif
        return _mm_or_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_maskz_mov_epi32(a.k, m);
            else
        #endif
        return _mm_and_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_xor_epi32(m, a.k, m, _mm_set1_epi32(-1));
            else
        #endif
        return _mm_xor_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_epi32(m, a.k, _mm_setzero_si128());
            else
        #endif
        return _mm_andnot_si128(a.m, m);
    }

    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        return _mm_slli_epi32(m, (int) Imm);
    }

    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        return std::is_signed_v<Value> ? _mm_srai_epi32(m, (int) Imm)
                                       : _mm_srli_epi32(m, (int) Imm);
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm_sll_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        return std::is_signed_v<Value>
                   ? _mm_sra_epi32(m, _mm_set1_epi64x((long long) k))
                   : _mm_srl_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        #if defined(ENOKI_X86_AVX2)
            return _mm_sllv_epi32(m, k.m);
        #else
            Derived out;
            ENOKI_TRACK_SCALAR("sl");
            for (size_t i = 0; i < Size; ++i)
                out.coeff(i) = coeff(i) << (size_t) k.coeff(i);
            return out;
        #endif
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        #if defined(ENOKI_X86_AVX2)
            return std::is_signed_v<Value> ? _mm_srav_epi32(m, k.m)
                                           : _mm_srlv_epi32(m, k.m);
        #else
            Derived out;
            ENOKI_TRACK_SCALAR("sr");
            for (size_t i = 0; i < Size; ++i)
                out.coeff(i) = coeff(i) >> (size_t) k.coeff(i);
            return out;
        #endif
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Imm> ENOKI_INLINE Derived rol_() const { return _mm_rol_epi32(m, (int) Imm); }
    template <size_t Imm> ENOKI_INLINE Derived ror_() const { return _mm_ror_epi32(m, (int) Imm); }
    ENOKI_INLINE Derived rol_(Ref k) const { return _mm_rolv_epi32(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm_rorv_epi32(m, k.m); }
#endif

    ENOKI_INLINE auto eq_(Ref a)  const {
        using Return = mask_t<Derived>;

        #if defined(ENOKI_X86_AVX512VL)
            return Return::from_k(_mm_cmpeq_epi32_mask(m, a.m));
        #else
            return Return(_mm_cmpeq_epi32(m, a.m));
        #endif
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(_mm_cmpneq_epi32_mask(m, a.m));
        #else
            return ~eq_(a);
        #endif
    }

    ENOKI_INLINE auto lt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm_cmpgt_epi32(a.m, m));
            } else {
                const __m128i offset = _mm_set1_epi32((int32_t) 0x80000000ul);
                return Return(_mm_cmpgt_epi32(_mm_sub_epi32(a.m, offset),
                                              _mm_sub_epi32(m, offset)));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                      ? _mm_cmplt_epi32_mask(m, a.m)
                                      : _mm_cmplt_epu32_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm_cmpgt_epi32(m, a.m));
            } else {
                const __m128i offset = _mm_set1_epi32((int32_t) 0x80000000ul);
                return Return(_mm_cmpgt_epi32(_mm_sub_epi32(m, offset),
                                              _mm_sub_epi32(a.m, offset)));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm_cmpgt_epi32_mask(m, a.m)
                                  : _mm_cmpgt_epu32_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto le_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm_cmple_epi32_mask(m, a.m)
                                           : _mm_cmple_epu32_mask(m, a.m));
        #else
            return ~gt_(a);
        #endif
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm_cmpge_epi32_mask(m, a.m)
                                           : _mm_cmpge_epu32_mask(m, a.m));
        #else
            return ~lt_(a);
        #endif
    }

    ENOKI_INLINE Derived min_(Ref a) const {
        return std::is_signed_v<Value> ? _mm_min_epi32(a.m, m)
                                       : _mm_min_epu32(a.m, m);
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        return std::is_signed_v<Value> ? _mm_max_epi32(a.m, m)
                                       : _mm_max_epu32(a.m, m);
    }

    ENOKI_INLINE Derived abs_() const {
        return std::is_signed_v<Value> ? _mm_abs_epi32(m) : m;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        #if !defined(ENOKI_X86_AVX512VL)
            return _mm_blendv_epi8(f.m, t.m, m.m);
        #else
            return _mm_mask_blend_epi32(m.k, f.m, t.m);
        #endif
    }

    template <int I0, int I1, int I2, int I3>
    ENOKI_INLINE Derived shuffle_() const {
        return _mm_shuffle_epi32(m, _MM_SHUFFLE(I3, I2, I1, I0));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        #if defined(ENOKI_X86_AVX)
            return _mm_castps_si128(_mm_permutevar_ps(_mm_castsi128_ps(m), index.m));
        #else
            return Base::shuffle_(index);
        #endif
    }

    ENOKI_INLINE Derived mulhi_(Ref a) const {
        Derived even, odd;
        if constexpr (std::is_signed_v<Value>) {
            even.m = _mm_srli_epi64(_mm_mul_epi32(m, a.m), 32);
            odd.m = _mm_mul_epi32(_mm_srli_epi64(m, 32), _mm_srli_epi64(a.m, 32));
        } else {
            even.m = _mm_srli_epi64(_mm_mul_epu32(m, a.m), 32);
            odd.m = _mm_mul_epu32(_mm_srli_epi64(m, 32), _mm_srli_epi64(a.m, 32));
        }

        #if defined(ENOKI_X86_AVX512VL)
            const mask_t<Derived> blend = mask_t<Derived>::from_k(0b0101);
        #else
            const mask_t<Derived> blend(Value(-1), Value(0), Value(-1), Value(0));
        #endif

        return select(blend, even, odd);
    }

#if defined(ENOKI_X86_AVX512CD) && defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived lzcnt_() const { return _mm_lzcnt_epi32(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(32) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op)                                     \
        ENOKI_INLINE Value name##_() const {                                  \
            __m128i t1 = _mm_shuffle_epi32(m, 0x4e);                          \
            __m128i t2 = _mm_##op##_epi32(m, t1);                             \
            t1 = _mm_shufflelo_epi16(t2, 0x4e);                               \
            t2 = _mm_##op##_epi32(t2, t1);                                    \
            return (Value) _mm_cvtsi128_si32(t2);                             \
        }

    #define ENOKI_HORIZONTAL_OP_SIGNED(name, op)                              \
        ENOKI_INLINE Value name##_() const {                                  \
            __m128i t1 = _mm_shuffle_epi32(m, 0x4e);                          \
            __m128i t2 = std::is_signed_v<Value> ? _mm_##op##_epi32(m, t1) :  \
                                                   _mm_##op##_epu32(m, t1);   \
            t1 = _mm_shufflelo_epi16(t2, 0x4e);                               \
            t2 = std::is_signed_v<Value> ? _mm_##op##_epi32(t2, t1) :         \
                                           _mm_##op##_epu32(t2, t1);          \
            return (Value) _mm_cvtsi128_si32(t2);                             \
        }

    ENOKI_HORIZONTAL_OP(hsum, add)
    ENOKI_HORIZONTAL_OP(hprod, mullo)
    ENOKI_HORIZONTAL_OP_SIGNED(hmin, min)
    ENOKI_HORIZONTAL_OP_SIGNED(hmax, max)

    #undef ENOKI_HORIZONTAL_OP
    #undef ENOKI_HORIZONTAL_OP_SIGNED

    ENOKI_INLINE bool all_()  const { return _mm_movemask_ps(_mm_castsi128_ps(m)) == 0xF;}
    ENOKI_INLINE bool any_()  const { return _mm_movemask_ps(_mm_castsi128_ps(m)) != 0x0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_ps(_mm_castsi128_ps(m)); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm_mask_mov_epi32(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm_mask_add_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm_mask_sub_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm_mask_mullo_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm_mask_or_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm_mask_and_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm_mask_xor_epi32(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        _mm_store_si128((__m128i *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_store_epi32(ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX2)
            _mm_maskstore_epi32((int *) ptr, mask.m, m);
        #else
            Base::store_(ptr, mask);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm_storeu_si128((__m128i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_storeu_epi32(ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX2)
            _mm_maskstore_epi32((int *) ptr, mask.m, m);
        #else
            Base::store_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 16 == 0);
        return _mm_load_si128((const __m128i *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }
    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_load_epi32(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX2)
            return _mm_maskload_epi32((const int *) ptr, mask.m);
        #else
            return Base::load_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm_loadu_si128((const __m128i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_loadu_epi32(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX2)
            return _mm_maskload_epi32((const int *) ptr, mask.m);
        #else
            return Base::load_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm_setzero_si128(); }

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm_mmask_i32gather_epi32(_mm_setzero_si128(), mask.k, index.m, (const int *) ptr, Stride);
            else
                return _mm256_mmask_i64gather_epi32(_mm_setzero_si128(), mask.k, index.m, (const int *) ptr, Stride);
        #else
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm_mask_i32gather_epi32(_mm_setzero_si128(), (const int *) ptr, index.m, mask.m, Stride);
            else
                return _mm256_mask_i64gather_epi32(_mm_setzero_si128(), (const int *) ptr, index.m, mask.m, Stride);
        #endif
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm_mask_i32scatter_epi32(ptr, mask.k, index.m, m, Stride);
        else
            _mm256_mask_i64scatter_epi32(ptr, mask.k, index.m, m, Stride);
    }
#endif

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        #if !defined(ENOKI_X86_AVX512VL)
            unsigned int k = (unsigned int) _mm_movemask_ps(_mm_castsi128_ps(mask.m));
            return coeff((size_t) (detail::tzcnt_scalar(k) & 3));
        #else
            return (Value) _mm_cvtsi128_si32(_mm_mask_compress_epi32(_mm_setzero_si128(), mask.k, m));
        #endif
    }

    template <typename T, typename Mask>
    ENOKI_INLINE size_t compress_(T *&ptr, const Mask &mask) const {
        #if !defined(ENOKI_X86_AVX512VL)
            unsigned int k = (unsigned int) _mm_movemask_ps(_mm_castsi128_ps(mask.m));

            /** Fancy LUT-based partitioning algorithm, see
                https://deplinenoise.files.wordpress.com/2015/03/gdc2015_afredriksson_simd.pdf */

            __m128i shuf = _mm_load_si128(((const __m128i *) detail::compress_lut_128) + k),
                    perm = _mm_shuffle_epi8(m, shuf);

            _mm_storeu_si128((__m128i *) ptr, perm);
        #else
            _mm_storeu_si128((__m128i *) ptr,
                _mm_mask_compress_epi32(_mm_setzero_si128(), mask.k, m));
            unsigned int k = (unsigned int) mask.k;
        #endif

        size_t kn = (size_t) _mm_popcnt_u32(k);
        ptr += kn;
        return kn;
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using SSE4.2 intrinsics (64 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<Value_, 2, IsMask_, Derived_, enable_if_int64_t<Value_>>
  : StaticArrayBase<Value_, 2, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 2, __m128i)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm_set1_epi64x((int64_t) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1) {
        alignas(16) Value data[2];
        data[0] = (Value) v0;
        data[1] = (Value) v1;
        m = _mm_load_si128((__m128i *) data);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_CONVERT(double) {
        if constexpr (std::is_signed_v<Value>)
            m = _mm_cvttpd_epi64(a.derived().m);
        else
            m = _mm_cvttpd_epu64(a.derived().m);
    }
#endif

    ENOKI_CONVERT(int64_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        int16_t ival;
        memcpy(&ival, a.derived().data(), 2);
        m = _mm_cvtepi8_epi64(
            _mm_cmpgt_epi8(_mm_cvtsi32_si128((int) ival), _mm_setzero_si128()));
    }

    ENOKI_REINTERPRET(float) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, float32[2] -> int64[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_castps_si128(_mm_setr_ps(v0, v0, v1, v1));
    }

    ENOKI_REINTERPRET(int32_t) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, int32[2] -> int64[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_setr_epi32(v0, v0, v1, v1);
    }

    ENOKI_REINTERPRET(uint32_t) {
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, uint32[2] -> int64[2])");
        auto v0 = a.derived().coeff(0), v1 = a.derived().coeff(1);
        m = _mm_setr_epi32((int32_t) v0, (int32_t) v0, (int32_t) v1,
                           (int32_t) v1);
    }

    ENOKI_REINTERPRET(double) : m(_mm_castpd_si128(a.derived().m)) { }
    ENOKI_REINTERPRET(int64_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2) {
        alignas(16) Value data[2];
        data[0] = (Value) a1.coeff(0);
        data[1] = (Value) a2.coeff(0);
        m = _mm_load_si128((__m128i *) data);
    }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(1)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm_add_epi64(m, a.m);   }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm_sub_epi64(m, a.m);   }
    ENOKI_INLINE Derived mul_(Ref a) const {
        #if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            return _mm_mullo_epi64(m, a.m);
        #else
            Derived result;
            ENOKI_TRACK_SCALAR("mul");
            for (size_t i = 0; i < Size; ++i)
                result.coeff(i) = coeff(i) * a.coeff(i);
            return result;
        #endif
    }
    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_epi64(m, a.k, _mm_set1_epi64x(-1));
            else
        #endif
        return _mm_or_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_maskz_mov_epi64(a.k, m);
            else
        #endif
        return _mm_and_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_xor_epi64(m, a.k, m, _mm_set1_epi64x(-1));
            else
        #endif
        return _mm_xor_si128(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm_mask_mov_epi64(m, a.k, _mm_setzero_si128());
            else
        #endif
        return _mm_andnot_si128(a.m, m);
    }

    template <size_t k> ENOKI_INLINE Derived sl_() const {
        return _mm_slli_epi64(m, (int) k);
    }

    template <size_t k> ENOKI_INLINE Derived sr_() const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_srai_epi64(m, (int) k);
            #else
                Derived out;
                ENOKI_TRACK_SCALAR("sr");
                for (size_t i = 0; i < Size; ++i)
                    out.coeff(i) = coeff(i) >> k;
                return out;
            #endif
        } else {
            return _mm_srli_epi64(m, (int) k);
        }
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm_sll_epi64(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_sra_epi64(m, _mm_set1_epi64x((long long) k));
            #else
                Derived out;
                ENOKI_TRACK_SCALAR("sr");
                for (size_t i = 0; i < Size; ++i)
                    out.coeff(i) = coeff(i) >> k;
                return out;
            #endif
        } else {
            return _mm_srl_epi64(m, _mm_set1_epi64x((long long) k));
        }
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        #if defined(ENOKI_X86_AVX2)
            return _mm_sllv_epi64(m, k.m);
        #else
            Derived out;
            ENOKI_TRACK_SCALAR("sl");
            for (size_t i = 0; i < Size; ++i)
                out.coeff(i) = coeff(i) << (unsigned int) k.coeff(i);
            return out;
        #endif
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_srav_epi64(m, k.m);
            #endif
        } else {
            #if defined(ENOKI_X86_AVX2)
                return _mm_srlv_epi64(m, k.m);
            #endif
        }
        Derived out;
        ENOKI_TRACK_SCALAR("sr");
        for (size_t i = 0; i < Size; ++i)
            out.coeff(i) = coeff(i) >> (unsigned int) k.coeff(i);
        return out;
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Imm> ENOKI_INLINE Derived rol_() const { return _mm_rol_epi64(m, (int) Imm); }
    template <size_t Imm> ENOKI_INLINE Derived ror_() const { return _mm_ror_epi64(m, (int) Imm); }
    ENOKI_INLINE Derived rol_(Ref k) const { return _mm_rolv_epi64(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm_rorv_epi64(m, k.m); }
#endif

    ENOKI_INLINE auto eq_(Ref a)  const {
        using Return = mask_t<Derived>;

        #if defined(ENOKI_X86_AVX512VL)
            return Return::from_k(_mm_cmpeq_epi64_mask(m, a.m));
        #else
            return Return(_mm_cmpeq_epi64(m, a.m));
        #endif
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(_mm_cmpneq_epi64_mask(m, a.m));
        #else
            return ~eq_(a);
        #endif
    }

    ENOKI_INLINE auto lt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm_cmpgt_epi64(a.m, m));
            } else {
                const __m128i offset =
                    _mm_set1_epi64x((long long) 0x8000000000000000ull);
                return Return(_mm_cmpgt_epi64(
                    _mm_sub_epi64(a.m, offset),
                    _mm_sub_epi64(m, offset)
                ));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm_cmplt_epi64_mask(m, a.m)
                                  : _mm_cmplt_epu64_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm_cmpgt_epi64(m, a.m));
            } else {
                const __m128i offset =
                    _mm_set1_epi64x((long long) 0x8000000000000000ull);
                return Return(_mm_cmpgt_epi64(
                    _mm_sub_epi64(m, offset),
                    _mm_sub_epi64(a.m, offset)
                ));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm_cmpgt_epi64_mask(m, a.m)
                                  : _mm_cmpgt_epu64_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto le_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm_cmple_epi64_mask(m, a.m)
                                           : _mm_cmple_epu64_mask(m, a.m));
        #else
            return ~gt_(a);
        #endif
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm_cmpge_epi64_mask(m, a.m)
                                           : _mm_cmpge_epu64_mask(m, a.m));
        #else
            return ~lt_(a);
        #endif
    }

    ENOKI_INLINE Derived min_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return std::is_signed_v<Value> ? _mm_min_epi64(a.m, m)
                                           : _mm_min_epu64(a.m, m);
        #else
            return select(derived() < a, derived(), a);
        #endif
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return std::is_signed_v<Value> ? _mm_max_epi64(a.m, m)
                                           : _mm_max_epu64(a.m, m);
        #else
            return select(derived() > a, derived(), a);
        #endif
    }

    ENOKI_INLINE Derived abs_() const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_abs_epi64(m);
            #else
                return select(derived() < zero<Derived>(),
                              ~derived() + Derived(Value(1)), derived());
            #endif
        } else {
            return m;
        }
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        #if !defined(ENOKI_X86_AVX512VL)
            return _mm_blendv_epi8(f.m, t.m, m.m);
        #else
            return _mm_mask_blend_epi64(m.k, f.m, t.m);
        #endif
    }

    ENOKI_INLINE Derived mulhi_(Ref a) const {
        ENOKI_TRACK_SCALAR("mulhi");
        return Derived(
            mulhi(coeff(0), a.coeff(0)),
            mulhi(coeff(1), a.coeff(1))
        );
    }

    template <int I0, int I1>
    ENOKI_INLINE Derived shuffle_() const {
        return _mm_shuffle_epi32(
            m, _MM_SHUFFLE(I1 * 2 + 1, I1 * 2, I0 * 2 + 1, I0 * 2));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        #if defined(ENOKI_X86_AVX)
            return _mm_castpd_si128(_mm_permutevar_pd(_mm_castsi128_pd(m), _mm_slli_epi64(index.m, 1)));
        #else
            return Base::shuffle_(index);
        #endif
    }

#if defined(ENOKI_X86_AVX512CD) && defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived lzcnt_() const { return _mm_lzcnt_epi64(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(64) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm_mask_mov_epi64(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm_mask_add_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm_mask_sub_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm_mask_mullo_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm_mask_or_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm_mask_and_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm_mask_xor_epi64(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op)                                     \
        ENOKI_INLINE Value name##_() const {                                  \
            Value t1 = Value(detail::mm_extract_epi64<1>(m));                 \
            Value t2 = Value(detail::mm_cvtsi128_si64(m));                    \
            return op;                                                        \
        }

    ENOKI_HORIZONTAL_OP(hsum,  t1 + t2)
    ENOKI_HORIZONTAL_OP(hprod, t1 * t2)
    ENOKI_HORIZONTAL_OP(hmin,  min(t1, t2))
    ENOKI_HORIZONTAL_OP(hmax,  max(t1, t2))

    #undef ENOKI_HORIZONTAL_OP

    ENOKI_INLINE bool all_()  const { return _mm_movemask_pd(_mm_castsi128_pd(m)) == 0x3;}
    ENOKI_INLINE bool any_()  const { return _mm_movemask_pd(_mm_castsi128_pd(m)) != 0x0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_pd(_mm_castsi128_pd(m)); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        _mm_store_si128((__m128i *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_store_epi64(ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX2)
            _mm_maskstore_epi64((long long *) ptr, mask.m, m);
        #else
            Base::store_(ptr, mask);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm_storeu_si128((__m128i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm_mask_storeu_epi64(ptr, mask.k, m);
        #elif defined(ENOKI_X86_AVX2)
            _mm_maskstore_epi64((long long *) ptr, mask.m, m);
        #else
            Base::store_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 16 == 0);
        return _mm_load_si128((const __m128i *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_load_epi64(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX2)
            return _mm_maskload_epi64((const long long *) ptr, mask.m);
        #else
            return Base::load_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm_loadu_si128((const __m128i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm_maskz_loadu_epi64(mask.k, ptr);
        #elif defined(ENOKI_X86_AVX2)
            return _mm_maskload_epi64((const long long *) ptr, mask.m);
        #else
            return Base::load_unaligned_(ptr, mask);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm_setzero_si128(); }

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            return Base::template gather_<Stride>(ptr, index, mask);
        } else {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm_mmask_i64gather_epi64(_mm_setzero_si128(), mask.k, index.m, (const long long *) ptr, Stride);
            #else
                return _mm_mask_i64gather_epi64(_mm_setzero_si128(), (const long long *) ptr, index.m, mask.m, Stride);
            #endif
        }
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            Base::template scatter_<Stride>(ptr, index, mask);
        else
            _mm_mask_i64scatter_epi64(ptr, mask.k, index.m, m, Stride);
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return (Value) detail::mm_cvtsi128_si64(_mm_mask_compress_epi64(_mm_setzero_si128(), mask.k, m));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        _mm_storeu_si128((__m128i *) ptr, _mm_mask_compress_epi64(_mm_setzero_si128(), mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }
#endif

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl for the n=3 case (single precision)
template <bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<float, 3, IsMask_, Derived_>
  : StaticArrayImpl<float, 4, IsMask_, Derived_> {
    using Base = StaticArrayImpl<float, 4, IsMask_, Derived_>;

    ENOKI_DECLARE_3D_ARRAY(StaticArrayImpl)

#if defined(ENOKI_X86_F16C)
    template <typename Derived2>
    ENOKI_INLINE StaticArrayImpl(
        const StaticArrayBase<half, 3, IsMask_, Derived2> &a) {
        uint16_t temp[4];
        memcpy(temp, a.derived().data(), sizeof(uint16_t) * 3);
        temp[3] = 0;
        m = _mm_cvtph_ps(_mm_loadl_epi64((const __m128i *) temp));
    }
#endif

    template <int I0, int I1, int I2>
    ENOKI_INLINE Derived shuffle_() const {
        return Base::template shuffle_<I0, I1, I2, 3>();
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op)                                     \
        ENOKI_INLINE Value name##_() const {                                  \
            __m128 t1 = _mm_movehl_ps(m, m);                                  \
            __m128 t2 = _mm_##op##_ss(m, t1);                                 \
            t1 = _mm_movehdup_ps(m);                                          \
            t1 = _mm_##op##_ss(t1, t2);                                       \
            return _mm_cvtss_f32(t1);                                         \
        }

    ENOKI_HORIZONTAL_OP(hsum, add)
    ENOKI_HORIZONTAL_OP(hprod, mul)
    ENOKI_HORIZONTAL_OP(hmin, min)
    ENOKI_HORIZONTAL_OP(hmax, max)

    #undef ENOKI_HORIZONTAL_OP

    ENOKI_INLINE Value dot_(Ref a) const {
        return _mm_cvtss_f32(_mm_dp_ps(m, a.m, 0b01110001));
    }

    ENOKI_INLINE bool all_()  const { return (_mm_movemask_ps(m) & 7) == 7; }
    ENOKI_INLINE bool any_()  const { return (_mm_movemask_ps(m) & 7) != 0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_ps(m) & 7; }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Loading/writing data (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    static ENOKI_INLINE auto mask_() {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k((__mmask8) 7);
        #else
            return mask_t<Derived>(_mm_castsi128_ps(_mm_setr_epi32(-1, -1, -1, 0)));
        #endif
    }

    using Base::load_;
    using Base::load_unaligned_;
    using Base::store_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        memcpy(ptr, &m, sizeof(Value) * 3);
    }
    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        store_(ptr);
    }
    static ENOKI_INLINE Derived load_(const void *ptr) {
        return Base::load_unaligned_(ptr);
    }
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        Derived result;
        memcpy(&result.m, ptr, sizeof(Value) * 3);
        return result;
    }

#if defined(ENOKI_X86_AVX)
    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        Base::store_(ptr, mask & mask_());
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        Base::store_unaligned_(ptr, mask & mask_());
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        return Base::load_(ptr, mask & mask_());
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return Base::load_unaligned_(ptr, mask & mask_());
    }
#endif

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        return Base::template gather_<Stride>(ptr, index, mask & mask_());
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        Base::template scatter_<Stride>(ptr, index, mask & mask_());
    }
#endif

    template <typename Mask>
    ENOKI_INLINE size_t compress_(float *&ptr, const Mask &mask) const {
        return Base::compress_(ptr, mask & mask_());
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl for the n=3 case (32 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(16)
    StaticArrayImpl<Value_, 3, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayImpl<Value_, 4, IsMask_, Derived_> {
    using Base = StaticArrayImpl<Value_, 4, IsMask_, Derived_>;

    ENOKI_DECLARE_3D_ARRAY(StaticArrayImpl)

    template <int I0, int I1, int I2>
    ENOKI_INLINE Derived shuffle_() const {
        return Base::template shuffle_<I0, I1, I2, 3>();
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    #define ENOKI_HORIZONTAL_OP(name, op)                                     \
        ENOKI_INLINE Value name##_() const {                                  \
            __m128i t1 = _mm_unpackhi_epi32(m, m);                            \
            __m128i t2 = _mm_##op##_epi32(m, t1);                             \
            t1 = _mm_shuffle_epi32(m, 1);                                     \
            t1 = _mm_##op##_epi32(t1, t2);                                    \
            return (Value) _mm_cvtsi128_si32(t1);                             \
        }

    #define ENOKI_HORIZONTAL_OP_SIGNED(name, op)                              \
        ENOKI_INLINE Value name##_() const {                                  \
            __m128i t2, t1 = _mm_unpackhi_epi32(m, m);                        \
            if constexpr (std::is_signed<Value>::value)                       \
                t2 = _mm_##op##_epi32(m, t1);                                 \
            else                                                              \
                t2 = _mm_##op##_epu32(m, t1);                                 \
            t1 = _mm_shuffle_epi32(m, 1);                                     \
            if constexpr (std::is_signed<Value>::value)                       \
                t1 = _mm_##op##_epi32(t1, t2);                                \
            else                                                              \
                t1 = _mm_##op##_epu32(t1, t2);                                \
            return (Value) _mm_cvtsi128_si32(t1);                             \
        }

    ENOKI_HORIZONTAL_OP(hsum, add)
    ENOKI_HORIZONTAL_OP(hprod, mullo)
    ENOKI_HORIZONTAL_OP_SIGNED(hmin, min)
    ENOKI_HORIZONTAL_OP_SIGNED(hmax, max)

    #undef ENOKI_HORIZONTAL_OP
    #undef ENOKI_HORIZONTAL_OP_SIGNED

    ENOKI_INLINE bool all_()  const { return (_mm_movemask_ps(_mm_castsi128_ps(m)) & 7) == 7;}
    ENOKI_INLINE bool any_()  const { return (_mm_movemask_ps(_mm_castsi128_ps(m)) & 7) != 0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm_movemask_ps(_mm_castsi128_ps(m)) & 7; }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Loading/writing data (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    static ENOKI_INLINE auto mask_() {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k((__mmask8) 7);
        #else
            return mask_t<Derived>(_mm_setr_epi32(-1, -1, -1, 0));
        #endif
    }

    using Base::load_;
    using Base::load_unaligned_;
    using Base::store_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        memcpy(ptr, &m, sizeof(Value) * 3);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        store_(ptr);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return Base::load_unaligned_(ptr);
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        Derived result;
        memcpy(&result.m, ptr, sizeof(Value) * 3);
        return result;
    }

#if defined(ENOKI_X86_AVX2)
    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        return Base::store_unaligned_(ptr, mask & mask_());
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        return Base::store_(ptr, mask & mask_());
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        return Base::load_(ptr, mask & mask_());
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return Base::load_unaligned_(ptr, mask & mask_());
    }
#endif

#if defined(ENOKI_X86_AVX2)
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        return Base::template gather_<Stride>(ptr, index, mask & mask_());
    }
#endif

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        Base::template scatter_<Stride>(ptr, index, mask & mask_());
    }
#endif

    template <typename T, typename Mask>
    ENOKI_INLINE size_t compress_(T *&ptr, const Mask &mask) const {
        return Base::compress_(ptr, mask & mask_());
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

#if defined(ENOKI_X86_AVX512VL)
template <typename Derived_>
ENOKI_DECLARE_KMASK(float, 4, Derived_, int)
template <typename Derived_>
ENOKI_DECLARE_KMASK(float, 3, Derived_, int)
template <typename Derived_>
ENOKI_DECLARE_KMASK(double, 2, Derived_, int)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 4, Derived_, enable_if_int32_t<Value_>)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 3, Derived_, enable_if_int32_t<Value_>)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 2, Derived_, enable_if_int64_t<Value_>)
#endif

NAMESPACE_END(enoki)
