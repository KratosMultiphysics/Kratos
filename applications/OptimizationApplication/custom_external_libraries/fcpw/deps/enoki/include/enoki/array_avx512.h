/*
    enoki/array_avx512.h -- Packed SIMD array (AVX512 specialization)

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
template <> struct is_native<float, 16> : std::true_type { } ;
template <> struct is_native<double, 8> : std::true_type { };
template <typename Value>    struct is_native<Value, 16, enable_if_int32_t<Value>> : std::true_type { };
template <typename Value>    struct is_native<Value, 8, enable_if_int64_t<Value>> : std::true_type { };
NAMESPACE_END(detail)

/// Partial overload of StaticArrayImpl using AVX512 intrinsics (single precision)
template <bool IsMask_, typename Derived_> struct alignas(64)
    StaticArrayImpl<float, 16, IsMask_, Derived_>
  : StaticArrayBase<float, 16, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(float, 16, __m512)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm512_set1_ps(value)) { }

    ENOKI_INLINE StaticArrayImpl(Value f0,  Value f1,  Value f2,  Value f3,
                                 Value f4,  Value f5,  Value f6,  Value f7,
                                 Value f8,  Value f9,  Value f10, Value f11,
                                 Value f12, Value f13, Value f14, Value f15)
        : m(_mm512_setr_ps(f0, f1, f2, f3, f4, f5, f6, f7, f8,
                           f9, f10, f11, f12, f13, f14, f15)) { }


    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(half)
        : m(_mm512_cvtph_ps(_mm256_loadu_si256((const __m256i *) a.derived().data()))) { }

    ENOKI_CONVERT(float) : m(a.derived().m) { }

    ENOKI_CONVERT(int32_t) : m(_mm512_cvtepi32_ps(a.derived().m)) { }

    ENOKI_CONVERT(uint32_t) : m(_mm512_cvtepu32_ps(a.derived().m)) { }

    ENOKI_CONVERT(double)
        : m(detail::concat(_mm512_cvtpd_ps(low(a).m),
                           _mm512_cvtpd_ps(high(a).m))) { }

#if defined(ENOKI_X86_AVX512DQ)
    ENOKI_CONVERT(int64_t)
        : m(detail::concat(_mm512_cvtepi64_ps(low(a).m),
                           _mm512_cvtepi64_ps(high(a).m))) { }

    ENOKI_CONVERT(uint64_t)
        : m(detail::concat(_mm512_cvtepu64_ps(low(a).m),
                           _mm512_cvtepu64_ps(high(a).m))) { }
#elif defined(ENOKI_X86_AVX512CD)
    /* Emulate uint64_t -> float conversion instead of falling
       back to scalar operations. This is quite a bit faster
       (>6x for unsigned, >5x for signed). */

    ENOKI_CONVERT(uint64_t) {
        using Int64 = int64_array_t<Derived2>;
        using Int32 = uint32_array_t<Derived2>;

        auto lz = lzcnt(a);
        auto shift = (63 - 23) - Int64(lz);
        auto abs_shift = abs(shift);
        auto nzero_mask = neq(a, 0ull);
        auto mant = select(shift > 0, a >> abs_shift, a << abs_shift);
        auto exp = sl<23>(uint64_t(127 + 63) - lz) & nzero_mask;
        auto comb = exp | (mant & 0x7fffffull);

        m = reinterpret_array<Derived>(Int32(comb)).m;
    }

    ENOKI_CONVERT(int64_t) {
        using Int32 = uint32_array_t<Derived2>;

        auto b = abs(a), lz = lzcnt(b);
        auto shift = (63 - 23) - lz;
        auto abs_shift = abs(shift);
        auto nzero_mask = neq(a, 0ll);
        auto mant = select(shift > 0, b >> abs_shift, b << abs_shift);
        auto sign = sr<32>(a) & 0x80000000ll;
        auto exp = sl<23>(int64_t(127 + 63) - lz) & nzero_mask;
        auto comb = exp | (mant & 0x7fffffll) | sign;

        m = reinterpret_array<Derived>(Int32(comb)).m;
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(float) : m(a.derived().m) { }

    ENOKI_REINTERPRET(int32_t) : m(_mm512_castsi512_ps(a.derived().m)) { }
    ENOKI_REINTERPRET(uint32_t) : m(_mm512_castsi512_ps(a.derived().m)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm512_castps512_ps256(m); }
    ENOKI_INLINE Array2 high_() const {
        #if defined(ENOKI_X86_AVX512DQ)
            return _mm512_extractf32x8_ps(m, 1);
        #else
            return _mm256_castpd_ps(_mm512_extractf64x4_pd(_mm512_castps_pd(m), 1));
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm512_add_ps(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm512_sub_ps(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm512_mul_ps(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const { return _mm512_div_ps(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_mask_mov_ps(m, a.k, _mm512_set1_ps(memcpy_cast<Value>(int32_t(-1))));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_or_ps(m, a.m);
            #else
                return _mm512_castsi512_ps(
                    _mm512_or_si512(_mm512_castps_si512(m), _mm512_castps_si512(a.m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_maskz_mov_ps(a.k, m);
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_and_ps(m, a.m);
            #else
                return _mm512_castsi512_ps(
                    _mm512_and_si512(_mm512_castps_si512(m), _mm512_castps_si512(a.m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_mask_mov_ps(m, a.k, _mm512_setzero_ps());
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_andnot_ps(a.m, m);
            #else
                return _mm512_castsi512_ps(
                    _mm512_andnot_si512(_mm512_castps_si512(a.m), _mm512_castps_si512(m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            const __m512 c = _mm512_set1_ps(memcpy_cast<Value>(int32_t(-1)));
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_mask_xor_ps(m, a.k, m, c);
            #else
                const __m512i v0 = _mm512_castps_si512(m);
                return _mm512_castsi512_ps(_mm512_mask_xor_epi32(v0, a.k, v0, c));
            #endif
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_xor_ps(m, a.m);
            #else
                return _mm512_castsi512_ps(
                    _mm512_xor_si512(_mm512_castps_si512(m), _mm512_castps_si512(a.m)));
            #endif
        }
    }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_LT_OQ));  }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_GT_OQ));  }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_LE_OQ));  }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_GE_OQ));  }
    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_EQ_OQ));  }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_ps_mask(m, a.m, _CMP_NEQ_UQ)); }

    ENOKI_INLINE Derived abs_() const { return andnot_(Derived(_mm512_set1_ps(-0.f))); }

    ENOKI_INLINE Derived min_(Ref b) const { return _mm512_min_ps(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return _mm512_max_ps(b.m, m); }
    ENOKI_INLINE Derived ceil_()     const { return _mm512_ceil_ps(m);     }
    ENOKI_INLINE Derived floor_()    const { return _mm512_floor_ps(m);    }
    ENOKI_INLINE Derived sqrt_()     const { return _mm512_sqrt_ps(m); }

    ENOKI_INLINE Derived round_() const {
        return _mm512_roundscale_ps(m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    ENOKI_INLINE Derived trunc_() const {
        return _mm512_roundscale_ps(m, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    }

    template <typename T>
    ENOKI_INLINE auto ceil2int_() const {
        if constexpr (sizeof(scalar_t<T>) == 4) {
            if constexpr (std::is_signed_v<scalar_t<T>>)
                return T(_mm512_cvt_roundps_epi32(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
            else
                return T(_mm512_cvt_roundps_epu32(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                using A = typename T::Array1;
                if constexpr (std::is_signed_v<scalar_t<T>>)
                    return T(
                        A(_mm512_cvt_roundps_epi64(low(derived()).m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC)),
                        A(_mm512_cvt_roundps_epi64(high(derived()).m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC))
                    );
                else
                    return T(
                        A(_mm512_cvt_roundps_epu64(low(derived()).m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC)),
                        A(_mm512_cvt_roundps_epu64(high(derived()).m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC))
                    );
            #else
                return Base::template ceil2int_<T>();
            #endif
        }
    }

    template <typename T>
    ENOKI_INLINE auto floor2int_() const {
        if constexpr (sizeof(scalar_t<T>) == 4) {
            if constexpr (std::is_signed_v<scalar_t<T>>)
                return T(_mm512_cvt_roundps_epi32(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
            else
                return T(_mm512_cvt_roundps_epu32(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                using A = typename T::Array1;

                if constexpr (std::is_signed_v<scalar_t<T>>)
                    return T(
                        A(_mm512_cvt_roundps_epi64(low(derived()).m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)),
                        A(_mm512_cvt_roundps_epi64(high(derived()).m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC))
                    );
                else
                    return T(
                        A(_mm512_cvt_roundps_epu64(low(derived()).m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC)),
                        A(_mm512_cvt_roundps_epu64(high(derived()).m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC))
                    );
            #else
                return Base::template floor2int_<T>();
            #endif
        }
    }

    ENOKI_INLINE Derived fmadd_   (Ref b, Ref c) const { return _mm512_fmadd_ps   (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsub_   (Ref b, Ref c) const { return _mm512_fmsub_ps   (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmadd_  (Ref b, Ref c) const { return _mm512_fnmadd_ps  (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmsub_  (Ref b, Ref c) const { return _mm512_fnmsub_ps  (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsubadd_(Ref b, Ref c) const { return _mm512_fmsubadd_ps(m, b.m, c.m); }
    ENOKI_INLINE Derived fmaddsub_(Ref b, Ref c) const { return _mm512_fmaddsub_ps(m, b.m, c.m); }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        return _mm512_mask_blend_ps(m.k, f.m, t.m);
    }

    template <size_t I0,  size_t I1,  size_t I2,  size_t I3,  size_t I4,
              size_t I5,  size_t I6,  size_t I7,  size_t I8,  size_t I9,
              size_t I10, size_t I11, size_t I12, size_t I13, size_t I14,
              size_t I15>
    ENOKI_INLINE Derived shuffle_() const {
        const __m512i idx =
            _mm512_setr_epi32(I0, I1, I2, I3, I4, I5, I6, I7, I8,
                              I9, I10, I11, I12, I13, I14, I15);
        return _mm512_permutexvar_ps(idx, m);
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return _mm512_permutexvar_ps(index.m, m);
    }

    ENOKI_INLINE Derived rcp_() const {
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28, use as is */
            return _mm512_rcp28_ps(m);
        #else
            /* Use best reciprocal approximation available on the current
               hardware and refine */
            __m512 r = _mm512_rcp14_ps(m); /* rel error < 2^-14 */

            /* Refine using one Newton-Raphson iteration */
            __m512 t0 = _mm512_add_ps(r, r),
                   t1 = _mm512_mul_ps(r, m);

            r = _mm512_fnmadd_ps(t1, r, t0);

            return _mm512_fixupimm_ps(r, m,
                _mm512_set1_epi32(0x0087A622), 0);
        #endif
    }

    ENOKI_INLINE Derived rsqrt_() const {
        #if defined(ENOKI_X86_AVX512ER)
            /* rel err < 2^28, use as is */
            return _mm512_rsqrt28_ps(m);
        #else
            __m512 r = _mm512_rsqrt14_ps(m); /* rel error < 2^-14 */

            /* Refine using one Newton-Raphson iteration */
            const __m512 c0 = _mm512_set1_ps(0.5f),
                         c1 = _mm512_set1_ps(3.0f);

            __m512 t0 = _mm512_mul_ps(r, c0),
                   t1 = _mm512_mul_ps(r, m);

            r = _mm512_mul_ps(_mm512_fnmadd_ps(t1, r, c1), t0);

            return _mm512_fixupimm_ps(r, m,
                _mm512_set1_epi32(0x0383A622), 0);
        #endif
    }

    ENOKI_INLINE Derived ldexp_(Ref arg) const { return _mm512_scalef_ps(m, arg.m); }

    ENOKI_INLINE std::pair<Derived, Derived> frexp_() const {
        return std::make_pair<Derived, Derived>(
            _mm512_getmant_ps(m, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src),
            _mm512_getexp_ps(m));
    }

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    //! @}
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_store_ps((Value *) ENOKI_ASSUME_ALIGNED(ptr, 64), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_mask_store_ps((Value *) ENOKI_ASSUME_ALIGNED(ptr, 64), mask.k, m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm512_storeu_ps((Value *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        _mm512_mask_storeu_ps((Value *) ptr, mask.k, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_load_ps((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_maskz_load_ps(mask.k, (const Value *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm512_loadu_ps((const Value *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return _mm512_maskz_loadu_ps(mask.k, (const Value *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return _mm512_setzero_ps(); }

#if defined(ENOKI_X86_AVX512PF)
    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *ptr, const Index &index, const Mask &mask) {
        constexpr auto Hint = Level == 1 ? _MM_HINT_T0 : _MM_HINT_T1;
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            if constexpr (Write)
                _mm512_mask_prefetch_i32scatter_ps((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i32gather_ps(index.m, mask.k, ptr, Stride, Hint);
        } else {
            if constexpr (Write) {
                _mm512_mask_prefetch_i64scatter_ps((void *) ptr, low(mask).k, low(index).m, Stride, Hint);
                _mm512_mask_prefetch_i64scatter_ps((void *) ptr, high(mask).k, high(index).m, Stride, Hint);
            } else {
                _mm512_mask_prefetch_i64gather_ps(low(index).m, low(mask).k, ptr, Stride, Hint);
                _mm512_mask_prefetch_i64gather_ps(high(index).m, high(mask).k, ptr, Stride, Hint);
            }
        }
    }
#endif

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            return _mm512_mask_i32gather_ps(_mm512_setzero_ps(), mask.k, index.m, (const float *) ptr, Stride);
        } else {
            return detail::concat(
                _mm512_mask_i64gather_ps(_mm256_setzero_ps(),  low(mask).k,  low(index).m, (const float *) ptr, Stride),
                _mm512_mask_i64gather_ps(_mm256_setzero_ps(), high(mask).k, high(index).m, (const float *) ptr, Stride));
        }
    }

    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            _mm512_mask_i32scatter_ps(ptr, mask.k, index.m, m, Stride);
        } else {
            _mm512_mask_i64scatter_ps(ptr, low(mask).k,   low(index).m,  low(derived()).m,  Stride);
            _mm512_mask_i64scatter_ps(ptr, high(mask).k, high(index).m, high(derived()).m, Stride);
        }
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return _mm_cvtss_f32(_mm512_castps512_ps128(_mm512_maskz_compress_ps(mask.k, m)));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(float *&ptr, const Mask &mask) const {
        _mm512_storeu_ps(ptr, _mm512_maskz_compress_ps(mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }

#if defined(ENOKI_X86_AVX512CD)
    template <size_t Stride, typename Index, typename Mask, typename Func, typename... Args>
    static ENOKI_INLINE void transform_(void *mem,
                                        Index index,
                                        const Mask &mask,
                                        const Func &func,
                                        const Args &... args) {
        Derived values = _mm512_mask_i32gather_ps(
            _mm512_undefined_ps(), mask.k, index.m, mem, (int) Stride);

        index.m = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), mask.k, index.m);

        __m512i conflicts = _mm512_conflict_epi32(index.m);
        __m512i perm_idx  = _mm512_sub_epi32(_mm512_set1_epi32(31), _mm512_lzcnt_epi32(conflicts));
        __mmask16 todo    = _mm512_mask_test_epi32_mask(mask.k, conflicts, _mm512_set1_epi32(-1));

        func(values, args...);

        ENOKI_NOUNROLL while (ENOKI_UNLIKELY(!_mm512_kortestz(todo, todo))) {
            __mmask16 cur = _mm512_mask_testn_epi32_mask(
                todo, conflicts, _mm512_broadcastmw_epi32(todo));
            values.m = _mm512_mask_permutexvar_ps(values.m, cur, perm_idx, values.m);

            __m512 backup(values.m);
            func(values, args...);

            values.m = _mm512_mask_mov_ps(backup, cur, values.m);
            todo = _mm512_kxor(todo, cur);
        }

        _mm512_mask_i32scatter_ps(mem, mask.k, index.m, values.m, (int) Stride);
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm512_mask_mov_ps(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm512_mask_add_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm512_mask_sub_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm512_mask_mul_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mdiv_   (const Derived &a, const Mask &mask) { m = _mm512_mask_div_ps(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_or_ps(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_ps(
                _mm512_or_si512(_mm512_castps_si512(m), mask.k,
                                _mm512_castps_si512(m), _mm512_castps_si512(a.m)));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_and_ps(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_ps(_mm512_and_si512(_mm512_castps_si512(m), mask.k,
                                                     _mm512_castps_si512(m),
                                                     _mm512_castps_si512(a.m)));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_xor_ps(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_ps(_mm512_xor_si512(_mm512_castps_si512(m), mask.k,
                                                     _mm512_castps_si512(m),
                                                     _mm512_castps_si512(a.m)));
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using AVX512 intrinsics (double precision)
template <bool IsMask_, typename Derived_> struct alignas(64)
    StaticArrayImpl<double, 8, IsMask_, Derived_>
  : StaticArrayBase<double, 8, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(double, 8, __m512d)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(const Value &value) : m(_mm512_set1_pd(value)) { }
    ENOKI_INLINE StaticArrayImpl(Value f0, Value f1, Value f2, Value f3,
                                 Value f4, Value f5, Value f6, Value f7)
        : m(_mm512_setr_pd(f0, f1, f2, f3, f4, f5, f6, f7)) { }


    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(half)
        : m(_mm512_cvtps_pd(
              _mm256_cvtph_ps(_mm_loadu_si128((const __m128i *) a.derived().data())))) { }

    ENOKI_CONVERT(float) : m(_mm512_cvtps_pd(a.derived().m)) { }

    ENOKI_CONVERT(double) : m(a.derived().m) { }

    ENOKI_CONVERT(int32_t) : m(_mm512_cvtepi32_pd(a.derived().m)) { }

    ENOKI_CONVERT(uint32_t) : m(_mm512_cvtepu32_pd(a.derived().m)) { }

#if defined(ENOKI_X86_AVX512DQ)
    ENOKI_CONVERT(int64_t)
        : m(_mm512_cvtepi64_pd(a.derived().m)) { }

    ENOKI_CONVERT(uint64_t)
        : m(_mm512_cvtepu64_pd(a.derived().m)) { }
#elif defined(ENOKI_X86_AVX512CD)
    /* Emulate uint64_t -> double conversion instead of falling
       back to scalar operations. This is quite a bit faster
       (>5.5x for unsigned, > for signed). */

    ENOKI_CONVERT(uint64_t) {
        using Int64 = int64_array_t<Derived2>;

        auto lz = lzcnt(a);
        auto shift = (63 - 52) - Int64(lz);
        auto abs_shift = abs(shift);
        auto nzero_mask = neq(a, 0ull);
        auto mant = select(shift > 0, a >> abs_shift, a << abs_shift);
        auto exp = sl<52>(uint64_t(1023 + 63) - lz) & nzero_mask;
        auto comb = exp | (mant & 0xfffffffffffffull);

        m = reinterpret_array<Derived>(comb).m;
    }

    ENOKI_CONVERT(int64_t) {
        auto b = abs(a), lz = lzcnt(b);
        auto shift = (63 - 52) - lz;
        auto abs_shift = abs(shift);
        auto nzero_mask = neq(a, 0ll);
        auto mant = select(shift > 0, b >> abs_shift, b << abs_shift);
        auto sign = a & 0x8000000000000000ull;
        auto exp = sl<52>(int64_t(1023 + 63) - lz) & nzero_mask;
        auto comb = exp | (mant & 0xfffffffffffffull) | sign;

        m = reinterpret_array<Derived>(comb).m;
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(double) : m(a.derived().m) { }

    ENOKI_REINTERPRET(int64_t) : m(_mm512_castsi512_pd(a.derived().m)) { }
    ENOKI_REINTERPRET(uint64_t) : m(_mm512_castsi512_pd(a.derived().m)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm512_castpd512_pd256(m); }
    ENOKI_INLINE Array2 high_() const { return _mm512_extractf64x4_pd(m, 1); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm512_add_pd(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm512_sub_pd(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm512_mul_pd(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const { return _mm512_div_pd(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_mask_mov_pd(m, a.k, _mm512_set1_pd(memcpy_cast<Value>(int64_t(-1))));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_or_pd(m, a.m);
            #else
                return _mm512_castsi512_pd(
                    _mm512_or_si512(_mm512_castpd_si512(m), _mm512_castpd_si512(a.m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_maskz_mov_pd(a.k, m);
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_and_pd(m, a.m);
            #else
                return _mm512_castsi512_pd(
                    _mm512_and_si512(_mm512_castpd_si512(m), _mm512_castpd_si512(a.m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            return _mm512_mask_mov_pd(m, a.k, _mm512_setzero_pd());
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_andnot_pd(a.m, m);
            #else
                return _mm512_castsi512_pd(
                    _mm512_andnot_si512(_mm512_castpd_si512(a.m), _mm512_castpd_si512(m)));
            #endif
        }
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        if constexpr (is_mask_v<T>) {
            const __m512 c = _mm512_set1_pd(memcpy_cast<Value>(int64_t(-1)));
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_mask_xor_pd(m, a.k, m, c);
            #else
                const __m512i v0 = _mm512_castpd_si512(m);
                return _mm512_castsi512_pd(_mm512_mask_xor_epi64(v0, a.k, v0, c));
            #endif
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                return _mm512_xor_pd(m, a.m);
            #else
                return _mm512_castsi512_pd(
                    _mm512_xor_si512(_mm512_castpd_si512(m), _mm512_castpd_si512(a.m)));
            #endif
        }
    }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_LT_OQ));  }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_GT_OQ));  }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_LE_OQ));  }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_GE_OQ));  }
    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_EQ_OQ));  }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_pd_mask(m, a.m, _CMP_NEQ_UQ)); }

    ENOKI_INLINE Derived abs_() const { return andnot_(Derived(_mm512_set1_pd(-0.0))); }

    ENOKI_INLINE Derived min_(Ref b) const { return _mm512_min_pd(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return _mm512_max_pd(b.m, m); }
    ENOKI_INLINE Derived ceil_()     const { return _mm512_ceil_pd(m);     }
    ENOKI_INLINE Derived floor_()    const { return _mm512_floor_pd(m);    }
    ENOKI_INLINE Derived sqrt_()     const { return _mm512_sqrt_pd(m); }

    template <typename T>
    ENOKI_INLINE auto ceil2int_() const {
        if constexpr (sizeof(scalar_t<T>) == 4) {
            if constexpr (std::is_signed_v<scalar_t<T>>)
                return T(_mm512_cvt_roundpd_epi32(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
            else
                return T(_mm512_cvt_roundpd_epu32(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                if constexpr (std::is_signed_v<scalar_t<T>>)
                    return T(_mm512_cvt_roundpd_epi64(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
                else
                    return T(_mm512_cvt_roundpd_epu64(m, _MM_FROUND_TO_POS_INF | _MM_FROUND_NO_EXC));
            #else
                return Base::template ceil2int_<T>();
            #endif
        }
    }

    template <typename T>
    ENOKI_INLINE auto floor2int_() const {
        if constexpr (sizeof(scalar_t<T>) == 4) {
            if constexpr (std::is_signed_v<scalar_t<T>>)
                return T(_mm512_cvt_roundpd_epi32(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
            else
                return T(_mm512_cvt_roundpd_epu32(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
        } else {
            #if defined(ENOKI_X86_AVX512DQ)
                if constexpr (std::is_signed_v<scalar_t<T>>)
                    return T(_mm512_cvt_roundpd_epi64(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
                else
                    return T(_mm512_cvt_roundpd_epu64(m, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));
            #else
                return Base::template floor2int_<T>();
            #endif
        }
    }

    ENOKI_INLINE Derived round_() const {
        return _mm512_roundscale_pd(m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    }

    ENOKI_INLINE Derived trunc_() const {
        return _mm512_roundscale_pd(m, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    }

    ENOKI_INLINE Derived fmadd_   (Ref b, Ref c) const { return _mm512_fmadd_pd   (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsub_   (Ref b, Ref c) const { return _mm512_fmsub_pd   (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmadd_  (Ref b, Ref c) const { return _mm512_fnmadd_pd  (m, b.m, c.m); }
    ENOKI_INLINE Derived fnmsub_  (Ref b, Ref c) const { return _mm512_fnmsub_pd  (m, b.m, c.m); }
    ENOKI_INLINE Derived fmsubadd_(Ref b, Ref c) const { return _mm512_fmsubadd_pd(m, b.m, c.m); }
    ENOKI_INLINE Derived fmaddsub_(Ref b, Ref c) const { return _mm512_fmaddsub_pd(m, b.m, c.m); }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        return _mm512_mask_blend_pd(m.k, f.m, t.m);
    }

    template <size_t I0, size_t I1, size_t I2, size_t I3, size_t I4, size_t I5,
              size_t I6, size_t I7>
    ENOKI_INLINE Derived shuffle_() const {
        const __m512i idx =
            _mm512_setr_epi64(I0, I1, I2, I3, I4, I5, I6, I7);
        return _mm512_permutexvar_pd(idx, m);
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return _mm512_permutexvar_pd(index.m, m);
    }

    ENOKI_INLINE Derived rcp_() const {
        /* Use best reciprocal approximation available on the current
           hardware and refine */
        __m512d r;

        #if defined(ENOKI_X86_AVX512ER)
            r = _mm512_rcp28_pd(m); /* rel err < 2^28 */
        #else
            r = _mm512_rcp14_pd(m); /* rel error < 2^-14 */
        #endif

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            __m512d t0 = _mm512_add_pd(r, r);
            __m512d t1 = _mm512_mul_pd(r, m);

            r = _mm512_fnmadd_pd(t1, r, t0);
        }

        return _mm512_fixupimm_pd(r, m,
            _mm512_set1_epi32(0x0087A622), 0);
    }

    ENOKI_INLINE Derived rsqrt_() const {
        /* Use best reciprocal square root approximation available
           on the current hardware and refine */
        __m512d r;
        #if defined(ENOKI_X86_AVX512ER)
            r = _mm512_rsqrt28_pd(m); /* rel err < 2^28 */
        #else
            r = _mm512_rsqrt14_pd(m); /* rel error < 2^-14 */
        #endif

        const __m512d c0 = _mm512_set1_pd(0.5),
                      c1 = _mm512_set1_pd(3.0);

        /* Refine using 1-2 Newton-Raphson iterations */
        ENOKI_UNROLL for (int i = 0; i < (has_avx512er ? 1 : 2); ++i) {
            __m512d t0 = _mm512_mul_pd(r, c0);
            __m512d t1 = _mm512_mul_pd(r, m);

            r = _mm512_mul_pd(_mm512_fnmadd_pd(t1, r, c1), t0);
        }

        return _mm512_fixupimm_pd(r, m,
            _mm512_set1_epi32(0x0383A622), 0);
    }


    ENOKI_INLINE Derived ldexp_(Ref arg) const { return _mm512_scalef_pd(m, arg.m); }

    ENOKI_INLINE std::pair<Derived, Derived> frexp_() const {
        return std::make_pair<Derived, Derived>(
            _mm512_getmant_pd(m, _MM_MANT_NORM_p5_1, _MM_MANT_SIGN_src),
            _mm512_getexp_pd(m));
    }

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    //! @}
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_store_pd((Value *) ENOKI_ASSUME_ALIGNED(ptr, 64), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_mask_store_pd((Value *) ENOKI_ASSUME_ALIGNED(ptr, 64), mask.k, m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm512_storeu_pd((Value *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        _mm512_mask_storeu_pd((Value *) ptr, mask.k, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_load_pd((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_maskz_load_pd(mask.k, (const Value *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm512_loadu_pd((const Value *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return _mm512_maskz_loadu_pd(mask.k, (const Value *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return _mm512_setzero_pd(); }

#if defined(ENOKI_X86_AVX512PF)
    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *ptr, const Index &index, const Mask &mask) {
        constexpr auto Hint = Level == 1 ? _MM_HINT_T0 : _MM_HINT_T1;
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            if (Write)
                _mm512_mask_prefetch_i32scatter_pd((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i32gather_pd(index.m, mask.k, ptr, Stride, Hint);
        } else {
            if (Write)
                _mm512_mask_prefetch_i64scatter_pd((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i64gather_pd(index.m, mask.k, ptr, Stride, Hint);
        }
    }
#endif

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            return _mm512_mask_i32gather_pd(_mm512_setzero_pd(), mask.k, index.m, (const double *) ptr, Stride);
        else
            return _mm512_mask_i64gather_pd(_mm512_setzero_pd(), mask.k, index.m, (const double *) ptr, Stride);
    }

    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm512_mask_i32scatter_pd(ptr, mask.k, index.m, m, Stride);
        else
            _mm512_mask_i64scatter_pd(ptr, mask.k, index.m, m, Stride);
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return _mm_cvtsd_f64(_mm512_castpd512_pd128(_mm512_maskz_compress_pd(mask.k, m)));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(double *&ptr, const Mask &mask) const {
        _mm512_storeu_pd(ptr, _mm512_maskz_compress_pd(mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }

#if defined(ENOKI_X86_AVX512CD)
    template <size_t Stride, typename Index, typename Mask, typename Func, typename... Args>
    static ENOKI_INLINE void transform_(void *mem,
                                        Index index,
                                        const Mask &mask,
                                        const Func &func,
                                        const Args &... args) {
        Derived values = _mm512_mask_i64gather_pd(
            _mm512_undefined_pd(), mask.k, index.m, mem, (int) Stride);

        index.m = _mm512_mask_mov_epi64(_mm512_set1_epi64(-1), mask.k, index.m);

        __m512i conflicts = _mm512_conflict_epi64(index.m);
        __m512i perm_idx  = _mm512_sub_epi64(_mm512_set1_epi64(63), _mm512_lzcnt_epi64(conflicts));
        __mmask8 todo     = _mm512_mask_test_epi64_mask(mask.k, conflicts, _mm512_set1_epi64(-1));

        func(values, args...);

        ENOKI_NOUNROLL while (ENOKI_UNLIKELY(todo)) {
            __mmask8 cur = _mm512_mask_testn_epi64_mask(
                todo, conflicts, _mm512_broadcastmb_epi64(todo));
            values.m = _mm512_mask_permutexvar_pd(values.m, cur, perm_idx, values.m);

            __m512d backup(values.m);
            func(values, args...);

            values.m = _mm512_mask_mov_pd(backup, cur, values.m);
            todo ^= cur;
        }

        _mm512_mask_i64scatter_pd(mem, mask.k, index.m, values.m, (int) Stride);
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm512_mask_mov_pd(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm512_mask_add_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm512_mask_sub_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm512_mask_mul_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mdiv_   (const Derived &a, const Mask &mask) { m = _mm512_mask_div_pd(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_or_pd(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_pd(_mm512_or_si512(_mm512_castpd_si512(m), mask.k,
                                                    _mm512_castpd_si512(m),
                                                    _mm512_castpd_si512(a.m)));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_and_pd(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_pd(_mm512_and_si512(_mm512_castpd_si512(m), mask.k,
                                                     _mm512_castpd_si512(m),
                                                     _mm512_castpd_si512(a.m)));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ)
            m = _mm512_mask_xor_pd(m, mask.k, m, a.m);
        #else
            m = _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(m), mask.k,
                                                     _mm512_castpd_si512(m),
                                                     _mm512_castpd_si512(a.m)));
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using AVX512 intrinsics (32 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(64)
    StaticArrayImpl<Value_, 16, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayBase<Value_, 16, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 16, __m512i)


    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm512_set1_epi32((int32_t) value)) { }

    ENOKI_INLINE StaticArrayImpl(Value f0,  Value f1,  Value f2,  Value f3,
                                 Value f4,  Value f5,  Value f6,  Value f7,
                                 Value f8,  Value f9,  Value f10, Value f11,
                                 Value f12, Value f13, Value f14, Value f15)
        : m(_mm512_setr_epi32(
              (int32_t) f0,  (int32_t) f1,  (int32_t) f2,  (int32_t) f3,
              (int32_t) f4,  (int32_t) f5,  (int32_t) f6,  (int32_t) f7,
              (int32_t) f8,  (int32_t) f9,  (int32_t) f10, (int32_t) f11,
              (int32_t) f12, (int32_t) f13, (int32_t) f14, (int32_t) f15)) { }

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(int32_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint32_t) : m(a.derived().m) { }

    ENOKI_CONVERT(float) {
        m = std::is_signed_v<Value> ? _mm512_cvttps_epi32(a.derived().m)
                                    : _mm512_cvttps_epu32(a.derived().m);
    }

    ENOKI_CONVERT(double) {
        m = std::is_signed_v<Value>
                ? detail::concat(_mm512_cvttpd_epi32(low(a).m),
                                 _mm512_cvttpd_epi32(high(a).m))
                : detail::concat(_mm512_cvttpd_epu32(low(a).m),
                                 _mm512_cvttpd_epu32(high(a).m));
    }

    ENOKI_CONVERT(int64_t)
        : m(detail::concat(_mm512_cvtepi64_epi32(low(a).m),
                           _mm512_cvtepi64_epi32(high(a).m))) { }

    ENOKI_CONVERT(uint64_t)
        : m(detail::concat(_mm512_cvtepi64_epi32(low(a).m),
                           _mm512_cvtepi64_epi32(high(a).m))) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(float) : m(_mm512_castps_si512(a.derived().m)) { }
    ENOKI_REINTERPRET(int32_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint32_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm512_castsi512_si256(m); }
    ENOKI_INLINE Array2 high_() const {
        #if defined(ENOKI_X86_AVX512DQ)
            return _mm512_extracti32x8_epi32(m, 1);
        #else
            return _mm512_extracti64x4_epi64(m, 1);
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm512_add_epi32(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm512_sub_epi32(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm512_mullo_epi32(m, a.m); }

    template <typename T>
    ENOKI_INLINE Derived or_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_mov_epi32(m, a.k, _mm512_set1_epi32(int32_t(-1)));
        else
            return _mm512_or_epi32(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived and_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_maskz_mov_epi32(a.k, m);
        else
            return _mm512_and_epi32(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived andnot_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_mov_epi32(m, a.k, _mm512_setzero_si512());
        else
            return _mm512_andnot_epi32(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived xor_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_xor_epi32(m, a.k, m, _mm512_set1_epi32(int32_t(-1)));
        else
            return _mm512_xor_epi32(m, a.m);
    }

    template <size_t k> ENOKI_INLINE Derived sl_() const {
        return _mm512_slli_epi32(m, (int) k);
    }

    template <size_t k> ENOKI_INLINE Derived sr_() const {
        return std::is_signed_v<Value> ? _mm512_srai_epi32(m, (int) k)
                                       : _mm512_srli_epi32(m, (int) k);
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm512_sll_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        return std::is_signed_v<Value>
                   ? _mm512_sra_epi32(m, _mm_set1_epi64x((long long) k))
                   : _mm512_srl_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        return _mm512_sllv_epi32(m, k.m);
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        return std::is_signed_v<Value> ? _mm512_srav_epi32(m, k.m)
                                       : _mm512_srlv_epi32(m, k.m);
    }

    ENOKI_INLINE Derived rol_(Ref k) const { return _mm512_rolv_epi32(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm512_rorv_epi32(m, k.m); }

    template <size_t Imm>
    ENOKI_INLINE Derived rol_() const { return _mm512_rol_epi32(m, (int) Imm); }

    template <size_t Imm>
    ENOKI_INLINE Derived ror_() const { return _mm512_ror_epi32(m, (int) Imm); }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_LT));  }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_GT));  }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_LE));  }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_GE));  }
    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_EQ));  }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi32_mask(m, a.m, _MM_CMPINT_NE)); }

    ENOKI_INLINE Derived min_(Ref a) const {
        return std::is_signed_v<Value> ? _mm512_min_epi32(a.m, m)
                                       : _mm512_min_epu32(a.m, m);
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        return std::is_signed_v<Value> ? _mm512_max_epi32(a.m, m)
                                       : _mm512_max_epu32(a.m, m);
    }

    ENOKI_INLINE Derived abs_() const {
        return std::is_signed_v<Value> ? _mm512_abs_epi32(m) : m;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        return _mm512_mask_blend_epi32(m.k, f.m, t.m);
    }

    template <size_t I0,  size_t I1,  size_t I2,  size_t I3,  size_t I4,
              size_t I5,  size_t I6,  size_t I7,  size_t I8,  size_t I9,
              size_t I10, size_t I11, size_t I12, size_t I13, size_t I14,
              size_t I15>
    ENOKI_INLINE Derived shuffle_() const {
        const __m512i idx =
            _mm512_setr_epi32(I0, I1, I2, I3, I4, I5, I6, I7, I8,
                              I9, I10, I11, I12, I13, I14, I15);
        return _mm512_permutexvar_epi32(idx, m);
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return _mm512_permutexvar_epi32(index.m, m);
    }

    ENOKI_INLINE Derived mulhi_(Ref a) const {
        auto blend = mask_t<Derived>::from_k(0b0101010101010101);
        Derived even, odd;

        if constexpr (std::is_signed_v<Value>) {
            even.m = _mm512_srli_epi64(_mm512_mul_epi32(m, a.m), 32);
            odd.m = _mm512_mul_epi32(_mm512_srli_epi64(m, 32),
                                     _mm512_srli_epi64(a.m, 32));
        } else {
            even.m = _mm512_srli_epi64(_mm512_mul_epu32(m, a.m), 32);
            odd.m = _mm512_mul_epu32(_mm512_srli_epi64(m, 32),
                                     _mm512_srli_epi64(a.m, 32));
        }

        return select(blend, even, odd);
    }

#if defined(ENOKI_X86_AVX512CD)
    ENOKI_INLINE Derived lzcnt_() const { return _mm512_lzcnt_epi32(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(32) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

#if defined(ENOKI_X86_AVX512VPOPCNTDQ)
    ENOKI_INLINE Derived popcnt_() const { return _mm512_popcnt_epi32(m); }
#endif

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    //! @}
    // -----------------------------------------------------------------------
    //

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_store_si512((__m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        assert((uintptr_t) ptr % 64 == 0);
        _mm512_mask_store_epi32((__m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64), mask.k, m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm512_storeu_si512((__m512i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        _mm512_mask_storeu_epi32((__m512i *) ptr, mask.k, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_load_si512((const __m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        assert((uintptr_t) ptr % 64 == 0);
        return _mm512_maskz_load_epi32(mask.k, (const __m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm512_loadu_si512((const __m512i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return _mm512_maskz_loadu_epi32(mask.k, (const __m512i *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return _mm512_setzero_si512(); }

#if defined(ENOKI_X86_AVX512PF)
    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *ptr, const Index &index, const Mask &mask) {
        constexpr auto Hint = Level == 1 ? _MM_HINT_T0 : _MM_HINT_T1;

        if constexpr (sizeof(scalar_t<Index>) == 4) {
            if (Write)
                _mm512_mask_prefetch_i32scatter_ps((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i32gather_ps(index.m, mask.k, ptr, Stride, Hint);
        } else {
            if (Write) {
                _mm512_mask_prefetch_i64scatter_ps((void *) ptr, low(mask).k, low(index).m, Stride, Hint);
                _mm512_mask_prefetch_i64scatter_ps((void *) ptr, high(mask).k, high(index).m, Stride, Hint);
            } else {
                _mm512_mask_prefetch_i64gather_ps(low(index).m, low(mask).k, ptr, Stride, Hint);
                _mm512_mask_prefetch_i64gather_ps(high(index).m, high(mask).k, ptr, Stride, Hint);
            }
        }
    }
#endif

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            return _mm512_mask_i32gather_epi32(_mm512_setzero_si512(), mask.k, index.m, (const float *) ptr, Stride);
        } else {
            return detail::concat(
                _mm512_mask_i64gather_epi32(_mm256_setzero_si256(),  low(mask).k,  low(index).m, (const float *) ptr, Stride),
                _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), high(mask).k, high(index).m, (const float *) ptr, Stride));
        }
    }

    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4) {
            _mm512_mask_i32scatter_epi32(ptr, mask.k, index.m, m, Stride);
        } else {
            _mm512_mask_i64scatter_epi32(ptr, low(mask).k,   low(index).m,  low(derived()).m,  Stride);
            _mm512_mask_i64scatter_epi32(ptr, high(mask).k, high(index).m, high(derived()).m, Stride);
        }
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return (Value) _mm_cvtsi128_si32(_mm512_castsi512_si128(_mm512_maskz_compress_epi32(mask.k, m)));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        _mm512_storeu_si512((__m512i *) ptr, _mm512_maskz_compress_epi32(mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }

#if defined(ENOKI_X86_AVX512CD)
    template <size_t Stride, typename Index, typename Mask, typename Func, typename... Args>
    static ENOKI_INLINE void transform_(void *mem,
                                        Index index,
                                        const Mask &mask,
                                        const Func &func,
                                        const Args &... args) {
        Derived values = _mm512_mask_i32gather_epi32(
            _mm512_undefined_epi32(), mask.k, index.m, mem, (int) Stride);

        index.m = _mm512_mask_mov_epi32(_mm512_set1_epi32(-1), mask.k, index.m);

        __m512i conflicts = _mm512_conflict_epi32(index.m);
        __m512i perm_idx  = _mm512_sub_epi32(_mm512_set1_epi32(31), _mm512_lzcnt_epi32(conflicts));
        __mmask16 todo    = _mm512_mask_test_epi32_mask(mask.k, conflicts, _mm512_set1_epi32(-1));

        func(values, args...);

        ENOKI_NOUNROLL while (ENOKI_UNLIKELY(!_mm512_kortestz(todo, todo))) {
            __mmask16 cur = _mm512_mask_testn_epi32_mask(
                todo, conflicts, _mm512_broadcastmw_epi32(todo));
            values.m = _mm512_mask_permutexvar_epi32(values.m, cur, perm_idx, values.m);

            __m512i backup(values.m);
            func(values, args...);

            values.m = _mm512_mask_mov_epi32(backup, cur, values.m);
            todo = _mm512_kxor(todo, cur);
        }

        _mm512_mask_i32scatter_epi32(mem, mask.k, index.m, values.m, (int) Stride);
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm512_mask_mov_epi32(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm512_mask_add_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm512_mask_sub_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm512_mask_mullo_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm512_mask_or_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm512_mask_and_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm512_mask_xor_epi32(m, mask.k, m, a.m); }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using AVX512 intrinsics (64 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(64)
    StaticArrayImpl<Value_, 8, IsMask_, Derived_, enable_if_int64_t<Value_>>
  : StaticArrayBase<Value_, 8, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 8, __m512i)


    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(const Value &value) : m(_mm512_set1_epi64((long long) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value f0, Value f1, Value f2, Value f3,
                                 Value f4, Value f5, Value f6, Value f7)
        : m(_mm512_setr_epi64((long long) f0, (long long) f1, (long long) f2,
                              (long long) f3, (long long) f4, (long long) f5,
                              (long long) f6, (long long) f7)) { }

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512DQ)
    ENOKI_CONVERT(float) {
        m = std::is_signed_v<Value> ? _mm512_cvttps_epi64(a.derived().m)
                                    : _mm512_cvttps_epu64(a.derived().m);
    }
#else
    /* Emulate float -> uint64 conversion instead of falling
       back to scalar operations. This is quite a bit faster (~4x!) */

    ENOKI_CONVERT(float) {
        using Int32  = int_array_t<Derived2>;
        using UInt32 = uint_array_t<Derived2>;
        using UInt64 = uint64_array_t<Derived2>;

        /* Shift out sign bit */
        auto b = reinterpret_array<UInt32>(a);
        b += b;

        auto mant = UInt64((b & 0xffffffu) | 0x1000000u);
        auto shift = (24 + 127) - Int32(sr<24>(b));
        auto abs_shift = UInt64(abs(shift));

        auto result = select(shift > 0, mant >> abs_shift, mant << abs_shift);

        if constexpr (std::is_signed_v<Value>)
            result[a < 0] = -result;

        m = result.m;
    }
#endif

    ENOKI_CONVERT(int32_t)
        : m(_mm512_cvtepi32_epi64(a.derived().m)) { }

    ENOKI_CONVERT(uint32_t)
        : m(_mm512_cvtepu32_epi64(a.derived().m)) { }

#if defined(ENOKI_X86_AVX512DQ)
    ENOKI_CONVERT(double) {
        m = std::is_signed_v<Value> ? _mm512_cvttpd_epi64(a.derived().m)
                                    : _mm512_cvttpd_epu64(a.derived().m);
    }
#else
    /* Emulate double -> uint64 conversion instead of falling
       back to scalar operations. This is quite a bit faster (>~11x!) */

    ENOKI_CONVERT(double) {
        using Int64  = int_array_t<Derived2>;
        using UInt64 = uint_array_t<Derived2>;

        /* Shift out sign bit */
        auto b = reinterpret_array<UInt64>(a);
        b += b;

        auto mant = (b & 0x1fffffffffffffull) | 0x20000000000000ull;
        auto shift = (53 + 1023) - Int64(sr<53>(b));
        auto abs_shift = UInt64(abs(shift));

        auto result = select(shift > 0, mant >> abs_shift, mant << abs_shift);

        if constexpr (std::is_signed_v<Value>)
            result[a < 0] = -result;

        m = result.m;
    }
#endif

    ENOKI_CONVERT(int64_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(double) : m(_mm512_castpd_si512(a.derived().m)) { }
    ENOKI_REINTERPRET(int64_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm512_castsi512_si256(m); }
    ENOKI_INLINE Array2 high_() const { return _mm512_extracti64x4_epi64(m, 1); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm512_add_epi64(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm512_sub_epi64(m, a.m); }

    ENOKI_INLINE Derived mul_(Ref a) const {
        #if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            return _mm512_mullo_epi64(m, a.m);
        #else
            __m512i h0    = _mm512_srli_epi64(m, 32);
            __m512i h1    = _mm512_srli_epi64(a.m, 32);
            __m512i low   = _mm512_mul_epu32(m, a.m);
            __m512i mix0  = _mm512_mul_epu32(m, h1);
            __m512i mix1  = _mm512_mul_epu32(h0, a.m);
            __m512i mix   = _mm512_add_epi64(mix0, mix1);
            __m512i mix_s = _mm512_slli_epi64(mix, 32);
            return  _mm512_add_epi64(mix_s, low);
        #endif
    }

    template <typename T>
    ENOKI_INLINE Derived or_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_mov_epi64(m, a.k, _mm512_set1_epi64(int64_t(-1)));
        else
            return _mm512_or_epi64(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived and_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_maskz_mov_epi64(a.k, m);
        else
            return _mm512_and_epi64(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived andnot_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_mov_epi64(m, a.k, _mm512_setzero_si512());
        else
            return _mm512_andnot_epi64(m, a.m);
    }

    template <typename T>
    ENOKI_INLINE Derived xor_ (const T &a) const {
        if constexpr (is_mask_v<T>)
            return _mm512_mask_xor_epi64(m, a.k, m, _mm512_set1_epi64(int64_t(-1)));
        else
            return _mm512_xor_epi64(m, a.m);
    }

    template <size_t k> ENOKI_INLINE Derived sl_() const {
        return _mm512_slli_epi64(m, (int) k);
    }

    template <size_t k> ENOKI_INLINE Derived sr_() const {
        return std::is_signed_v<Value> ? _mm512_srai_epi64(m, (int) k)
                                       : _mm512_srli_epi64(m, (int) k);
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm512_sll_epi64(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        return std::is_signed_v<Value>
                   ? _mm512_sra_epi64(m, _mm_set1_epi64x((long long) k))
                   : _mm512_srl_epi64(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        return _mm512_sllv_epi64(m, k.m);
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        return std::is_signed_v<Value> ? _mm512_srav_epi64(m, k.m)
                                       : _mm512_srlv_epi64(m, k.m);
    }

    ENOKI_INLINE Derived rol_(Ref k) const { return _mm512_rolv_epi64(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm512_rorv_epi64(m, k.m); }

    template <size_t Imm>
    ENOKI_INLINE Derived rol_() const { return _mm512_rol_epi64(m, (int) Imm); }

    template <size_t Imm>
    ENOKI_INLINE Derived ror_() const { return _mm512_ror_epi64(m, (int) Imm); }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_LT)); }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_GT)); }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_LE)); }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_GE)); }
    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_EQ)); }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>::from_k(_mm512_cmp_epi64_mask(m, a.m, _MM_CMPINT_NE)); }

    ENOKI_INLINE Derived min_(Ref a) const {
        return std::is_signed_v<Value> ? _mm512_min_epi64(a.m, m)
                                       : _mm512_min_epu64(a.m, m);
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        return std::is_signed_v<Value> ? _mm512_max_epi64(a.m, m)
                                       : _mm512_max_epu64(a.m, m);
    }

    ENOKI_INLINE Derived abs_() const {
        return std::is_signed_v<Value> ? _mm512_abs_epi64(m) : m;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        return _mm512_mask_blend_epi64(m.k, f.m, t.m);
    }

    template <size_t I0, size_t I1, size_t I2, size_t I3, size_t I4, size_t I5,
              size_t I6, size_t I7>
    ENOKI_INLINE Derived shuffle_() const {
        const __m512i idx =
            _mm512_setr_epi64(I0, I1, I2, I3, I4, I5, I6, I7);
        return _mm512_permutexvar_epi64(idx, m);
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return _mm512_permutexvar_epi64(index.m, m);
    }

    ENOKI_INLINE Derived mulhi_(Ref b) const {
        if (std::is_unsigned_v<Value>) {
            const __m512i low_bits = _mm512_set1_epi64(0xffffffffu);
            __m512i al = m, bl = b.m;
            __m512i ah = _mm512_srli_epi64(al, 32);
            __m512i bh = _mm512_srli_epi64(bl, 32);

            // 4x unsigned 32x32->64 bit multiplication
            __m512i albl = _mm512_mul_epu32(al, bl);
            __m512i albh = _mm512_mul_epu32(al, bh);
            __m512i ahbl = _mm512_mul_epu32(ah, bl);
            __m512i ahbh = _mm512_mul_epu32(ah, bh);

            // Calculate a possible carry from the low bits of the multiplication.
            __m512i carry = _mm512_add_epi64(
                _mm512_srli_epi64(albl, 32),
                _mm512_add_epi64(_mm512_and_epi64(albh, low_bits),
                                 _mm512_and_epi64(ahbl, low_bits)));

            __m512i s0 = _mm512_add_epi64(ahbh, _mm512_srli_epi64(carry, 32));
            __m512i s1 = _mm512_add_epi64(_mm512_srli_epi64(albh, 32),
                                          _mm512_srli_epi64(ahbl, 32));

            return _mm512_add_epi64(s0, s1);
        } else {
            const Derived mask(0xffffffff);
            const Derived a = derived();
            Derived ah = sr<32>(a), bh = sr<32>(b),
                    al = a & mask, bl = b & mask;

            Derived albl_hi = _mm512_srli_epi64(_mm512_mul_epu32(m, b.m), 32);

            Derived t = ah * bl + albl_hi;
            Derived w1 = al * bh + (t & mask);

            return ah * bh + sr<32>(t) + sr<32>(w1);
        }
    }

#if defined(ENOKI_X86_AVX512CD)
    ENOKI_INLINE Derived lzcnt_() const { return _mm512_lzcnt_epi64(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(64) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

#if defined(ENOKI_X86_AVX512VPOPCNTDQ)
    ENOKI_INLINE Derived popcnt_() const { return _mm512_popcnt_epi64(m); }
#endif

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    //! @}
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uint64_t) ptr % 64 == 0);
        _mm512_store_si512((__m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        assert((uint64_t) ptr % 64 == 0);
        _mm512_mask_store_epi64((__m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64), mask.k, m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm512_storeu_si512((__m512i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        _mm512_mask_storeu_epi64((__m512i *) ptr, mask.k, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uint64_t) ptr % 64 == 0);
        return _mm512_load_si512((const __m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        assert((uint64_t) ptr % 64 == 0);
        return _mm512_maskz_load_epi64(mask.k, (const __m512i *) ENOKI_ASSUME_ALIGNED(ptr, 64));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm512_loadu_si512((const __m512i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        return _mm512_maskz_loadu_epi64(mask.k, (const __m512i *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return _mm512_setzero_si512(); }

#if defined(ENOKI_X86_AVX512PF)
    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *ptr, const Index &index, const Mask &mask) {
        constexpr auto Hint = Level == 1 ? _MM_HINT_T0 : _MM_HINT_T1;

        if constexpr (sizeof(scalar_t<Index>) == 4) {
            if constexpr (Write)
                _mm512_mask_prefetch_i32scatter_pd((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i32gather_pd(index.m, mask.k, ptr, Stride, Hint);
        } else {
            if constexpr (Write)
                _mm512_mask_prefetch_i64scatter_pd((void *) ptr, mask.k, index.m, Stride, Hint);
            else
                _mm512_mask_prefetch_i64gather_pd(index.m, mask.k, ptr, Stride, Hint);
        }
    }
#endif

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            return _mm512_mask_i32gather_epi64(_mm512_setzero_si512(), mask.k, index.m, (const float *) ptr, Stride);
        else
            return _mm512_mask_i64gather_epi64(_mm512_setzero_si512(), mask.k, index.m, (const float *) ptr, Stride);
    }


    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm512_mask_i32scatter_epi64(ptr, mask.k, index.m, m, Stride);
        else
            _mm512_mask_i64scatter_epi64(ptr, mask.k, index.m, m, Stride);
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        return (Value) _mm_cvtsi128_si64(_mm512_castsi512_si128(_mm512_maskz_compress_epi64(mask.k, m)));
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        _mm512_storeu_si512((__m512i *) ptr, _mm512_maskz_compress_epi64(mask.k, m));
        size_t kn = (size_t) _mm_popcnt_u32(mask.k);
        ptr += kn;
        return kn;
    }

#if defined(ENOKI_X86_AVX512CD)
    template <size_t Stride, typename Index, typename Mask, typename Func, typename... Args>
    static ENOKI_INLINE void transform_(void *mem,
                                        Index index,
                                        const Mask &mask,
                                        const Func &func,
                                        const Args &... args) {
        Derived values = _mm512_mask_i64gather_epi64(
            _mm512_undefined_epi32(), mask.k, index.m, mem, (int) Stride);

        index.m = _mm512_mask_mov_epi64(_mm512_set1_epi64(-1), mask.k, index.m);

        __m512i conflicts = _mm512_conflict_epi64(index.m);
        __m512i perm_idx  = _mm512_sub_epi64(_mm512_set1_epi64(63), _mm512_lzcnt_epi64(conflicts));
        __mmask8 todo     = _mm512_mask_test_epi64_mask(mask.k, conflicts, _mm512_set1_epi64(-1));

        func(values, args...);

        ENOKI_NOUNROLL while (ENOKI_UNLIKELY(todo)) {
            __mmask8 cur = _mm512_mask_testn_epi64_mask(
                todo, conflicts, _mm512_broadcastmb_epi64(todo));
            values.m = _mm512_mask_permutexvar_epi64(values.m, cur, perm_idx, values.m);

            __m512i backup(values.m);
            func(values, args...);

            values.m = _mm512_mask_mov_epi64(backup, cur, values.m);
            todo ^= cur;
        }

        _mm512_mask_i64scatter_epi64(mem, mask.k, index.m, values.m, (int) Stride);
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm512_mask_mov_epi64(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm512_mask_add_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm512_mask_sub_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            m = _mm512_mask_mullo_epi64(m, mask.k, m, a.m);
        #else
            m = select(mask, a * derived(), derived()).m;
        #endif
    }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm512_mask_or_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm512_mask_and_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm512_mask_xor_epi64(m, mask.k, m, a.m); }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

template <typename Derived_>
ENOKI_DECLARE_KMASK(float, 16, Derived_, int)
template <typename Derived_>
ENOKI_DECLARE_KMASK(double, 8, Derived_, int)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 16, Derived_, enable_if_int32_t<Value_>)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 8, Derived_, enable_if_int64_t<Value_>)

NAMESPACE_END(enoki)
