/*
    enoki/array_avx2.h -- Packed SIMD array (AVX2 specialization)

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
template <typename Value>    struct is_native<Value, 8, enable_if_int32_t<Value>> : std::true_type { };
template <typename Value>    struct is_native<Value, 4, enable_if_int64_t<Value>> : std::true_type { };
template <typename Value>    struct is_native<Value, 3, enable_if_int64_t<Value>> : std::true_type { };
NAMESPACE_END(detail)

/// Partial overload of StaticArrayImpl using AVX intrinsics (32 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(32)
    StaticArrayImpl<Value_, 8, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayBase<Value_, 8, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 8, __m256i)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(_mm256_set1_epi32((int32_t) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3,
                                 Value v4, Value v5, Value v6, Value v7)
        : m(_mm256_setr_epi32((int32_t) v0, (int32_t) v1, (int32_t) v2, (int32_t) v3,
                              (int32_t) v4, (int32_t) v5, (int32_t) v6, (int32_t) v7)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(float) {
        if constexpr (std::is_signed_v<Value>) {
            m = _mm256_cvttps_epi32(a.derived().m);
        } else {
            #if defined(ENOKI_X86_AVX512VL)
                m = _mm256_cvttps_epu32(a.derived().m);
            #else
                constexpr uint32_t limit = 1u << 31;
                const __m256  limit_f = _mm256_set1_ps((float) limit);
                const __m256i limit_i = _mm256_set1_epi32((int) limit);

                __m256 v = a.derived().m;

                __m256i mask =
                    _mm256_castps_si256(_mm256_cmp_ps(v, limit_f, _CMP_GE_OQ));

                __m256i b2 = _mm256_add_epi32(
                    _mm256_cvttps_epi32(_mm256_sub_ps(v, limit_f)), limit_i);

                __m256i b1 = _mm256_cvttps_epi32(v);

                m = _mm256_blendv_epi8(b1, b2, mask);
            #endif
        }
    }

    ENOKI_CONVERT(int32_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint32_t) : m(a.derived().m) { }

    ENOKI_CONVERT(double) {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512F)
                m = _mm512_cvttpd_epi32(a.derived().m);
            #else
                m = detail::concat(_mm256_cvttpd_epi32(low(a).m),
                                   _mm256_cvttpd_epi32(high(a).m));
            #endif
        } else {
            #if defined(ENOKI_X86_AVX512F)
                m = _mm512_cvttpd_epu32(a.derived().m);
            #else
                ENOKI_TRACK_SCALAR("Constructor (converting, double[8] -> [u]int32[8])");
                for (size_t i = 0; i < Size; ++i)
                    coeff(i) = Value(a.derived().coeff(i));
            #endif
        }
    }

    ENOKI_CONVERT(int64_t) {
        #if defined(ENOKI_X86_AVX512F)
            m = _mm512_cvtepi64_epi32(a.derived().m);
        #else
            m = detail::mm512_cvtepi64_epi32(low(a).m, high(a).m);
        #endif
    }

    ENOKI_CONVERT(uint64_t) {
        #if defined(ENOKI_X86_AVX512F)
            m = _mm512_cvtepi64_epi32(a.derived().m);
        #else
            m = detail::mm512_cvtepi64_epi32(low(a).m, high(a).m);
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        uint64_t ival;
        memcpy(&ival, a.derived().data(), 8);
        __m128i value = _mm_cmpgt_epi8(detail::mm_cvtsi64_si128((long long) ival),
                                       _mm_setzero_si128());
        m = _mm256_cvtepi8_epi32(value);
    }

    ENOKI_REINTERPRET(float) : m(_mm256_castps_si256(a.derived().m)) { }
    ENOKI_REINTERPRET(int32_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint32_t) : m(a.derived().m) { }

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_REINTERPRET(double) : m(_mm256_movm_epi32(a.derived().k)) { }
    ENOKI_REINTERPRET(int64_t) : m(_mm256_movm_epi32(a.derived().k)) { }
    ENOKI_REINTERPRET(uint64_t) : m(_mm256_movm_epi32(a.derived().k)) { }
#elif defined(ENOKI_X86_AVX512F)
    ENOKI_REINTERPRET(double)
        : m(_mm512_castsi512_si256(_mm512_maskz_mov_epi32(
              (__mmask16) a.derived().k, _mm512_set1_epi32(int32_t(-1))))) { }
    ENOKI_REINTERPRET(int64_t)
        : m(_mm512_castsi512_si256(_mm512_maskz_mov_epi32(
              (__mmask16) a.derived().k, _mm512_set1_epi32(int32_t(-1))))) { }
    ENOKI_REINTERPRET(uint64_t)
        : m(_mm512_castsi512_si256(_mm512_maskz_mov_epi32(
              (__mmask16) a.derived().k, _mm512_set1_epi32(int32_t(-1))))) { }
#else
    ENOKI_REINTERPRET(double)
        : m(detail::mm512_cvtepi64_epi32(_mm256_castpd_si256(low(a).m),
                                         _mm256_castpd_si256(high(a).m))) { }
    ENOKI_REINTERPRET(int64_t)
        : m(detail::mm512_cvtepi64_epi32(low(a).m, high(a).m)) { }
    ENOKI_REINTERPRET(uint64_t)
        : m(detail::mm512_cvtepi64_epi32(low(a).m, high(a).m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm256_castsi256_si128(m); }
    ENOKI_INLINE Array2 high_() const { return _mm256_extractf128_si256(m, 1); }

    //! @}
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm256_add_epi32(m, a.m);   }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm256_sub_epi32(m, a.m);   }
    ENOKI_INLINE Derived mul_(Ref a) const { return _mm256_mullo_epi32(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_mov_epi32(m, a.k, _mm256_set1_epi32(-1));
            else
        #endif
        return _mm256_or_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_maskz_mov_epi32(a.k, m);
            else
        #endif
        return _mm256_and_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_xor_epi32(m, a.k, m, _mm256_set1_epi32(-1));
            else
        #endif
        return _mm256_xor_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_mov_epi32(m, a.k, _mm256_setzero_si256());
            else
        #endif
        return _mm256_andnot_si256(a.m, m);
    }

    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        return _mm256_slli_epi32(m, (int) Imm);
    }

    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        return std::is_signed_v<Value> ? _mm256_srai_epi32(m, (int) Imm)
                                       : _mm256_srli_epi32(m, (int) Imm);
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm256_sll_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        return std::is_signed_v<Value>
                   ? _mm256_sra_epi32(m, _mm_set1_epi64x((long long) k))
                   : _mm256_srl_epi32(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        return _mm256_sllv_epi32(m, k.m);
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        return std::is_signed_v<Value> ? _mm256_srav_epi32(m, k.m)
                                       : _mm256_srlv_epi32(m, k.m);
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Imm> ENOKI_INLINE Derived rol_() const { return _mm256_rol_epi32(m, (int) Imm); }
    template <size_t Imm> ENOKI_INLINE Derived ror_() const { return _mm256_ror_epi32(m, (int) Imm); }
    ENOKI_INLINE Derived rol_(Ref k) const { return _mm256_rolv_epi32(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm256_rorv_epi32(m, k.m); }
#endif

    ENOKI_INLINE auto eq_(Ref a)  const {
        using Return = mask_t<Derived>;

        #if defined(ENOKI_X86_AVX512VL)
            return Return::from_k(_mm256_cmpeq_epi32_mask(m, a.m));
        #else
            return Return(_mm256_cmpeq_epi32(m, a.m));
        #endif
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(_mm256_cmpneq_epi32_mask(m, a.m));
        #else
            return ~eq_(a);
        #endif
    }

    ENOKI_INLINE auto lt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm256_cmpgt_epi32(a.m, m));
            } else {
                const __m256i offset = _mm256_set1_epi32((int32_t) 0x80000000ul);
                return Return(_mm256_cmpgt_epi32(_mm256_sub_epi32(a.m, offset),
                                                 _mm256_sub_epi32(m, offset)));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm256_cmplt_epi32_mask(m, a.m)
                                  : _mm256_cmplt_epu32_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm256_cmpgt_epi32(m, a.m));
            } else {
                const __m256i offset = _mm256_set1_epi32((int32_t) 0x80000000ul);
                return Return(_mm256_cmpgt_epi32(_mm256_sub_epi32(m, offset),
                                                 _mm256_sub_epi32(a.m, offset)));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm256_cmpgt_epi32_mask(m, a.m)
                                  : _mm256_cmpgt_epu32_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto le_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm256_cmple_epi32_mask(m, a.m)
                                           : _mm256_cmple_epu32_mask(m, a.m));
        #else
            return ~gt_(a);
        #endif
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm256_cmpge_epi32_mask(m, a.m)
                                           : _mm256_cmpge_epu32_mask(m, a.m));
        #else
            return ~lt_(a);
        #endif
    }

    ENOKI_INLINE Derived min_(Ref a) const {
        return std::is_signed_v<Value> ? _mm256_min_epi32(a.m, m)
                                       : _mm256_min_epu32(a.m, m);
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        return std::is_signed_v<Value> ? _mm256_max_epi32(a.m, m)
                                       : _mm256_max_epu32(a.m, m);
    }

    ENOKI_INLINE Derived abs_() const {
        return std::is_signed_v<Value> ? _mm256_abs_epi32(m) : m;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        #if !defined(ENOKI_X86_AVX512VL)
            return _mm256_blendv_epi8(f.m, t.m, m.m);
        #else
            return _mm256_mask_blend_epi32(m.k, f.m, t.m);
        #endif
    }

    template <int I0, int I1, int I2, int I3, int I4, int I5, int I6, int I7>
    ENOKI_INLINE Derived shuffle_() const {
        return _mm256_permutevar8x32_epi32(m,
            _mm256_setr_epi32(I0, I1, I2, I3, I4, I5, I6, I7));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return _mm256_permutevar8x32_epi32(m, index.m);
    }

    ENOKI_INLINE Derived mulhi_(Ref a) const {
        Derived even, odd;

        if constexpr (std::is_signed_v<Value>) {
            even.m = _mm256_srli_epi64(_mm256_mul_epi32(m, a.m), 32);
            odd.m = _mm256_mul_epi32(_mm256_srli_epi64(m, 32), _mm256_srli_epi64(a.m, 32));
        } else {
            even.m = _mm256_srli_epi64(_mm256_mul_epu32(m, a.m), 32);
            odd.m = _mm256_mul_epu32(_mm256_srli_epi64(m, 32), _mm256_srli_epi64(a.m, 32));
        }

        #if defined(ENOKI_X86_AVX512VL)
            const mask_t<Derived> blend = mask_t<Derived>::from_k(0b01010101);
        #else
            const mask_t<Derived> blend(Value(-1), Value(0), Value(-1), Value(0),
                                        Value(-1), Value(0), Value(-1), Value(0));
        #endif

        return select(blend, even, odd);
    }

#if defined(ENOKI_X86_AVX512CD) && defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived lzcnt_() const { return _mm256_lzcnt_epi32(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(32) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    ENOKI_INLINE bool all_() const { return _mm256_movemask_ps(_mm256_castsi256_ps(m)) == 0xFF; }
    ENOKI_INLINE bool any_() const { return _mm256_movemask_ps(_mm256_castsi256_ps(m)) != 0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm256_movemask_ps(_mm256_castsi256_ps(m)); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm256_mask_mov_epi32(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm256_mask_add_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm256_mask_sub_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm256_mask_mullo_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm256_mask_or_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm256_mask_and_epi32(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm256_mask_xor_epi32(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 32 == 0);
        _mm256_store_si256((__m256i *) ENOKI_ASSUME_ALIGNED(ptr, 32), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm256_mask_store_epi32((Value *) ptr, mask.k, m);
        #else
            _mm256_maskstore_epi32((int *) ptr, mask.m, m);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm256_storeu_si256((__m256i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm256_mask_storeu_epi32((Value *) ptr, mask.k, m);
        #else
            _mm256_maskstore_epi32((int *) ptr, mask.m, m);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 32 == 0);
        return _mm256_load_si256((const __m256i *) ENOKI_ASSUME_ALIGNED(ptr, 32));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm256_maskz_load_epi32(mask.k, ptr);
        #else
            return _mm256_maskload_epi32((const int *) ptr, mask.m);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm256_loadu_si256((const __m256i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm256_maskz_loadu_epi32(mask.k, ptr);
        #else
            return _mm256_maskload_epi32((const int *) ptr, mask.m);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm256_setzero_si256(); }

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL) && defined(ENOKI_X86_AVX512DQ)
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm256_mmask_i32gather_epi32(_mm256_setzero_si256(), mask.k, index.m, (const int *) ptr, Stride);
            else
                return _mm512_mask_i64gather_epi32(_mm256_setzero_si256(), mask.k, index.m, (const int *) ptr, Stride);
        #else
            if constexpr (sizeof(scalar_t<Index>) == 4) {
                return _mm256_mask_i32gather_epi32(
                    _mm256_setzero_si256(), (const int *) ptr, index.m, mask.m, Stride);
            } else {
                return Derived(
                    _mm256_mask_i64gather_epi32(_mm_setzero_si128(), (const int *) ptr, low(index).m, low(mask).m, Stride),
                    _mm256_mask_i64gather_epi32(_mm_setzero_si128(), (const int *) ptr, high(index).m, high(mask).m, Stride)
                );
            }
        #endif
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm256_mask_i32scatter_epi32(ptr, mask.k, index.m, m, Stride);
        else
            _mm512_mask_i64scatter_epi32(ptr, mask.k, index.m, m, Stride);
    }
#endif

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        #if !defined(ENOKI_X86_AVX512VL)
            unsigned int k = (unsigned int) _mm256_movemask_ps(_mm256_castsi256_ps(mask.m));
            return coeff((size_t) (detail::tzcnt_scalar(k) & 7));
        #else
            return (Value) _mm_cvtsi128_si32(_mm256_castsi256_si128(
                _mm256_mask_compress_epi32(_mm256_setzero_si256(), mask.k, m)));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            size_t kn = (size_t) _mm_popcnt_u32(mask.k);
            _mm256_storeu_si256((__m256i *) ptr, _mm256_maskz_compress_epi32(mask.k, m));
            ptr += kn;
            return kn;
        #elif defined(ENOKI_X86_64) // requires _pdep_u64
            /** Clever BMI2-based partitioning algorithm by Christoph Diegelmann
                see https://goo.gl/o3ysMN for context */

            unsigned int k = (unsigned int) _mm256_movemask_epi8(mask.m);
            uint32_t wanted_indices = _pext_u32(0x76543210, k);
            uint64_t expanded_indices = _pdep_u64((uint64_t) wanted_indices,
                                                  0x0F0F0F0F0F0F0F0Full);
            size_t kn = (size_t) (_mm_popcnt_u32(k) >> 2);

            __m128i bytevec = detail::mm_cvtsi64_si128((long long) expanded_indices);
            __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);
            __m256i perm = _mm256_permutevar8x32_epi32(m, shufmask);

            _mm256_storeu_si256((__m256i *) ptr, perm);
            ptr += kn;
            return kn;
        #else
            return Base::compress_(ptr, mask);
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl using AVX intrinsics (64 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(32)
    StaticArrayImpl<Value_, 4, IsMask_, Derived_, enable_if_int64_t<Value_>>
  : StaticArrayBase<Value_, 4, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 4, __m256i)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(const Value &value)
        : m(_mm256_set1_epi64x((long long) value)) { }

    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3)
        : m(_mm256_setr_epi64x((long long) v0, (long long) v1,
                               (long long) v2, (long long) v3)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
    ENOKI_CONVERT(float) {
        m = std::is_signed_v<Value> ? _mm256_cvttps_epi64(a.derived().m)
                                    : _mm256_cvttps_epu64(a.derived().m);
    }

    ENOKI_CONVERT(double) {
        m = std::is_signed_v<Value> ? _mm256_cvttpd_epi64(a.derived().m)
                                    : _mm256_cvttpd_epu64(a.derived().m);
    }
#endif
    ENOKI_CONVERT(int32_t)  : m(_mm256_cvtepi32_epi64(a.derived().m)) { }
    ENOKI_CONVERT(uint32_t) : m(_mm256_cvtepu32_epi64(a.derived().m)) { }

    ENOKI_CONVERT(int64_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(bool) {
        int ival;
        memcpy(&ival, a.derived().data(), 4);
        m = _mm256_cvtepi8_epi64(
            _mm_cmpgt_epi8(_mm_cvtsi32_si128(ival), _mm_setzero_si128()));
    }

    ENOKI_REINTERPRET(float)
        : m(_mm256_cvtepi32_epi64(_mm_castps_si128(a.derived().m))) { }
    ENOKI_REINTERPRET(int32_t) : m(_mm256_cvtepi32_epi64(a.derived().m)) { }
    ENOKI_REINTERPRET(uint32_t) : m(_mm256_cvtepi32_epi64(a.derived().m)) { }

    ENOKI_REINTERPRET(double) : m(_mm256_castpd_si256(a.derived().m)) { }
    ENOKI_REINTERPRET(int64_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint64_t) : m(a.derived().m) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m(detail::concat(a1.m, a2.m)) { }

    ENOKI_INLINE Array1 low_()  const { return _mm256_castsi256_si128(m); }
    ENOKI_INLINE Array2 high_() const { return _mm256_extractf128_si256(m, 1); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return _mm256_add_epi64(m, a.m);   }
    ENOKI_INLINE Derived sub_(Ref a) const { return _mm256_sub_epi64(m, a.m);   }

    ENOKI_INLINE Derived mul_(Ref a) const {
        #if defined(ENOKI_X86_AVX512DQ) && defined(ENOKI_X86_AVX512VL)
            return _mm256_mullo_epi64(m, a.m);
        #else
            __m256i h0    = _mm256_srli_epi64(m, 32);
            __m256i h1    = _mm256_srli_epi64(a.m, 32);
            __m256i low   = _mm256_mul_epu32(m, a.m);
            __m256i mix0  = _mm256_mul_epu32(m, h1);
            __m256i mix1  = _mm256_mul_epu32(h0, a.m);
            __m256i mix   = _mm256_add_epi64(mix0, mix1);
            __m256i mix_s = _mm256_slli_epi64(mix, 32);
            return _mm256_add_epi64(mix_s, low);
        #endif
    }

    ENOKI_INLINE Derived mulhi_(Ref b) const {
        if constexpr (std::is_unsigned_v<Value>) {
            const __m256i low_bits = _mm256_set1_epi64x(0xffffffffu);
            __m256i al = m, bl = b.m;

            __m256i ah = _mm256_srli_epi64(al, 32);
            __m256i bh = _mm256_srli_epi64(bl, 32);

            // 4x unsigned 32x32->64 bit multiplication
            __m256i albl = _mm256_mul_epu32(al, bl);
            __m256i albh = _mm256_mul_epu32(al, bh);
            __m256i ahbl = _mm256_mul_epu32(ah, bl);
            __m256i ahbh = _mm256_mul_epu32(ah, bh);

            // Calculate a possible carry from the low bits of the multiplication.
            __m256i carry = _mm256_add_epi64(
                _mm256_srli_epi64(albl, 32),
                _mm256_add_epi64(_mm256_and_si256(albh, low_bits),
                                 _mm256_and_si256(ahbl, low_bits)));

            __m256i s0 = _mm256_add_epi64(ahbh, _mm256_srli_epi64(carry, 32));
            __m256i s1 = _mm256_add_epi64(_mm256_srli_epi64(albh, 32),
                                          _mm256_srli_epi64(ahbl, 32));

            return _mm256_add_epi64(s0, s1);

        } else {
            const Derived mask(0xffffffff);
            const Derived a = derived();
            Derived ah = sr<32>(a), bh = sr<32>(b),
                    al = a & mask, bl = b & mask;

            Derived albl_hi = _mm256_srli_epi64(_mm256_mul_epu32(m, b.m), 32);

            Derived t = ah * bl + albl_hi;
            Derived w1 = al * bh + (t & mask);

            return ah * bh + sr<32>(t) + sr<32>(w1);
        }
    }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_mov_epi64(m, a.k, _mm256_set1_epi64x(-1));
            else
        #endif
        return _mm256_or_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_maskz_mov_epi64(a.k, m);
            else
        #endif
        return _mm256_and_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_xor_epi64(m, a.k, m, _mm256_set1_epi64x(-1));
            else
        #endif
        return _mm256_xor_si256(m, a.m);
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (is_mask_v<T>)
                return _mm256_mask_mov_epi64(m, a.k, _mm256_setzero_si256());
            else
        #endif
        return _mm256_andnot_si256(a.m, m);
    }

    template <size_t k> ENOKI_INLINE Derived sl_() const {
        return _mm256_slli_epi64(m, (int) k);
    }

    template <size_t k> ENOKI_INLINE Derived sr_() const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm256_srai_epi64(m, (int) k);
            #else
                const __m256i offset = _mm256_set1_epi64x((long long) 0x8000000000000000ull);
                __m256i s1 = _mm256_srli_epi64(_mm256_add_epi64(m, offset), (int) k);
                __m256i s2 = _mm256_srli_epi64(offset, (int) k);
                return _mm256_sub_epi64(s1, s2);
            #endif
        } else {
            return _mm256_srli_epi64(m, (int) k);
        }
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return _mm256_sll_epi64(m, _mm_set1_epi64x((long long) k));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm256_sra_epi64(m, _mm_set1_epi64x((long long) k));
            #else
                const __m256i offset = _mm256_set1_epi64x((long long) 0x8000000000000000ull);
                __m128i s0 = _mm_set1_epi64x((long long) k);
                __m256i s1 = _mm256_srl_epi64(_mm256_add_epi64(m, offset), s0);
                __m256i s2 = _mm256_srl_epi64(offset, s0);
                return _mm256_sub_epi64(s1, s2);
            #endif
        } else {
            return _mm256_srl_epi64(m, _mm_set1_epi64x((long long) k));
        }
    }

    ENOKI_INLINE Derived sl_(Ref k) const {
        return _mm256_sllv_epi64(m, k.m);
    }

    ENOKI_INLINE Derived sr_(Ref k) const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm256_srav_epi64(m, k.m);
            #else
                const __m256i offset = _mm256_set1_epi64x((long long) 0x8000000000000000ull);
                __m256i s1 = _mm256_srlv_epi64(_mm256_add_epi64(m, offset), k.m);
                __m256i s2 = _mm256_srlv_epi64(offset, k.m);
                return _mm256_sub_epi64(s1, s2);
            #endif
        } else {
            return _mm256_srlv_epi64(m, k.m);
        }
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Imm> ENOKI_INLINE Derived rol_() const { return _mm256_rol_epi64(m, (int) Imm); }
    template <size_t Imm> ENOKI_INLINE Derived ror_() const { return _mm256_ror_epi64(m, (int) Imm); }
    ENOKI_INLINE Derived rol_(Ref k) const { return _mm256_rolv_epi64(m, k.m); }
    ENOKI_INLINE Derived ror_(Ref k) const { return _mm256_rorv_epi64(m, k.m); }
#endif

    ENOKI_INLINE auto eq_(Ref a)  const {
        using Return = mask_t<Derived>;

        #if defined(ENOKI_X86_AVX512VL)
            return Return::from_k(_mm256_cmpeq_epi64_mask(m, a.m));
        #else
            return Return(_mm256_cmpeq_epi64(m, a.m));
        #endif
    }

    ENOKI_INLINE auto neq_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(_mm256_cmpneq_epi64_mask(m, a.m));
        #else
            return ~eq_(a);
        #endif
    }

    ENOKI_INLINE auto lt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm256_cmpgt_epi64(a.m, m));
            } else {
                const __m256i offset =
                    _mm256_set1_epi64x((long long) 0x8000000000000000ull);
                return Return(_mm256_cmpgt_epi64(
                    _mm256_sub_epi64(a.m, offset),
                    _mm256_sub_epi64(m, offset)
                ));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm256_cmplt_epi64_mask(m, a.m)
                                  : _mm256_cmplt_epu64_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        using Return = mask_t<Derived>;

        #if !defined(ENOKI_X86_AVX512VL)
            if constexpr (std::is_signed_v<Value>) {
                return Return(_mm256_cmpgt_epi64(m, a.m));
            } else {
                const __m256i offset =
                    _mm256_set1_epi64x((long long) 0x8000000000000000ull);
                return Return(_mm256_cmpgt_epi64(
                    _mm256_sub_epi64(m, offset),
                    _mm256_sub_epi64(a.m, offset)
                ));
            }
        #else
            return Return::from_k(std::is_signed_v<Value>
                                  ? _mm256_cmpgt_epi64_mask(m, a.m)
                                  : _mm256_cmpgt_epu64_mask(m, a.m));
        #endif
    }

    ENOKI_INLINE auto le_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm256_cmple_epi64_mask(m, a.m)
                                           : _mm256_cmple_epu64_mask(m, a.m));
        #else
            return ~gt_(a);
        #endif
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return mask_t<Derived>::from_k(std::is_signed_v<Value>
                                           ? _mm256_cmpge_epi64_mask(m, a.m)
                                           : _mm256_cmpge_epu64_mask(m, a.m));
        #else
            return ~lt_(a);
        #endif
    }

    ENOKI_INLINE Derived min_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return std::is_signed_v<Value> ? _mm256_min_epi64(a.m, m)
                                           : _mm256_min_epu64(a.m, m);
        #else
            return select(derived() < a, derived(), a);
        #endif
    }

    ENOKI_INLINE Derived max_(Ref a) const {
        #if defined(ENOKI_X86_AVX512VL)
            return std::is_signed_v<Value> ? _mm256_max_epi64(a.m, m)
                                           : _mm256_max_epu64(a.m, m);
        #else
            return select(derived() > a, derived(), a);
        #endif
    }

    ENOKI_INLINE Derived abs_() const {
        if constexpr (std::is_signed_v<Value>) {
            #if defined(ENOKI_X86_AVX512VL)
                return _mm256_abs_epi64(m);
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
            return _mm256_blendv_epi8(f.m, t.m, m.m);
        #else
            return _mm256_mask_blend_epi64(m.k, f.m, t.m);
        #endif
    }

    template <int I0, int I1, int I2, int I3>
    ENOKI_INLINE Derived shuffle_() const {
        return _mm256_permute4x64_epi64(m, _MM_SHUFFLE(I3, I2, I1, I0));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }

#if defined(ENOKI_X86_AVX512CD) && defined(ENOKI_X86_AVX512VL)
    ENOKI_INLINE Derived lzcnt_() const { return _mm256_lzcnt_epi64(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(64) - lzcnt(~derived() & (derived() - Value(1))); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Masked versions of key operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_X86_AVX512VL)
    template <typename Mask>
    ENOKI_INLINE void massign_(const Derived &a, const Mask &mask) { m = _mm256_mask_mov_epi64(m, mask.k, a.m); }
    template <typename Mask>
    ENOKI_INLINE void madd_   (const Derived &a, const Mask &mask) { m = _mm256_mask_add_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void msub_   (const Derived &a, const Mask &mask) { m = _mm256_mask_sub_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mmul_   (const Derived &a, const Mask &mask) { m = _mm256_mask_mullo_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mor_    (const Derived &a, const Mask &mask) { m = _mm256_mask_or_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mand_   (const Derived &a, const Mask &mask) { m = _mm256_mask_and_epi64(m, mask.k, m, a.m); }
    template <typename Mask>
    ENOKI_INLINE void mxor_   (const Derived &a, const Mask &mask) { m = _mm256_mask_xor_epi64(m, mask.k, m, a.m); }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------
    //
    ENOKI_INLINE Value hsum_()  const { return hsum(low_() + high_()); }
    ENOKI_INLINE Value hprod_() const { return hprod(low_() * high_()); }
    ENOKI_INLINE Value hmin_()  const { return hmin(min(low_(), high_())); }
    ENOKI_INLINE Value hmax_()  const { return hmax(max(low_(), high_())); }

    ENOKI_INLINE bool all_() const { return _mm256_movemask_pd(_mm256_castsi256_pd(m)) == 0xF; }
    ENOKI_INLINE bool any_() const { return _mm256_movemask_pd(_mm256_castsi256_pd(m)) != 0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm256_movemask_pd(_mm256_castsi256_pd(m)); }
    ENOKI_INLINE size_t count_() const { return (size_t) _mm_popcnt_u32(bitmask_()); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 32 == 0);
        _mm256_store_si256((__m256i *) ENOKI_ASSUME_ALIGNED(ptr, 32), m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm256_mask_store_epi64(ptr, mask.k, m);
        #else
            _mm256_maskstore_epi64((long long *) ptr, mask.m, m);
        #endif
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        _mm256_storeu_si256((__m256i *) ptr, m);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            _mm256_mask_storeu_epi64(ptr, mask.k, m);
        #else
            _mm256_maskstore_epi64((long long *) ptr, mask.m, m);
        #endif
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 32 == 0);
        return _mm256_load_si256((const __m256i *) ENOKI_ASSUME_ALIGNED(ptr, 32));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm256_maskz_load_epi64(mask.k, ptr);
        #else
            return _mm256_maskload_epi64((const long long *) ptr, mask.m);
        #endif
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return _mm256_loadu_si256((const __m256i *) ptr);
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *ptr, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            return _mm256_maskz_loadu_epi64(mask.k, ptr);
        #else
            return _mm256_maskload_epi64((const long long *) ptr, mask.m);
        #endif
    }

    static ENOKI_INLINE Derived zero_() { return _mm256_setzero_si256(); }

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        #if defined(ENOKI_X86_AVX512VL)
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm256_mmask_i32gather_epi64(_mm256_setzero_si256(), mask.k, index.m, (const long long *) ptr, Stride);
            else
                return _mm256_mmask_i64gather_epi64(_mm256_setzero_si256(), mask.k, index.m, (const long long *) ptr, Stride);
        #else
            if constexpr (sizeof(scalar_t<Index>) == 4)
                return _mm256_mask_i32gather_epi64(_mm256_setzero_si256(), (const long long *) ptr, index.m, mask.m, Stride);
            else
                return _mm256_mask_i64gather_epi64(_mm256_setzero_si256(), (const long long *) ptr, index.m, mask.m, Stride);
        #endif
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        if constexpr (sizeof(scalar_t<Index>) == 4)
            _mm256_mask_i32scatter_epi64(ptr, mask.k, index.m, m, Stride);
        else
            _mm256_mask_i64scatter_epi64(ptr, mask.k, index.m, m, Stride);
    }
#endif

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            return (Value) detail::mm_cvtsi128_si64(_mm256_castsi256_si128(
                _mm256_mask_compress_epi64(_mm256_setzero_si256(), mask.k, m)));
        #else
            unsigned int k =
                (unsigned int) _mm256_movemask_pd(_mm256_castsi256_pd(mask.m));
            return coeff((size_t) (tzcnt(k) & 3));
        #endif
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        #if defined(ENOKI_X86_AVX512VL)
            size_t kn = (size_t) _mm_popcnt_u32(mask.k);
            _mm256_storeu_si256((__m256i *) ptr, _mm256_maskz_compress_epi64(mask.k, m));
            ptr += kn;
            return kn;
        #elif defined(ENOKI_X86_64) // requires _pdep_u64
            /** Clever BMI2-based partitioning algorithm by Christoph Diegelmann
                see https://goo.gl/o3ysMN for context */

            unsigned int k = (unsigned int) _mm256_movemask_epi8(mask.m);
            uint32_t wanted_indices = _pext_u32(0x76543210, k);
            uint64_t expanded_indices = _pdep_u64((uint64_t) wanted_indices,
                                                  0x0F0F0F0F0F0F0F0Full);
            size_t kn = (size_t) (_mm_popcnt_u32(k) >> 3);

            __m128i bytevec = detail::mm_cvtsi64_si128((long long) expanded_indices);
            __m256i shufmask = _mm256_cvtepu8_epi32(bytevec);

            __m256i perm = _mm256_permutevar8x32_epi32(m, shufmask);

            _mm256_storeu_si256((__m256i *) ptr, perm);
            ptr += kn;
            return kn;
        #else
            return Base::compress_(ptr, mask);
        #endif
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

/// Partial overload of StaticArrayImpl for the n=3 case (64 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct alignas(32)
    StaticArrayImpl<Value_, 3, IsMask_, Derived_, enable_if_int64_t<Value_>>
  : StaticArrayImpl<Value_, 4, IsMask_, Derived_> {
    using Base = StaticArrayImpl<Value_, 4, IsMask_, Derived_>;

    ENOKI_DECLARE_3D_ARRAY(StaticArrayImpl)

    template <int I0, int I1, int I2>
    ENOKI_INLINE Derived shuffle_() const {
        return Base::template shuffle_<I0, I1, I2, 3>();
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &idx) const {
        return Base::shuffle_(idx);
    }

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_() const {
        Value result = coeff(0);
        for (size_t i = 1; i < 3; ++i)
            result += coeff(i);
        return result;
    }

    ENOKI_INLINE Value hprod_() const {
        Value result = coeff(0);
        for (size_t i = 1; i < 3; ++i)
            result *= coeff(i);
        return result;
    }

    ENOKI_INLINE Value hmin_() const {
        Value result = coeff(0);
        for (size_t i = 1; i < 3; ++i)
            result = std::min(result, coeff(i));
        return result;
    }

    ENOKI_INLINE Value hmax_() const {
        Value result = coeff(0);
        for (size_t i = 1; i < 3; ++i)
            result = std::max(result, coeff(i));
        return result;
    }

    ENOKI_INLINE bool all_()  const { return (_mm256_movemask_pd(_mm256_castsi256_pd(m)) & 7) == 7;}
    ENOKI_INLINE bool any_()  const { return (_mm256_movemask_pd(_mm256_castsi256_pd(m)) & 7) != 0; }

    ENOKI_INLINE uint32_t bitmask_() const { return (uint32_t) _mm256_movemask_pd(_mm256_castsi256_pd(m)) & 7; }
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
            return mask_t<Derived>(_mm256_setr_epi64x(
                (int64_t) -1, (int64_t) -1, (int64_t) -1, (int64_t) 0));
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

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        return Base::template gather_<Stride>(ptr, index, mask & mask_());
    }

#if defined(ENOKI_X86_AVX512VL)
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        Base::template scatter_<Stride>(ptr, index, mask & mask_());
    }
#endif

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Value_ *&ptr, const Mask &mask) const {
        return Base::compress_(ptr, mask & mask_());
    }

    //! @}
    // -----------------------------------------------------------------------
} ENOKI_MAY_ALIAS;

#if defined(ENOKI_X86_AVX512VL)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 8, Derived_, enable_if_int32_t<Value_>)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 4, Derived_, enable_if_int64_t<Value_>)
template <typename Value_, typename Derived_>
ENOKI_DECLARE_KMASK(Value_, 3, Derived_, enable_if_int64_t<Value_>)
#endif

NAMESPACE_END(enoki)
