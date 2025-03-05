/*
    enoki/array_kmask.h -- Abstraction around AVX512 'k' mask registers

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

NAMESPACE_BEGIN(enoki)

/// SFINAE macro for constructors that reinterpret another type
#define ENOKI_REINTERPRET_KMASK(Value)                                         \
    template <typename Value2, typename Derived2, bool IsMask2,                \
              enable_if_t<detail::is_same_v<Value2, Value>> = 0>               \
    ENOKI_INLINE KMaskBase(                                                    \
        const StaticArrayBase<Value2, Size, IsMask2, Derived2> &a,             \
        detail::reinterpret_flag)

#define ENOKI_REINTERPRET_KMASK_SIZE(Value, Size)                              \
    template <typename Value2, typename Derived2, bool IsMask2,                \
              enable_if_t<detail::is_same_v<Value2, Value>> = 0>               \
    ENOKI_INLINE KMaskBase(                                                    \
        const StaticArrayBase<Value2, Size, IsMask2, Derived2> &a,             \
        detail::reinterpret_flag)

template <typename Value_, size_t Size_> struct KMask;

template <typename Value_, size_t Size_, typename Derived_>
struct KMaskBase : StaticArrayBase<Value_, Size_, true, Derived_> {
    using Register = std::conditional_t<(Size_ > 8), __mmask16, __mmask8>;
    using Derived = Derived_;
    using Base = StaticArrayBase<Value_, Size_, true, Derived_>;
    using Base::Size;
    using Base::derived;
    static constexpr bool IsNative = true;
    static constexpr bool IsKMask = true;
    static constexpr Register BitMask = Register((1 << Size_) - 1);

    ENOKI_ARRAY_DEFAULTS(KMaskBase)

#if defined(NDEBUG)
    KMaskBase() = default;
#else
    KMaskBase() : k(BitMask) { }
#endif

    template <typename Array, enable_if_t<std::is_same_v<Register, typename Array::Derived::Register>> = 0>
    ENOKI_INLINE KMaskBase(const Array &other, detail::reinterpret_flag) : k(other.derived().k) { }

    template <typename T, enable_if_t<std::is_same_v<bool, T> || std::is_same_v<int, T>> = 0>
    ENOKI_INLINE KMaskBase(const T &b, detail::reinterpret_flag)
        : k(bool(b) ? BitMask : Register(0)) { }

    ENOKI_REINTERPRET_KMASK(bool) {
        __m128i value;
        if constexpr (Size == 16)
            value = _mm_loadu_si128((__m128i *) a.derived().data());
        else if constexpr (Size == 8)
            value = _mm_loadl_epi64((const __m128i *) a.derived().data());
        else if constexpr (Size == 4 || Size == 3)
            value = _mm_cvtsi32_si128(*((const int *) a.derived().data()));
        else if constexpr (Size == 2)
            value = _mm_cvtsi32_si128((int) *((const short *) a.derived().data()));
        else
            static_assert(detail::false_v<Value2>, "Unsupported number of elements");

#if defined(ENOKI_X86_AVX512VL) && defined(ENOKI_X86_AVX512BW)
        k = (Register) _mm_test_epi8_mask(value, _mm_set1_epi8((char) 0xFF));
#else
        k = (Register) _mm512_test_epi32_mask(_mm512_cvtepi8_epi32(value),
                                              _mm512_set1_epi8((char) 0xFF));
#endif
    }

#if !defined(ENOKI_X86_AVX512VL)
    ENOKI_REINTERPRET_KMASK_SIZE(float, 8)    : k((Register) _mm256_movemask_ps(a.derived().m)) { }
    ENOKI_REINTERPRET_KMASK_SIZE(int32_t, 8)  : k((Register) _mm256_movemask_ps(_mm256_castsi256_ps(a.derived().m))) { }
    ENOKI_REINTERPRET_KMASK_SIZE(uint32_t, 8) : k((Register) _mm256_movemask_ps(_mm256_castsi256_ps(a.derived().m))) { }
#endif

    ENOKI_REINTERPRET_KMASK_SIZE(double, 16)   { k = _mm512_kunpackb(high(a).k, low(a).k); }
    ENOKI_REINTERPRET_KMASK_SIZE(int64_t, 16)  { k = _mm512_kunpackb(high(a).k, low(a).k); }
    ENOKI_REINTERPRET_KMASK_SIZE(uint64_t, 16) { k = _mm512_kunpackb(high(a).k, low(a).k); }

    template <typename T> ENOKI_INLINE static Derived from_k(const T &k) {
        Derived result;
        result.k = (Register) k;
        return result;
    }

    ENOKI_INLINE Derived eq_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kxnor(k, a.k));
        else
            return Derived::from_k(~(k ^ a.k));
    }

    ENOKI_INLINE Derived neq_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kxor(k, a.k));
        else
            return Derived::from_k(k ^ a.k);
    }

    ENOKI_INLINE Derived or_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kor(k, a.k));
        else
            return Derived::from_k(k | a.k);
    }

    ENOKI_INLINE Derived and_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kand(k, a.k));
        else
            return Derived::from_k(k & a.k);
    }

    ENOKI_INLINE Derived andnot_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kandn(a.k, k));
        else
            return Derived::from_k(k & ~a.k);
    }

    ENOKI_INLINE Derived xor_(const Derived &a) const {
        if constexpr (Size_ == 16) /* Use intrinsic if possible */
            return Derived::from_k(_mm512_kxor(k, a.k));
        else
            return Derived::from_k(k ^ a.k);
    }

    ENOKI_INLINE Derived not_() const {
        if constexpr (Size_ == 16)
            return Derived::from_k(_mm512_knot(k));
        else
            return Derived::from_k(~k);
    }

    static ENOKI_INLINE Derived select_(const Derived &m, const Derived &t, const Derived &f) {
        if constexpr (Size_ == 16)
            return Derived::from_k(_mm512_kor(_mm512_kand (m.k, t.k),
                                              _mm512_kandn(m.k, f.k)));
        else
            return Derived::from_k((m.k & t.k) | (~m.k & f.k));
    }

    ENOKI_INLINE bool all_() const {
        if constexpr (Size_ == 16)
            return _mm512_kortestc(k, k);
        else if constexpr (Size_ == 8)
            return k == BitMask;
        else
            return (k & BitMask) == BitMask;
    }

    ENOKI_INLINE bool any_() const {
        if constexpr (Size_ == 16)
            return !_mm512_kortestz(k, k);
        else if constexpr (Size_ == 8)
            return k != 0;
        else
            return (k & BitMask) != 0;
    }

    ENOKI_INLINE uint32_t bitmask_() const {
        if constexpr (Size_ == 8 || Size_ == 16)
            return (uint32_t) k;
        else
            return (uint32_t) (k & BitMask);
    }

    ENOKI_INLINE size_t count_() const {
        return (size_t) _mm_popcnt_u32(bitmask_());
    }

    ENOKI_INLINE bool bit_(size_t i) const {
        return (k & ((Register) 1 << i)) != 0;
    }

    ENOKI_INLINE void set_bit_(size_t i, bool value) {
        k = (Register) (k ^ ((-value ^ k) & ((Register) 1 << i)));
    }

    ENOKI_INLINE auto coeff(size_t i) const {
        return MaskBit<const Derived &>(derived(), i);
    }

    ENOKI_INLINE auto coeff(size_t i) {
        return MaskBit<Derived &>(derived(), i);
    }

    static Derived zero_() { return Derived::from_k(0); }

    template <typename Return = KMask<Value_, Size_ / 2>>
    ENOKI_INLINE Return low_() const {
        if constexpr (Size == 16)
            return Return::from_k(__mmask8(k));
        else
            return Return::from_k(Return::BitMask & k);
    }

    template <typename Return = KMask<Value_, Size_ / 2>>
    ENOKI_INLINE Return high_()  const {
        return Return::from_k(k >> (Size_ / 2));
    }

    ENOKI_INLINE void store_(void *ptr) const {
        store_unaligned_(ptr);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        memcpy(ptr, &k, sizeof(Register));
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return load_unaligned_(ptr);
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        Derived result;
        memcpy(&result.k, ptr, sizeof(Register));
        return result;
    }

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index_, const Mask &mask) {
        using UInt32 = Array<uint32_t, Size>;

        UInt32 index_32 = UInt32(index_),
               index, offset;

        if (Size == 2) {
            index  = sr<1>(index_32);
            offset = Index(1) << (index_32 & (uint32_t) 0x1);
        } else if (Size == 4) {
            index  = sr<2>(index_32);
            offset = Index(1) << (index_32 & (uint32_t) 0x3);
        } else {
            index  = sr<3>(index_32);
            offset = Index(1) << (index_32 & (uint32_t) 0x7);
        }

#if 0
        const uint8_t *in = (const uint8_t *) ptr;
        Register bit = 1, accum = 0;
        for (size_t i = 0; i < Size; ++i) {
            if ((bool) mask.coeff(i) && (in[index.coeff(i)] & offset.coeff(i)) != 0)
                accum |= bit;
            bit <<= 1;
        }
        return Derived::from_k(accum);
#else
        return Derived(neq(gather<UInt32, 1>(ptr, index, mask) & offset, (uint32_t) 0));
#endif
    }

    template <typename Array, enable_if_t<std::is_same_v<Register, typename Array::Derived::Register>> = 0>
    ENOKI_INLINE Derived& operator=(const Array &other) {
        k = other.derived().k;
        return derived();
    }

    template <typename T, enable_if_t<std::is_same_v<bool, T> || std::is_same_v<int, T>> = 0>
    ENOKI_INLINE Derived& operator=(const T &b) {
        k = bool(b) ? BitMask : Register(0);
        return derived();
    }

    Register k;
};

template <typename Value_, size_t Size_>
struct KMask : KMaskBase<Value_, Size_, KMask<Value_, Size_>> {
    using Base = KMaskBase<Value_, Size_, KMask<Value_, Size_>>;

    ENOKI_ARRAY_IMPORT(Base, KMask)
};

#define ENOKI_DECLARE_KMASK(Type, Size, Derived, SFINAE)                       \
    struct StaticArrayImpl<Type, Size, true, Derived, SFINAE>                  \
        : KMaskBase<Type, Size, Derived> {                                     \
        using Base = KMaskBase<Type, Size, Derived>;                           \
        ENOKI_ARRAY_DEFAULTS(StaticArrayImpl)                                  \
        using Base::Base;                                                      \
        using Base::operator=;                                                 \
    };

NAMESPACE_END(enoki)
