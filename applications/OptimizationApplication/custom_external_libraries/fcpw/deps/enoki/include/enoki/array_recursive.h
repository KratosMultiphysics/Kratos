/*
    enoki/array_recursive.h -- Template specialization that recursively
    instantiates Array instances with smaller sizes when the requested packet
    float array size is not directly supported by the processor's SIMD
    instructions

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_generic.h>

NAMESPACE_BEGIN(enoki)

template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayImpl<Value_, Size_, IsMask_, Derived_,
                       enable_if_t<detail::array_config<Value_, Size_>::use_recursive_impl>>
    : StaticArrayBase<Value_, Size_, IsMask_, Derived_> {

    using Base = StaticArrayBase<Value_, Size_, IsMask_, Derived_>;

    ENOKI_ARRAY_IMPORT_BASIC(Base, StaticArrayImpl)

    using typename Base::Array1;
    using typename Base::Array2;
    using Base::Size1;
    using Base::Size2;
    using Ref = const Derived &;
    static constexpr bool IsRecursive = true;

    StaticArrayImpl() = default;

    /// Initialize all entries with a constant
    ENOKI_INLINE StaticArrayImpl(const Value &value) : a1(value), a2(value) { }

    /// Initialize from a list of component values
    template <typename... Ts, enable_if_t<sizeof...(Ts) == Size &&
        std::conjunction_v<detail::is_constructible<Value, Ts>...>> = 0>
    ENOKI_INLINE StaticArrayImpl(Ts... args) {
        alignas(alignof(Array1)) Value storage[Size] = { (Value) args... };
        a1 = load<Array1>(storage);
        a2 = load<Array2>(storage + Size1);
    }

    /// Construct from the two sub-components
    template <typename T1, typename T2,
              enable_if_t<T1::Size == Size1 && T2::Size == Size2> = 0>
    ENOKI_INLINE StaticArrayImpl(const T1 &a1, const T2 &a2)
        : a1(a1), a2(a2) { }

    /// Cast another array
    template <size_t Size2, typename Value2,
              typename Derived2, enable_if_t<Derived2::Size == Size_> = 0>
    ENOKI_INLINE StaticArrayImpl(
        const StaticArrayBase<Value2, Size2, IsMask_, Derived2> &a)
        : a1(low(a)), a2(high(a)) { }

    /// Reinterpret another array
    template <typename Value2, size_t Size2,
              bool IsMask2, typename Derived2, enable_if_t<Derived2::Size == Size_> = 0>
    ENOKI_INLINE StaticArrayImpl(
        const StaticArrayBase<Value2, Size2, IsMask2, Derived2> &a,
        detail::reinterpret_flag)
        : a1(low (a), detail::reinterpret_flag()),
          a2(high(a), detail::reinterpret_flag()) { }

    /// Reinterpret another array (masks)
    template <bool M = IsMask_, enable_if_t<M> = 0>
    ENOKI_INLINE StaticArrayImpl(bool b, detail::reinterpret_flag)
        : a1(b, detail::reinterpret_flag()),
          a2(b, detail::reinterpret_flag()) { }

    template <bool M = IsMask_, enable_if_t<!M> = 0>
    ENOKI_INLINE StaticArrayImpl &operator=(Value_ v) {
        *this = StaticArrayImpl(v);
        return *this;
    }

    template <bool M = IsMask_, enable_if_t<M> = 0>
    ENOKI_INLINE StaticArrayImpl &operator=(bool v) {
        *this = StaticArrayImpl(v, detail::reinterpret_flag());
        return *this;
    }

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return Derived(a1 + a.a1, a2 + a.a2); }
    ENOKI_INLINE Derived sub_(Ref a) const { return Derived(a1 - a.a1, a2 - a.a2); }
    ENOKI_INLINE Derived mul_(Ref a) const { return Derived(a1 * a.a1, a2 * a.a2); }
    ENOKI_INLINE Derived div_(Ref a) const { return Derived(a1 / a.a1, a2 / a.a2); }
    ENOKI_INLINE Derived mod_(Ref a) const { return Derived(a1 % a.a1, a2 % a.a2); }

    ENOKI_INLINE Derived mulhi_(Ref a) const {
        return Derived(mulhi(a1, a.a1), mulhi(a2, a.a2));
    }

    ENOKI_INLINE Derived fmod_(Ref a) const {
        return Derived(fmod(a1, a.a1), fmod(a2, a.a2));
    }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>(a1 <  a.a1, a2 <  a.a2); }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>(a1 >  a.a1, a2 >  a.a2); }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>(a1 <= a.a1, a2 <= a.a2); }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>(a1 >= a.a1, a2 >= a.a2); }
    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>(eq(a1, a.a1), eq(a2, a.a2)); }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>(neq(a1, a.a1), neq(a2, a.a2)); }

    ENOKI_INLINE Derived min_(Ref a) const { return Derived(min(a1, a.a1), min(a2, a.a2)); }
    ENOKI_INLINE Derived max_(Ref a) const { return Derived(max(a1, a.a1), max(a2, a.a2)); }
    ENOKI_INLINE Derived abs_() const { return Derived(abs(a1), abs(a2)); }
    ENOKI_INLINE Derived ceil_() const { return Derived(ceil(a1), ceil(a2)); }
    ENOKI_INLINE Derived floor_() const { return Derived(floor(a1), floor(a2)); }
    ENOKI_INLINE Derived sqrt_() const { return Derived(sqrt(a1), sqrt(a2)); }
    ENOKI_INLINE Derived round_() const { return Derived(round(a1), round(a2)); }
    ENOKI_INLINE Derived trunc_() const { return Derived(trunc(a1), trunc(a2)); }
    ENOKI_INLINE Derived rcp_() const { return Derived(rcp(a1), rcp(a2)); }
    ENOKI_INLINE Derived rsqrt_() const { return Derived(rsqrt(a1), rsqrt(a2)); }
    ENOKI_INLINE Derived not_() const { return Derived(~a1, ~a2); }
    ENOKI_INLINE Derived neg_() const { return Derived(-a1, -a2); }

    ENOKI_INLINE Derived fmadd_(Ref b, Ref c) const {
        return Derived(fmadd(a1, b.a1, c.a1), fmadd(a2, b.a2, c.a2));
    }

    ENOKI_INLINE Derived fnmadd_(Ref b, Ref c) const {
        return Derived(fnmadd(a1, b.a1, c.a1), fnmadd(a2, b.a2, c.a2));
    }

    ENOKI_INLINE Derived fmsub_(Ref b, Ref c) const {
        return Derived(fmsub(a1, b.a1, c.a1), fmsub(a2, b.a2, c.a2));
    }

    ENOKI_INLINE Derived fnmsub_(Ref b, Ref c) const {
        return Derived(fnmsub(a1, b.a1, c.a1), fnmsub(a2, b.a2, c.a2));
    }

    ENOKI_INLINE Derived fmaddsub_(Ref b, Ref c) const {
        return Derived(fmaddsub(a1, b.a1, c.a1), fmaddsub(a2, b.a2, c.a2));
    }

    ENOKI_INLINE Derived fmsubadd_(Ref b, Ref c) const {
        return Derived(fmsubadd(a1, b.a1, c.a1), fmsubadd(a2, b.a2, c.a2));
    }

    template <typename T> ENOKI_INLINE Derived or_(const T &a) const {
        return Derived(a1 | low(a), a2 | high(a));
    }

    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const {
        return Derived(andnot(a1, low(a)), andnot(a2, high(a)));
    }

    template <typename T> ENOKI_INLINE Derived and_(const T &a) const {
        return Derived(a1 & low(a), a2 & high(a));
    }

    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const {
        return Derived(a1 ^ low(a), a2 ^ high(a));
    }

    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        return Derived(sl<Imm>(a1), sl<Imm>(a2));
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return Derived(a1 << k, a2 << k);
    }

    ENOKI_INLINE Derived sl_(Ref a) const {
        return Derived(a1 << a.a1, a2 << a.a2);
    }

    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        return Derived(sr<Imm>(a1), sr<Imm>(a2));
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        return Derived(a1 >> k, a2 >> k);
    }

    ENOKI_INLINE Derived sr_(Ref a) const {
        return Derived(a1 >> a.a1, a2 >> a.a2);
    }

    template <size_t Imm> ENOKI_INLINE Derived rol_() const {
        return Derived(rol<Imm>(a1), rol<Imm>(a2));
    }

    template <size_t Imm> ENOKI_INLINE Derived ror_() const {
        return Derived(ror<Imm>(a1), ror<Imm>(a2));
    }

    ENOKI_INLINE Derived rol_(Ref arg) const {
        return Derived(rol(a1, arg.a1), rol(a2, arg.a2));
    }

    ENOKI_INLINE Derived ror_(Ref arg) const {
        return Derived(ror(a1, arg.a1), ror(a2, arg.a2));
    }

    template <typename Mask>
    static ENOKI_INLINE Derived select_(const Mask &m, Ref t, Ref f) {
        return Derived(select(m.a1, t.a1, f.a1),
                       select(m.a2, t.a2, f.a2));
    }

    template <size_t Imm> ENOKI_INLINE Derived ror_array_() const {
        if constexpr (Size1 == Size2) {
            static_assert(
                Imm <= Size1 && Imm <= Size2,
                "ror_array(): Refusing to rotate a recursively defined array by an "
                "amount that is larger than the recursive array sizes.");
            const mask_t<Array1> mask = arange<Array1>() >= Scalar(Imm);

            Array1 a1_r = ror_array<Imm>(a1);
            Array2 a2_r = ror_array<Imm>(a2);

            return Derived(
                select(mask, a1_r, a2_r),
                select(mask, a2_r, a1_r)
            );
        } else {
            return Base::template ror_array_<Imm>();
        }
    }

    template <size_t Imm> ENOKI_INLINE Derived rol_array_() const {
        if constexpr (Size1 == Size2) {
            static_assert(
                Imm <= Size1 && Imm <= Size2,
                "rol_array(): Refusing to rotate a recursively defined array "
                "by an amount that is larger than the recursive array sizes.");
            const mask_t<Array1> mask = arange<Array1>() < Scalar(Size1 - Imm);

            Array1 a1_r = rol_array<Imm>(a1);
            Array2 a2_r = rol_array<Imm>(a2);

            return Derived(
                select(mask, a1_r, a2_r),
                select(mask, a2_r, a1_r)
            );
        } else {
            return Base::template rol_array_<Imm>();
        }
    }

    Derived ldexp_(Ref a) const {
        return Derived(ldexp(a1, a.a1), ldexp(a2, a.a2));
    }

    std::pair<Derived, Derived> frexp_() const {
        auto r1 = frexp(a1);
        auto r2 = frexp(a2);
        return std::make_pair<Derived, Derived>(
            Derived(r1.first, r2.first),
            Derived(r1.second, r2.second)
        );
    }


    template <typename T>
    ENOKI_INLINE auto ceil2int_() const {
        return T(ceil2int<typename T::Array1>(a1),
                 ceil2int<typename T::Array2>(a2));
    }

    template <typename T>
    ENOKI_INLINE auto floor2int_() const {
        return T(floor2int<typename T::Array1>(a1),
                 floor2int<typename T::Array2>(a2));
    }

    Derived lzcnt_() const  { return Derived(lzcnt(a1),  lzcnt(a2));  }
    Derived tzcnt_() const  { return Derived(tzcnt(a1),  tzcnt(a2));  }
    Derived popcnt_() const { return Derived(popcnt(a1), popcnt(a2)); }

    template<size_t... Is, size_t ... Is2>
    static constexpr auto split_(std::index_sequence<Is...>,
                                 std::index_sequence<Is2...>) {
        constexpr std::array<size_t, sizeof...(Is)> out { Is ... };
        return std::make_pair(std::index_sequence<out[Is2]...>(),
                              std::index_sequence<out[Is2 + Size1]...>());
    }

    template <size_t... Indices> ENOKI_INLINE Derived shuffle_() const {
        if constexpr (Size1 != Size2) {
            return Base::template shuffle_<Indices...>();
        } else {
            constexpr auto indices = split_(std::index_sequence<Indices...>(),
                                            std::make_index_sequence<Size1>());
            return shuffle_impl_(indices.first, indices.second);
        }
    }

    template <size_t... Indices1, typename T= size_t, size_t... Indices2>
    ENOKI_INLINE Derived shuffle_impl_(std::index_sequence<Indices1...>,
                                       std::index_sequence<Indices2...>) const {
        using Int = int_array_t<Array1>;
        Array1 a1l = a1.template shuffle_<(size_t) std::min(Size1 - 1, Indices1)...>(),
               a1h = a2.template shuffle_<(size_t) std::max((ssize_t) 0, (ssize_t) Indices1 - (ssize_t) Size1)...>(),
               a1f = select(Int(Indices1...) < Int(Size1), a1l, a1h);

        Array2 a2l = a1.template shuffle_<std::min(Size1 - 1, Indices2)...>(),
               a2h = a2.template shuffle_<(size_t) std::max((ssize_t) 0, (ssize_t) Indices2 - (ssize_t) Size1)...>(),
               a2f = select(Int(Indices2...) < Int(Size1), a2l, a2h);

        return Derived(a1f, a2f);
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        if constexpr (Size1 != Size2) {
            return Base::shuffle_(index);
        } else {
            auto il = low(index), ih = high(index);

            decltype(il) size = scalar_t<Index>(Size1);

            Array1 a1l = a1.shuffle_(il),
                   a1h = a2.shuffle_(il - size),
                   a1f = select(il < size, a1l, a1h);

            Array2 a2l = a1.shuffle_(ih),
                   a2h = a2.shuffle_(ih - size),
                   a2f = select(ih < size, a2l, a2h);

            return Derived(a1f, a2f);
        }
    }

    #define ENOKI_MASKED_OPERATOR(name)                                        \
        template <typename Mask>                                               \
        ENOKI_INLINE void m##name##_(Ref value, const Mask &mask) {            \
            a1.m##name##_(low(value), low(mask));                              \
            a2.m##name##_(high(value), high(mask));                            \
        }

    ENOKI_MASKED_OPERATOR(assign)
    ENOKI_MASKED_OPERATOR(add)
    ENOKI_MASKED_OPERATOR(sub)
    ENOKI_MASKED_OPERATOR(mul)
    ENOKI_MASKED_OPERATOR(div)
    ENOKI_MASKED_OPERATOR(and)
    ENOKI_MASKED_OPERATOR(or)
    ENOKI_MASKED_OPERATOR(xor)

    #undef ENOKI_MASKED_OPERATOR

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_() const {
        if constexpr (Size1 == Size2)
            return hsum(a1 + a2);
        else
            return hsum(a1) + hsum(a2);
    }

    ENOKI_INLINE Value hprod_() const {
        if constexpr (Size1 == Size2)
            return hprod(a1 * a2);
        else
            return hprod(a1) * hprod(a2);
    }

    ENOKI_INLINE Value hmin_() const {
        if constexpr (Size1 == Size2)
            return hmin(min(a1, a2));
        else
            return min(hmin(a1), hmin(a2));
    }

    ENOKI_INLINE Value hmax_() const {
        if constexpr (Size1 == Size2)
            return hmax(max(a1, a2));
        else
            return max(hmax(a1), hmax(a2));
    }

    ENOKI_INLINE Value dot_(Ref a) const {
        if constexpr (Size1 == Size2)
            return hsum(fmadd(a1, a.a1, a2 * a.a2));
        else
            return dot(a1, a.a1) + dot(a2, a.a2);
    }

    ENOKI_INLINE bool all_() const {
        if constexpr (Size1 == Size2)
            return all(a1 & a2);
        else
            return all(a1) && all(a2);
    }

    ENOKI_INLINE bool any_() const {
        if constexpr (Size1 == Size2)
            return any(a1 | a2);
        else
            return any(a1) || any(a2);
    }

    ENOKI_INLINE size_t count_() const { return count(a1) + count(a2); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    ENOKI_INLINE void store_(void *mem) const {
        store((uint8_t *) mem, a1);
        store((uint8_t *) mem + sizeof(Array1), a2);
    }

    template <typename Mask>
    ENOKI_INLINE void store_(void *mem, const Mask &mask) const {
        store((uint8_t *) mem, a1, low(mask));
        store((uint8_t *) mem + sizeof(Array1), a2, high(mask));
    }

    ENOKI_INLINE void store_unaligned_(void *mem) const {
        store_unaligned((uint8_t *) mem, a1);
        store_unaligned((uint8_t *) mem + sizeof(Array1), a2);
    }

    template <typename Mask>
    ENOKI_INLINE void store_unaligned_(void *mem, const Mask &mask) const {
        store_unaligned((uint8_t *) mem, a1, low(mask));
        store_unaligned((uint8_t *) mem + sizeof(Array1), a2, high(mask));
    }

    static ENOKI_INLINE Derived load_(const void *mem) {
        return Derived(
            load<Array1>((uint8_t *) mem),
            load<Array2>((uint8_t *) mem + sizeof(Array1))
        );
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *mem, const Mask &mask) {
        return Derived(
            load<Array1>((uint8_t *) mem, low(mask)),
            load<Array2>((uint8_t *) mem + sizeof(Array1), high(mask))
        );
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *a) {
        return Derived(
            load_unaligned<Array1>((uint8_t *) a),
            load_unaligned<Array2>((uint8_t *) a + sizeof(Array1))
        );
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *a, const Mask &mask) {
        return Derived(
            load_unaligned<Array1>((uint8_t *) a, low(mask)),
            load_unaligned<Array2>((uint8_t *) a + sizeof(Array1), high(mask))
        );
    }

    static ENOKI_INLINE Derived zero_() {
        return Derived(zero<Array1>(), zero<Array2>());
    }

    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *ptr, const Index &index, const Mask &mask) {
        prefetch<Array1, Write, Level, Stride>(ptr, low(index), low(mask));
        prefetch<Array2, Write, Level, Stride>(ptr, high(index), high(mask));
    }

    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *ptr, const Index &index, const Mask &mask) {
        return Derived(
            gather<Array1, Stride>(ptr, low(index), low(mask)),
            gather<Array2, Stride>(ptr, high(index), high(mask))
        );
    }

    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *ptr, const Index &index, const Mask &mask) const {
        scatter<Stride>(ptr, a1, low(index), low(mask));
        scatter<Stride>(ptr, a2, high(index), high(mask));
    }

    template <size_t Stride, typename Index, typename Func, typename... Args, typename Mask>
    static ENOKI_INLINE void transform_(void *ptr, const Index &index, const Mask &,
                                        const Func &func, const Args &... args) {
        transform<Array1, Stride>(ptr, low(index),  func, low(args)...);
        transform<Array2, Stride>(ptr, high(index), func, high(args)...);
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        if constexpr (Size1 == Size2) {
            return extract(select(low(mask), a1, a2), low(mask) | high(mask));
        } else {
            if (ENOKI_LIKELY(any(low(mask))))
                return extract(a1, low(mask));
            else
                return extract(a2, high(mask));
        }
    }

    template <typename T, typename Mask>
    ENOKI_INLINE size_t compress_(T *&ptr, const Mask &mask) const {
        size_t r0 = compress(ptr, a1, low(mask));
        size_t r1 = compress(ptr, a2, high(mask));
        return r0 + r1;
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Component access
    // -----------------------------------------------------------------------

    ENOKI_INLINE const Array1& low_()  const { return a1; }
    ENOKI_INLINE const Array2& high_() const { return a2; }

    ENOKI_INLINE decltype(auto) coeff(size_t i) const {
        if constexpr (Size1 == Size2)
            return ((i < Size1) ? a1 : a2).coeff(i % Size1);
        else
            return (i < Size1) ? a1.coeff(i) : a2.coeff(i - Size1);
    }

    ENOKI_INLINE decltype(auto) coeff(size_t i) {
        if constexpr (Size1 == Size2)
            return ((i < Size1) ? a1 : a2).coeff(i % Size1);
        else
            return (i < Size1) ? a1.coeff(i) : a2.coeff(i - Size1);
    }

    //! @}
    // -----------------------------------------------------------------------

    Array1 a1;
    Array2 a2;
};

NAMESPACE_END(enoki)
