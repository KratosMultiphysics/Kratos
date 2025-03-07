/*
    enoki/array.h -- Main header file for the Enoki array class and
    various template specializations

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable: 4146) // warning C4146: unary minus operator applied to unsigned type, result still unsigned
#  pragma warning(disable: 4554) // warning C4554: '>>': check operator precedence for possible error; use parentheses to clarify precedence
#  pragma warning(disable: 4702) // warning C4702: unreachable code
#  pragma warning(disable: 4522) // warning C4522: multiple assignment operators specified
#  pragma warning(disable: 4310) // warning C4310: cast truncates constant value
#  pragma warning(disable: 4127) // warning C4127: conditional expression is constant
#elif defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include <enoki/array_generic.h>

#include <enoki/array_math.h>

#if defined(ENOKI_ARM_NEON) || defined(ENOKI_X86_SSE42)
#  include <enoki/array_recursive.h>
#endif

#if defined(ENOKI_X86_AVX512F)
#  include <enoki/array_kmask.h>
#endif

#if defined(ENOKI_X86_SSE42)
#  include <enoki/array_sse42.h>
#endif

#if defined(ENOKI_X86_AVX)
#  include <enoki/array_avx.h>
#endif

#if defined(ENOKI_X86_AVX2)
#  include <enoki/array_avx2.h>
#endif

#if defined(ENOKI_X86_AVX512F)
#  include <enoki/array_avx512.h>
#endif

#if defined(ENOKI_ARM_NEON)
#  include <enoki/array_neon.h>
#endif

#include <enoki/array_idiv.h>
#include <enoki/array_call.h>
#include <enoki/array_enum.h>
#include <enoki/array_utils.h>
#include <enoki/array_macro.h>

#include <enoki/half.h>

NAMESPACE_BEGIN(enoki)

template <typename Value_, size_t Size_>
struct Array : StaticArrayImpl<Value_, Size_, false, Array<Value_, Size_>> {

    using Base = StaticArrayImpl<Value_, Size_, false, Array<Value_, Size_>>;

    using ArrayType = Array;
    using MaskType = Mask<Value_, Size_>;

    /// Type alias for creating a similar-shaped array over a different type
    template <typename T> using ReplaceValue = Array<T, Size_>;

    ENOKI_ARRAY_IMPORT(Base, Array)
};

template <typename Value_, size_t Size_>
struct Mask : StaticArrayImpl<Value_, Size_, true, Mask<Value_, Size_>> {

    using Base = StaticArrayImpl<Value_, Size_, true, Mask<Value_, Size_>>;

    using ArrayType = Array<Value_, Size_>;
    using MaskType = Mask;

    /// Type alias for creating a similar-shaped array over a different type
    template <typename T> using ReplaceValue = Mask<T, Size_>;

    Mask() = default;

    template <typename T> Mask(T &&value)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }

    template <typename T> Mask(T &&value, detail::reinterpret_flag)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }

    /// Construct from sub-arrays
    template <typename T1, typename T2, typename T = Mask, enable_if_t<
              array_depth_v<T1> == array_depth_v<T> && array_size_v<T1> == Base::Size1 &&
              array_depth_v<T2> == array_depth_v<T> && array_size_v<T2> == Base::Size2 &&
              Base::Size2 != 0> = 0>
    Mask(const T1 &a1, const T2 &a2)
        : Base(a1, a2) { }

    template <typename... Ts,
        enable_if_t<(sizeof...(Ts) == Base::Size || sizeof...(Ts) == Base::ActualSize) && Size_ != 1 &&
                    std::conjunction_v<detail::is_not_reinterpret_flag<Ts>...>> = 0>
    Mask(Ts&&... ts) : Base(std::forward<Ts>(ts)...) { }

    ENOKI_ARRAY_IMPORT_BASIC(Base, Mask)
    using Base::operator=;
};

template <typename Value_, size_t Size_>
struct Packet : StaticArrayImpl<Value_, Size_, false, Packet<Value_, Size_>> {

    using Base = StaticArrayImpl<Value_, Size_, false, Packet<Value_, Size_>>;

    using ArrayType = Packet;
    using MaskType = PacketMask<Value_, Size_>;

    static constexpr bool BroadcastPreferOuter = false;

    /// Type alias for creating a similar-shaped array over a different type
    template <typename T> using ReplaceValue = Packet<T, Size_>;

    ENOKI_ARRAY_IMPORT(Base, Packet)
};

template <typename Value_, size_t Size_>
struct PacketMask : StaticArrayImpl<Value_, Size_, true, PacketMask<Value_, Size_>> {

    using Base = StaticArrayImpl<Value_, Size_, true, PacketMask<Value_, Size_>>;

    static constexpr bool BroadcastPreferOuter = false;

    using ArrayType = Packet<Value_, Size_>;
    using MaskType = PacketMask;

    /// Type alias for creating a similar-shaped array over a different type
    template <typename T> using ReplaceValue = PacketMask<T, Size_>;

    PacketMask() = default;

    template <typename T> PacketMask(T &&value)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }

    template <typename T> PacketMask(T &&value, detail::reinterpret_flag)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }

    /// Construct from sub-arrays
    template <typename T1, typename T2, typename T = PacketMask, enable_if_t<
              array_depth_v<T1> == array_depth_v<T> && array_size_v<T1> == Base::Size1 &&
              array_depth_v<T2> == array_depth_v<T> && array_size_v<T2> == Base::Size2 &&
              Base::Size2 != 0> = 0>
    PacketMask(const T1 &a1, const T2 &a2)
        : Base(a1, a2) { }

    template <typename... Ts,
        enable_if_t<(sizeof...(Ts) == Base::Size || sizeof...(Ts) == Base::ActualSize) && Size_ != 1 &&
                    std::conjunction_v<detail::is_not_reinterpret_flag<Ts>...>> = 0>
    PacketMask(Ts&&... ts) : Base(std::forward<Ts>(ts)...) { }

    ENOKI_ARRAY_IMPORT_BASIC(Base, PacketMask)
    using Base::operator=;
};

NAMESPACE_END(enoki)

#if defined(_MSC_VER)
#  pragma warning(pop)
#elif defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop
#endif
