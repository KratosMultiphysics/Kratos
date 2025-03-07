/*
    enoki/morton.h -- Morton/Z-order curve encoding and decoding routines

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>
    Includes contributions by Sebastien Speierer

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>

#if defined(_MSC_VER)
#  pragma warning (push)
#  pragma warning (disable: 4310) // cast truncates constant value
#endif

NAMESPACE_BEGIN(enoki)
NAMESPACE_BEGIN(detail)

/// Generate bit masks for the functions \ref scatter_bits() and \ref gather_bits()
template <typename Value> constexpr Value morton_magic(size_t dim, size_t level) {
    size_t n_bits = sizeof(Value) * 8;
    size_t max_block_size = n_bits / dim;
    size_t block_size = std::min(size_t(1) << (level - 1), max_block_size);
    size_t count = 0;

    Value mask  = Value(1) << (n_bits - 1),
         value = Value(0);

    for (size_t i = 0; i < n_bits; ++i) {
        value >>= 1;

        if (count < max_block_size && (i / block_size) % dim == 0) {
            count++;
            value |= mask;
        }
    }

    return value;
}

/// Bit scatter function. \c Dimension defines the final distance between two output bits
template <size_t, typename Value, size_t Level, enable_if_t<Level == 0> = 0>
ENOKI_INLINE Value scatter_bits(Value x) { return x; }

template <size_t Dimension, typename Value,
          size_t Level = clog2i(sizeof(Value) * 8),
          enable_if_t<Level != 0 && (!(has_avx2 && has_x86_64) || !std::is_integral_v<Value>)> = 0>
ENOKI_INLINE Value scatter_bits(Value x) {
    using Scalar = scalar_t<Value>;

    constexpr Scalar magic = morton_magic<Scalar>(Dimension, Level);
    constexpr size_t shift_maybe = (1 << (Level - 1)) * (Dimension - 1);
    constexpr size_t shift = (shift_maybe < sizeof(Scalar) * 8) ? shift_maybe : 0;

    if constexpr (shift != 0)
        x |= sl<shift>(x);

    x &= magic;

    return scatter_bits<Dimension, Value, Level - 1>(x);
}

template <size_t, typename Value, size_t Level,
          enable_if_t<Level == 0> = 0>
ENOKI_INLINE Value gather_bits(Value x) { return x; }

/// Bit gather function. \c Dimension defines the final distance between two input bits
template <size_t Dimension, typename Value,
          size_t Level = clog2i(sizeof(Value) * 8),
          enable_if_t<Level != 0 && (!(has_avx2 && has_x86_64) || !std::is_integral_v<Value>)> = 0>
ENOKI_INLINE Value gather_bits(Value x) {
    using Scalar = scalar_t<Value>;

    constexpr size_t ilevel = clog2i(sizeof(Value) * 8) - Level + 1;
    constexpr Scalar magic = morton_magic<Scalar>(Dimension, ilevel);
    constexpr size_t shift_maybe = (1 << (ilevel - 1)) * (Dimension - 1);
    constexpr size_t shift = (shift_maybe < sizeof(Scalar) * 8) ? shift_maybe : 0;

    x &= magic;

    if constexpr (shift != 0)
        x |= sr<shift>(x);

    return gather_bits<Dimension, Value, Level - 1>(x);
}

#if defined(ENOKI_X86_AVX2) && defined(ENOKI_X86_64)
template <size_t Dimension, typename Value,
          enable_if_t<std::is_integral_v<Value>> = 0>
ENOKI_INLINE Value scatter_bits(Value x) {
    constexpr Value magic = morton_magic<Value>(Dimension, 1);
    if constexpr (sizeof(Value) <= 4)
        return Value(_pdep_u32((uint32_t) x, (uint32_t) magic));
    else
        return Value(_pdep_u64((uint64_t) x, (uint64_t) magic));
}

template <size_t Dimension, typename Value,
          enable_if_t<std::is_integral_v<Value>> = 0>
ENOKI_INLINE Value gather_bits(Value x) {
    constexpr Value magic = morton_magic<Value>(Dimension, 1);
    if constexpr (sizeof(Value) <= 4)
        return Value(_pext_u32((uint32_t) x, (uint32_t) magic));
    else
        return Value(_pext_u64((uint64_t) x, (uint64_t) magic));
}
#endif

template <typename Array, size_t Index,
          enable_if_t<Index == 0> = 0>
ENOKI_INLINE void morton_decode_helper(value_t<Array> value, Array &out) {
    out.coeff(0) = gather_bits<Array::Size>(value);
}

template <typename Array, size_t Index = array_size_v<Array> - 1,
          enable_if_t<Index != 0> = 0>
ENOKI_INLINE void morton_decode_helper(value_t<Array> value, Array &out) {
    out.coeff(Index) = gather_bits<Array::Size>(sr<Index>(value));
    morton_decode_helper<Array, Index - 1>(value, out);
}

NAMESPACE_END(detail)

/// Convert a N-dimensional integer array into the Morton/Z-order curve encoding
template <typename Array, size_t Index, typename Return = value_t<Array>,
          enable_if_t<Index == 0> = 0>
ENOKI_INLINE Return morton_encode(Array a) {
    return detail::scatter_bits<Array::Size>(a.coeff(0));
}

/// Convert a N-dimensional integer array into the Morton/Z-order curve encoding
template <typename Array, size_t Index = array_size_v<Array> - 1,
          typename Return = value_t<Array>, enable_if_t<Index != 0> = 0>
ENOKI_INLINE Return morton_encode(Array a) {
    static_assert(std::is_unsigned_v<scalar_t<Array>>, "morton_encode() requires unsigned arguments");
    return sl<Index>(detail::scatter_bits<Array::Size>(a.coeff(Index))) |
           morton_encode<Array, Index - 1>(a);
}

/// Convert Morton/Z-order curve encoding into a N-dimensional integer array
template <typename Array, typename Value = value_t<Array>>
ENOKI_INLINE Array morton_decode(Value value) {
    static_assert(std::is_unsigned_v<scalar_t<Array>>, "morton_decode() requires unsigned arguments");
    Array result;
    detail::morton_decode_helper(value, result);
    return result;
}

NAMESPACE_END(enoki)

#if defined(_MSC_VER)
#  pragma warning (pop)
#endif
