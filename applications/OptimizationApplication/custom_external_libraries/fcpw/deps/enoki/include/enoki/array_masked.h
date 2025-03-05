/*
    enoki/array_masked.h -- Helper classes for masked assignments and
    in-place operators

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

NAMESPACE_BEGIN(enoki)

// -----------------------------------------------------------------------
//! @{ \name Masked array helper classes
// -----------------------------------------------------------------------

NAMESPACE_BEGIN(detail)

template <typename T> struct MaskedValue {
    MaskedValue(T &d, bool m) : d(d), m(m) { }

    template <typename T2> ENOKI_INLINE void operator =(const T2 &value) { if (m) d = value; }
    template <typename T2> ENOKI_INLINE void operator+=(const T2 &value) { if (m) d += value; }
    template <typename T2> ENOKI_INLINE void operator-=(const T2 &value) { if (m) d -= value; }
    template <typename T2> ENOKI_INLINE void operator*=(const T2 &value) { if (m) d *= value; }
    template <typename T2> ENOKI_INLINE void operator/=(const T2 &value) { if (m) d /= value; }
    template <typename T2> ENOKI_INLINE void operator|=(const T2 &value) { if (m) d |= value; }
    template <typename T2> ENOKI_INLINE void operator&=(const T2 &value) { if (m) d &= value; }
    template <typename T2> ENOKI_INLINE void operator^=(const T2 &value) { if (m) d ^= value; }

    T &d;
    bool m;
};

template <typename T> struct MaskedArray : ArrayBase<value_t<T>, MaskedArray<T>> {
    using Mask     = mask_t<T>;
    using Scalar   = MaskedValue<scalar_t<T>>;
    using MaskType = MaskedArray<Mask>;
    using Value    = std::conditional_t<is_scalar_v<value_t<T>>,
                                     MaskedValue<value_t<T>>,
                                     MaskedArray<value_t<T>>>;
    using UnderlyingValue = value_t<T>;
    static constexpr size_t Size = array_size_v<T>;
    static constexpr bool IsMaskedArray = true;

    MaskedArray(T &d, const Mask &m) : d(d), m(m) { }

    template <typename T2> ENOKI_INLINE void operator =(const T2 &value) { d.massign_(value, m); }
    template <typename T2> ENOKI_INLINE void operator+=(const T2 &value) { d.madd_(value, m); }
    template <typename T2> ENOKI_INLINE void operator-=(const T2 &value) { d.msub_(value, m); }
    template <typename T2> ENOKI_INLINE void operator*=(const T2 &value) { d.mmul_(value, m); }
    template <typename T2> ENOKI_INLINE void operator/=(const T2 &value) { d.mdiv_(value, m); }
    template <typename T2> ENOKI_INLINE void operator|=(const T2 &value) { d.mor_(value, m); }
    template <typename T2> ENOKI_INLINE void operator&=(const T2 &value) { d.mand_(value, m); }
    template <typename T2> ENOKI_INLINE void operator^=(const T2 &value) { d.mxor_(value, m); }

    /// Type alias for a similar-shaped array over a different type
    template <typename T2> using ReplaceValue = MaskedArray<typename T::template ReplaceValue<T2>>;

    T &d;
    Mask m;
};

NAMESPACE_END(detail)

template <typename Value_, size_t Size_>
struct Array<detail::MaskedArray<Value_>, Size_>
    : detail::MaskedArray<Array<Value_, Size_>> {
    using Base = detail::MaskedArray<Array<Value_, Size_>>;
    using Base::Base;
    using Base::operator=;
    Array(const Base &b) : Base(b) { }
};

template <typename T, typename Mask>
ENOKI_INLINE auto masked(T &value, const Mask &mask) {
    if constexpr (std::is_same_v<Mask, bool>)
        return detail::MaskedValue<T>{ value, mask };
    else
        return struct_support_t<T>::masked(value, mask);
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
