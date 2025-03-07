/*
    enoki/array_call.h -- Enoki arrays of pointers, support for
    array (virtual) method calls

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

NAMESPACE_BEGIN(enoki)

template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayImpl<Value_, Size_, IsMask_, Derived_,
                       enable_if_t<detail::array_config<Value_, Size_>::use_enum_impl>>
    : StaticArrayImpl<std::underlying_type_t<Value_>, Size_, IsMask_, Derived_> {

    using UnderlyingType = std::underlying_type_t<Value_>;

    using Base = StaticArrayImpl<UnderlyingType, Size_, IsMask_, Derived_>;

    ENOKI_ARRAY_DEFAULTS(StaticArrayImpl)
    using Base::derived;

    using Value = std::conditional_t<IsMask_, typename Base::Value, Value_>;
    using Scalar = std::conditional_t<IsMask_, typename Base::Scalar, Value_>;

    StaticArrayImpl() = default;
    StaticArrayImpl(Value value) : Base(UnderlyingType(value)) { }

    template <typename T, enable_if_t<!std::is_enum_v<T>> = 0>
    StaticArrayImpl(const T &b) : Base(b) { }

    template <typename T, enable_if_t<!is_array_v<T>> = 0>
    StaticArrayImpl(const T &v1, const T &v2) : Base(v1, v2) { }

    template <typename T>
    StaticArrayImpl(const T &b, detail::reinterpret_flag)
        : Base(b, detail::reinterpret_flag()) { }

    template <typename T1, typename T2, typename T = StaticArrayImpl, enable_if_t<
              array_depth_v<T1> == array_depth_v<T> && array_size_v<T1> == Base::Size1 &&
              array_depth_v<T2> == array_depth_v<T> && array_size_v<T2> == Base::Size2 &&
              Base::Size2 != 0> = 0>
    StaticArrayImpl(const T1 &a1, const T2 &a2)
        : Base(a1, a2) { }

    ENOKI_INLINE decltype(auto) coeff(size_t i) const {
        using Coeff = decltype(Base::coeff(i));
        if constexpr (std::is_same_v<Coeff, const typename Base::Value &>)
            return (const Value &) Base::coeff(i);
        else
            return Base::coeff(i);
    }

    ENOKI_INLINE decltype(auto) coeff(size_t i) {
        using Coeff = decltype(Base::coeff(i));
        if constexpr (std::is_same_v<Coeff, typename Base::Value &>)
            return (Value &) Base::coeff(i);
        else
            return Base::coeff(i);
    }

    template <typename T, typename Mask>
    ENOKI_INLINE size_t compress_(T *&ptr, const Mask &mask) const {
        return Base::compress_((UnderlyingType *&) ptr, mask);
    }

    template <typename T> Derived_& operator=(T&& t) {
        ENOKI_MARK_USED(t);
        if constexpr (std::is_same_v<T, std::nullptr_t>)
            return (Derived_ &) Base::operator=(UnderlyingType(0));
        else if constexpr (std::is_convertible_v<T, Value>)
            return (Derived_ &) Base::operator=(UnderlyingType(t));
        else
            return (Derived_ &) Base::operator=(std::forward<T>(t));
    }
};

NAMESPACE_END(enoki)
