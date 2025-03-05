/*
    enoki/stl.h -- vectorization support for STL pairs, tuples, and arrays

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>

NAMESPACE_BEGIN(enoki)

template <typename Arg0, typename Arg1> struct struct_support<std::pair<Arg0, Arg1>> {
    static constexpr bool IsDynamic =
        enoki::is_dynamic_v<Arg0> || enoki::is_dynamic_v<Arg1>;
    using Dynamic = std::pair<enoki::make_dynamic_t<Arg0>, enoki::make_dynamic_t<Arg1>>;
    using Value = std::pair<Arg0, Arg1>;

    static ENOKI_INLINE size_t slices(const Value &value) {
        return enoki::slices(value.first);
    }

    static ENOKI_INLINE size_t packets(const Value &value) {
        return enoki::packets(value.first);
    }

    static ENOKI_INLINE void set_slices(Value &value, size_t size) {
        enoki::set_slices(value.first, size);
        enoki::set_slices(value.second, size);
    }

    template <typename T2>
    static ENOKI_INLINE auto packet(T2 &&value, size_t i) {
        return std::pair<decltype(enoki::packet(value.first, i)),
                         decltype(enoki::packet(value.second, i))>(
            enoki::packet(value.first, i), enoki::packet(value.second, i));
    }

    template <typename T2>
    static ENOKI_INLINE auto slice(T2 &&value, size_t i) {
        return std::pair<decltype(enoki::slice(value.first, i)),
                         decltype(enoki::slice(value.second, i))>(
            enoki::slice(value.first, i), enoki::slice(value.second, i));
    }

    template <typename T2>
    static ENOKI_INLINE auto slice_ptr(T2 &&value, size_t i) {
        return std::pair<decltype(enoki::slice_ptr(value.first, i)),
                         decltype(enoki::slice_ptr(value.second, i))>(
            enoki::slice_ptr(value.first, i), enoki::slice_ptr(value.second, i));
    }

    template <typename T2>
    static ENOKI_INLINE auto ref_wrap(T2 &&value) {
        return std::pair<decltype(enoki::ref_wrap(value.first)),
                         decltype(enoki::ref_wrap(value.second))>(
            enoki::ref_wrap(value.first), enoki::ref_wrap(value.second));
    }

    template <typename T2, typename Mask>
    static ENOKI_INLINE auto masked(T2 &&value, const Mask &mask) {
        return std::pair<decltype(enoki::masked(value.first, mask)),
                         decltype(enoki::masked(value.second, mask))>(
            enoki::masked(value.first, mask), enoki::masked(value.second, mask));
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE void scatter(T2 &dst, const Value &value, const Index &index, const Mask &mask) {
        enoki::scatter(dst.first, value.first, index, mask);
        enoki::scatter(dst.second, value.second, index, mask);
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE Value gather(const T2 &src, const Index &index, const Mask &mask) {
        return Value(
            enoki::gather<Arg0>(src.first, index, mask),
            enoki::gather<Arg1>(src.second, index, mask)
        );
    }

    static ENOKI_INLINE Value zero(size_t size) {
        return Value(enoki::zero<Arg0>(size), enoki::zero<Arg1>(size));
    }

    static ENOKI_INLINE Value empty(size_t size) {
        return Value(enoki::empty<Arg0>(size), enoki::empty<Arg1>(size));
    }
};

template <typename... Args> struct struct_support<std::tuple<Args...>> {
    static constexpr bool IsDynamic = std::disjunction_v<enoki::is_dynamic<Args>...>;
    using Dynamic = std::tuple<enoki::make_dynamic_t<Args>...>;
    using Value = std::tuple<Args...>;

    static ENOKI_INLINE size_t slices(const Value &value) {
        return enoki::slices(std::get<0>(value));
    }

    static ENOKI_INLINE size_t packets(const Value &value) {
        return enoki::packets(std::get<0>(value));
    }

    static ENOKI_INLINE void set_slices(Value &value, size_t size) {
        set_slices(value, size, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2>
    static ENOKI_INLINE auto packet(T2 &&value, size_t i) {
        return packet(std::forward<T2>(value), i, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice(T2 &&value, size_t i) {
        return slice(std::forward<T2>(value), i, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice_ptr(T2 &&value, size_t i) {
        return slice_ptr(std::forward<T2>(value), i, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2>
    static ENOKI_INLINE auto ref_wrap(T2 &&value) {
        return ref_wrap(std::forward<T2>(value), std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2, typename Mask>
    static ENOKI_INLINE auto masked(T2 &&value, const Mask &mask) {
        return masked(value, mask, std::make_index_sequence<sizeof...(Args)>());
    }

    static ENOKI_INLINE Value zero(size_t size) {
        return Value(enoki::zero<Args>(size)...);
    }

    static ENOKI_INLINE Value empty(size_t size) {
        return Value(enoki::empty<Args>(size)...);
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE void scatter(T2 &dst, const Value &value, const Index &index, const Mask &mask) {
        scatter(dst, value, index, mask, std::make_index_sequence<sizeof...(Args)>());
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE Value gather(const T2 &src, const Index &index, const Mask &mask) {
        return gather(src, index, mask, std::make_index_sequence<sizeof...(Args)>());
    }
private:
    template <size_t... Index>
    static ENOKI_INLINE void set_slices(Value &value, size_t i, std::index_sequence<Index...>) {
        bool unused[] = { (enoki::set_slices(std::get<Index>(value), i), false)..., false };
        (void) unused;
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto packet(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::tuple<decltype(enoki::packet(std::get<Index>(value), i))...>(
            enoki::packet(std::get<Index>(value), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::tuple<decltype(enoki::slice(std::get<Index>(value), i))...>(
            enoki::slice(std::get<Index>(value), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice_ptr(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::tuple<decltype(enoki::slice_ptr(std::get<Index>(value), i))...>(
            enoki::slice_ptr(std::get<Index>(value), i)...);
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto ref_wrap(T2 &&value, std::index_sequence<Index...>) {
        return std::tuple<decltype(enoki::ref_wrap(std::get<Index>(value)))...>(
            enoki::ref_wrap(std::get<Index>(value))...);
    }

    template <typename T2, typename Mask, size_t... Index>
    static ENOKI_INLINE auto masked(T2 &&value, const Mask &mask, std::index_sequence<Index...>) {
        return std::tuple<decltype(enoki::masked(std::get<Index>(value), mask))...>(
            enoki::masked(std::get<Index>(value), mask)...);
    }

    template <typename T2, typename Index, typename Mask, size_t... Is>
    static ENOKI_INLINE void scatter(T2 &dst, const Value &value, const Index &index, const Mask &mask, std::index_sequence<Is...>) {
        bool unused[] = { (enoki::scatter(std::get<Is>(dst),
                                          std::get<Is>(value), index, mask), false)..., false };
        ENOKI_MARK_USED(unused);
    }

    template <typename T2, typename Index, typename Mask, size_t... Is>
    static ENOKI_INLINE Value gather(const T2 &src, const Index &index, const Mask &mask, std::index_sequence<Is...>) {
        return Value(
            enoki::gather<std::tuple_element_t<Is, Value>>(std::get<Is>(src), index, mask)...
        );
    }
};

template <typename T, size_t Size> struct struct_support<std::array<T, Size>> {
    static constexpr bool IsDynamic = enoki::is_dynamic_v<T>;
    using Dynamic = std::array<enoki::make_dynamic_t<T>, Size>;
    using Value = std::array<T, Size>;

    static ENOKI_INLINE size_t slices(const Value &value) {
        return enoki::slices(value[0]);
    }

    static ENOKI_INLINE size_t packets(const Value &value) {
        return enoki::packets(value[0]);
    }

    static ENOKI_INLINE void set_slices(Value &value, size_t size) {
        for (size_t i = 0; i < Size; ++i)
            enoki::set_slices(value[i], size);
    }

    template <typename T2>
    static ENOKI_INLINE auto packet(T2 &&value, size_t i) {
        return packet(std::forward<T2>(value), i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice(T2 &&value, size_t i) {
        return slice(std::forward<T2>(value), i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto slice_ptr(T2 &&value, size_t i) {
        return slice_ptr(std::forward<T2>(value), i, std::make_index_sequence<Size>());
    }

    template <typename T2>
    static ENOKI_INLINE auto ref_wrap(T2 &&value) {
        return ref_wrap(std::forward<T2>(value), std::make_index_sequence<Size>());
    }

    template <typename T2, typename Mask>
    static ENOKI_INLINE auto masked(T2 &value, const Mask &mask) {
        return masked(value, mask, std::make_index_sequence<Size>());
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE void scatter(T2 &dst, const Value &value, const Index &index, const Mask &mask) {
        scatter(dst, value, index, mask, std::make_index_sequence<Size>());
    }

    template <typename T2, typename Index, typename Mask>
    static ENOKI_INLINE Value gather(const T2 &src, const Index &index, const Mask &mask) {
        return gather(src, index, mask, std::make_index_sequence<Size>());
    }

    static ENOKI_INLINE auto zero(size_t size) {
        return zero(size, std::make_index_sequence<Size>());
    }

    static ENOKI_INLINE auto empty(size_t size) {
        return empty(size, std::make_index_sequence<Size>());
    }
private:
    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto packet(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::array<decltype(enoki::packet(value[0], i)), Size>{{
            enoki::packet(value[Index], i)...}};
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::array<decltype(enoki::slice(value[0], i)), Size>{{
            enoki::slice(value[Index], i)...}};
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto slice_ptr(T2 &&value, size_t i, std::index_sequence<Index...>) {
        return std::array<decltype(enoki::slice_ptr(value[0], i)), Size>{{
            enoki::slice_ptr(value[Index], i)...}};
    }

    template <typename T2, size_t... Index>
    static ENOKI_INLINE auto ref_wrap(T2 &&value, std::index_sequence<Index...>) {
        return std::array<decltype(enoki::ref_wrap(value[0])), Size>{{
            enoki::ref_wrap(value[Index])...}};
    }

    template <typename T2, typename Mask, size_t... Index>
    static ENOKI_INLINE auto masked(T2 &value, const Mask &mask, std::index_sequence<Index...>) {
        return std::array<decltype(enoki::masked(value[0], mask)), Size>{{
            enoki::masked(value[Index], mask)...}};
    }

    template <size_t... Index>
    static ENOKI_INLINE auto zero(size_t size, std::index_sequence<Index...>) {
        return Value{{ zero<T>(Index, size)... }};
    }

    template <size_t... Index>
    static ENOKI_INLINE auto empty(size_t size, std::index_sequence<Index...>) {
        return Value{{ empty<T>(Index, size)... }};
    }

    template <typename T2, typename Index, typename Mask, size_t... Is>
    static ENOKI_INLINE void scatter(T2 &dst, const Value &value, const Index &index, const Mask &mask, std::index_sequence<Is...>) {
        bool unused[] = { (enoki::scatter(dst[Is], value[Is], index, mask), false)..., false };
        ENOKI_MARK_USED(unused);
    }

    template <typename T2, typename Index, typename Mask, size_t... Is>
    static ENOKI_INLINE Value gather(const T2 &src, const Index &index, const Mask &mask, std::index_sequence<Is...>) {
        return Value{
            enoki::gather<T>(src[Is], index, mask)...
        };
    }
};

NAMESPACE_END(enoki)
