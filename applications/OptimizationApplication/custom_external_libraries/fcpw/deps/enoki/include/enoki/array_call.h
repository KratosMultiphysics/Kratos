/*
    enoki/array_call.h -- Enoki arrays of pointers, support for
    array (virtual) method calls

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_generic.h>

NAMESPACE_BEGIN(enoki)

template <typename Class, typename Storage> struct call_support {
    call_support(const Storage &) { }
};

template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayImpl<Value_, Size_, IsMask_, Derived_,
                       enable_if_t<detail::array_config<Value_, Size_>::use_pointer_impl>>
    : StaticArrayImpl<uintptr_t, Size_, IsMask_, Derived_> {

    using UnderlyingType = std::uintptr_t;

    using Base = StaticArrayImpl<UnderlyingType, Size_, IsMask_, Derived_>;

    ENOKI_ARRAY_DEFAULTS(StaticArrayImpl)
    using Base::derived;

    using Value = std::conditional_t<IsMask_, typename Base::Value, Value_>;
    using Scalar = std::conditional_t<IsMask_, typename Base::Scalar, Value_>;

    StaticArrayImpl() = default;
    StaticArrayImpl(Value value) : Base(UnderlyingType(value)) { }
    StaticArrayImpl(std::nullptr_t) : Base(UnderlyingType(0)) { }

    template <typename T, enable_if_t<!std::is_pointer_v<T>> = 0>
    StaticArrayImpl(const T &b) : Base(b) { }

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

    auto operator->() const {
        using BaseType = std::decay_t<std::remove_pointer_t<scalar_t<Derived_>>>;
        return call_support<BaseType, Derived_>(derived());
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

NAMESPACE_BEGIN(detail)
template <typename, template <typename...> typename T, typename... Args>
struct is_callable : std::false_type {};
template <template <typename...> typename T, typename... Args>
struct is_callable<std::void_t<T<Args...>>, T, Args...> : std::true_type { };
template <template <typename...> typename T, typename... Args>
constexpr bool is_callable_v = is_callable<void, T, Args...>::value;

template <typename Guide, typename Result, typename = int> struct vectorize_result {
    using type = Result;
};

template <typename Guide, typename Result> struct vectorize_result<Guide, Result, enable_if_t<is_scalar_v<Result>>> {
    using type = replace_scalar_t<array_t<Guide>, Result, false>;
};

template <typename T, typename Perm>
decltype(auto) gather_helper(T&& v, const Perm &perm) {
    ENOKI_MARK_USED(perm);
    using DT = std::decay_t<T>;
    if constexpr (!is_cuda_array_v<DT> && !std::is_class_v<DT>)
        return v;
    else
        return gather<std::decay_t<DT>, 0, true, true>(v, perm);
}

template <typename Storage_> struct call_support_base {
    using Storage = Storage_;
    using InstancePtr = value_t<Storage_>;
    using Mask = mask_t<Storage_>;
    call_support_base(const Storage &self) : self(self) { }
    const Storage &self;

    template <typename Func, typename InputMask,
              typename Tuple, size_t ... Indices>
    ENOKI_INLINE auto dispatch(Func func, InputMask mask_, Tuple tuple,
                               std::index_sequence<Indices...>) const {
        Mask mask = Mask(mask_) & neq(self, nullptr);

        using FuncResult = decltype(func(
            std::declval<InstancePtr>(),
            mask,
            std::get<Indices>(tuple)...
        ));

        if constexpr (!std::is_void_v<FuncResult>) {
            using Result = typename vectorize_result<Mask, FuncResult>::type;
            Result result = zero<Result>(self.size());

            if constexpr (!is_cuda_array_v<Storage>) {
                while (any(mask)) {
                    InstancePtr value      = extract(self, mask);
                    Mask active            = mask & eq(self, value);
                    mask                   = andnot(mask, active);
                    masked(result, active) = func(value, active, std::get<Indices>(tuple)...);
                }
            } else {
                auto partitioned = partition(self);

                if (partitioned.size() == 1 && partitioned[0].first != nullptr) {
                    result = func(partitioned[0].first, true,
                                  std::get<Indices>(tuple)...);
                } else {
                    for (auto [value, permutation] : partitioned) {
                        if (value == nullptr)
                            continue;

                        Result temp = func(value, gather_helper(mask, permutation),
                                           gather_helper(std::get<Indices>(tuple),
                                                         permutation)...);

                        scatter<0, true, true>(result, temp, permutation);
                    }
                }
            }

            return result;
        } else {
            if constexpr (!is_cuda_array_v<Storage>) {
                while (any(mask)) {
                    InstancePtr value = extract(self, mask);
                    Mask active       = mask & eq(self, value);
                    mask              = andnot(mask, active);
                    func(value, active, std::get<Indices>(tuple)...);
                }
            } else {
                auto partitioned = partition(self);

                if (partitioned.size() == 1 && partitioned[0].first != nullptr) {
                    func(partitioned[0].first, true, std::get<Indices>(tuple)...);
                } else {
                    for (auto [value, permutation] : partitioned) {
                        if (value == nullptr)
                            continue;

                        func(value, gather_helper(mask, permutation),
                             gather_helper(std::get<Indices>(tuple),
                                           permutation)...);
                    }
                }
            }
        }
    }
};

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-value"
#endif

template <typename... Ts>
inline constexpr bool last_of(Ts... values) { return (false, ..., values); }

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

NAMESPACE_END(detail)

#define ENOKI_CALL_SUPPORT_FRIEND()                                            \
    template <typename, typename> friend struct enoki::call_support;

#define ENOKI_CALL_SUPPORT_BEGIN(Class_)                                       \
    namespace enoki {                                                          \
    template <typename Storage>                                                \
    struct call_support<Class_, Storage> : detail::call_support_base<Storage> {\
        using Base = detail::call_support_base<Storage>;                       \
        using Base::Base;                                                      \
        using typename Base::Mask;                                             \
        using Class = Class_;                                                  \
        using typename Base::InstancePtr;                                      \
        using Base::self;                                                      \
        auto operator-> () { return this; }

#define ENOKI_CALL_SUPPORT_TEMPLATE_BEGIN(Class_)                              \
    namespace enoki {                                                          \
    template <typename Storage, typename... Ts>                                \
    struct call_support<Class_<Ts...>, Storage>                                \
        : detail::call_support_base<Storage> {                                 \
        using Base = detail::call_support_base<Storage>;                       \
        using Base::Base;                                                      \
        using typename Base::Mask;                                             \
        using Class = Class_<Ts...>;                                           \
        using typename Base::InstancePtr;                                      \
        using Base::self;                                                      \
        auto operator-> () { return this; }

#define ENOKI_CALL_SUPPORT_METHOD(func)                                        \
private:                                                                       \
    template <typename... Args>                                                \
    using __##func##_t =                                                       \
        decltype(std::declval<InstancePtr>()->func(std::declval<Args>()...));  \
                                                                               \
public:                                                                        \
    template <typename... Args> auto func(Args&&... args) const {              \
        auto lambda = [](InstancePtr instance, const Mask &mask,               \
                         auto &&... a) ENOKI_INLINE_LAMBDA {                   \
            ENOKI_MARK_USED(mask);                                             \
            /* Does the method accept a mask argument? If so, provide. */      \
            if constexpr (detail::is_callable_v<__##func##_t, decltype(a)...,  \
                                                Mask>)                         \
                return instance->func(a..., mask);                             \
            else                                                               \
                return instance->func(a...);                                   \
        };                                                                     \
        /* Was a mask provided to this function? If not, set to all ones. */   \
        auto args_tuple = std::tie(args...);                                   \
        if constexpr (detail::last_of(is_mask_v<Args>...)) {                   \
            return Base::dispatch(                                             \
                lambda, std::get<sizeof...(Args) - 1>(args_tuple), args_tuple, \
                std::make_index_sequence<sizeof...(Args) - 1>());              \
        } else {                                                               \
            return Base::dispatch(                                             \
                lambda, true, args_tuple,                                      \
                std::make_index_sequence<sizeof...(Args)>());                  \
        }                                                                      \
    }

#define ENOKI_CALL_SUPPORT_GETTER_TYPE(name, field, type)                      \
    template <                                                                 \
        typename Field = decltype(Class::field),                               \
        typename Return = replace_scalar_t<Storage, type, false>>              \
    Return name(Mask mask = true) const {                                      \
        using IntType = replace_scalar_t<Storage, std::uintptr_t, false>;      \
        auto offset =                                                          \
           IntType(self) + (std::uintptr_t) &(((Class *) nullptr)->field);     \
        mask &= neq(self, nullptr);                                            \
        return gather<Return, 1>(nullptr, offset, mask);                       \
    }

#define ENOKI_CALL_SUPPORT_GETTER(name, field)                                 \
    ENOKI_CALL_SUPPORT_GETTER_TYPE(name, field, Field)

#define ENOKI_CALL_SUPPORT_END(Name)                                           \
        };                                                                     \
    }

#define ENOKI_CALL_SUPPORT_TEMPLATE_END(Name)                                  \
    ENOKI_CALL_SUPPORT_END(Name)

NAMESPACE_END(enoki)
