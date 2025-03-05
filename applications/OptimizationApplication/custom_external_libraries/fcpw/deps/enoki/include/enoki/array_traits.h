/*
    enoki/array_traits.h -- Type traits for Enoki arrays

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include "fwd.h"
#include <cstdint>
#include <cmath>
#include <cassert>
#include <array>
#include <limits>
#include <iostream>
#include <string>
#include <stdexcept>
#include <tuple>
#include <memory>

NAMESPACE_BEGIN(enoki)

// -----------------------------------------------------------------------
//! @{ \name General type traits (not specific to Enoki arrays)
// -----------------------------------------------------------------------

/// Convenience wrapper around std::enable_if
template <bool B> using enable_if_t = std::enable_if_t<B, int>;

constexpr size_t Dynamic = (size_t) -1;

namespace detail {
    /// Identity function for types
    template <typename T, typename...> struct identity {
        using type = T;
    };

    template <template <typename...> typename B, typename T>
    struct is_base_of_impl {
    private:
        template <typename... Ts>
        static constexpr std::true_type test(const B<Ts...> *);
        static constexpr std::false_type test(...);

    public:
        using type = decltype(test(std::declval<T *>()));
    };

    template <typename, template <typename...> typename Op, typename... Ts>
    struct detector : std::false_type { };

    template <template <typename...> typename Op, typename... Ts>
    struct detector<std::void_t<Op<Ts...>>, Op, Ts...>
        : std::true_type { };

    template <typename... > constexpr bool false_v = false;
}

template <typename... Ts> using identity_t = typename detail::identity<Ts...>::type;

template <template<typename ...> class Op, class... Args>
constexpr bool is_detected_v = detail::detector<void, Op, Args...>::value;

/// Check if 'T' is a subtype of a given template 'B'
template <template <typename...> typename B, typename T>
using is_base_of = typename detail::is_base_of_impl<B, T>::type;

template <template <typename...> typename B, typename T>
constexpr bool is_base_of_v = is_base_of<B, T>::value;

/// Check if T is an integer of a given size (supports both 'int' and 'long' family)
template <typename T> using is_int8 = std::bool_constant<std::is_integral_v<T> && sizeof(T) == 1>;
template <typename T> constexpr bool is_int8_v = is_int8<T>::value;

template <typename T> using is_int16 = std::bool_constant<std::is_integral_v<T> && sizeof(T) == 2>;
template <typename T> constexpr bool is_int16_v = is_int16<T>::value;

template <typename T> using is_int32 = std::bool_constant<std::is_integral_v<T> && sizeof(T) == 4>;
template <typename T> constexpr bool is_int32_v = is_int32<T>::value;

template <typename T> using is_int64 = std::bool_constant<std::is_integral_v<T> && sizeof(T) == 8>;
template <typename T> constexpr bool is_int64_v = is_int64<T>::value;

template <typename T> constexpr bool is_float_v = std::is_same_v<T, float>;
template <typename T> constexpr bool is_double_v = std::is_same_v<T, double>;

template <typename T> using is_std_float = std::bool_constant<is_float_v<T> || is_double_v<T>>;
template <typename T> constexpr bool is_std_float_v = is_std_float<T>::value;

template <typename T> using is_std_int = std::bool_constant<is_int32_v<T> || is_int64_v<T>>;
template <typename T> constexpr bool is_std_int_v = is_std_int<T>::value;

template <typename T> using is_std_type = std::bool_constant<is_std_int_v<T> || is_std_float_v<T>>;
template <typename T> constexpr bool is_std_type_v = is_std_type<T>::value;

template <typename T> using enable_if_int32_t = enable_if_t<is_int32_v<T>>;
template <typename T> using enable_if_int64_t = enable_if_t<is_int64_v<T>>;
template <typename T> using enable_if_std_int_v = enable_if_t<is_std_int_v<T>>;
template <typename T> using enable_if_std_float_v = enable_if_t<is_std_float_v<T>>;
template <typename T> using enable_if_std_type_v = enable_if_t<is_std_type_v<T>>;

template <typename T> constexpr bool is_scalar_v = std::is_scalar_v<std::decay_t<T>>;

namespace detail {
    /// Value equivalence between arithmetic type to work around subtle issues between 'long' vs 'long long' on OSX
    template <typename T0, typename T1>
    struct is_same {
        static constexpr bool value =
            sizeof(T0) == sizeof(T1) &&
            std::is_floating_point_v<T0> == std::is_floating_point_v<T1> &&
            std::is_signed_v<T0> == std::is_signed_v<T1> &&
            std::is_arithmetic_v<T0> == std::is_arithmetic_v<T1>;
    };

    template <typename T0, typename T1>
    static constexpr bool is_same_v = is_same<T0, T1>::value;

    template <typename T> using has_size = std::enable_if_t<std::decay_t<T>::Size != Dynamic>;
    template <typename T> constexpr bool has_size_v = is_detected_v<has_size, T>;

    template <typename T> using is_masked_array = std::enable_if_t<T::IsMaskedArray>;
    template <typename T> constexpr bool is_masked_array_v = is_detected_v<is_masked_array, T>;
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Type traits for Enoki arrays
// -----------------------------------------------------------------------

/// Is 'T' an Enoki array? (any variant)
template <typename T> using is_array = is_base_of<ArrayBase, std::decay_t<T>>;
template <typename T> constexpr bool is_array_v = is_array<T>::value;
template <typename T> using enable_if_array_t = enable_if_t<is_array_v<T>>;
template <typename T> using enable_if_not_array_t = enable_if_t<!is_array_v<T>>;

template <typename... Ts> using is_array_any = std::disjunction<is_array<Ts>...>;
template <typename... Ts> constexpr bool is_array_any_v = is_array_any<Ts...>::value;
template <typename... Ts> using enable_if_array_any_t = enable_if_t<is_array_any_v<Ts...>>;

template <typename T> using is_static_array = std::bool_constant<is_array_v<T> && detail::has_size_v<T>>;
template <typename T> constexpr bool is_static_array_v = is_static_array<T>::value;
template <typename T> using enable_if_static_array_t = enable_if_t<is_static_array_v<T>>;

template <typename T> using is_dynamic_array = std::bool_constant<is_array_v<T> && !detail::has_size_v<T>>;
template <typename T> constexpr bool is_dynamic_array_v = is_dynamic_array<T>::value;
template <typename T> using enable_if_dynamic_array_t = enable_if_t<is_dynamic_array_v<T>>;

namespace detail {
    template <typename T, typename = int> struct value {
        using type = std::decay_t<T>;
    };

    template <typename T, typename = int> struct packet_ {
        using type = std::decay_t<T>;
    };

    template <typename T> struct value<T, enable_if_array_t<T>> {
        using type = typename std::decay_t<T>::Derived::Value;
    };

    template <typename T>
    struct packet_<
        T, enable_if_t<is_array_v<T> && !detail::is_masked_array_v<T>>> {
        using type = typename std::decay_t<T>::Derived::Value;
    };

    template <typename T>
    struct packet_<
        T, enable_if_t<is_array_v<T> && detail::is_masked_array_v<T>>> {
        using type = typename std::decay_t<T>::Derived::UnderlyingValue;
    };
}

/// Type trait to access the value type of an array
template <typename T> using value_t = typename detail::value<T>::type;

/// Is 'T' an Enoki mask or a boolean?
template <typename T, typename = int> struct is_mask {
    static constexpr bool value = std::is_same_v<std::decay_t<T>, bool>;
};

template <typename T> struct is_mask<MaskBit<T>> {
    static constexpr bool value = true;
};

template <typename T> struct is_mask<T, enable_if_array_t<T>> {
    static constexpr bool value = std::decay_t<T>::Derived::IsMask;
};

template <typename T> constexpr bool is_mask_v = is_mask<T>::value;
template <typename T> using enable_if_mask_t = enable_if_t<is_mask_v<T>>;
template <typename T> using enable_if_not_mask_t = enable_if_t<!is_mask_v<T>>;

/// Is 'T' implemented using a recursive implementation?
template <typename T, typename = int> struct is_recursive_array {
    static constexpr bool value = false;
};

template <typename T> struct is_recursive_array<T, enable_if_array_t<T>> {
    static constexpr bool value = std::decay_t<T>::Derived::IsRecursive;
};

template <typename T> constexpr bool is_recursive_array_v = is_recursive_array<T>::value;
template <typename T> using enable_if_recursive_t = enable_if_t<is_recursive_array_v<T>>;

/// Does this array compute derivatives using automatic differentiation?
template <typename T, typename = int> struct is_diff_array {
    static constexpr bool value = false;
};

template <typename T> struct is_diff_array<T, enable_if_array_t<T>> {
    static constexpr bool value = std::decay_t<T>::Derived::IsDiff;
};

template <typename T> constexpr bool is_diff_array_v = is_diff_array<T>::value;
template <typename T> using enable_if_diff_array_t = enable_if_t<is_diff_array_v<T>>;

/// Does this array reside on the GPU (via CUDA)?
template <typename T, typename = int> struct is_cuda_array {
    static constexpr bool value = false;
};

template <typename T> struct is_cuda_array<T, enable_if_array_t<T>> {
    static constexpr bool value = std::decay_t<T>::Derived::IsCUDA;
};

template <typename T> constexpr bool is_cuda_array_v = is_cuda_array<T>::value;
template <typename T> using enable_if_cuda_t = enable_if_t<is_cuda_array_v<T>>;

/// Determine the depth of a nested Enoki array (scalars evaluate to zero)
template <typename T, typename = int> struct array_depth {
    static constexpr size_t value = 0;
};

template <typename T> struct array_depth<T, enable_if_array_t<T>> {
    static constexpr size_t value = std::decay_t<T>::Derived::Depth;
};

template <typename T> constexpr size_t array_depth_v = array_depth<T>::value;

/// Determine the size of a nested Enoki array (scalars evaluate to one)
template <typename T, typename = int> struct array_size {
    static constexpr size_t value = 1;
};

template <typename T> struct array_size<T, enable_if_static_array_t<T>> {
    static constexpr size_t value = std::decay_t<T>::Derived::Size;
};

template <typename T> struct array_size<T, enable_if_dynamic_array_t<T>> {
    static constexpr size_t value = Dynamic;
};

template <typename T> constexpr size_t array_size_v = array_size<T>::value;

namespace detail {
    template <typename T, size_t>
    struct prepend_index { };

    template <size_t... Index, size_t Value>
    struct prepend_index<std::index_sequence<Index...>, Value> {
        using type = std::index_sequence<Value, Index...>;
    };

    template <typename T, size_t Value>
    using prepend_index_t = typename prepend_index<T, Value>::type;
}

/// Determine the shape of an array
template <typename T, typename = int> struct array_shape {
    using type = std::index_sequence<>;
};

template <typename T>
using array_shape_t = typename array_shape<T>::type;

template <typename T> struct array_shape<T, enable_if_array_t<T>> {
    using type = detail::prepend_index_t<array_shape_t<value_t<T>>, array_size_v<T>>;
};

namespace detail {
    template <typename T, typename = int> struct scalar {
        using type = std::decay_t<T>;
    };

    template <typename T> struct scalar<T, enable_if_array_t<T>> {
        using type = typename std::decay_t<T>::Derived::Scalar;
    };

    template <typename T> using packet_t = typename detail::packet_<T>::type;
}

/// Type trait to access the base scalar type underlying a potentially nested array
template <typename T> using scalar_t = typename detail::scalar<T>::type;

struct BitRef;

namespace detail {
    /// Copy modifier flags (const/pointer/lvalue/rvalue reference from 'S' to 'T')
    template <typename S, typename T> struct copy_flags {
    private:
        using R = std::remove_reference_t<S>;
        using T1 = std::conditional_t<std::is_const_v<R>, std::add_const_t<T>, T>;
        using T2 = std::conditional_t<std::is_pointer_v<S>,
                                      std::add_pointer_t<T1>, T1>;
        using T3 = std::conditional_t<std::is_lvalue_reference_v<S>,
                                      std::add_lvalue_reference_t<T2>, T2>;
        using T4 = std::conditional_t<std::is_rvalue_reference_v<S>,
                                      std::add_rvalue_reference_t<T3>, T3>;

    public:
        using type = T4;
    };

    template <typename S, typename T>
    using copy_flags_t = typename detail::copy_flags<S, T>::type;

    template <typename T, bool CopyFlags, typename = int> struct mask {
        using type = bool;
    };

    template <typename T, bool CopyFlags> struct mask<T&, CopyFlags, enable_if_t<is_scalar_v<T>>> {
        using type = BitRef;
    };

    template <typename T, bool CopyFlags> struct mask<T, CopyFlags, enable_if_array_t<T>> {
    private:
        using Mask = copy_flags_t<T, typename std::decay_t<T>::Derived::MaskType>;
    public:
        using type = std::conditional_t<CopyFlags, detail::copy_flags_t<T, Mask>, Mask>;
    };

    template <typename T, bool CopyFlags, typename = int> struct array { };

    template <typename T, bool CopyFlags> struct array<T, CopyFlags, enable_if_array_t<T>> {
    private:
        using Array = copy_flags_t<T, typename std::decay_t<T>::Derived::ArrayType>;
    public:
        using type = std::conditional_t<CopyFlags, detail::copy_flags_t<T, Array>, Array>;
    };
}

/// Type trait to access the mask type underlying an array
template <typename T, bool CopyFlags = true> using mask_t = typename detail::mask<T, CopyFlags>::type;

/// Type trait to access the array type underlying a mask
template <typename T, bool CopyFlags = true> using array_t = typename detail::array<T, CopyFlags>::type;

/// Extract the most deeply nested Enoki array type from a list of arguments
template <typename... Args> struct deepest_array;
template <> struct deepest_array<> { using type = void; };

template <typename Arg, typename... Args> struct deepest_array<Arg, Args...> {
private:
    using T0 = Arg;
    using T1 = typename deepest_array<Args...>::type;

    // Give precedence to dynamic arrays
    static constexpr size_t D0 = array_depth_v<T0>;
    static constexpr size_t D1 = array_depth_v<T1>;

public:
    using type = std::conditional_t<(D1 > D0 || D0 == 0), T1, T0>;
};

template <typename... Args> using deepest_array_t = typename deepest_array<Args...>::type;

namespace detail {
    template <typename... Ts> struct expr;
}

/// Type trait to compute the type of an arithmetic expression involving Ts...
template <typename... Ts> using expr_t = typename detail::expr<Ts...>::type;

namespace detail {
    /// Type trait to compute the result of a unary expression
    template <typename Array, typename T> struct expr_1;

    template <typename T> struct expr_1<T, T> {
    private:
        using Td        = std::decay_t<T>;
        using Entry     = value_t<T>;
        using EntryExpr = expr_t<Entry>;

    public:
        using type = std::conditional_t<
            std::is_same_v<Entry, EntryExpr>,
            Td, typename Td::Derived::template ReplaceValue<EntryExpr>
        >;
    };

    template <typename T>
    struct expr_1<void, T> { using type = std::decay_t<T>; };

    /// Type trait to compute the result of a n-ary expression involving types (T, Ts...)
    template <typename Array, typename T, typename... Ts>
    struct expr_n {
    private:
        using Value = expr_t<detail::packet_t<T>, detail::packet_t<Ts>...>;
    public:
        using type  = typename std::decay_t<Array>::Derived::template ReplaceValue<Value>;
    };

    template <typename T, typename... Ts>
    struct expr_n<void, T, Ts...> {
        using type = decltype(std::declval<T>() + std::declval<expr_t<Ts...>>());
    };

    template <typename T1, typename T2> struct expr_n<void, T1*, T2*> { using type = std::common_type_t<T1*, T2*>; };
    template <typename T> struct expr_n<void, T*, std::nullptr_t> { using type = T*; };
    template <typename T> struct expr_n<void, T*, unsigned long long> { using type = T*; };
    template <typename T> struct expr_n<void, T*, unsigned long> { using type = T*; };
    template <typename T> struct expr_n<void, std::nullptr_t, T*> { using type = T*; };
    template <typename T, typename T2> struct expr_n<void, T, enoki::divisor_ext<T2>> { using type = T2; };
    template <typename T, typename T2> struct expr_n<void, T, enoki::divisor<T2>> { using type = T2; };
    template <> struct expr_n<void, bool, bool> { using type = bool; };

    /// Type trait to compute the result of arbitrary expressions
    template <typename... Ts> struct expr    : detail::expr_n<deepest_array_t<Ts...>, Ts...> { };
    template <typename T>     struct expr<T> : detail::expr_1<deepest_array_t<T>,     T>     { };
}

namespace detail {
    template <typename T, typename = int> struct array_broadcast_outer {
        static constexpr bool value = true;
    };

    template <typename T> struct array_broadcast_outer<T, enable_if_array_t<T>> {
        static constexpr bool value = std::decay_t<T>::Derived::BroadcastPreferOuter;
    };

    template <typename T> constexpr bool array_broadcast_outer_v = array_broadcast_outer<T>::value;

    /// Convenience class to choose an arithmetic type based on its size and flavor
    template <size_t Size> struct type_chooser { };

    template <> struct type_chooser<1> {
        using Int = int8_t;
        using UInt = uint8_t;
    };

    template <> struct type_chooser<2> {
        using Int = int16_t;
        using UInt = uint16_t;
        using Float = half;
    };

    template <> struct type_chooser<4> {
        using Int = int32_t;
        using UInt = uint32_t;
        using Float = float;
    };

    template <> struct type_chooser<8> {
        using Int = int64_t;
        using UInt = uint64_t;
        using Float = double;
    };
}

/// Replace the base scalar type of a (potentially nested) array
template <typename T, typename Value, bool CopyFlags = true, typename = int>
struct replace_scalar { };

template <typename T, typename Value, bool CopyFlags = true>
using replace_scalar_t = typename replace_scalar<T, Value, CopyFlags>::type;

template <typename T, typename Value, bool CopyFlags> struct replace_scalar<T, Value, CopyFlags, enable_if_not_array_t<T>> {
    using type = std::conditional_t<CopyFlags, detail::copy_flags_t<T, Value>, Value>;
};

template <typename T, typename Value, bool CopyFlags> struct replace_scalar<T, Value, CopyFlags, enable_if_array_t<T>> {
private:
    using Entry = replace_scalar_t<detail::packet_t<T>, Value, CopyFlags>;
    using Array = typename std::decay_t<T>::Derived::template ReplaceValue<Entry>;
public:
    using type = std::conditional_t<CopyFlags, detail::copy_flags_t<T, Array>, Array>;
};

/// Integer-based version of a given array class
template <typename T, bool CopyFlags = true>
using int_array_t = replace_scalar_t<T, typename detail::type_chooser<sizeof(scalar_t<T>)>::Int, CopyFlags>;

/// Unsigned integer-based version of a given array class
template <typename T, bool CopyFlags = true>
using uint_array_t = replace_scalar_t<T, typename detail::type_chooser<sizeof(scalar_t<T>)>::UInt, CopyFlags>;

/// Floating point-based version of a given array class
template <typename T, bool CopyFlags = true>
using float_array_t = replace_scalar_t<T, typename detail::type_chooser<sizeof(scalar_t<T>)>::Float, CopyFlags>;


template <typename T, bool CopyFlags = true> using int32_array_t   = replace_scalar_t<T, int32_t, CopyFlags>;
template <typename T, bool CopyFlags = true> using uint32_array_t  = replace_scalar_t<T, uint32_t, CopyFlags>;
template <typename T, bool CopyFlags = true> using int64_array_t   = replace_scalar_t<T, int64_t, CopyFlags>;
template <typename T, bool CopyFlags = true> using uint64_array_t  = replace_scalar_t<T, uint64_t, CopyFlags>;
template <typename T, bool CopyFlags = true> using float16_array_t = replace_scalar_t<T, half, CopyFlags>;
template <typename T, bool CopyFlags = true> using float32_array_t = replace_scalar_t<T, float, CopyFlags>;
template <typename T, bool CopyFlags = true> using float64_array_t = replace_scalar_t<T, double, CopyFlags>;
template <typename T, bool CopyFlags = true> using bool_array_t    = replace_scalar_t<T, bool, CopyFlags>;
template <typename T, bool CopyFlags = true> using size_array_t    = replace_scalar_t<T, size_t, CopyFlags>;
template <typename T, bool CopyFlags = true> using ssize_array_t   = replace_scalar_t<T, ssize_t, CopyFlags>;

//! @}
// -----------------------------------------------------------------------

template <typename T> using struct_support_t = struct_support<std::decay_t<T>>;

// -----------------------------------------------------------------------
//! @{ \name Type enumeration
// -----------------------------------------------------------------------

enum class EnokiType { Invalid = 0, Int8, UInt8, Int16, UInt16,
                       Int32, UInt32, Int64, UInt64, Float16,
                       Float32, Float64, Bool, Pointer };

template <typename T, typename = int> struct enoki_type {
    static constexpr EnokiType value = EnokiType::Invalid;
};

template <typename T> struct enoki_type<T, enable_if_t<is_int8_v<T>>> {
    static constexpr EnokiType value =
        std::is_signed_v<T> ? EnokiType::Int8 : EnokiType::UInt8;
};

template <typename T> struct enoki_type<T, enable_if_t<is_int16_v<T>>> {
    static constexpr EnokiType value =
        std::is_signed_v<T> ? EnokiType::Int16 : EnokiType::UInt16;
};

template <typename T> struct enoki_type<T, enable_if_t<is_int32_v<T>>> {
    static constexpr EnokiType value =
        std::is_signed_v<T> ? EnokiType::Int32 : EnokiType::UInt32;
};

template <typename T> struct enoki_type<T, enable_if_t<is_int64_v<T>>> {
    static constexpr EnokiType value =
        std::is_signed_v<T> ? EnokiType::Int64 : EnokiType::UInt64;
};

template <typename T> struct enoki_type<T, enable_if_t<std::is_enum_v<T>>> {
    static constexpr EnokiType value = enoki_type<std::underlying_type_t<T>>::value;
};

template <> struct enoki_type<half> {
    static constexpr EnokiType value = EnokiType::Float16;
};

template <> struct enoki_type<float> {
    static constexpr EnokiType value = EnokiType::Float32;
};

template <> struct enoki_type<double> {
    static constexpr EnokiType value = EnokiType::Float64;
};

template <> struct enoki_type<bool> {
    static constexpr EnokiType value = EnokiType::Bool;
};

template <typename T> struct enoki_type<T *> {
    static constexpr EnokiType value = EnokiType::Pointer;
};

template <typename T> constexpr EnokiType enoki_type_v = enoki_type<T>::value;

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Type trait to inspect the return/argument types of functions
// -----------------------------------------------------------------------

template <typename T, typename SFINAE = void> struct function_traits { };

// Vanilla function
template <typename R, typename... A> struct function_traits<R(*)(A...)> {
    using Args = std::tuple<A...>;
    using Return = R;
};

// Method
template <typename C, typename R, typename... A> struct function_traits<R(C::*)(A...)> {
    using Class = C;
    using Args = std::tuple<A...>;
    using Return = R;
};

// Method (const)
template <typename C, typename R, typename... A> struct function_traits<R(C::*)(A...) const> {
    using Class = C;
    using Args = std::tuple<A...>;
    using Return = R;
};

// Lambda function -- strip lambda closure and delegate back to ``function_traits``
template <typename F>
struct function_traits<
    F, std::enable_if_t<std::is_member_function_pointer_v<decltype(
           &std::remove_reference_t<F>::operator())>>>
    : function_traits<decltype(&std::remove_reference_t<F>::operator())> { };

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
