/*
    enoki/array_generic.h -- Generic array implementation that forwards
    all operations to the underlying data type (usually without making use of
    hardware vectorization)

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_static.h>
#include <functional>

NAMESPACE_BEGIN(nanogui)
template <typename Value, size_t Size> struct Array;
NAMESPACE_END(nanogui)

NAMESPACE_BEGIN(enoki)

namespace detail {
    template <typename StorageType, typename T>
    using is_constructible = std::bool_constant<
         std::is_constructible_v<StorageType, T> &&
        !std::is_same_v<std::decay_t<T>, reinterpret_flag>>;

    template <typename T>
    using is_not_reinterpret_flag = std::bool_constant<
        !std::is_same_v<std::decay_t<T>, reinterpret_flag>>;

    template <typename Source, typename Target>
    constexpr bool broadcast =
        !is_static_array_v<Source> || array_size_v<Source> != Target::Size ||
        !(array_depth_v<Source> == array_depth_v<Target> ||
          (array_depth_v<Source> < array_depth_v<Target> &&
           detail::array_broadcast_outer_v<Source>));

    template <typename Value, size_t Size, typename = int>
    struct is_native {
        static constexpr bool value = false;
    };

    template <typename Value, size_t Size>
    constexpr bool is_native_v = is_native<Value, Size>::value;

    /**
     * \brief The class StaticArrayImpl has several different implementations.
     * This class specifies which one to use.
     */
    template <typename Value, size_t Size>
    struct array_config {
        /// Use SSE/AVX/NEON implementation
        static constexpr bool use_native_impl =
            is_native_v<Value, Size>;

        /// Reduce to several recursive operations
        static constexpr bool use_recursive_impl =
            !use_native_impl &&
            is_std_type_v<Value> &&
            has_vectorization &&
            Size > 3;

        /// Special case for arrays of enumerations
        static constexpr bool use_enum_impl =
            std::is_enum_v<Value>;

        /// Special case for arrays of pointers of classes
        static constexpr bool use_pointer_impl =
             std::is_pointer_v<Value> &&
            !std::is_arithmetic_v<std::remove_pointer_t<Value>>;

        /// Catch-all for anything that wasn't matched so far
        static constexpr bool use_generic_impl =
            !use_native_impl &&
            !use_recursive_impl &&
            !use_enum_impl &&
            !use_pointer_impl;
    };

    template <typename T>
    using has_bitmask = decltype(std::declval<T>().bitmask_());
    template <typename T>
    constexpr bool has_bitmask_v = is_detected_v<has_bitmask, T>;
}

/// Macro to initialize uninitialized floating point arrays with 1 bits (NaN/-1) in debug mode
#if defined(NDEBUG)
#define ENOKI_TRIVIAL_CONSTRUCTOR(Value)                                       \
    template <typename T = Value,                                              \
              enable_if_t<std::is_default_constructible_v<T>> = 0>             \
    ENOKI_INLINE StaticArrayImpl() { }
#else
#define ENOKI_TRIVIAL_CONSTRUCTOR(Value)                                       \
    template <typename T = Value, enable_if_t<std::is_scalar_v<T>> = 0>        \
    ENOKI_INLINE StaticArrayImpl()                                             \
        : StaticArrayImpl(memcpy_cast<T>(int_array_t<T>(-1))) { }              \
    template <typename T = Value,                                              \
              enable_if_t<!std::is_scalar_v<T> &&                              \
                           std::is_default_constructible_v<T>> = 0>            \
    ENOKI_INLINE StaticArrayImpl() {}
#endif


/// SFINAE macro for constructors that convert from another type
#define ENOKI_CONVERT(Value)                                                   \
    template <typename Value2, typename Derived2,                               \
              enable_if_t<detail::is_same_v<Value2, Value>> = 0>               \
    ENOKI_INLINE StaticArrayImpl(                                              \
        const StaticArrayBase<Value2, Size, IsMask_, Derived2> &a)

/// SFINAE macro for constructors that reinterpret another type
#define ENOKI_REINTERPRET(Value)                                               \
    template <typename Value2, typename Derived2, bool IsMask2,                \
              enable_if_t<detail::is_same_v<Value2, Value>> = 0>               \
    ENOKI_INLINE StaticArrayImpl(                                              \
        const StaticArrayBase<Value2, Size, IsMask2, Derived2> &a,             \
        detail::reinterpret_flag)

#define ENOKI_ARRAY_DEFAULTS(Array)                                            \
    Array(const Array &) = default;                                            \
    Array(Array &&) = default;                                                 \
    Array &operator=(const Array &) = default;                                 \
    Array &operator=(Array &&) = default;

/// Import the essentials when declaring an array subclass
#define ENOKI_ARRAY_IMPORT_BASIC(Base, Array)                                  \
    ENOKI_ARRAY_DEFAULTS(Array)                                                \
    using typename Base::Derived;                                              \
    using typename Base::Value;                                                \
    using typename Base::Scalar;                                               \
    using Base::Size;                                                          \
    using Base::derived;                                                       \

/// Import the essentials when declaring an array subclass (+constructor/assignment op)
#define ENOKI_ARRAY_IMPORT(Base, Array)                                        \
    ENOKI_ARRAY_IMPORT_BASIC(Base, Array)                                      \
    using Base::Base;                                                          \
    using Base::operator=;


/// Internal macro for native StaticArrayImpl overloads (SSE, AVX, ..)
#define ENOKI_NATIVE_ARRAY(Value_, Size_, Register_)                           \
    using Base =                                                               \
        StaticArrayBase<Value_, Size_, IsMask_, Derived_>;                     \
    ENOKI_ARRAY_IMPORT_BASIC(Base, StaticArrayImpl)                            \
    using typename Base::Array1;                                               \
    using typename Base::Array2;                                               \
    using Base::ActualSize;                                                    \
    using Ref = const Derived &;                                               \
    using Register = Register_;                                                \
    static constexpr bool IsNative = true;                                     \
    Register m;                                                                \
    ENOKI_TRIVIAL_CONSTRUCTOR(Value_)                                          \
    ENOKI_INLINE StaticArrayImpl(Register value) : m(value) {}                 \
    ENOKI_INLINE StaticArrayImpl(Register value, detail::reinterpret_flag)     \
        : m(value) { }                                                         \
    ENOKI_INLINE StaticArrayImpl(bool b, detail::reinterpret_flag)             \
        : StaticArrayImpl(b ? memcpy_cast<Value_>(int_array_t<Value>(-1))      \
                            : memcpy_cast<Value_>(int_array_t<Value>(0))) { }  \
    template <typename Value2, size_t Size2, typename Derived2,                \
              enable_if_t<is_scalar_v<Value2>> = 0>                            \
    ENOKI_INLINE StaticArrayImpl(                                              \
        const StaticArrayBase<Value2, Size2, IsMask_, Derived2> &a)            \
        : Base(a) { }                                                          \
    ENOKI_INLINE StaticArrayImpl &operator=(const Derived &v) {                \
        m = v.m;                                                               \
        return *this;                                                          \
    }                                                                          \
    template <typename T> ENOKI_INLINE StaticArrayImpl &operator=(const T &v) {\
        return operator=(Derived(v)); return *this;                            \
    }                                                                          \
    ENOKI_INLINE Value& raw_coeff_(size_t i) {                                 \
        union Data {                                                           \
            Register value;                                                    \
            Value data[Size_];                                                 \
        };                                                                     \
        return ((Data *) &m)->data[i];                                         \
    }                                                                          \
    ENOKI_INLINE const Value& raw_coeff_(size_t i) const {                     \
        union Data {                                                           \
            Register value;                                                    \
            Value data[Size_];                                                 \
        };                                                                     \
        return ((const Data *) &m)->data[i];                                   \
    }                                                                          \
    ENOKI_INLINE decltype(auto) coeff(size_t i) {                              \
        if constexpr (Derived::IsMask)                                         \
            return MaskBit<Derived &>(derived(), i);                           \
        else                                                                   \
            return raw_coeff_(i);                                              \
    }                                                                          \
    ENOKI_INLINE decltype(auto) coeff(size_t i) const {                        \
        if constexpr (Derived::IsMask)                                         \
            return MaskBit<const Derived &>(derived(), i);                     \
        else                                                                   \
            return raw_coeff_(i);                                              \
    }                                                                          \
    ENOKI_INLINE bool bit_(size_t i) const {                                   \
        return detail::convert_mask(raw_coeff_(i));                            \
    }                                                                          \
    ENOKI_INLINE void set_bit_(size_t i, bool value) {                         \
        raw_coeff_(i) = reinterpret_array<Value>(value);                       \
    }

/// Internal macro for native StaticArrayImpl overloads -- 3D special case
#define ENOKI_DECLARE_3D_ARRAY(Array)                                          \
    ENOKI_ARRAY_DEFAULTS(Array)                                                \
    using typename Base::Value;                                                \
    using typename Base::Derived;                                              \
    using typename Base::Ref;                                                  \
    using Base::m;                                                             \
    using Base::coeff;                                                         \
    static constexpr size_t Size = 3;                                          \
    Array() = default;                                                         \
    ENOKI_INLINE Array(Value v) : Base(v) { }                                  \
    ENOKI_INLINE Array(Value f1, Value f2, Value f3)                           \
        : Base(f1, f2, f3, (Value) 0) { }                                      \
    ENOKI_INLINE Array(Value f1, Value f2, Value f3, Value f4)                 \
        : Base(f1, f2, f3, f4) { }                                             \
    ENOKI_INLINE Array(typename Base::Register r) : Base(r) { }                \
    ENOKI_INLINE Array(typename Base::Register r, detail::reinterpret_flag)    \
        : Base(r, detail::reinterpret_flag()) { }                              \
    ENOKI_INLINE Array(bool b, detail::reinterpret_flag)                       \
        : Base(b, detail::reinterpret_flag()) { }                              \
    template <typename Value2, typename Derived2>                              \
    ENOKI_INLINE Array(const StaticArrayBase<Value2, 4, IsMask_, Derived2> &a) \
        : Base(a) { }                                                          \
    template <typename Value2, bool IsMask2, typename Derived2>                \
    ENOKI_INLINE Array(const StaticArrayBase<Value2, 4, IsMask2, Derived2> &a, \
                       detail::reinterpret_flag)                               \
        : Base(a, detail::reinterpret_flag()) { }                              \
    template <typename Value2, typename Derived2>                              \
    ENOKI_INLINE Array(const StaticArrayBase<Value2, 3, IsMask_, Derived2>&a) {\
        ENOKI_TRACK_SCALAR("Constructor (conversion, 3D case)");               \
        Base::operator=(Derived(Value(a.derived().coeff(0)),                   \
                                Value(a.derived().coeff(1)),                   \
                                Value(a.derived().coeff(2))));                 \
    }                                                                          \
    template <typename Value2, typename Derived2, bool IsMask2>                \
    ENOKI_INLINE Array(const StaticArrayBase<Value2, 3, IsMask2, Derived2> &a, \
                       detail::reinterpret_flag) {                             \
        ENOKI_TRACK_SCALAR("Constructor (reinterpreting, 3D case)");           \
        Base::operator=(                                                       \
            Derived(reinterpret_array<Value>(a.derived().coeff(0)),            \
                    reinterpret_array<Value>(a.derived().coeff(1)),            \
                    reinterpret_array<Value>(a.derived().coeff(2))));          \
    }                                                                          \
    template <typename T> Array &operator=(T &&value) {                        \
        return (Array&) Base::operator=(Derived(value));                       \
    }


template <typename Value_, size_t Size_, bool IsMask_, typename Derived_, typename = int>
struct StaticArrayImpl;

template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayImpl<
    Value_, Size_, IsMask_, Derived_,
    enable_if_t<detail::array_config<Value_, Size_>::use_generic_impl>>
    : StaticArrayBase<std::conditional_t<IsMask_, mask_t<Value_>, Value_>,
                      Size_, IsMask_, Derived_> {

    using Base =
        StaticArrayBase<std::conditional_t<IsMask_, mask_t<Value_>, Value_>,
                        Size_, IsMask_, Derived_>;

    using typename Base::Derived;
    using typename Base::Value;
    using typename Base::Scalar;
    using typename Base::Array1;
    using typename Base::Array2;

    using Base::Size;
    using Base::derived;

    using StorageType =
        std::conditional_t<std::is_reference_v<Value> && Size_ != 0,
                           std::reference_wrapper<std::remove_reference_t<Value>>,
                           std::remove_reference_t<Value>>;

    using Ref = std::remove_reference_t<Value> &;
    using ConstRef = const std::remove_reference_t<Value> &;

    StaticArrayImpl(const StaticArrayImpl &) = default;
    StaticArrayImpl(StaticArrayImpl &&) = default;

    /// Trivial constructor
    ENOKI_TRIVIAL_CONSTRUCTOR(Value)

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable:4244) // warning C4244: 'argument': conversion from 'int' to 'Value_', possible loss of data
#  pragma warning(disable:4554) // warning C4554: '>>': check operator precedence for possible error; use parentheses to clarify precedence
#  pragma warning(disable:4702) // warning C4702: unreachable code
#elif defined(__GNUC__)
// Don't be so noisy about sign conversion in constructor
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#  pragma GCC diagnostic ignored "-Wdouble-promotion"
#  pragma GCC diagnostic ignored "-Wunused-value"
#endif

    template <typename Src>
    using cast_t = std::conditional_t<
        std::is_scalar_v<Value> ||
            !std::is_same_v<std::decay_t<Value>, std::decay_t<Src>>,
        expr_t<Value>,
        std::conditional_t<std::is_reference_v<Src>, Src, Src &&>>;

    /// Construct from component values
    template <typename... Ts, enable_if_t<sizeof...(Ts) == Size_ && Size_ != 1 &&
              std::conjunction_v<detail::is_constructible<StorageType, Ts>...>> = 0>
    ENOKI_INLINE StaticArrayImpl(Ts&&... ts)
        : m_data{{ cast_t<Ts>(ts)... }} {
        ENOKI_CHKSCALAR("Constructor (component values)");
    }

    /// Construct from a scalar or another array
    template <typename T, typename ST = StorageType,
              enable_if_t<!std::is_default_constructible_v<ST>> = 0>
    ENOKI_INLINE StaticArrayImpl(T &&value)
        : StaticArrayImpl(std::forward<T>(value),
                          std::make_index_sequence<Derived::Size>()) { }

    template <typename T, typename ST = StorageType,
              enable_if_t<!std::is_default_constructible_v<ST>> = 0>
    ENOKI_INLINE StaticArrayImpl(T &&value, detail::reinterpret_flag)
        : StaticArrayImpl(std::forward<T>(value),
                          std::make_index_sequence<Derived::Size>()) { }

    /// Construct from a scalar or another array (potential optimizations)
    template <typename T, typename ST = StorageType,
              enable_if_t<std::is_default_constructible_v<ST>> = 0>
    ENOKI_INLINE StaticArrayImpl(T &&value) {
        if constexpr (Derived::IsMask) {
            derived() = Derived(value, detail::reinterpret_flag());
        } else if constexpr (is_recursive_array_v<T> &&
                             array_depth_v<T> == array_depth_v<Derived>) {
            derived() = Derived(Array1(low(value)), Array2(high(value)));
        } else {
            assign_(std::forward<T>(value),
                    std::make_index_sequence<Derived::Size>());
        }
    }

    /// Reinterpret another array (potential optimizations)
    template <typename T, typename ST = StorageType,
              enable_if_t<std::is_default_constructible_v<ST>> = 0>
    ENOKI_INLINE StaticArrayImpl(T&& value, detail::reinterpret_flag) {
        if constexpr (is_recursive_array_v<T> &&
                      array_depth_v<T> == array_depth_v<Derived>) {
            derived() = Derived(reinterpret_array<Array1>(low(value)),
                                reinterpret_array<Array2>(high(value)));
        } else {
            assign_(std::forward<T>(value), detail::reinterpret_flag(),
                    std::make_index_sequence<Derived::Size>());
        }
    }

    template <typename T> ENOKI_INLINE StaticArrayImpl &operator=(T &&value) {
        assign_(std::forward<T>(value),
                std::make_index_sequence<Derived::Size>());
        return *this;
    }

    StaticArrayImpl& operator=(const StaticArrayImpl& value) {
        assign_(value, std::make_index_sequence<Derived::Size>());
        return *this;
    }

    StaticArrayImpl& operator=(StaticArrayImpl& value) {
        assign_(value, std::make_index_sequence<Derived::Size>());
        return *this;
    }

    StaticArrayImpl& operator=(StaticArrayImpl&& value) {
        assign_(std::move(value), std::make_index_sequence<Derived::Size>());
        return *this;
    }

    /// Construct from sub-arrays
    template <typename T1, typename T2, typename T = StaticArrayImpl, enable_if_t<
              array_depth_v<T1> == array_depth_v<T> && array_size_v<T1> == Base::Size1 &&
              array_depth_v<T2> == array_depth_v<T> && array_size_v<T2> == Base::Size2 &&
              Base::Size2 != 0> = 0>
    StaticArrayImpl(const T1 &a1, const T2 &a2)
        : StaticArrayImpl(a1, a2, std::make_index_sequence<Base::Size1>(),
                                  std::make_index_sequence<Base::Size2>()) { }

private:
    template <typename T, size_t... Is, enable_if_t<!detail::broadcast<T, Derived>> = 0>
    ENOKI_INLINE StaticArrayImpl(T&& value, std::index_sequence<Is...>)
        : m_data{{ cast_t<decltype(value.coeff(0))>(value.coeff(Is))... }} {
        ENOKI_CHKSCALAR("Copy constructor");
    }

    template <typename T, enable_if_t<detail::broadcast<T, Derived>> = 0, size_t... Is>
    ENOKI_INLINE StaticArrayImpl(T&& value, std::index_sequence<Is...>)
        : m_data{{ (Is, value)... }} {
        ENOKI_CHKSCALAR("Copy constructor (broadcast)");
    }

    template <typename T1, typename T2, size_t... Index1, size_t... Index2>
    ENOKI_INLINE StaticArrayImpl(const T1 &a1, const T2 &a2,
                                 std::index_sequence<Index1...>,
                                 std::index_sequence<Index2...>)
        : m_data{{ a1.coeff(Index1)..., a2.coeff(Index2)... }} {
        ENOKI_CHKSCALAR("Copy constructor (from 2 components)");
    }

    template <typename T, size_t... Is>
    ENOKI_INLINE void assign_(T&& value, std::index_sequence<Is...>) {
        if constexpr (std::is_same_v<array_shape_t<T>, array_shape_t<Derived>> &&
                      std::is_same_v<Value, half>) {
            #if defined(ENOKI_X86_F16C)
                using Value2 = value_t<T>;

                if constexpr (std::is_same_v<Value2, double>) {
                    derived() = float32_array_t<T, false>(value);
                    return;
                } else if constexpr (std::is_same_v<Value2, float>) {
                    if constexpr (Size == 4) {
                        long long result = detail::mm_cvtsi128_si64(_mm_cvtps_ph(
                            value.derived().m, _MM_FROUND_CUR_DIRECTION));
                        memcpy(m_data.data(), &result, sizeof(long long));
                        return;
                    } else if constexpr (Size == 8) {
                        __m128i result = _mm256_cvtps_ph(value.derived().m,
                                                         _MM_FROUND_CUR_DIRECTION);
                        _mm_storeu_si128((__m128i *) m_data.data(), result);
                        return;
                    }
                    #if defined(ENOKI_X86_AVX512F)
                        if constexpr (Size == 16) {
                            __m256i result = _mm512_cvtps_ph(value.derived().m,
                                                             _MM_FROUND_CUR_DIRECTION);
                            _mm256_storeu_si256((__m256i *) m_data.data(), result);
                            return;
                        }
                    #endif
                }
            #endif
        }

        constexpr bool Move = !std::is_lvalue_reference_v<T> && !is_scalar_v<Value> &&
                               std::is_same_v<value_t<T>, value_t<Derived>>;
        ENOKI_MARK_USED(Move);

        if constexpr (std::is_same_v<std::decay_t<T>, nanogui::Array<Value, Size>>) {
            for (size_t i = 0; i < Size; ++i)
                coeff(i) = value[i];
        } else if constexpr (detail::broadcast<T, Derived>) {
            auto s = static_cast<cast_t<T>>(value);
            bool unused[] = { (coeff(Is) = s, false)..., false };
            (void) unused; (void) s;
        } else {
            if constexpr (Move) {
                bool unused[] = { (coeff(Is) = std::move(value.derived().coeff(Is)), false)..., false };
                (void) unused;
            } else {
                using Src = decltype(value.derived().coeff(0));
                bool unused[] = { (coeff(Is) = cast_t<Src>(value.derived().coeff(Is)), false)..., false };
                (void) unused;
            }
        }
    }

    template <typename T, size_t... Is>
    ENOKI_INLINE void assign_(T&& value, detail::reinterpret_flag, std::index_sequence<Is...>) {
        if constexpr (std::is_same_v<array_shape_t<T>, array_shape_t<Derived>> &&
                      std::is_same_v<Value, bool> && detail::has_bitmask_v<T>) {
            #if defined(ENOKI_X86_AVX512VL)
                if constexpr (Size == 16) {
                    _mm_storeu_si128((__m128i *) data(),
                                     _mm_maskz_set1_epi8((__mmask16) value.bitmask_(), (char) 1));
                    return;
                } else if constexpr (Size == 8) {
                    uint64_t result = (uint64_t) detail::mm_cvtsi128_si64(
                        _mm_maskz_set1_epi8((__mmask8) value.bitmask_(), (char) 1));
                    memcpy(data(), &result, sizeof(uint64_t));
                    return;
                } else if constexpr (Size == 4) {
                    uint32_t result = (uint32_t) _mm_cvtsi128_si32(
                        _mm_maskz_set1_epi8((__mmask8) value.bitmask_(), (char) 1));
                    memcpy(data(), &result, sizeof(uint32_t));
                    return;
                }
            #elif defined(ENOKI_X86_AVX2) && defined(ENOKI_X86_64)
                uint32_t k = value.bitmask_();
                if constexpr (Size == 16) {
                    uint64_t low = (uint64_t) _pdep_u64(k,      0x0101010101010101ull);
                    uint64_t hi  = (uint64_t) _pdep_u64(k >> 8, 0x0101010101010101ull);
                    memcpy((uint8_t *) data(), &low, sizeof(uint64_t));
                    memcpy((uint8_t *) data() + sizeof(uint64_t), &hi, sizeof(uint64_t));
                    return;
                } else if constexpr (Size == 8) {
                    uint64_t result = (uint64_t) _pdep_u64(k, 0x0101010101010101ull);
                    memcpy(data(), &result, sizeof(uint64_t));
                    return;
                } else if constexpr (Size == 4) {
                    uint32_t result = (uint32_t) _pdep_u32(k, 0x01010101ull);
                    memcpy(data(), &result, sizeof(uint32_t));
                    return;
                }
            #endif
        }

        if constexpr(detail::broadcast<T, Derived>) {
            bool unused[] = { (coeff(Is) = reinterpret_array<Value>(value), false)..., false };
            (void) unused;
        } else {
            bool unused[] = { (coeff(Is) = reinterpret_array<Value>(value.coeff(Is)), false)..., false };
            (void) unused;
        }
    }

#if defined(_MSC_VER)
#  pragma warning(pop)
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

public:
    /// Return the size in bytes
    size_t nbytes() const {
        if constexpr (is_dynamic_v<Value>) {
            size_t result = 0;
            for (size_t i = 0; i < Derived::Size; ++i)
                result += coeff(i).nbytes();
            return result;
        } else {
            return Base::nbytes();
        }
    }

    /// Arithmetic NOT operation
    ENOKI_INLINE Derived not_() const {
        Derived result;
        ENOKI_CHKSCALAR("not");
        for (size_t i = 0; i < Derived::Size; ++i) {
            if constexpr (IsMask_)
                (Value &) result.coeff(i) = !(Value) derived().coeff(i);
            else
                (Value &) result.coeff(i) = ~(Value) derived().coeff(i);
        }
        return result;
    }

    /// Arithmetic unary negation operation
    ENOKI_INLINE Derived neg_() const {
        Derived result;
        ENOKI_CHKSCALAR("neg");
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = - (Value) derived().coeff(i);
        return result;
    }

    /// Array indexing operator
    ENOKI_INLINE Ref coeff(size_t i) {
        ENOKI_CHKSCALAR("coeff");
        return m_data[i];
    }

    /// Array indexing operator (const)
    ENOKI_INLINE ConstRef coeff(size_t i) const {
        ENOKI_CHKSCALAR("coeff");
        return m_data[i];
    }

    /// Recursive array indexing operator (const)
    template <typename... Args, enable_if_t<(sizeof...(Args) >= 1)> = 0>
    ENOKI_INLINE decltype(auto) coeff(size_t i0, Args... other) const {
        return coeff(i0).coeff(size_t(other)...);
    }

    /// Recursive array indexing operator
    template <typename... Args, enable_if_t<(sizeof...(Args) >= 1)> = 0>
    ENOKI_INLINE decltype(auto) coeff(size_t i0, Args... other) {
        return coeff(i0).coeff(size_t(other)...);
    }

    StorageType *data() { return m_data.data(); }
    const StorageType *data() const { return m_data.data(); }

private:
    std::array<StorageType, Size> m_data;
};

struct BitRef {
private:
    struct BitWrapper {
        virtual bool get() = 0;
        virtual void set(bool value) = 0;
        virtual ~BitWrapper() = default;
    };

    std::unique_ptr<BitWrapper> accessor;
public:
    BitRef(bool &b) {
        struct BoolWrapper : BitWrapper {
            BoolWrapper(bool& data) : data(data) { }
            bool get() override { return data; }
            void set(bool value) override { data = value; }
            bool &data;
        };
        accessor = std::make_unique<BoolWrapper>(b);
    }

    template <typename T>
    BitRef(MaskBit<T> b) {
        struct MaskBitWrapper : BitWrapper {
            MaskBitWrapper(MaskBit<T> data) : data(data) { }
            bool get() override { return (bool) data; }
            void set(bool value) override { data = value; }
            MaskBit<T> data;
        };
        accessor = std::make_unique<MaskBitWrapper>(b);
    }

    operator bool() const { return accessor->get(); }
    BitRef& operator=(bool value) { accessor->set(value); return *this; }
};

NAMESPACE_END(enoki)
