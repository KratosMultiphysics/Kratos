/*
    enoki/array_macro.h -- Code generation macros for custom data structures

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

// The main idea of this macro is borrowed from https://github.com/swansontec/map-macro
// (C) William Swanson, Paul Fultz
#define ENOKI_EVAL_0(...) __VA_ARGS__
#define ENOKI_EVAL_1(...) ENOKI_EVAL_0(ENOKI_EVAL_0(ENOKI_EVAL_0(__VA_ARGS__)))
#define ENOKI_EVAL_2(...) ENOKI_EVAL_1(ENOKI_EVAL_1(ENOKI_EVAL_1(__VA_ARGS__)))
#define ENOKI_EVAL_3(...) ENOKI_EVAL_2(ENOKI_EVAL_2(ENOKI_EVAL_2(__VA_ARGS__)))
#define ENOKI_EVAL_4(...) ENOKI_EVAL_3(ENOKI_EVAL_3(ENOKI_EVAL_3(__VA_ARGS__)))
#define ENOKI_EVAL(...)   ENOKI_EVAL_4(ENOKI_EVAL_4(ENOKI_EVAL_4(__VA_ARGS__)))
#define ENOKI_MAP_END(...)
#define ENOKI_MAP_OUT
#define ENOKI_MAP_COMMA ,
#define ENOKI_MAP_GET_END() 0, ENOKI_MAP_END
#define ENOKI_MAP_NEXT_0(test, next, ...) next ENOKI_MAP_OUT
#define ENOKI_MAP_NEXT_1(test, next) ENOKI_MAP_NEXT_0(test, next, 0)
#define ENOKI_MAP_NEXT(test, next) ENOKI_MAP_NEXT_1(ENOKI_MAP_GET_END test, next)
#define ENOKI_EXTRACT_0(next, ...) next

#if defined(_MSC_VER) // MSVC is not as eager to expand macros, hence this workaround
#define ENOKI_MAP_EXPR_NEXT_1(test, next) \
    ENOKI_EVAL_0(ENOKI_MAP_NEXT_0(test, ENOKI_MAP_COMMA next, 0))
#define ENOKI_MAP_STMT_NEXT_1(test, next) \
    ENOKI_EVAL_0(ENOKI_MAP_NEXT_0(test, next, 0))
#else
#define ENOKI_MAP_EXPR_NEXT_1(test, next) \
    ENOKI_MAP_NEXT_0(test, ENOKI_MAP_COMMA next, 0)
#define ENOKI_MAP_STMT_NEXT_1(test, next) \
    ENOKI_MAP_NEXT_0(test, next, 0)
#endif

#define ENOKI_MAP_EXPR_NEXT(test, next) \
    ENOKI_MAP_EXPR_NEXT_1 (ENOKI_MAP_GET_END test, next)
#define ENOKI_MAP_STMT_NEXT(test, next) \
    ENOKI_MAP_STMT_NEXT_1 (ENOKI_MAP_GET_END test, next)

#define ENOKI_MAP_TEMPLATE_FWD_0(x, peek, ...) \
    typename T##x ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_TEMPLATE_FWD_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_TEMPLATE_FWD_1(x, peek, ...) \
    typename T##x ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_TEMPLATE_FWD_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_DECL_FWD_0(x, peek, ...) \
    T##x &&x ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_DECL_FWD_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_DECL_FWD_1(x, peek, ...) \
    T##x &&x ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_DECL_FWD_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_BASE_FWD_0(x, peek, ...) \
    std::forward<T##x>(x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_BASE_FWD_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_BASE_FWD_1(x, peek, ...) \
    std::forward<T##x>(x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_BASE_FWD_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_FWD_0(x, peek, ...) \
    x(std::forward<T##x>(x)) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_FWD_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_FWD_1(x, peek, ...) \
    x(std::forward<T##x>(x)) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_FWD_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_COPY_0(x, peek, ...) \
    x(x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_COPY_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_COPY_1(x, peek, ...) \
    x(x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_COPY_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_COPY_V_0(v, x, peek, ...) \
    x(v.x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_COPY_V_1)(v, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_COPY_V_1(v, x, peek, ...) \
    x(v.x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_COPY_V_0)(v, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_MOVE_V_0(v, x, peek, ...) \
    x(std::move(v.x)) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_MOVE_V_1)(v, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_MOVE_V_1(v, x, peek, ...) \
    x(std::move(v.x)) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_MOVE_V_0)(v, peek, __VA_ARGS__)

#define ENOKI_MAP_STMT_ASSIGN_0(v, x, peek, ...)                               \
    this->x = v.x;                                                             \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_MAP_STMT_ASSIGN_1)(v, peek, __VA_ARGS__)
#define ENOKI_MAP_STMT_ASSIGN_1(v, x, peek, ...)                               \
    this->x = v.x;                                                             \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_MAP_STMT_ASSIGN_0)(v, peek, __VA_ARGS__)

#define ENOKI_MAP_STMT_MOVE_0(v, x, peek, ...)                                 \
    this->x = std::move(v.x);                                                  \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_MAP_STMT_MOVE_1)(v, peek, __VA_ARGS__)
#define ENOKI_MAP_STMT_MOVE_1(v, x, peek, ...)                                 \
    this->x = std::move(v.x);                                                  \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_MAP_STMT_MOVE_0)(v, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_F1_0(f, v, x, peek, ...) \
    f(v.x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F1_1)(f, v, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_F1_1(f, v, x, peek, ...) \
    f(v.x) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F1_0)(f, v, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_F2_0(f, v, t, x, peek, ...) \
    f(v.x, t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F2_1)(f, v, t, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_F2_1(f, v, t, x, peek, ...) \
    f(v.x, t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F2_0)(f, v, t, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_F3_0(f, m, v, t, x, peek, ...) \
    f(m.x, v.x, t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F3_1)(f, m, v, t, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_F3_1(f, m, v, t, x, peek, ...) \
    f(m.x, v.x, t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_F3_0)(f, m, v, t, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_T2_0(f, t, x, peek, ...) \
    f<decltype(Value::x)>(t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_T2_1)(f, t, peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_T2_1(f, t, x, peek, ...) \
    f<decltype(Value::x)>(t) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_T2_0)(f, t, peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_GATHER_0(x, peek, ...) \
    enoki::gather<decltype(Value::x)>(src.x, index, mask) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_GATHER_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_GATHER_1(x, peek, ...) \
    enoki::gather<decltype(Value::x)>(src.x, index, mask) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_GATHER_0)(peek, __VA_ARGS__)

#define ENOKI_MAP_EXPR_SCATTER_0(x, peek, ...) \
    enoki::scatter(dst.x, value.x, index, mask) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_SCATTER_1)(peek, __VA_ARGS__)
#define ENOKI_MAP_EXPR_SCATTER_1(x, peek, ...) \
    enoki::scatter(dst.x, value.x, index, mask) ENOKI_MAP_EXPR_NEXT(peek, ENOKI_MAP_EXPR_SCATTER_0)(peek, __VA_ARGS__)

#define ENOKI_USING_MEMBERS_0(base, x, peek, ...)                               \
    using base::x;                                                           \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_MEMBERS_1)(base, peek, __VA_ARGS__)
#define ENOKI_USING_MEMBERS_1(base, x, peek, ...)                               \
    using base::x;                                                           \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_MEMBERS_0)(base, peek, __VA_ARGS__)
#define ENOKI_USING_MEMBERS_2(base, peek, ...) \
    ENOKI_EVAL(ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_MEMBERS_0)(base, peek, __VA_ARGS__))

#define ENOKI_USING_TYPES_0(base, x, peek, ...)                               \
    using x = typename base::x;                                             \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_TYPES_1)(base, peek, __VA_ARGS__)
#define ENOKI_USING_TYPES_1(base, x, peek, ...)                               \
    using x = typename base::x;                                             \
    ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_TYPES_0)(base, peek, __VA_ARGS__)
#define ENOKI_USING_TYPES_2(base, peek, ...) \
    ENOKI_EVAL(ENOKI_MAP_STMT_NEXT(peek, ENOKI_USING_TYPES_0)(base, peek, __VA_ARGS__))

// ENOKI_MAP_TEMPLATE_FWD(a1, a2, ...) expands to typename Ta1, typename Ta2, ...
#define ENOKI_MAP_TEMPLATE_FWD(...) \
    ENOKI_EVAL(ENOKI_MAP_TEMPLATE_FWD_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_DECL_FWD(a1, a2, ...) expands to Ta1 &&a1, Ta2&& a2...
#define ENOKI_MAP_EXPR_DECL_FWD(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_DECL_FWD_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_BASE_FWD(a1, a2, ...) expands to std::forward<Ta1>(a1), std::std::forward<Ta2>(a2), ...
#define ENOKI_MAP_EXPR_BASE_FWD(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_BASE_FWD_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_FWD(a1, a2, ...) expands to a1(std::forward<Ta1>(a1)), a2(std::std::forward<Ta2>(a2)), ...
#define ENOKI_MAP_EXPR_FWD(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_FWD_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_COPY(a1, a2, ...) expands to a1(a1), a2(a2), ...
#define ENOKI_MAP_EXPR_COPY(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_COPY_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_COPY_V(v, a1, a2, ...) expands to a1(v.a1), a2(v.a2), ...
#define ENOKI_MAP_EXPR_COPY_V(v, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_COPY_V_0(v, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_MOVE_V(v, a1, a2, ...) expands to a1(std::move(v.a1)), a2(std::move(v.a2)), ...
#define ENOKI_MAP_EXPR_MOVE_V(v, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_MOVE_V_0(v, __VA_ARGS__, (), 0))

// ENOKI_MAP_STMT_ASSIGN(v, a1, a2, ...) expands to this->a1 = v.a1; ..
#define ENOKI_MAP_STMT_ASSIGN(v, ...) \
    ENOKI_EVAL(ENOKI_MAP_STMT_ASSIGN_0(v, __VA_ARGS__, (), 0))

// ENOKI_MAP_STMT_MOVE(v, a1, a2, ...) expands to this->a1 = std::move(v.a1); ..
#define ENOKI_MAP_STMT_MOVE(v, ...) \
    ENOKI_EVAL(ENOKI_MAP_STMT_MOVE_0(v, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_F1(f, v, a1, a2, ...) expands to f(v.a1), f(v.a2), ...
#define ENOKI_MAP_EXPR_F1(f, v, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_F1_0(f, v, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_F2(f, v, t, a1, a2, ...) expands to f(v.a1, t), f(v.a2, t), ...
#define ENOKI_MAP_EXPR_F2(f, v, t, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_F2_0(f, v, t, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_T2(f, v, t, a1, a2, ...) expands to f<decltype(Value::a1)>(t), f<decltype(Value::a2>>(t), ...
#define ENOKI_MAP_EXPR_T2(f, v, t, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_T2_0(f, v, t, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_F3(f, m, v, t, a1, a2, ...) expands to f(m.a1, v.a1, t), f(m.a2, v.a2, t), ...
#define ENOKI_MAP_EXPR_F3(f, v, t, ...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_F3_0(f, v, t, __VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_GATHER(a1, a2, ...) expands to enoki::gather<decltype(Value::a1)>(src.a1, index, mask), ..
#define ENOKI_MAP_EXPR_GATHER(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_GATHER_0(__VA_ARGS__, (), 0))

// ENOKI_MAP_EXPR_SCATTER(a1, a2, ...) expands to enoki::scatter(dst.a1, src.a1, index, mask), ..
#define ENOKI_MAP_EXPR_SCATTER(...) \
    ENOKI_EVAL(ENOKI_MAP_EXPR_SCATTER_0(__VA_ARGS__, (), 0))

// ENOKI_USING_TYPES(base, a1, a2, ...) expands to using a1 = typename base::a1; using a2 = typename base::a2; ...
#define ENOKI_USING_TYPES(...) \
    ENOKI_EVAL_0(ENOKI_USING_TYPES_2(__VA_ARGS__, (), 0))

// ENOKI_USING_MEMBERS(base, a1, a2, ...) expands to using base::a1; using base::a2; ...
#define ENOKI_USING_MEMBERS(...) \
    ENOKI_EVAL_0(ENOKI_USING_MEMBERS_2(__VA_ARGS__, (), 0))


#define ENOKI_STRUCT(Struct, ...)                                              \
    Struct() = default;                                                        \
    template <ENOKI_MAP_TEMPLATE_FWD(__VA_ARGS__)>                             \
    ENOKI_INLINE Struct(ENOKI_MAP_EXPR_DECL_FWD(__VA_ARGS__))                  \
        : ENOKI_MAP_EXPR_FWD(__VA_ARGS__) { }                                  \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct(const Struct<Args...> &value)                          \
        : ENOKI_MAP_EXPR_COPY_V(value, __VA_ARGS__) { }                        \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct(Struct<Args...> &&value)                               \
        : ENOKI_MAP_EXPR_MOVE_V(value, __VA_ARGS__) { }                        \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct &operator=(const Struct<Args...> &value) {             \
        ENOKI_MAP_STMT_ASSIGN(value, __VA_ARGS__)                              \
        return *this;                                                          \
    }                                                                          \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct &operator=(Struct<Args...> &&value) {                  \
        ENOKI_MAP_STMT_MOVE(value, __VA_ARGS__)                                \
        return *this;                                                          \
    }

#define ENOKI_BASE_FIELDS(...) __VA_ARGS__
#define ENOKI_DERIVED_FIELDS(...) __VA_ARGS__

#define ENOKI_DERIVED_STRUCT(Struct, Base, BaseFields, StructFields)           \
    Struct() = default;                                                        \
    template <ENOKI_MAP_TEMPLATE_FWD(BaseFields),                              \
              ENOKI_MAP_TEMPLATE_FWD(StructFields)>                            \
    ENOKI_INLINE Struct(ENOKI_MAP_EXPR_DECL_FWD(BaseFields),                   \
           ENOKI_MAP_EXPR_DECL_FWD(StructFields))                              \
        : Base(ENOKI_MAP_EXPR_BASE_FWD(BaseFields)),                           \
          ENOKI_MAP_EXPR_FWD(StructFields) { }                                 \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct(const Struct<Args...> &value)                          \
        : Base(value), ENOKI_MAP_EXPR_COPY_V(value, StructFields) { }          \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct(Struct<Args...> &&value)                               \
        : Base(std::move(value)),                                              \
          ENOKI_MAP_EXPR_MOVE_V(value, StructFields) { }                       \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct &operator=(const Struct<Args...> &value) {             \
        Base::operator=(value);                                                \
        ENOKI_MAP_STMT_ASSIGN(value, StructFields)                             \
        return *this;                                                          \
    }                                                                          \
    template <typename... Args>                                                \
    ENOKI_INLINE Struct &operator=(Struct<Args...> &&value) {                  \
        Base::operator=(std::move(value));                                     \
        ENOKI_MAP_STMT_MOVE(value, StructFields)                               \
        return *this;                                                          \
    }                                                                          \
    template <typename Mask, enoki::enable_if_mask_t<Mask> = 0>                \
    auto operator[](const Mask &m) { return masked(*this, m); }                \


#define ENOKI_STRUCT_SUPPORT(Struct, ...)                                      \
    NAMESPACE_BEGIN(enoki)                                                     \
    template <typename... Args> struct struct_support<Struct<Args...>> {       \
        static constexpr bool IsDynamic =                                      \
            std::disjunction_v<enoki::is_dynamic<Args>...>;                    \
        using Dynamic = Struct<enoki::make_dynamic_t<Args>...>;                \
        using Value = Struct<Args...>;                                         \
        template <typename T, typename Arg>                                    \
        using ArgType =                                                        \
            std::conditional_t<std::is_const_v<std::remove_reference_t<T>>,    \
                               const Arg &, Arg &>;                            \
        static ENOKI_INLINE size_t packets(const Value &value) {               \
            return enoki::packets(                                             \
                value.ENOKI_EVAL_0(ENOKI_EXTRACT_0(__VA_ARGS__)));             \
        }                                                                      \
        static ENOKI_INLINE size_t slices(const Value &value) {                \
            return enoki::slices(                                              \
                value.ENOKI_EVAL_0(ENOKI_EXTRACT_0(__VA_ARGS__)));             \
        }                                                                      \
        static void set_slices(Value &value, size_t size) {                    \
            ENOKI_MAP_EXPR_F2(enoki::set_slices, value, size, __VA_ARGS__);    \
        }                                                                      \
        template <typename Mem, typename Mask>                                 \
        static ENOKI_INLINE size_t compress(Mem &mem, const Value &value,      \
                                            const Mask &mask) {                \
            return ENOKI_MAP_EXPR_F3(enoki::compress, mem, value,              \
                                     mask, __VA_ARGS__);                       \
        }                                                                      \
        template <typename Src, typename Index, typename Mask>                 \
        static ENOKI_INLINE Value gather(Src &src, const Index &index,         \
                                         const Mask &mask) {                   \
            return Value(ENOKI_MAP_EXPR_GATHER(__VA_ARGS__));                  \
        }                                                                      \
        template <typename Dst, typename Index, typename Mask>                 \
        static void scatter(Dst &dst, const Value &value, const Index &index,  \
                            const Mask &mask) {                                \
            ENOKI_MAP_EXPR_SCATTER(__VA_ARGS__);                               \
        }                                                                      \
        template <typename T>                                                  \
        static ENOKI_INLINE auto slice(T &&value, size_t index) {              \
            using Value = Struct<decltype(enoki::slice(std::declval<           \
                ArgType<T, Args>>(), index))...>;                              \
            return Value(ENOKI_MAP_EXPR_F2(enoki::slice, value, index,         \
                                           __VA_ARGS__));                      \
        }                                                                      \
        template <typename T>                                                  \
        static ENOKI_INLINE auto slice_ptr(T &&value, size_t index) {          \
            using Value = Struct<decltype(enoki::slice_ptr(std::declval<       \
                ArgType<T, Args>>(), index))...>;                              \
            return Value(ENOKI_MAP_EXPR_F2(enoki::slice_ptr, value, index,     \
                                           __VA_ARGS__));                      \
        }                                                                      \
        template <typename T>                                                  \
        static ENOKI_INLINE auto packet(T &&value, size_t index) {             \
            using Value = Struct<decltype(enoki::packet(std::declval<          \
                ArgType<T, Args>>(), index))...>;                              \
            return Value(ENOKI_MAP_EXPR_F2(enoki::packet, value, index,        \
                                           __VA_ARGS__));                      \
        }                                                                      \
        template <typename T> static ENOKI_INLINE auto ref_wrap(T &&value) {   \
            using Value = Struct<decltype(enoki::ref_wrap(std::declval<        \
                ArgType<T, Args>>()))...>;                                     \
            return Value(ENOKI_MAP_EXPR_F1(enoki::ref_wrap, value,             \
                                           __VA_ARGS__));                      \
        }                                                                      \
        template <typename T> static ENOKI_INLINE auto detach(T &&value) {     \
            using Value = Struct<decltype(enoki::detach(std::declval<          \
                ArgType<T, Args>>()))...>;                                     \
            return Value(ENOKI_MAP_EXPR_F1(enoki::detach, value,               \
                                           __VA_ARGS__));                      \
        }                                                                      \
        template <typename T, typename M> static ENOKI_INLINE                  \
        auto masked(T& value, const M & mask) {                                \
            using Value = Struct<decltype(enoki::masked(                       \
                        std::declval<Args &>(), mask))...>;                    \
            return Value(ENOKI_MAP_EXPR_F2(enoki::masked,                      \
                                           value, mask, __VA_ARGS__) );        \
        }                                                                      \
        static ENOKI_INLINE auto zero(size_t size) {                           \
            return Value(ENOKI_EVAL_0(                                         \
                ENOKI_MAP_EXPR_T2(enoki::zero, size, __VA_ARGS__)));           \
        }                                                                      \
        static ENOKI_INLINE auto empty(size_t size) {                          \
            return Value(ENOKI_EVAL_0(                                         \
                ENOKI_MAP_EXPR_T2(enoki::empty, size, __VA_ARGS__)));          \
        }                                                                      \
    };                                                                         \
    NAMESPACE_END(enoki)

#define ENOKI_PINNED_OPERATOR_NEW(Type)                                        \
    void *operator new(size_t size) {                                          \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            return enoki::cuda_host_malloc(size);                              \
        else                                                                   \
            return ::operator new(size);                                       \
    }                                                                          \
    void *operator new(size_t size, std::align_val_t align) {                  \
        ENOKI_MARK_USED(align);                                                \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            return enoki::cuda_host_malloc(size);                              \
        else                                                                   \
            return ::operator new(size, align);                                \
    }                                                                          \
    void *operator new[](size_t size) {                                        \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            return enoki::cuda_host_malloc(size);                              \
        else                                                                   \
            return ::operator new[](size);                                     \
    }                                                                          \
                                                                               \
    void *operator new[](size_t size, std::align_val_t align) {                \
        ENOKI_MARK_USED(align);                                                \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            return enoki::cuda_host_malloc(size);                              \
        else                                                                   \
            return ::operator new[](size, align);                              \
    }                                                                          \
                                                                               \
    void operator delete(void *ptr) {                                          \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            enoki::cuda_host_free(ptr);                                        \
        else                                                                   \
            return ::operator delete(ptr);                                     \
    }                                                                          \
                                                                               \
    void operator delete(void *ptr, std::align_val_t align) {                  \
        ENOKI_MARK_USED(align);                                                \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            enoki::cuda_host_free(ptr);                                        \
        else                                                                   \
            return ::operator delete(ptr, align);                              \
    }                                                                          \
                                                                               \
    void operator delete[](void *ptr) {                                        \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            enoki::cuda_host_free(ptr);                                        \
        else                                                                   \
            return ::operator delete[](ptr);                                   \
    }                                                                          \
                                                                               \
    void operator delete[](void *ptr, std::align_val_t align) {                \
        ENOKI_MARK_USED(align);                                                \
        if constexpr (enoki::is_cuda_array_v<Type>)                            \
            enoki::cuda_host_free(ptr);                                        \
        else                                                                   \
            return ::operator delete[](ptr, align);                            \
    }

