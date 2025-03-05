/*
    enoki/array_router.h -- Helper functions which route function calls
    in the enoki namespace to the intended recipients

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include "array_traits.h"
#include "array_fallbacks.h"

NAMESPACE_BEGIN(enoki)

/// Define an unary operation
#define ENOKI_ROUTE_UNARY(name, func)                                          \
    template <typename T, enable_if_array_t<T> = 0>                            \
    ENOKI_INLINE auto name(const T &a) {                                       \
        return eval(a).func##_();                                              \
    }

/// Define an unary operation with an immediate argument (e.g. sr<5>(...))
#define ENOKI_ROUTE_UNARY_IMM(name, func)                                      \
    template <size_t Imm, typename T, enable_if_array_t<T> = 0>                \
    ENOKI_INLINE auto name(const T &a) {                                       \
        return eval(a).template func##_<Imm>(); /* Forward to array */         \
    }

/// Define an unary operation with a fallback expression for scalar arguments
#define ENOKI_ROUTE_UNARY_SCALAR(name, func, expr)                             \
    template <typename T> ENOKI_INLINE auto name(const T &a) {                 \
        if constexpr (!is_array_v<T>)                                          \
            return expr; /* Scalar fallback implementation */                  \
        else                                                                   \
            return eval(a).func##_(); /* Forward to array */                   \
    }

/// Define an unary operation with an immediate argument and a scalar fallback
#define ENOKI_ROUTE_UNARY_SCALAR_IMM(name, func, expr)                         \
    template <size_t Imm, typename T> ENOKI_INLINE auto name(const T &a) {     \
        if constexpr (!is_array_v<T>)                                          \
            return expr; /* Scalar fallback implementation */                  \
        else                                                                   \
            return eval(a).template func##_<Imm>(); /* Forward to array */     \
    }

/// Define a binary operation
#define ENOKI_ROUTE_BINARY(name, func)                                         \
    template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>     \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2) {                       \
        using E = expr_t<T1, T2>;                                              \
        if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)          \
            return a1.derived().func##_(a2.derived());                         \
        else                                                                   \
            return name(static_cast<const E &>(a1),                            \
                        static_cast<const E &>(a2));                           \
    }

/// Define a binary operation for bit operations
#define ENOKI_ROUTE_BINARY_BITOP(name, func)                                   \
    template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>     \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2) {                       \
        using E = expr_t<T1, T2>;                                              \
        if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)          \
            return a1.derived().func##_(a2.derived());                         \
        else if constexpr (is_mask_v<T2> && !is_array_v<T2>)                   \
            return a1.derived().func##_((const mask_t<T1> &) a2);              \
        else if constexpr (is_array_v<T2>) {                                   \
            if constexpr (std::decay_t<T2>::IsMask)                            \
                return a1.derived().func##_((const mask_t<T1> &) a2.derived());\
            else                                                               \
                return name(static_cast<const E &>(a1),                        \
                            static_cast<const E &>(a2));                       \
        } else {                                                               \
            return name(static_cast<const E &>(a1),                            \
                        static_cast<const E &>(a2));                           \
        }                                                                      \
    }

/// Define a binary operation (but only restrict to cases where 'cond' is true)
#define ENOKI_ROUTE_BINARY_COND(name, func, cond)                              \
    template <typename T1, typename T2,                                        \
              enable_if_t<cond> = 0,                                           \
              enable_if_array_any_t<T1, T2> = 0>                               \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2) {                       \
        using E = expr_t<T1, T2>;                                              \
        if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)          \
            return a1.derived().func##_(a2.derived());                         \
        else                                                                   \
            return name(static_cast<const E &>(a1),                            \
                        static_cast<const E &>(a2));                           \
    }

#define ENOKI_ROUTE_BINARY_SHIFT(name, func)                                   \
    template <typename T1, typename T2,                                        \
              enable_if_t<std::is_arithmetic_v<scalar_t<T1>>> = 0,             \
              enable_if_array_any_t<T1, T2> = 0>                               \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2) {                       \
        using E = expr_t<T1, T2>;                                              \
        if constexpr (std::is_integral_v<T2>)                                  \
            return eval(a1).func##_((size_t) a2);                              \
        else if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)     \
            return a1.derived().func##_(a2.derived());                         \
        else                                                                   \
            return name(static_cast<const E &>(a1),                            \
                        static_cast<const E &>(a2));                           \
    }

/// Define a binary operation with a fallback expression for scalar arguments
#define ENOKI_ROUTE_BINARY_SCALAR(name, func, expr)                            \
    template <typename T1, typename T2>                                        \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2) {                       \
        using E = expr_t<T1, T2>;                                              \
        if constexpr (is_array_any_v<T1, T2>) {                                \
            if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)      \
                return a1.derived().func##_(a2.derived());                     \
            else                                                               \
                return name(static_cast<const E &>(a1),                        \
                            static_cast<const E &>(a2));                       \
        } else {                                                               \
            return expr;                                                       \
        }                                                                      \
    }

/// Define a ternary operation
#define ENOKI_ROUTE_TERNARY_SCALAR(name, func, expr)                           \
    template <typename T1, typename T2, typename T3>                           \
    ENOKI_INLINE auto name(const T1 &a1, const T2 &a2, const T3 &a3) {         \
        using E = expr_t<T1, T2, T3>;                                          \
        if constexpr (is_array_any_v<T1, T2, T3>) {                            \
            if constexpr (std::is_same_v<T1, E> &&                             \
                          std::is_same_v<T2, E> &&                             \
                          std::is_same_v<T3, E>)                               \
                return a1.derived().func##_(a2.derived(), a3.derived());       \
            else                                                               \
                return name(static_cast<const E &>(a1),                        \
                            static_cast<const E &>(a2),                        \
                            static_cast<const E &>(a3));                       \
        } else {                                                               \
            return expr;                                                       \
        }                                                                      \
    }

/// Macro for compound assignment operators (operator+=, etc.)
#define ENOKI_ROUTE_COMPOUND_OPERATOR(op)                                      \
    template <typename T1, enable_if_t<is_array_v<T1> &&                       \
                                      !std::is_const_v<T1>> = 0, typename T2>  \
    ENOKI_INLINE T1 &operator op##=(T1 &a1, const T2 &a2) {                    \
        a1 = a1 op a2;                                                         \
        return a1;                                                             \
    }

template <typename T, enable_if_array_t<T> = 0>
ENOKI_INLINE decltype(auto) eval(const T& x) {
    if constexpr (std::is_same_v<std::decay_t<T>, expr_t<T>>)
        return x.derived();
    else
        return expr_t<T>(x);
}

ENOKI_ROUTE_UNARY(operator-, neg)
ENOKI_ROUTE_UNARY(operator~, not)
ENOKI_ROUTE_UNARY(operator!, not)

ENOKI_ROUTE_BINARY_COND(operator+, add, !std::is_pointer_v<scalar_t<T1>> && !std::is_pointer_v<scalar_t<T2>>)
ENOKI_ROUTE_BINARY_COND(operator-, sub, !std::is_pointer_v<scalar_t<T1>> && !std::is_pointer_v<scalar_t<T2>>)
ENOKI_ROUTE_BINARY(operator*, mul)

ENOKI_ROUTE_BINARY_SHIFT(operator<<, sl)
ENOKI_ROUTE_BINARY_SHIFT(operator>>, sr)

ENOKI_ROUTE_UNARY_SCALAR_IMM(sl, sl, a << Imm)
ENOKI_ROUTE_UNARY_SCALAR_IMM(sr, sr, a >> Imm)

ENOKI_ROUTE_BINARY_BITOP(operator&,  and)
ENOKI_ROUTE_BINARY_BITOP(operator&&, and)
ENOKI_ROUTE_BINARY_BITOP(operator|,  or)
ENOKI_ROUTE_BINARY_BITOP(operator||, or)
ENOKI_ROUTE_BINARY_BITOP(operator^,  xor)
ENOKI_ROUTE_BINARY_SCALAR(andnot, andnot, a1 & !a2)

ENOKI_ROUTE_BINARY(operator<,  lt)
ENOKI_ROUTE_BINARY(operator<=, le)
ENOKI_ROUTE_BINARY(operator>,  gt)
ENOKI_ROUTE_BINARY(operator>=, ge)

ENOKI_ROUTE_BINARY_SCALAR(eq,  eq,  a1 == a2)
ENOKI_ROUTE_BINARY_SCALAR(neq, neq, a1 != a2)

ENOKI_ROUTE_COMPOUND_OPERATOR(+)
ENOKI_ROUTE_COMPOUND_OPERATOR(-)
ENOKI_ROUTE_COMPOUND_OPERATOR(*)
ENOKI_ROUTE_COMPOUND_OPERATOR(/)
ENOKI_ROUTE_COMPOUND_OPERATOR(^)
ENOKI_ROUTE_COMPOUND_OPERATOR(|)
ENOKI_ROUTE_COMPOUND_OPERATOR(&)
ENOKI_ROUTE_COMPOUND_OPERATOR(<<)
ENOKI_ROUTE_COMPOUND_OPERATOR(>>)

ENOKI_ROUTE_BINARY_SCALAR(max,   max,  (std::decay_t<E>) std::max((E) a1, (E) a2))
ENOKI_ROUTE_BINARY_SCALAR(min,   min,  (std::decay_t<E>) std::min((E) a1, (E) a2))

ENOKI_ROUTE_BINARY_SCALAR(dot,   dot,   (E) a1 * (E) a2)

ENOKI_ROUTE_BINARY_SCALAR(mulhi, mulhi, detail::mulhi_scalar(a1, a2))

ENOKI_ROUTE_UNARY_SCALAR(abs, abs, detail::abs_scalar(a))

ENOKI_ROUTE_TERNARY_SCALAR(fmadd,  fmadd,  detail::fmadd_scalar((E)  a1, (E) a2, (E)  a3))
ENOKI_ROUTE_TERNARY_SCALAR(fmsub,  fmsub,  detail::fmadd_scalar((E)  a1, (E) a2, (E) -a3))
ENOKI_ROUTE_TERNARY_SCALAR(fnmadd, fnmadd, detail::fmadd_scalar((E) -a1, (E) a2, (E)  a3))
ENOKI_ROUTE_TERNARY_SCALAR(fnmsub, fnmsub, detail::fmadd_scalar((E) -a1, (E) a2, (E) -a3))
ENOKI_ROUTE_TERNARY_SCALAR(fmaddsub, fmaddsub, fmsub(a1, a2, a3))
ENOKI_ROUTE_TERNARY_SCALAR(fmsubadd, fmsubadd, fmadd(a1, a2, a3))

ENOKI_ROUTE_UNARY_SCALAR(rcp, rcp, 1 / a)
ENOKI_ROUTE_UNARY_SCALAR(rsqrt, rsqrt, detail::rsqrt_scalar(a))

ENOKI_ROUTE_UNARY_SCALAR(popcnt, popcnt, detail::popcnt_scalar(a))
ENOKI_ROUTE_UNARY_SCALAR(lzcnt, lzcnt, detail::lzcnt_scalar(a))
ENOKI_ROUTE_UNARY_SCALAR(tzcnt, tzcnt, detail::tzcnt_scalar(a))

ENOKI_ROUTE_UNARY_SCALAR(all,   all,   (bool) a)
ENOKI_ROUTE_UNARY_SCALAR(any,   any,   (bool) a)
ENOKI_ROUTE_UNARY_SCALAR(count, count, (size_t) ((bool) a ? 1 : 0))
ENOKI_ROUTE_UNARY_SCALAR(reverse, reverse, a)
ENOKI_ROUTE_UNARY_SCALAR(psum,  psum,  a)
ENOKI_ROUTE_UNARY_SCALAR(hsum,  hsum,  a)
ENOKI_ROUTE_UNARY_SCALAR(hprod, hprod, a)
ENOKI_ROUTE_UNARY_SCALAR(hmin,  hmin,  a)
ENOKI_ROUTE_UNARY_SCALAR(hmax,  hmax,  a)
ENOKI_ROUTE_UNARY_SCALAR(hmean, hmean,  a)

ENOKI_ROUTE_UNARY_SCALAR(all_inner,   all_inner,   (bool) a)
ENOKI_ROUTE_UNARY_SCALAR(any_inner,   any_inner,   (bool) a)
ENOKI_ROUTE_UNARY_SCALAR(count_inner, count_inner, (size_t) ((bool) a ? 1 : 0))
ENOKI_ROUTE_UNARY_SCALAR(psum_inner,  psum_inner,  a)
ENOKI_ROUTE_UNARY_SCALAR(hsum_inner,  hsum_inner,  a)
ENOKI_ROUTE_UNARY_SCALAR(hprod_inner, hprod_inner, a)
ENOKI_ROUTE_UNARY_SCALAR(hmin_inner,  hmin_inner,  a)
ENOKI_ROUTE_UNARY_SCALAR(hmax_inner,  hmax_inner,  a)
ENOKI_ROUTE_UNARY_SCALAR(hmean_inner, hmean_inner,  a)

ENOKI_ROUTE_UNARY_SCALAR(sqrt,  sqrt,  std::sqrt(a))
ENOKI_ROUTE_UNARY_SCALAR(floor, floor, std::floor(a))
ENOKI_ROUTE_UNARY_SCALAR(ceil,  ceil,  std::ceil(a))
ENOKI_ROUTE_UNARY_SCALAR(round, round, std::rint(a))
ENOKI_ROUTE_UNARY_SCALAR(trunc, trunc, std::trunc(a))

ENOKI_ROUTE_UNARY_IMM(rol_array, rol_array)
ENOKI_ROUTE_UNARY_IMM(ror_array, ror_array)

template <typename T> auto none(const T &value) {
    return !any(value);
}

template <typename T> auto none_inner(const T &value) {
    return !any_inner(value);
}

/// Floating point division
template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0,
          enable_if_t<std::is_floating_point_v<scalar_t<expr_t<T1, T2>>>> = 0>
ENOKI_INLINE auto operator/(const T1 &a1, const T2 &a2) {
    using E = expr_t<T1, T2>;
    using T = expr_t<scalar_t<T1>, T2>;

    if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)
        return a1.derived().div_(a2.derived());
    else if constexpr (array_depth_v<T1> > array_depth_v<T2>)
        return static_cast<const E &>(a1) * // reciprocal approximation
               rcp((const T &) a2);
    else
        return operator/(static_cast<const E &>(a1),
                         static_cast<const E &>(a2));
}

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0,
          enable_if_t<!std::is_floating_point_v<scalar_t<expr_t<T1, T2>>> &&
                       is_array_v<T2>> = 0>
ENOKI_INLINE auto operator/(const T1 &a1, const T2 &a2) {
    using E = expr_t<T1, T2>;

    if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)
        return a1.derived().div_(a2.derived());
    else
        return operator/(static_cast<const E &>(a1),
                         static_cast<const E &>(a2));
}

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0,
          enable_if_t<!std::is_floating_point_v<scalar_t<expr_t<T1, T2>>> &&
                       is_array_v<T2>> = 0>
ENOKI_INLINE auto operator%(const T1 &a1, const T2 &a2) {
    using E = expr_t<T1, T2>;

    if constexpr (std::is_same_v<T1, E> && std::is_same_v<T2, E>)
        return a1.derived().mod_(a2.derived());
    else
        return operator%(static_cast<const E &>(a1),
                         static_cast<const E &>(a2));
}

/// Shuffle the entries of an array
template <size_t... Is, typename T>
ENOKI_INLINE auto shuffle(const T &a) {
    if constexpr (is_array_v<T>) {
        return eval(a).template shuffle_<Is...>();
    } else {
        static_assert(sizeof...(Is) == 1 && (... && (Is == 0)), "Shuffle argument out of bounds!");
        return a;
    }
}

template <typename Array, typename Index,
          enable_if_t<is_array_v<Array> && is_array_v<Index> && std::is_integral_v<scalar_t<Index>>> = 0>
ENOKI_INLINE Array shuffle(const Array &a, const Index &idx) {
    if constexpr (Index::Depth > Array::Depth) {
        Array result;
        for (size_t i = 0; i < Array::Size; ++i)
            result.coeff(i) = shuffle(a.derived().coeff(i), idx);
        return result;
    } else {
        return a.derived().shuffle_((int_array_t<Array> &) idx);
    }
}

//// Compute the square of the given value
template <typename T> ENOKI_INLINE auto sqr(const T &value) {
    return value * value;
}

//// Convert radians to degrees
template <typename T> ENOKI_INLINE auto rad_to_deg(const T &a) {
    return a * scalar_t<T>(180 / M_PI);
}

/// Convert degrees to radians
template <typename T> ENOKI_INLINE auto deg_to_rad(const T &a) {
    return a * scalar_t<T>(M_PI / 180);
}

template <typename T> ENOKI_INLINE auto sign_mask() {
    using Scalar = scalar_t<T>;
    using UInt = uint_array_t<Scalar>;
    return memcpy_cast<Scalar>(UInt(1) << (sizeof(UInt) * 8 - 1));
}

template <typename T, typename Expr = expr_t<T>>
ENOKI_INLINE Expr sign(const T &a) {
    using Scalar = scalar_t<T>;

    if constexpr (array_depth_v<Expr> >= 2) {
        Expr result;
        for (size_t i = 0; i < Expr::Size; ++i)
            result.coeff(i) = sign(a.coeff(i));
        return result;
    } else if constexpr (!std::is_signed_v<Scalar>) {
        return Expr(Scalar(1));
    } else if constexpr (!std::is_floating_point_v<Scalar> || is_diff_array_v<Expr>) {
        return select(a < Scalar(0), Expr(Scalar(-1)), Expr(Scalar(1)));
    } else if constexpr (is_scalar_v<Expr>) {
        return std::copysign(Scalar(1), a);
    } else {
        return (sign_mask<T>() & a) | Expr(Scalar(1));
    }
}

template <typename T1, typename T2, typename Expr = expr_t<T1, T2>>
ENOKI_INLINE Expr copysign(const T1 &a1, const T2 &a2) {
    using Scalar1 = scalar_t<T1>;
    using Scalar2 = scalar_t<T2>;

    static_assert(std::is_same_v<Scalar1, Scalar2> || !std::is_signed_v<Scalar1>,
                  "copysign(): Incompatible input arguments!");

    if constexpr (!std::is_same_v<T1, Expr> || !std::is_same_v<T2, Expr>) {
        return copysign((const Expr &) a1, (const Expr &) a2);
    } else if constexpr (array_depth_v<Expr> >= 2) {
        Expr result;
        for (size_t i = 0; i < Expr::Size; ++i)
            result.coeff(i) = copysign(a1.coeff(i), a2.coeff(i));
        return result;
    } else if constexpr (!std::is_floating_point_v<Scalar1>) {
        return select((a1 ^ a2) < Scalar1(0), a1, -a1);
    } else if constexpr (is_scalar_v<Expr>) {
        return std::copysign(a1, a2);
    } else if constexpr (is_diff_array_v<Expr>) {
        return abs(a1) * sign(a2);
    } else {
        return abs(a1) | (sign_mask<Expr>() & a2);
    }
}

template <typename T1, typename T2, typename Expr = expr_t<T1, T2>>
ENOKI_INLINE Expr copysign_neg(const T1 &a1, const T2 &a2) {
    using Scalar1 = scalar_t<T1>;
    using Scalar2 = scalar_t<T2>;

    static_assert(std::is_same_v<Scalar1, Scalar2> || !std::is_signed_v<Scalar1>,
                  "copysign_neg(): Incompatible input arguments!");

    if constexpr (!std::is_same_v<T1, Expr> || !std::is_same_v<T2, Expr>) {
        return copysign_neg((const Expr &) a1, (const Expr &) a2);
    } else if constexpr (array_depth_v<Expr> >= 2) {
        Expr result;
        for (size_t i = 0; i < Expr::Size; ++i)
            result.coeff(i) = copysign_neg(a1.coeff(i), a2.coeff(i));
        return result;
    } else if constexpr (!std::is_floating_point_v<Scalar1>) {
        return select((a1 ^ a2) < Scalar1(0), -a1, a1);
    } else if constexpr (is_scalar_v<Expr>) {
        return std::copysign(a1, -a2);
    } else if constexpr (is_diff_array_v<Expr>) {
        return abs(a1) * -sign(a2);
    } else {
        return abs(a1) | andnot(sign_mask<Expr>(), a2);
    }
}

template <typename T1, typename T2, typename Expr = expr_t<T1, T2>>
ENOKI_INLINE Expr mulsign(const T1 &a1, const T2 &a2) {
    using Scalar1 = scalar_t<T1>;
    using Scalar2 = scalar_t<T2>;

    static_assert(std::is_same_v<Scalar1, Scalar2> || !std::is_signed_v<Scalar1>,
                  "mulsign(): Incompatible input arguments!");

    if constexpr (!std::is_same_v<T1, Expr> || !std::is_same_v<T2, Expr>) {
        return mulsign((const Expr &) a1, (const Expr &) a2);
    } else if constexpr (array_depth_v<Expr> >= 2) {
        Expr result;
        for (size_t i = 0; i < Expr::Size; ++i)
            result.coeff(i) = mulsign(a1.coeff(i), a2.coeff(i));
        return result;
    } else if constexpr (!std::is_floating_point_v<Scalar1>) {
        return select(a2 < Scalar1(0), -a1, a1);
    } else if constexpr (is_scalar_v<Expr>) {
        return a1 * std::copysign(Scalar1(1), a2);
    } else if constexpr (is_diff_array_v<Expr>) {
        return a1 * sign(a2);
    } else {
        return a1 ^ (sign_mask<Expr>() & a2);
    }
}

template <typename T1, typename T2, typename Expr = expr_t<T1, T2>>
ENOKI_INLINE Expr mulsign_neg(const T1 &a1, const T2 &a2) {
    using Scalar1 = scalar_t<T1>;
    using Scalar2 = scalar_t<T2>;

    static_assert(std::is_same_v<Scalar1, Scalar2> || !std::is_signed_v<Scalar1>,
                  "mulsign_neg(): Incompatible input arguments!");

    if constexpr (!std::is_same_v<T1, Expr> || !std::is_same_v<T2, Expr>) {
        return mulsign_neg((const Expr &) a1, (const Expr &) a2);
    } else if constexpr (array_depth_v<Expr> >= 2) {
        Expr result;
        for (size_t i = 0; i < Expr::Size; ++i)
            result.coeff(i) = mulsign_neg(a1.coeff(i), a2.coeff(i));
        return result;
    } else if constexpr (!std::is_floating_point_v<Scalar1>) {
        return select(a2 < Scalar1(0), a1, -a1);
    } else if constexpr (is_scalar_v<Expr>) {
        return a1 * std::copysign(Scalar1(1), -a2);
    } else if constexpr (is_diff_array_v<Expr>) {
        return a1 * -sign(a2);
    } else {
        return a1 ^ andnot(sign_mask<Expr>(), a2);
    }
}

template <typename M, typename T, typename F>
ENOKI_INLINE auto select(const M &m, const T &t, const F &f) {
    using E = expr_t<T, F>;

    if constexpr (!is_array_v<E>)
        return (bool) m ? (E) t : (E) f;
    else if constexpr (std::is_same_v<M, mask_t<E>> &&
                       std::is_same_v<T, E> &&
                       std::is_same_v<F, E>)
        return E::select_(m.derived(), t.derived(), f.derived());
    else
        return select((const mask_t<E> &) m, (const E &) t, (const E &) f);
}

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
ENOKI_INLINE bool operator==(const T1 &a1, const T2 &a2) {
    return all_nested(eq(a1, a2));
}

template <typename T1, typename T2, enable_if_array_any_t<T1, T2> = 0>
ENOKI_INLINE bool operator!=(const T1 &a1, const T2 &a2) {
    return any_nested(neq(a1, a2));
}

namespace detail {
    template <typename T>
    using has_ror = decltype(std::declval<T>().template ror_<0>());
    template <typename T>
    constexpr bool has_ror_v = is_detected_v<has_ror, T>;
}

/// Bit-level rotate left (with immediate offset value)
template <size_t Imm, typename T>
ENOKI_INLINE auto rol(const T &a) {
    constexpr size_t Mask = 8 * sizeof(scalar_t<T>) - 1u;
    using UInt = uint_array_t<T>;

    if constexpr (detail::has_ror_v<T>)
        return a.template rol_<Imm>();
    else
        return sl<Imm & Mask>(a) | T(sr<((~Imm + 1u) & Mask)>(UInt(a)));
}

/// Bit-level rotate right (with immediate offset value)
template <typename T1, typename T2>
ENOKI_INLINE auto rol(const T1 &a1, const T2 &a2) {
    if constexpr (detail::has_ror_v<T1>) {
        return a1.rol_(a2);
    } else {
        using U1 = uint_array_t<T1>;
        using U2 = uint_array_t<T2>;
        using Expr = expr_t<T1, T2>;
        constexpr scalar_t<U2> Mask = 8 * sizeof(scalar_t<Expr>) - 1u;

        U1 u1 = (U1) a1; U2 u2 = (U2) a2;
        return Expr((u1 << u2) | (u1 >> ((~u2 + 1u) & Mask)));
    }
}

/// Bit-level rotate right (with scalar or array offset value)
template <size_t Imm, typename T>
ENOKI_INLINE T ror(const T &a) {
    constexpr size_t Mask = 8 * sizeof(scalar_t<T>) - 1u;
    using UInt = uint_array_t<T>;

    if constexpr (detail::has_ror_v<T>)
        return a.template ror_<Imm>();
    else
        return T(sr<Imm & Mask>(UInt(a))) | sl<((~Imm + 1u) & Mask)>(a);
}

/// Bit-level rotate right (with scalar or array offset value)
template <typename T1, typename T2>
ENOKI_INLINE auto ror(const T1 &a1, const T2 &a2) {
    if constexpr (detail::has_ror_v<T1>) {
        return a1.ror_(a2);
    } else {
        using U1 = uint_array_t<T1>;
        using U2 = uint_array_t<T2>;
        using Expr = expr_t<T1, T2>;
        constexpr scalar_t<U2> Mask = 8 * sizeof(scalar_t<Expr>) - 1u;

        U1 u1 = (U1) a1; U2 u2 = (U2) a2;
        return Expr((u1 >> u2) | (u1 << ((~u2 + 1u) & Mask)));
    }
}

/// Fast implementation for computing the base 2 log of an integer.
template <typename T> ENOKI_INLINE auto log2i(T value) {
    return scalar_t<T>(sizeof(scalar_t<T>) * 8 - 1) - lzcnt(value);
}

template <typename T> struct MaskBit {
    MaskBit(T &mask, size_t index) : mask(mask), index(index) { }
    operator bool() const { return mask.bit_(index); }
    MaskBit &operator=(bool b) { mask.set_bit_(index, b); return *this; }
private:
    T mask;
    size_t index;
};

template <typename Target, typename Source>
ENOKI_INLINE Target reinterpret_array(const Source &src) {
    if constexpr (std::is_same_v<Source, Target>) {
        return src;
    } else if constexpr (std::is_constructible_v<Target, const Source &, detail::reinterpret_flag>) {
        return Target(src, detail::reinterpret_flag());
    } else if constexpr (is_scalar_v<Source> && is_scalar_v<Target>) {
        if constexpr (sizeof(Source) == sizeof(Target)) {
            return memcpy_cast<Target>(src);
        } else {
            using SrcInt = int_array_t<Source>;
            using TrgInt = int_array_t<Target>;
            if constexpr (std::is_same_v<Target, bool>)
                return memcpy_cast<SrcInt>(src) != 0 ? true : false;
            else
                return memcpy_cast<Target>(memcpy_cast<SrcInt>(src) != 0 ? TrgInt(-1) : TrgInt(0));
        }
    } else {
        static_assert(detail::false_v<Source, Target>, "reinterpret_array(): don't know what to do!");
    }
}

template <typename Target, typename T>
ENOKI_INLINE Target reinterpret_array(const MaskBit<T> &src) {
    return reinterpret_array<Target>((bool) src);
}

/// Element-wise test for NaN values
template <typename T>
ENOKI_INLINE auto isnan(const T &a) { return !eq(a, a); }

/// Element-wise test for +/- infinity
template <typename T>
ENOKI_INLINE auto isinf(const T &a) {
    return eq(abs(a), std::numeric_limits<scalar_t<T>>::infinity());
}

/// Element-wise test for finiteness
template <typename T>
ENOKI_INLINE auto isfinite(const T &a) {
    return abs(a) < std::numeric_limits<scalar_t<T>>::infinity();
}

/// Extract the low elements from an array of even size
template <typename Array, enable_if_t<(Array::Size > 1 && Array::Size != -1)> = 0>
auto low(const Array &a) { return a.derived().low_(); }

/// Extract the high elements from an array of even size
template <typename Array, enable_if_t<(Array::Size > 1 && Array::Size != -1)> = 0>
auto high(const Array &a) { return a.derived().high_(); }

template <typename T, typename Arg>
T floor2int(const Arg &a) {
    if constexpr (is_array_v<Arg>)
        return a.template floor2int_<T>();
    else
        return detail::floor2int_scalar<T>(a);
}

template <typename T, typename Arg>
T ceil2int(const Arg &a) {
    if constexpr (is_array_v<Arg>)
        return a.template ceil2int_<T>();
    else
        return detail::ceil2int_scalar<T>(a);
}

// -----------------------------------------------------------------------
//! @{ \name Miscellaneous routines for vector spaces
// -----------------------------------------------------------------------

template <typename T1, typename T2>
ENOKI_INLINE auto abs_dot(const T1 &a1, const T2 &a2) {
    return abs(dot(a1, a2));
}

template <typename T> ENOKI_INLINE auto norm(const T &v) {
    return sqrt(dot(v, v));
}

template <typename T> ENOKI_INLINE auto squared_norm(const T &v) {
    return dot(v, v);
}

template <typename T> ENOKI_INLINE auto normalize(const T &v) {
    return v * rsqrt(squared_norm(v));
}

template <typename T, enable_if_t<is_dynamic_array_v<T>> = 0>
ENOKI_INLINE auto partition(const T &v) {
    return v.partition_();
}

template <typename T1, typename T2,
          enable_if_t<array_size_v<T1> == 3 &&
                      array_size_v<T2> == 3> = 0>
ENOKI_INLINE auto cross(const T1 &v1, const T2 &v2) {
#if defined(ENOKI_ARM_32) || defined(ENOKI_ARM_64)
    return fnmadd(
        shuffle<2, 0, 1>(v1), shuffle<1, 2, 0>(v2),
        shuffle<1, 2, 0>(v1) * shuffle<2, 0, 1>(v2)
    );
#else
    return fmsub(shuffle<1, 2, 0>(v1),  shuffle<2, 0, 1>(v2),
                 shuffle<2, 0, 1>(v1) * shuffle<1, 2, 0>(v2));
#endif
}

template <typename T> decltype(auto) detach(T &value) {
    if constexpr (is_array_v<T>) {
        if constexpr (!is_diff_array_v<T>)
            return value;
        else if constexpr (array_depth_v<T> == 1)
            return value.value_();
        else
            return struct_support_t<T>::detach(value);
    } else {
        return struct_support_t<T>::detach(value);
    }
}

template <typename T> decltype(auto) gradient(T &&value) {
    if constexpr (is_array_v<T>) {
        if constexpr (!is_diff_array_v<T>)
            return value;
        else if constexpr (array_depth_v<T> == 1)
            return value.gradient_();
        else
            return struct_support_t<T>::gradient(value);
    } else {
        return struct_support_t<T>::gradient(value);
    }
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Initialization, loading/writing data
// -----------------------------------------------------------------------

template <typename T> ENOKI_INLINE T zero(size_t size = 1);
template <typename T> ENOKI_INLINE T empty(size_t size = 1);

/// Construct an index sequence, i.e. 0, 1, 2, ..
template <typename Array, enable_if_dynamic_array_t<Array> = 0>
ENOKI_INLINE Array arange(size_t end = 1) {
    return Array::arange_(0, (ssize_t) end, 1);
}

template <typename Array, enable_if_static_array_t<Array> = 0>
ENOKI_INLINE Array arange(size_t end = Array::Size) {
    assert(end == Array::Size);
    (void) end;
    return Array::arange_(0, (ssize_t) Array::Size, 1);
}

template <typename Arg, enable_if_not_array_t<Arg> = 0>
ENOKI_INLINE Arg arange(size_t end = 1) {
    assert(end == 1);
    (void) end;
    return Arg(0);
}

template <typename T>
ENOKI_INLINE T arange(ssize_t start, ssize_t end, ssize_t step = 1) {
    if constexpr (is_static_array_v<T>) {
        assert(end - start == (ssize_t) T::Size * step);
        return T::arange_(start, end, step);
    } else if constexpr (is_dynamic_array_v<T>) {
        return T::arange_(start, end, step);
    } else {
        assert(end - start == step);
        (void) end;
        (void) step;
        return T(start);
    }
}

/// Construct an array that linearly interpolates from min..max
template <typename Array, enable_if_dynamic_array_t<Array> = 0>
ENOKI_INLINE Array linspace(scalar_t<Array> min, scalar_t<Array> max, size_t size = 1) {
    return Array::linspace_(min, max, size);
}

template <typename Array, enable_if_static_array_t<Array> = 0>
ENOKI_INLINE Array linspace(scalar_t<Array> min, scalar_t<Array> max, size_t size = Array::Size) {
    assert(size == Array::Size);
    (void) size;
    return Array::linspace_(min, max);
}

/// Construct an array that linearly interpolates from min..max (scalar fallback)
template <typename Arg, enable_if_not_array_t<Arg> = 0>
ENOKI_INLINE Arg linspace(scalar_t<Arg> min, scalar_t<Arg>, size_t size = 1) {
    assert(size == 1);
    (void) size;
    return min;
}

template <typename Outer, typename Inner,
          typename Return = replace_scalar_t<Outer, Inner>>
ENOKI_INLINE Return full(const Inner &inner, size_t size = 1) {
    ENOKI_MARK_USED(size);
    if constexpr (std::is_scalar_v<Return>)
        return inner;
    else
        return Return::full_(inner, size);
}

/// Load an array from aligned memory
template <typename T> ENOKI_INLINE T load(const void *mem) {
    if constexpr (is_array_v<T>) {
        return T::load_(mem);
    } else {
        assert((uintptr_t) mem % alignof(T) == 0);
        return *static_cast<const T *>(mem);
    }
}

/// Load an array from aligned memory (masked)
template <typename T> ENOKI_INLINE T load(const void *mem, const mask_t<T> &mask) {
    if constexpr (is_array_v<T>) {
        return T::load_(mem, mask);
    } else {
        if (mask) {
            assert((uintptr_t) mem % alignof(T) == 0);
            return *static_cast<const T *>(mem);
        } else {
            return T(0);
        }
    }
}

/// Load an array from unaligned memory
template <typename T> ENOKI_INLINE T load_unaligned(const void *mem) {
    if constexpr (is_array_v<T>)
        return T::load_unaligned_(mem);
    else
        return *static_cast<const T *>(mem);
}

/// Load an array from unaligned memory (masked)
template <typename T> ENOKI_INLINE T load_unaligned(const void *mem, const mask_t<T> &mask) {
    if constexpr (is_array_v<T>)
        return T::load_unaligned_(mem, mask);
    else
        return mask ? *static_cast<const T *>(mem) : T(0);
}

/// Store an array to aligned memory
template <typename T> ENOKI_INLINE void store(void *mem, const T &value) {
    if constexpr (is_array_v<T>) {
        value.store_(mem);
    } else {
        assert((uintptr_t) mem % alignof(T) == 0);
        *static_cast<T *>(mem) = value;
    }
}

/// Store an array to aligned memory (masked)
template <typename T> ENOKI_INLINE void store(void *mem, const T &value, const mask_t<T> &mask) {
    if constexpr (is_array_v<T>) {
        value.store_(mem, mask);
    } else {
        if (mask) {
            assert((uintptr_t) mem % alignof(T) == 0);
            *static_cast<T *>(mem) = value;
        }
    }
}

/// Store an array to unaligned memory
template <typename T> ENOKI_INLINE void store_unaligned(void *mem, const T &value) {
    if constexpr (is_array_v<T>)
        value.store_unaligned_(mem);
    else
        *static_cast<T *>(mem) = value;
}

/// Store an array to unaligned memory (masked)
template <typename T> ENOKI_INLINE void store_unaligned(void *mem, const T &value, const mask_t<T> &mask) {
    if constexpr (is_array_v<T>)
        value.store_unaligned_(mem, mask);
    else if (mask)
        *static_cast<T *>(mem) = value;
}

template <typename T1, typename T2,
          enable_if_array_any_t<T1, T2> = 0> auto concat(const T1 &a1, const T2 &a2) {
    static_assert(std::is_same_v<scalar_t<T1>, scalar_t<T2>>,
                  "concat(): Scalar types must be identical");

    constexpr size_t Depth1 = array_depth_v<T1>,
                     Depth2 = array_depth_v<T2>,
                     Depth = std::max(Depth1, Depth2),
                     Size1 = array_size_v<T1>,
                     Size2 = array_size_v<T2>,
                     Size  = Size1 + Size2;

    using Value = expr_t<value_t<T1>, value_t<T2>>;
    using Result = Array<Value, Size>;
    if constexpr (Result::Size1 == Size1 && Result::Size2 == Size2 &&
                  Depth1 == 1 && Depth2 == 1) {
        return Result(a1, a2);
    } else if constexpr (Depth1 == 1 && Depth2 == 0 && T1::ActualSize == Size) {
        Result result(a1);
        #if defined(ENOKI_X86_SSE42)
            if constexpr (std::is_same_v<value_t<T1>, float>)
                result.m = _mm_insert_ps(result.m, _mm_set_ss(a2), 0b00110000);
            else
        #endif
        result.coeff(Size1) = a2;
        return result;
    } else {
        Result result;
        if constexpr (Depth1 == Depth) {
            for (size_t i = 0; i < Size1; ++i)
                result.coeff(i) = a1.derived().coeff(i);
        } else {
            result.coeff(0) = a1;
        }
        if constexpr (Depth2 == Depth) {
            for (size_t i = 0; i < Size2; ++i)
                result.coeff(i + Size1) = a2.derived().coeff(i);
        } else {
            result.coeff(Size1) = a2;
        }
        return result;
    }
}

namespace detail {
    template <typename Return, size_t Offset, typename T, size_t... Index>
    static ENOKI_INLINE Return extract(const T &a, std::index_sequence<Index...>) {
        return Return(a.coeff(Index + Offset)...);
    }
}

template <size_t Size, typename T,
          typename Return = Array<value_t<T>, Size>>
ENOKI_INLINE Return head(const T &a) {
    if constexpr (T::ActualSize == Return::ActualSize) {
        return a;
    } else if constexpr (T::Size1 == Size) {
        return low(a);
    } else {
        static_assert(Size <= array_size_v<T>, "Array size mismatch");
        return detail::extract<Return, 0>(a, std::make_index_sequence<Size>());
    }
}

template <size_t Size, typename T,
          typename Return = Array<value_t<T>, Size>>
ENOKI_INLINE Return tail(const T &a) {
    if constexpr (T::Size == Return::Size) {
        return a;
    } else if constexpr (T::Size2 == Size) {
        return high(a);
    } else {
        static_assert(Size <= array_size_v<T>, "Array size mismatch");
        return detail::extract<Return, T::Size - Size>(a, std::make_index_sequence<Size>());
    }
}

/// Masked extraction operation
template <typename Array, typename Mask>
ENOKI_INLINE auto extract(const Array &value, const Mask &mask) {
    if constexpr (is_array_v<Array>)
        return (value_t<Array>) value.extract_(mask);
    else
        return value;
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name CUDA-specific forward declarations
// -----------------------------------------------------------------------

/* Documentation in 'cuda.h' */
extern ENOKI_IMPORT void cuda_trace_printf(const char *, uint32_t, uint32_t*);
extern ENOKI_IMPORT void cuda_var_mark_dirty(uint32_t);
extern ENOKI_IMPORT void cuda_eval(bool log_assembly = false);
extern ENOKI_IMPORT void cuda_sync();
extern ENOKI_IMPORT void cuda_set_scatter_gather_operand(uint32_t index, bool gather = false);
extern ENOKI_IMPORT void cuda_set_log_level(uint32_t);
extern ENOKI_IMPORT uint32_t cuda_log_level();

/// Fancy templated 'printf', which extracts the indices of Enoki arrays
template <typename... Args> void cuda_printf(const char *fmt, const Args&... args) {
    uint32_t indices[] = { args.index()..., 0 };
    cuda_trace_printf(fmt, (uint32_t) sizeof...(Args), indices);
}

template <typename T, enable_if_t<!is_diff_array_v<T> && !is_cuda_array_v<T>> = 0>
ENOKI_INLINE void set_label(T&, const char *) { }


//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Scatter/gather/prefetch operations
// -----------------------------------------------------------------------

NAMESPACE_BEGIN(detail)

template <typename Array, bool Packed, size_t Mult1 = 0, size_t Mult2 = 1, typename Guide = Array,
          typename Func, typename Index1, typename Index2, typename Mask>
ENOKI_INLINE decltype(auto) do_recursive(const Func &func, const Index1 &offset1, const Index2 &offset2, const Mask &mask) {
    if constexpr (array_depth_v<Index1> + array_depth_v<Index2> != array_depth_v<Array>) {
        using NewIndex      = enoki::Array<scalar_t<Index1>, Guide::Size>;
        using CombinedIndex = replace_scalar_t<Index2, NewIndex>;

        constexpr size_t Size = (Packed || (array_depth_v<Index1> + array_depth_v<Index2> + 1 != array_depth_v<Array>)) ?
            Guide::Size : enoki::Array<scalar_t<Guide>, Guide::Size>::ActualSize;  /* Deal with n=3 special case */

        CombinedIndex combined_offset =
            CombinedIndex(offset2 * scalar_t<Index1>(Size)) +
            full<Index2>(arange<NewIndex>());

        return do_recursive<Array, Packed, Mult1, Mult2 * Size, value_t<Guide>>(
            func, offset1, combined_offset, mask);
    } else {
        using CombinedIndex = replace_scalar_t<Index2, Index1>;

        CombinedIndex combined_offset =
            CombinedIndex(offset2) +
            enoki::full<Index2>(offset1) * scalar_t<Index1>(Mult1 == 0 ? Mult2 : Mult1);

        return func(combined_offset, full<Index2>(mask));
    }
}

template <typename T> constexpr size_t fix_stride(size_t Stride) {
    if (has_avx2) {
       if (Stride % 8 == 0)      return 8;
       else if (Stride % 4 == 0) return 4;
       else                      return 1;
    }
    return Stride;
}

NAMESPACE_END(detail)

/// Masked prefetch operation
template <typename Array, bool Write = false, size_t Level = 2, size_t Stride_ = 0, bool Packed = true,
          typename Index, typename Mask = mask_t<replace_scalar_t<Index, scalar_t<Array>>>>
ENOKI_INLINE void prefetch(const void *mem, const Index &index, const identity_t<Mask> &mask = true) {
    static_assert(is_std_int_v<scalar_t<Index>>, "prefetch(): expected a signed 32/64-bit integer as 'index' argument!");
    constexpr size_t ScalarSize = sizeof(scalar_t<Array>);

    if constexpr (!is_array_v<Array> && !is_array_v<Index>) {
        /* Scalar case */
        #if defined(ENOKI_X86_SSE42)
            if (mask) {
                constexpr size_t Stride = (Stride_ != 0) ? Stride_ : ScalarSize;
                const uint8_t *ptr = (const uint8_t *) mem + index * Index(Stride);
                constexpr auto Hint = Level == 1 ? _MM_HINT_T0 : _MM_HINT_T1;
                _mm_prefetch((char *) ptr, Hint);
            }
        #else
            (void) mem; (void) index; (void) mask;
        #endif
    } else if constexpr (std::is_same_v<array_shape_t<Array>, array_shape_t<Index>>) {
        /* Forward to the array-specific implementation */
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ : ScalarSize,
                         Stride2 = detail::fix_stride<Array>(Stride);
        Index index2 = Stride != Stride2 ? index * scalar_t<Index>(Stride / Stride2) : index;
        Array::template prefetch_<Write, Level, Stride2>(mem, index2, mask);
    } else if constexpr (array_depth_v<Array> > array_depth_v<Index>) {
        /* Dimension mismatch, reduce to a sequence of gather operations */
        static_assert((Stride_ / ScalarSize) * ScalarSize == Stride_,
                      "Stride must be divisible by sizeof(Scalar)");
        return detail::do_recursive<Array, Packed, Stride_ / ScalarSize>(
            [mem](const auto &index2, const auto &mask2) ENOKI_INLINE_LAMBDA {
                constexpr size_t ScalarSize2 = sizeof(scalar_t<Array>); // needed for MSVC
                prefetch<Array, Write, Level, ScalarSize2>(mem, index2, mask2);
            },
            index, scalar_t<Index>(0), mask);
    } else {
        static_assert(detail::false_v<Array>, "prefetch(): don't know what to do with the input arguments!");
    }
}

/// Masked gather operation
template <typename Array, size_t Stride_ = 0, bool Packed = true, bool Masked = true,
          typename Index, typename Mask = mask_t<replace_scalar_t<Index, scalar_t<Array>>>>
ENOKI_INLINE Array gather(const void *mem, const Index &index, const identity_t<Mask> &mask) {
    static_assert(is_std_int_v<scalar_t<Index>>, "gather(): expected a signed 32/64-bit integer as 'index' argument!");
    constexpr size_t ScalarSize = sizeof(scalar_t<Array>);

    if constexpr (!is_array_v<Array> && !is_array_v<Index>) {
        /* Scalar case */
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ : ScalarSize;
        const Array *ptr = (const Array *) ((const uint8_t *) mem + index * Index(Stride));
        return mask ? *ptr : Array(0);
    } else if constexpr (std::is_same_v<array_shape_t<Array>, array_shape_t<Index>>) {
        /* Forward to the array-specific implementation */
        constexpr size_t Stride  = (Stride_ != 0) ? Stride_ : ScalarSize,
                         Stride2 = detail::fix_stride<Array>(Stride);
        Index index2 = Stride != Stride2 ? index * scalar_t<Index>(Stride / Stride2) : index;
        return Array::template gather_<Stride2>(mem, index2, mask);
    } else if constexpr (array_depth_v<Array> == 1 && array_depth_v<Index> == 0) {
        /* Turn into a load */
        ENOKI_MARK_USED(mask);
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ :
            (Packed ? (sizeof(value_t<Array>) * array_size_v<Array>) : (sizeof(Array)));
        if constexpr (Masked)
            return load_unaligned<Array>((uint8_t *) mem + Stride * (size_t) index, mask);
        else
            return load_unaligned<Array>((uint8_t *) mem + Stride * (size_t) index);
    } else if constexpr (array_depth_v<Array> > array_depth_v<Index>) {
        /* Dimension mismatch, reduce to a sequence of gather operations */
        static_assert((Stride_ / ScalarSize) * ScalarSize == Stride_,
                      "Stride must be divisible by sizeof(Scalar)");
        return detail::do_recursive<Array, Packed, Stride_ / ScalarSize>(
            [mem](const auto &index2, const auto &mask2) ENOKI_INLINE_LAMBDA {
                constexpr size_t ScalarSize2 = sizeof(scalar_t<Array>); // needed for MSVC
                return gather<Array, ScalarSize2>(mem, index2, mask2);
            },
            index, scalar_t<Index>(0), mask);
    } else {
        static_assert(detail::false_v<Array>, "gather(): don't know what to do with the input arguments!");
    }
}

/// Masked scatter operation
template <size_t Stride_ = 0, bool Packed = true, bool Masked = true, typename Array, typename Index,
          typename Mask = mask_t<replace_scalar_t<Index, scalar_t<Array>>>>
ENOKI_INLINE void scatter(void *mem, const Array &value, const Index &index, const identity_t<Mask> &mask) {
    static_assert(is_std_int_v<scalar_t<Index>>, "scatter(): expected a signed 32/64-bit integer as 'index' argument!");
    constexpr size_t ScalarSize = sizeof(scalar_t<Array>);

    if constexpr (!is_array_v<Array> && !is_array_v<Index>) {
        /* Scalar case */
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ : ScalarSize;
        Array *ptr = (Array *) ((uint8_t *) mem + index * Index(Stride));
        if (mask)
            *ptr = value;
    } else if constexpr (std::is_same_v<array_shape_t<Array>, array_shape_t<Index>>) {
        /* Forward to the array-specific implementation */
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ : ScalarSize,
                         Stride2 = detail::fix_stride<Array>(Stride);
        Index index2 = Stride != Stride2 ? index * scalar_t<Index>(Stride / Stride2) : index;
        value.template scatter_<Stride2>(mem, index2, mask);
    } else if constexpr (array_depth_v<Array> == 1 && array_depth_v<Index> == 0) {
        /* Turn into a store */
        ENOKI_MARK_USED(mask);
        constexpr size_t Stride = (Stride_ != 0) ? Stride_ :
            (Packed ? (sizeof(value_t<Array>) * array_size_v<Array>) : (sizeof(Array)));
        if constexpr (Masked)
            return store_unaligned((uint8_t *) mem + Stride * (size_t) index, value, mask);
        else
            return store_unaligned((uint8_t *) mem + Stride * (size_t) index, value);
    } else if constexpr (array_depth_v<Array> > array_depth_v<Index>) {
        /* Dimension mismatch, reduce to a sequence of gather operations */
        static_assert((Stride_ / ScalarSize) * ScalarSize == Stride_,
                      "Stride must be divisible by sizeof(Scalar)");
        detail::do_recursive<Array, Packed, Stride_ / ScalarSize>(
            [mem, &value](const auto &index2, const auto &mask2) ENOKI_INLINE_LAMBDA {
                constexpr size_t ScalarSize2 = sizeof(scalar_t<Array>); // needed for MSVC
                scatter<ScalarSize2, Masked>(mem, value, index2, mask2);
            },
            index, scalar_t<Index>(0), mask);
    } else {
        static_assert(detail::false_v<Array>, "scatter(): don't know what to do with the input arguments!");
    }
}

template <typename Array, size_t Stride = 0, bool Packed = true, typename Index>
ENOKI_INLINE Array gather(const void *mem, const Index &index) {
    return gather<Array, Stride, Packed, false>(mem, index, true);
}

template <size_t Stride_ = 0, bool Packed = true, typename Array, typename Index>
ENOKI_INLINE void scatter(void *mem, const Array &value, const Index &index) {
    scatter<Stride_, Packed, false>(mem, value, index, true);
}

#if defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wunused-value"
#endif

/// Conflict-free modification operation
template <typename Arg, size_t Stride = sizeof(scalar_t<Arg>),
          typename Func, typename Index, typename... Args>
void transform(void *mem, const Index &index, Func &&func, Args&&... args) {
    static_assert(is_std_int_v<scalar_t<Index>>,
                  "transform(): index argument must be a 32/64-bit integer array!");
    if constexpr (is_array_v<Arg>) {
        using Int = int_array_t<Arg>;
        if constexpr ((false, ..., is_mask_v<Args>))
            Arg::template transform_<Stride>(mem, (const Int &) index, (..., args), func, args...);
        else
            Arg::template transform_<Stride>(mem, (const Int &) index, mask_t<Arg>(true),
                                             func, args..., mask_t<Arg>(true));
    } else {
        Arg& ref = *(Arg *) ((uint8_t *) mem + index * Index(Stride));
        if constexpr ((false, ..., is_mask_v<Args>)) {
            if ((..., args))
                func(ref, args...);
        } else {
            func(ref, args..., true);
        }
    }
}

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

/// Conflict-free scatter-add update
template <size_t Stride_ = 0, typename Arg, typename Index>
ENOKI_INLINE void scatter_add(void *mem, const Arg &value, const Index &index, mask_t<Arg> mask = true) {
    static_assert(is_std_int_v<scalar_t<Index>>,
                  "scatter_add(): index argument must be a 32/64-bit integer array!");
    constexpr size_t Stride = Stride_ == 0 ? sizeof(scalar_t<Arg>) : Stride_;

    if constexpr (is_array_v<Arg>) {
        value.template scatter_add_<Stride>(mem, index, mask);
    } else {
        Arg& ref = *(Arg *) ((uint8_t *) mem + index * Index(Stride));
        if (mask)
            ref += value;
    }
}

/// Prefetch operations with an array source
template <typename Array, bool Write = false, size_t Level = 2, size_t Stride = 0,
          bool Packed = true, typename Source, typename... Args,
          enable_if_t<array_depth_v<Source> == 1> = 0>
ENOKI_INLINE void prefetch(const Source &source, const Args &... args) {
    prefetch<Array, Write, Level, Stride, Packed>(source.data(), args...);
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Nested horizontal reduction operators
// -----------------------------------------------------------------------

template <typename T> auto hsum_nested(const T &a) {
    if constexpr (array_depth_v<T> == 1)
        return hsum(a);
    else if constexpr (is_array_v<T>)
        return hsum_nested(hsum(a));
    else
        return a;
}

template <typename T> auto hprod_nested(const T &a) {
    if constexpr (array_depth_v<T> == 1)
        return hprod(a);
    else if constexpr (is_array_v<T>)
        return hprod_nested(hprod(a));
    else
        return a;
}

template <typename T> auto hmin_nested(const T &a) {
    if constexpr (array_depth_v<T> == 1)
        return hmin(a);
    else if constexpr (is_array_v<T>)
        return hmin_nested(hmin(a));
    else
        return a;
}

template <typename T> auto hmax_nested(const T &a) {
    if constexpr (array_depth_v<T> == 1)
        return hmax(a);
    else if constexpr (is_array_v<T>)
        return hmax_nested(hmax(a));
    else
        return a;
}

template <typename T> auto hmean_nested(const T &a) {
    if constexpr (array_depth_v<T> == 1)
        return hmean(a);
    else if constexpr (is_array_v<T>)
        return hmean_nested(hmean(a));
    else
        return a;
}

template <typename T> auto count_nested(const T &a) {
    if constexpr (is_array_v<T>)
        return hsum_nested(count(a));
    else
        return count(a);
}

template <typename T> auto any_nested(const T &a) {
    if constexpr (is_array_v<T>)
        return any_nested(any(a));
    else
        return any(a);
}

template <typename T> auto all_nested(const T &a) {
    if constexpr (is_array_v<T>)
        return all_nested(all(a));
    else
        return all(a);
}

template <typename T> auto none_nested(const T &a) {
    return !any_nested(a);
}

/// Convert an array with 1 entry into a scalar or throw an error
template <typename T> scalar_t<T> scalar_cast(const T &v) {
    static_assert(array_depth_v<T> <= 1);
    if constexpr (is_array_v<T>) {
        if (v.size() != 1)
            throw std::runtime_error("scalar_cast(): array should be of size 1!");
        return v.coeff(0);
    } else {
        return v;
    }
}

template <typename T1, typename T2>
bool allclose(const T1 &a, const T2 &b, float rtol = 1e-5f, float atol = 1e-8f,
              bool equal_nan = false) {
    auto cond = abs(a - b) <= abs(b) * rtol + atol;

    if constexpr (std::is_floating_point_v<scalar_t<T1>> &&
                  std::is_floating_point_v<scalar_t<T2>>) {
        if (equal_nan)
            cond |= isnan(a) & isnan(b);
    }

    return all_nested(cond);
}

//! @}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
//! @{ \name Reduction operators that return a default argument when
//           invoked using CUDA arrays
// -----------------------------------------------------------------------

template <bool Default, typename T> auto any_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return any(value);
    }
}

template <bool Default, typename T> auto any_nested_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return any_nested(value);
    }
}

template <bool Default, typename T> auto none_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return none(value);
    }
}

template <bool Default, typename T> auto none_nested_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return none_nested(value);
    }
}

template <bool Default, typename T> auto all_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return all(value);
    }
}

template <bool Default, typename T> auto all_nested_or(const T &value) {
    if constexpr (is_cuda_array_v<T>) {
        ENOKI_MARK_USED(value);
        return Default;
    } else {
        return all_nested(value);
    }
}

//! @}
// -----------------------------------------------------------------------

#undef ENOKI_ROUTE_UNARY
#undef ENOKI_ROUTE_UNARY_IMM
#undef ENOKI_ROUTE_UNARY_SCALAR
#undef ENOKI_ROUTE_UNARY_SCALAR_IMM
#undef ENOKI_ROUTE_BINARY
#undef ENOKI_ROUTE_BINARY_BITOP
#undef ENOKI_ROUTE_BINARY_COND
#undef ENOKI_ROUTE_BINARY_SHIFT
#undef ENOKI_ROUTE_BINARY_SCALAR
#undef ENOKI_ROUTE_TERNARY
#undef ENOKI_ROUTE_COMPOUND_OPERATOR

NAMESPACE_END(enoki)
