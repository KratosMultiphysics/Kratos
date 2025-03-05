/*
    enoki/array_idiv.h -- fast precomputed integer division by constants based
    on libdivide (https://github.com/ridiculousfish/libdivide)

    Copyright (C) 2010 ridiculous_fish

    This software is provided 'as-is', without any express or implied
    warranty.  In no event will the authors be held liable for any damages
    arising from the use of this software.

    Permission is granted to anyone to use this software for any purpose,
    including commercial applications, and to alter it and redistribute it
    freely, subject to the following restrictions:

    1. The origin of this software must not be misrepresented; you must not
       claim that you wrote the original software. If you use this software
       in a product, an acknowledgment in the product documentation would be
       appreciated but is not required.
    2. Altered source versions must be plainly marked as such, and must not be
       misrepresented as being the original software.
    3. This notice may not be removed or altered from any source distribution.

    libdivide@ridiculousfish.com

*/

#pragma once

#include <enoki/array_generic.h>

NAMESPACE_BEGIN(enoki)

NAMESPACE_BEGIN(detail)

// -----------------------------------------------------------------------
//! @{ \name Precomputation for division by integer constants
// -----------------------------------------------------------------------

template <bool UseIntrinsic = false>
std::pair<uint32_t, uint32_t> div_wide(uint32_t u1, uint32_t u0, uint32_t v) {
#if defined(__GNUC__) && (defined(ENOKI_X86_32) || defined(ENOKI_X86_64))
    if constexpr (UseIntrinsic) {
        uint32_t res, rem;
        __asm__("divl %[v]"
                : "=a"(res), "=d"(rem)
                : [v] "r"(v), "a"(u0), "d"(u1));
        return { res, rem };
    }
#endif

    uint64_t u = (((uint64_t) u1) << 32) | u0;

    return { (uint32_t) (u / v),
             (uint32_t) (u % v) };
}

template <bool UseIntrinsic = false>
std::pair<uint64_t, uint64_t> div_wide(uint64_t u1, uint64_t u0, uint64_t d) {
#if defined(__GNUC__) && defined(ENOKI_X86_64)
    if constexpr (UseIntrinsic) {
        uint64_t res, rem;
        __asm__("divq %[v]"
                : "=a"(res), "=d"(rem)
                : [v]"r"(d), "a"(u0), "d"(u1));
        return { res, rem };
    }
#endif

#if defined(__SIZEOF_INT128__)
    __uint128_t n = (((__uint128_t) u1) << 64) | u0;
    return {
        (uint64_t) (n / d),
        (uint64_t) (n % d)
    };
#else
    // Code taken from Hacker's Delight:
    // http://www.hackersdelight.org/HDcode/divlu.c.
    // License permits inclusion here per:
    // http://www.hackersdelight.org/permissions.htm

    const uint64_t b = (1ULL << 32); // Number base (16 bits).
    uint64_t un1, un0,  // Norm. dividend LSD's.
    vn1, vn0,           // Norm. divisor digits.
    q1, q0,             // Quotient digits.
    un64, un21, un10,   // Dividend digit pairs.
    rhat;               // A remainder.
    int s;              // Shift amount for norm.

    if (u1 >= d) // overflow
        return { (uint64_t) -1, (uint64_t) -1 };

    // count leading zeros
    s = (int) (63 - log2i(d)); // 0 <= s <= 63.
    if (s > 0) {
        d = d << s;         // Normalize divisor.
        un64 = (u1 << s) | ((u0 >> (64 - s)) & uint64_t(-s >> 31));
        un10 = u0 << s;     // Shift dividend left.
    } else {
        // Avoid undefined behavior.
        un64 = u1 | u0;
        un10 = u0;
    }

    vn1 = d >> 32;            // Break divisor up into
    vn0 = d & 0xFFFFFFFF;     // two 32-bit digits.

    un1 = un10 >> 32;         // Break right half of
    un0 = un10 & 0xFFFFFFFF;  // dividend into two digits.

    q1 = un64/vn1;            // Compute the first
    rhat = un64 - q1*vn1;     // quotient digit, q1.

again1:
    if (q1 >= b || q1*vn0 > b*rhat + un1) {
        q1 = q1 - 1;
        rhat = rhat + vn1;
        if (rhat < b)
            goto again1;
    }

    un21 = un64*b + un1 - q1*d;  // Multiply and subtract.

    q0 = un21/vn1;            // Compute the second
    rhat = un21 - q0*vn1;     // quotient digit, q0.

again2:
    if (q0 >= b || q0 * vn0 > b * rhat + un0) {
        q0 = q0 - 1;
        rhat = rhat + vn1;
        if (rhat < b)
            goto again2;
    }

    return {
        q1*b + q0,
        (un21*b + un0 - q0*d) >> s
    };
#endif
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(detail)

#if defined(_MSC_VER)
#  pragma pack(push)
#  pragma pack(1)
#endif

template <typename T, bool UseIntrinsic>
struct divisor<T, UseIntrinsic, enable_if_t<std::is_unsigned_v<T>>> {
    T multiplier;
    uint8_t shift;

    divisor() = default;

    ENOKI_INLINE divisor(T d) {
        /* Division by +/-1 is not supported by the
           precomputation-based approach */
        assert(d != 1);
        shift = (uint8_t) log2i(d);

        if ((d & (d - 1)) == 0) {
            /* Power of two */
            multiplier = 0;
            shift--;
        } else {
            /* General case */
            auto [m, rem] =
                detail::div_wide<UseIntrinsic>(T(1) << shift, T(0), d);
            multiplier = m * 2 + 1;
            assert(rem > 0 && rem < d);

            T rem2 = rem * 2;
            if (rem2 >= d || rem2 < rem)
                multiplier += 1;
        }
    }

    template <typename T2>
    ENOKI_INLINE auto operator()(const T2 &value) const {
        using Expr = decltype(value + value);
        auto q = mulhi(Expr(multiplier), value);
        auto t = sr<1>(value - q) + q;
        return t >> shift;
    }
} ENOKI_PACK;

template <typename T, bool UseIntrinsic>
struct divisor<T, UseIntrinsic, enable_if_t<std::is_signed_v<T>>> {
    using U = std::make_unsigned_t<T>;

    T multiplier;
    uint8_t shift;

    divisor() = default;

    ENOKI_INLINE divisor(T d) {
        /* Division by +/-1 is not supported by the
           precomputation-based approach */
        assert(d != 1 && d != -1);

        U ad = d < 0 ? (U) -d : (U) d;
        shift = (uint8_t) log2i(ad);

        if ((ad & (ad - 1)) == 0) {
            /* Power of two */
            multiplier = 0;
        } else {
            /* General case */
            auto [m, rem] =
                detail::div_wide<UseIntrinsic>(U(1) << (shift - 1), U(0), ad);
            multiplier = T(m * 2 + 1);

            U rem2 = rem * 2;
            if (rem2 >= ad || rem2 < rem)
                multiplier += 1;
        }
        if (d < 0)
            shift |= 0x80;
    }

    template <typename T2>
    ENOKI_INLINE auto operator()(const T2 &value) const {
        using Expr = decltype(value + value);
        uint8_t shift_ = shift & 0x3f;
        Expr sign(int8_t(shift) >> 7);

        auto q = mulhi(Expr(multiplier), value) + value;
        auto q_sign = sr<sizeof(T) * 8 - 1>(q);
        q += q_sign & ((T(1) << shift_) - (multiplier == 0 ? 1 : 0));

        return ((q >> shift_) ^ sign) - sign;
    }
} ENOKI_PACK;

/// Stores *both* the original divisor + magic number
template <typename T> struct divisor_ext : divisor<T> {
    T value;
    ENOKI_INLINE divisor_ext(T value) : divisor<T>(value), value(value) { }
} ENOKI_PACK;

#if defined(_MSC_VER)
#  pragma pack(pop)
#endif

template <typename T, enable_if_t<std::is_integral_v<scalar_t<T>>> = 0>
ENOKI_INLINE auto operator/(const T &a, const divisor<scalar_t<T>> &div) {
    return div(a);
}

template <typename T, enable_if_t<std::is_integral_v<scalar_t<T>>> = 0>
ENOKI_INLINE auto operator%(const T &a, const divisor_ext<scalar_t<T>> &div) {
    return a - div(a) * div.value;
}

// -----------------------------------------------------------------------
//! @{ \name Arithmetic operations for pointer arrays
// -----------------------------------------------------------------------

template <typename T1, typename T2,
          typename S1 = scalar_t<T1>, typename S2 = scalar_t<T2>,
          enable_if_t<std::is_pointer_v<S1> || std::is_pointer_v<S2>> = 0,
          enable_if_array_any_t<T1, T2> = 0>
ENOKI_INLINE auto operator-(const T1 &a1_, const T2 &a2_) {
    using Int = std::conditional_t<sizeof(void *) == 8, int64_t, int32_t>;
    using T1i = replace_scalar_t<T1, Int, false>;
    using T2i = replace_scalar_t<T2, Int, false>;
    using Ti  = expr_t<T1i, T2i>;
    using T   = expr_t<T1, T2>;

    constexpr Int InstanceSize    = sizeof(std::remove_pointer_t<scalar_t<T1>>),
                  LogInstanceSize = detail::clog2i(InstanceSize);

    constexpr bool PointerDiff = std::is_pointer_v<S1> &&
                                 std::is_pointer_v<S2>;

    using Ret = std::conditional_t<PointerDiff, Ti, T>;
    Ti a1 = Ti((T1i) a1_),
       a2 = Ti((T2i) a2_);

    if constexpr (InstanceSize == 1) {
        return Ret(a1.sub_(a2));
    } else if constexpr ((1 << LogInstanceSize) == InstanceSize) {
        if constexpr (PointerDiff)
            return Ret(a1.sub_(a2).template sr_<LogInstanceSize>());
        else
            return Ret(a1.sub_(a2.template sl_<LogInstanceSize>()));
    } else {
        if constexpr (PointerDiff)
            return Ret(a1.sub_(a2) / InstanceSize);
        else
            return Ret(a1.sub_(a2 * InstanceSize));
    }
}


template <typename T1, typename T2,
          typename S1 = scalar_t<T1>, typename S2 = scalar_t<T2>,
          enable_if_t<std::is_pointer_v<S1> && !std::is_pointer_v<S2>> = 0,
          enable_if_array_any_t<T1, T2> = 0>
ENOKI_INLINE auto operator+(const T1 &a1_, const T2 &a2_) {
    using Int = std::conditional_t<sizeof(void *) == 8, int64_t, int32_t>;
    using T1i = replace_scalar_t<T1, Int, false>;
    using T2i = replace_scalar_t<T2, Int, false>;
    using Ti  = expr_t<T1i, T2i>;
    using Ret = expr_t<T1, T2>;

    constexpr Int InstanceSize    = sizeof(std::remove_pointer_t<scalar_t<T1>>),
                  LogInstanceSize = detail::clog2i(InstanceSize);

    Ti a1 = Ti((T1i) a1_),
       a2 = Ti((T2i) a2_);

    if constexpr (InstanceSize == 1)
        return Ret(a1.add_(a2));
    if constexpr ((1 << LogInstanceSize) == InstanceSize)
        return Ret(a1.add_(a2.template sl_<LogInstanceSize>()));
    else
        return Ret(a1.add_(a2 * InstanceSize));
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
