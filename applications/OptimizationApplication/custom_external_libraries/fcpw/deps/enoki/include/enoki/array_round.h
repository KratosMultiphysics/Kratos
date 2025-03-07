/*
    enoki/array_round.h -- Fallback for nonstandard rounding modes

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using ENOKI instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_generic.h>

NAMESPACE_BEGIN(enoki)

#if defined(ENOKI_X86_64) || defined(ENOKI_X86_32)
/// RAII wrapper that saves and restores the FP Control/Status Register
template <RoundingMode Mode> struct set_rounding_mode {
    set_rounding_mode() : value(_mm_getcsr()) {
        unsigned int csr = value & ~(unsigned int) _MM_ROUND_MASK;
        switch (Mode) {
            case RoundingMode::Nearest: csr |= _MM_ROUND_NEAREST; break;
            case RoundingMode::Down: csr |= _MM_ROUND_DOWN; break;
            case RoundingMode::Up: csr |= _MM_ROUND_UP; break;
            case RoundingMode::Zero: csr |= _MM_ROUND_TOWARD_ZERO; break;
        }
        _mm_setcsr(csr);
    }

    ~set_rounding_mode() {
        _mm_setcsr(value);
    }

    unsigned int value;
};
#else
template <RoundingMode Mode> struct set_rounding_mode {
    // Don't know how to change rounding modes on this platform :(
};
#endif

template <typename Value_, size_t Size_, bool Approx_, RoundingMode Mode_, bool IsMask_, typename Derived_>
struct StaticArrayImpl<Value_, Size_, Approx_, Mode_, IsMask_, Derived_,
                       enable_if_t<detail::array_config<Value_, Size_, Mode_>::use_rounding_fallback_impl>>
    : StaticArrayImpl<Value_, Size_, Approx_, RoundingMode::Default, IsMask_, Derived_> {

    using Base = StaticArrayImpl<Value_, Size_, Approx_, RoundingMode::Default, IsMask_, Derived_>;
    using Derived = Derived_;

    using Base::derived;

    /// Rounding mode of arithmetic operations
    static constexpr RoundingMode Mode = Mode_;

    template <typename Arg, enable_if_t<std::is_same_v<value_t<Arg>, Value_>> = 0>
    ENOKI_INLINE StaticArrayImpl(Arg&& arg) : Base(std::forward<Arg>(arg)) { }

    template <typename... Args>
    ENOKI_INLINE StaticArrayImpl(Args&&... args) : Base(std::forward<Args>(args)...) { }

    template <typename Arg, enable_if_t<!std::is_same_v<value_t<Arg>, Value_>> = 0>
    ENOKI_NOINLINE StaticArrayImpl(Arg&& arg) {
        set_rounding_mode<Mode_> mode; (void) mode;
        using Base2 = std::conditional_t<IsMask_,
            Array<Value_, Size_, Approx_, RoundingMode::Default>,
            Packet<Value_, Size_, Approx_, RoundingMode::Default>>;
        Base::operator=(Base2(std::forward<Arg>(arg)));
    }

    template <typename Arg, enable_if_t<std::is_same_v<value_t<Arg>, Value_>> = 0>
    ENOKI_NOINLINE Derived& operator=(Arg&& arg) {
        Base::operator=(std::forward<Arg>(arg));
        return derived();
    }

    template <typename Arg, enable_if_t<!std::is_same_v<value_t<Arg>, Value_>> = 0>
    ENOKI_NOINLINE Derived& operator=(Arg&& arg) {
        set_rounding_mode<Mode_> mode; (void) mode;
        using Base2 = std::conditional_t<IsMask_,
            Array<Value_, Size_, Approx_, RoundingMode::Default>,
            Packet<Value_, Size_, Approx_, RoundingMode::Default>>;
        Base::operator=(Base2(std::forward<Arg>(arg)));
        return derived();
    }

    ENOKI_NOINLINE Derived add_(const Derived &a) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::add_(a);
    }

    ENOKI_NOINLINE Derived sub_(const Derived &a) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::sub_(a);
    }

    ENOKI_NOINLINE Derived mul_(const Derived &a) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::mul_(a);
    }

    ENOKI_NOINLINE Derived div_(const Derived &a) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::div_(a);
    }

    ENOKI_NOINLINE Derived sqrt_() const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::sqrt_();
    }

    ENOKI_NOINLINE Derived fmadd_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fmadd_(b, c);
    }

    ENOKI_NOINLINE Derived fmsub_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fmsub_(b, c);
    }

    ENOKI_NOINLINE Derived fnmadd_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fnmadd_(b, c);
    }

    ENOKI_NOINLINE Derived fnmsub_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fnmsub_(b, c);
    }

    ENOKI_NOINLINE Derived fmsubadd_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fmsubadd_(b, c);
    }

    ENOKI_NOINLINE Derived fmaddsub_(const Derived &b, const Derived &c) const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::fmaddsub_(b, c);
    }

    ENOKI_NOINLINE Value_ hsum() const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::hsum_();
    }

    ENOKI_NOINLINE Value_ hprod() const {
        set_rounding_mode<Mode_> mode; (void) mode;
        return Base::hprod_();
    }
};

NAMESPACE_END(enoki)
