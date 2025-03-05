/*
    enoki/half.h -- minimal half precision number type

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array_traits.h>

NAMESPACE_BEGIN(enoki)
struct half;
NAMESPACE_END(enoki)

NAMESPACE_BEGIN(std)
template<> struct is_floating_point<enoki::half> : true_type { };
template<> struct is_arithmetic<enoki::half> : true_type { };
template<> struct is_signed<enoki::half> : true_type { };
NAMESPACE_END(std)

NAMESPACE_BEGIN(enoki)
struct half {
    uint16_t value;

    half()
    #if !defined(NDEBUG)
        : value(0x7FFF) /* Initialize with NaN */
    #endif
    { }

    #define ENOKI_IF_SCALAR template <typename Value, enable_if_t<std::is_arithmetic_v<Value>> = 0>

    ENOKI_IF_SCALAR half(Value val) : value(float32_to_float16(float(val))) { }

    half operator+(half h) const { return half(float(*this) + float(h)); }
    half operator-(half h) const { return half(float(*this) - float(h)); }
    half operator*(half h) const { return half(float(*this) * float(h)); }
    half operator/(half h) const { return half(float(*this) / float(h)); }

    half operator-() const { return half(-float(*this)); }

    ENOKI_IF_SCALAR friend half operator+(Value val, half h) { return half(val) + h; }
    ENOKI_IF_SCALAR friend half operator-(Value val, half h) { return half(val) - h; }
    ENOKI_IF_SCALAR friend half operator*(Value val, half h) { return half(val) * h; }
    ENOKI_IF_SCALAR friend half operator/(Value val, half h) { return half(val) / h; }

    half& operator+=(half h) { return operator=(*this + h); }
    half& operator-=(half h) { return operator=(*this - h); }
    half& operator*=(half h) { return operator=(*this * h); }
    half& operator/=(half h) { return operator=(*this / h); }

    bool operator==(half h) const { return float(*this) == float(h); }
    bool operator!=(half h) const { return float(*this) != float(h); }
    bool operator<(half h) const  { return float(*this) < float(h); }
    bool operator>(half h) const  { return float(*this) > float(h); }
    bool operator<=(half h) const { return float(*this) <= float(h); }
    bool operator>=(half h) const { return float(*this) >= float(h); }

    ENOKI_IF_SCALAR operator Value() const { return Value(float16_to_float32(value)); }

    static half from_binary(uint16_t value) { half h; h.value = value; return h; }

    friend std::ostream &operator<<(std::ostream &os, const half &h) {
        os << float(h);
        return os;
    }

    #undef ENOKI_IF_SCALAR
private:
    /*
       Value float32<->float16 conversion code by Paul A. Tessier (@Phernost)
       Used with permission by the author, who released this code into the public domain
     */
    union Bits {
        float f;
        int32_t si;
        uint32_t ui;
    };

    static constexpr int const shift = 13;
    static constexpr int const shiftSign = 16;

    static constexpr int32_t const infN = 0x7F800000;  // flt32 infinity
    static constexpr int32_t const maxN = 0x477FE000;  // max flt16 normal as a flt32
    static constexpr int32_t const minN = 0x38800000;  // min flt16 normal as a flt32
    static constexpr int32_t const signN = (int32_t) 0x80000000; // flt32 sign bit

    static constexpr int32_t const infC = infN >> shift;
    static constexpr int32_t const nanN = (infC + 1) << shift; // minimum flt16 nan as a flt32
    static constexpr int32_t const maxC = maxN >> shift;
    static constexpr int32_t const minC = minN >> shift;
    static constexpr int32_t const signC = signN >> shiftSign; // flt16 sign bit

    static constexpr int32_t const mulN = 0x52000000; // (1 << 23) / minN
    static constexpr int32_t const mulC = 0x33800000; // minN / (1 << (23 - shift))

    static constexpr int32_t const subC = 0x003FF; // max flt32 subnormal down shifted
    static constexpr int32_t const norC = 0x00400; // min flt32 normal down shifted

    static constexpr int32_t const maxD = infC - maxC - 1;
    static constexpr int32_t const minD = minC - subC - 1;

public:
    static uint16_t float32_to_float16(float value) {
        #if defined(ENOKI_X86_F16C)
            return (uint16_t) _mm_cvtsi128_si32(
                _mm_cvtps_ph(_mm_set_ss(value), _MM_FROUND_CUR_DIRECTION));
        #elif defined(ENOKI_ARM_NEON)
            return memcpy_cast<uint16_t>((__fp16) value);
        #else
            Bits v, s;
            v.f = value;
            uint32_t sign = (uint32_t) (v.si & signN);
            v.si ^= sign;
            sign >>= shiftSign; // logical shift
            s.si = mulN;
            s.si = (int32_t) (s.f * v.f); // correct subnormals
            v.si ^= (s.si ^ v.si) & -(minN > v.si);
            v.si ^= (infN ^ v.si) & -((infN > v.si) & (v.si > maxN));
            v.si ^= (nanN ^ v.si) & -((nanN > v.si) & (v.si > infN));
            v.ui >>= shift; // logical shift
            v.si ^= ((v.si - maxD) ^ v.si) & -(v.si > maxC);
            v.si ^= ((v.si - minD) ^ v.si) & -(v.si > subC);
            return (uint16_t) (v.ui | sign);
        #endif
    }

    static float float16_to_float32(uint16_t value) {
        #if defined(ENOKI_X86_F16C)
            return _mm_cvtss_f32(_mm_cvtph_ps(_mm_cvtsi32_si128((int32_t) value)));
        #elif defined(ENOKI_ARM_NEON)
            return (float) memcpy_cast<__fp16>(value);
        #else
            Bits v;
            v.ui = value;
            int32_t sign = v.si & signC;
            v.si ^= sign;
            sign <<= shiftSign;
            v.si ^= ((v.si + minD) ^ v.si) & -(v.si > subC);
            v.si ^= ((v.si + maxD) ^ v.si) & -(v.si > maxC);
            Bits s;
            s.si = mulC;
            s.f *= float(v.si);
            int32_t mask = -(norC > v.si);
            v.si <<= shift;
            v.si ^= (s.si ^ v.si) & mask;
            v.si |= sign;
            return v.f;
        #endif
    }
};

NAMESPACE_END(enoki)

NAMESPACE_BEGIN(std)

template<> struct numeric_limits<enoki::half> {
    static constexpr bool is_signed = true;
    static constexpr bool is_exact = false;
    static constexpr bool is_modulo = false;
    static constexpr bool is_iec559 = true;
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr int digits = 11;
    static constexpr int digits10 = 3;
    static constexpr int max_digits10 = 5;
    static constexpr int radix = 2;
    static constexpr int min_exponent = -13;
    static constexpr int min_exponent10 = -4;
    static constexpr int max_exponent = 16;
    static constexpr int max_exponent10 = 4;
    static constexpr float_denorm_style has_denorm = denorm_present;
    static constexpr float_round_style round_style = round_indeterminate;
    static enoki::half min() noexcept { return enoki::half::from_binary(0x0400); }
    static enoki::half lowest() noexcept { return enoki::half::from_binary(0xFBFF); }
    static enoki::half max() noexcept { return enoki::half::from_binary(0x7BFF); }
    static enoki::half epsilon() noexcept { return enoki::half::from_binary(0x1400); }
    static enoki::half round_error() noexcept { return enoki::half::from_binary(0x3C00); }
    static enoki::half infinity() noexcept { return enoki::half::from_binary(0x7C00); }
    static enoki::half quiet_NaN() noexcept { return enoki::half::from_binary(0x7FFF); }
    static enoki::half signaling_NaN() noexcept { return enoki::half::from_binary(0x7DFF); }
    static enoki::half denorm_min() noexcept { return enoki::half::from_binary(0x0001); }
};

NAMESPACE_END(std)

