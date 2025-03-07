/*
    enoki/array_neon.h -- Packed SIMD array (ARM NEON specialization)

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyrighe (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include "array_generic.h"

NAMESPACE_BEGIN(enoki)
NAMESPACE_BEGIN(detail)

template <> struct is_native<float, 4> : std::true_type { };
template <> struct is_native<float, 3> : std::true_type { };
template <typename Value> struct is_native<Value, 4, enable_if_int32_t<Value>> : std::true_type { };
template <typename Value> struct is_native<Value, 3, enable_if_int32_t<Value>> : std::true_type { };

#if defined(ENOKI_ARM_64)
    template <> struct is_native<double, 2> : std::true_type { };
    template <typename Value>    struct is_native<Value, 2, enable_if_int64_t<Value>> : std::true_type { };
#endif

static constexpr uint64_t arm_shuffle_helper_(int i) {
    if (i == 0)
        return 0x03020100;
    else if (i == 1)
        return 0x07060504;
    else if (i == 2)
        return 0x0B0A0908;
    else
        return 0x0F0E0D0C;
}

NAMESPACE_END(detail)

ENOKI_INLINE uint64x2_t vmvnq_u64(uint64x2_t a) {
    return vreinterpretq_u64_u32(vmvnq_u32(vreinterpretq_u32_u64(a)));
}

ENOKI_INLINE int64x2_t vmvnq_s64(int64x2_t a) {
    return vreinterpretq_s64_s32(vmvnq_s32(vreinterpretq_s32_s64(a)));
}

/// Partial overload of StaticArrayImpl using ARM NEON intrinsics (single precision)
template <bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<float, 4, IsMask_, Derived_>
  : StaticArrayBase<float, 4, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(float, 4, float32x4_t)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(vdupq_n_f32(value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3)
        : m{v0, v1, v2, v3} { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(float) : m(a.derived().m) { }
    ENOKI_CONVERT(int32_t) : m(vcvtq_f32_s32(vreinterpretq_s32_u32(a.derived().m))) { }
    ENOKI_CONVERT(uint32_t) : m(vcvtq_f32_u32(a.derived().m)) { }
    ENOKI_CONVERT(half) : m(vcvt_f32_f16(vld1_f16((const __fp16 *) a.data()))) { }
#if defined(ENOKI_ARM_64)
    ENOKI_CONVERT(double) : m(vcvtx_high_f32_f64(vcvtx_f32_f64(low(a).m), high(a).m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

#define ENOKI_REINTERPRET_BOOL(type, target)                   \
    ENOKI_REINTERPRET(type) {                                  \
        m = vreinterpretq_##target##_u32(uint32x4_t {          \
            reinterpret_array<uint32_t>(a.derived().coeff(0)), \
            reinterpret_array<uint32_t>(a.derived().coeff(1)), \
            reinterpret_array<uint32_t>(a.derived().coeff(2)), \
            reinterpret_array<uint32_t>(a.derived().coeff(3))  \
        });                                                    \
    }

    ENOKI_REINTERPRET(float) : m(a.derived().m) { }
    ENOKI_REINTERPRET(int32_t) : m(vreinterpretq_f32_u32(a.derived().m)) { }
    ENOKI_REINTERPRET(uint32_t) : m(vreinterpretq_f32_u32(a.derived().m)) { }
#if defined(ENOKI_ARM_64)
    ENOKI_REINTERPRET(int64_t) : m(vreinterpretq_f32_u32(vcombine_u32(vmovn_u64(low(a).m), vmovn_u64(high(a).m)))) { }
    ENOKI_REINTERPRET(uint64_t) : m(vreinterpretq_f32_u32(vcombine_u32(vmovn_u64(low(a).m), vmovn_u64(high(a).m)))) { }
    ENOKI_REINTERPRET(double) : m(vreinterpretq_f32_u32(vcombine_u32(
        vmovn_u64(vreinterpretq_u64_f64(low(a).m)),
        vmovn_u64(vreinterpretq_u64_f64(high(a).m))))) { }
#else
    ENOKI_REINTERPRET_BOOL(int64_t, f32)
    ENOKI_REINTERPRET_BOOL(uint64_t, f32)
    ENOKI_REINTERPRET_BOOL(double, f32)
#endif

    ENOKI_REINTERPRET_BOOL(bool, f32)

#undef ENOKI_REINTERPRET_BOOL

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m{a1.coeff(0), a1.coeff(1), a2.coeff(0), a2.coeff(1)} { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0), coeff(1)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(2), coeff(3)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return vaddq_f32(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return vsubq_f32(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return vmulq_f32(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const {
        #if defined(ENOKI_ARM_64)
            return vdivq_f32(m, a.m);
        #else
            return Base::div_(a);
        #endif
    }

#if defined(ENOKI_ARM_FMA)
    ENOKI_INLINE Derived fmadd_(Ref b, Ref c) const { return vfmaq_f32(c.m, m, b.m); }
    ENOKI_INLINE Derived fnmadd_(Ref b, Ref c) const { return vfmsq_f32(c.m, m, b.m); }
    ENOKI_INLINE Derived fmsub_(Ref b, Ref c) const { return vfmaq_f32(vnegq_f32(c.m), m, b.m); }
    ENOKI_INLINE Derived fnmsub_(Ref b, Ref c) const { return vfmsq_f32(vnegq_f32(c.m), m, b.m); }
#else
    ENOKI_INLINE Derived fmadd_(Ref b, Ref c) const { return vmlaq_f32(c.m, m, b.m); }
    ENOKI_INLINE Derived fnmadd_(Ref b, Ref c) const { return vmlsq_f32(c.m, m, b.m); }
    ENOKI_INLINE Derived fmsub_(Ref b, Ref c) const { return vmlaq_f32(vnegq_f32(c.m), m, b.m); }
    ENOKI_INLINE Derived fnmsub_(Ref b, Ref c) const { return vmlsq_f32(vnegq_f32(c.m), m, b.m); }
#endif

    template <typename T> ENOKI_INLINE Derived or_ (const T &a) const { return vreinterpretq_f32_s32(vorrq_s32(vreinterpretq_s32_f32(m), vreinterpretq_s32_f32(a.m))); }
    template <typename T> ENOKI_INLINE Derived and_(const T &a) const { return vreinterpretq_f32_s32(vandq_s32(vreinterpretq_s32_f32(m), vreinterpretq_s32_f32(a.m))); }
    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const { return vreinterpretq_f32_s32(vbicq_s32(vreinterpretq_s32_f32(m), vreinterpretq_s32_f32(a.m))); }
    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const { return vreinterpretq_f32_s32(veorq_s32(vreinterpretq_s32_f32(m), vreinterpretq_s32_f32(a.m))); }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f32_u32(vcltq_f32(m, a.m))); }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f32_u32(vcgtq_f32(m, a.m))); }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f32_u32(vcleq_f32(m, a.m))); }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f32_u32(vcgeq_f32(m, a.m))); }

    ENOKI_INLINE auto eq_ (Ref a) const {
        if constexpr (!IsMask_)
            return mask_t<Derived>(vreinterpretq_f32_u32(vceqq_f32(m, a.m)));
        else
            return mask_t<Derived>(vceqq_u32(vreinterpretq_f32_u32(m), vreinterpretq_f32_u32(a.m)));
    }

    ENOKI_INLINE auto neq_ (Ref a) const {
        if constexpr (!IsMask_)
            return mask_t<Derived>(vreinterpretq_f32_u32(vmvnq_u32(vceqq_f32(m, a.m))));
        else
            return mask_t<Derived>(vmvnq_u32(vceqq_u32(vreinterpretq_f32_u32(m), vreinterpretq_f32_u32(a.m))));
    }

    ENOKI_INLINE Derived abs_()      const { return vabsq_f32(m); }
    ENOKI_INLINE Derived neg_()      const { return vnegq_f32(m); }
    ENOKI_INLINE Derived not_()      const { return vreinterpretq_f32_s32(vmvnq_s32(vreinterpretq_s32_f32(m))); }

    ENOKI_INLINE Derived min_(Ref b) const { return vminq_f32(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return vmaxq_f32(b.m, m); }

#if defined(ENOKI_ARM_64)
    ENOKI_INLINE Derived round_()    const { return vrndnq_f32(m);     }
    ENOKI_INLINE Derived floor_()    const { return vrndmq_f32(m);     }
    ENOKI_INLINE Derived ceil_()     const { return vrndpq_f32(m);     }
#endif

    ENOKI_INLINE Derived sqrt_() const {
        #if defined(ENOKI_ARM_64)
            return vsqrtq_f32(m);
        #else
            const float32x4_t inf = vdupq_n_f32(std::numeric_limits<float>::infinity());
            float32x4_t r = vrsqrteq_f32(m);
            uint32x4_t inf_or_zero = vorrq_u32(vceqq_f32(r, inf), vceqq_f32(m, inf));
            r = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(r, r), m));
            r = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(r, r), m));
            r = vmulq_f32(r, m);
            return vbslq_f32(inf_or_zero, m, r);
        #endif
    }

    ENOKI_INLINE Derived rcp_() const {
        float32x4_t r = vrecpeq_f32(m);
        r = vmulq_f32(r, vrecpsq_f32(r, m));
        r = vmulq_f32(r, vrecpsq_f32(r, m));
        return r;
    }

    ENOKI_INLINE Derived rsqrt_() const {
        float32x4_t r = vrsqrteq_f32(m);
        r = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(r, r), m));
        r = vmulq_f32(r, vrsqrtsq_f32(vmulq_f32(r, r), m));
        return r;
    }

    template <typename Mask_>
    static ENOKI_INLINE Derived select_(const Mask_ &m, Ref t, Ref f) {
        return vbslq_f32(vreinterpretq_u32_f32(m.m), t.m, f.m);
    }

    template <int I0, int I1, int I2, int I3>
    ENOKI_INLINE Derived shuffle_() const {
        /// Based on https://stackoverflow.com/a/32537433/1130282
        switch (I3 + I2*10 + I1*100 + I0*1000) {
            case 0123: return m;
            case 0000: return vdupq_lane_f32(vget_low_f32(m), 0);
            case 1111: return vdupq_lane_f32(vget_low_f32(m), 1);
            case 2222: return vdupq_lane_f32(vget_high_f32(m), 0);
            case 3333: return vdupq_lane_f32(vget_high_f32(m), 1);
            case 1032: return vrev64q_f32(m);
            case 0101: { float32x2_t vt = vget_low_f32(m); return vcombine_f32(vt, vt); }
            case 2323: { float32x2_t vt = vget_high_f32(m); return vcombine_f32(vt, vt); }
            case 1010: { float32x2_t vt = vrev64_f32(vget_low_f32(m)); return vcombine_f32(vt, vt); }
            case 3232: { float32x2_t vt = vrev64_f32(vget_high_f32(m)); return vcombine_f32(vt, vt); }
            case 0132: return vcombine_f32(vget_low_f32(m), vrev64_f32(vget_high_f32(m)));
            case 1023: return vcombine_f32(vrev64_f32(vget_low_f32(m)), vget_high_f32(m));
            case 2310: return vcombine_f32(vget_high_f32(m), vrev64_f32(vget_low_f32(m)));
            case 3201: return vcombine_f32(vrev64_f32(vget_high_f32(m)), vget_low_f32(m));
            case 3210: return vcombine_f32(vrev64_f32(vget_high_f32(m)), vrev64_f32(vget_low_f32(m)));
#if defined(ENOKI_ARM_64)
            case 0022: return vtrn1q_f32(m, m);
            case 1133: return vtrn2q_f32(m, m);
            case 0011: return vzip1q_f32(m, m);
            case 2233: return vzip2q_f32(m, m);
            case 0202: return vuzp1q_f32(m, m);
            case 1313: return vuzp2q_f32(m, m);
#endif
            case 1230: return vextq_f32(m, m, 1);
            case 2301: return vextq_f32(m, m, 2);
            case 3012: return vextq_f32(m, m, 3);

            default: {
                constexpr uint64_t prec0 = detail::arm_shuffle_helper_(I0) |
                                          (detail::arm_shuffle_helper_(I1) << 32);
                constexpr uint64_t prec1 = detail::arm_shuffle_helper_(I2) |
                                          (detail::arm_shuffle_helper_(I3) << 32);

                uint8x8x2_t tbl;
                tbl.val[0] = vreinterpret_u8_f32(vget_low_f32(m));
                tbl.val[1] = vreinterpret_u8_f32(vget_high_f32(m));

                uint8x8_t idx1 = vreinterpret_u8_u32(vcreate_u32(prec0));
                uint8x8_t idx2 = vreinterpret_u8_u32(vcreate_u32(prec1));

                float32x2_t l = vreinterpret_f32_u8(vtbl2_u8(tbl, idx1));
                float32x2_t h = vreinterpret_f32_u8(vtbl2_u8(tbl, idx2));

                return vcombine_f32(l, h);
            }
        }
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_ARM_64)
    ENOKI_INLINE Value hmax_() const { return vmaxvq_f32(m); }
    ENOKI_INLINE Value hmin_() const { return vminvq_f32(m); }
    ENOKI_INLINE Value hsum_() const { return vaddvq_f32(m); }

    bool all_() const {
        if constexpr (Derived::Size == 4)
            return vmaxvq_s32(vreinterpretq_s32_f32(m)) < 0;
        else
            return Base::all_();
    }

    bool any_() const {
        if constexpr (Derived::Size == 4)
            return vminvq_s32(vreinterpretq_s32_f32(m)) < 0;
        else
            return Base::any_();
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        vst1q_f32((Value *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        vst1q_f32((Value *) ptr, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 16 == 0);
        return vld1q_f32((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return vld1q_f32((const Value *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return vdupq_n_f32(0.f); }

    //! @}
    // -----------------------------------------------------------------------
};

#if defined(ENOKI_ARM_64)
/// Partial overload of StaticArrayImpl using ARM NEON intrinsics (double precision)
template <bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<double, 2, IsMask_, Derived_>
  : StaticArrayBase<double, 2, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(double, 2, float64x2_t)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(vdupq_n_f64(value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1) : m{v0, v1} { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(double) : m(a.derived().m) { }
    ENOKI_CONVERT(int64_t) : m(vcvtq_f64_s64(vreinterpretq_s64_u64(a.derived().m))) { }
    ENOKI_CONVERT(uint64_t) : m(vcvtq_f64_u64(a.derived().m)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(double) : m(a.derived().m) { }
    ENOKI_REINTERPRET(int64_t) : m(vreinterpretq_f64_u64(a.derived().m)) { }
    ENOKI_REINTERPRET(uint64_t) : m(vreinterpretq_f64_u64(a.derived().m)) { }
    ENOKI_REINTERPRET(bool) {
        m = vreinterpretq_f64_u64(uint64x2_t {
            reinterpret_array<uint64_t>(a.derived().coeff(0)),
            reinterpret_array<uint64_t>(a.derived().coeff(1))
        });
    }
    ENOKI_REINTERPRET(float) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_f64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    ENOKI_REINTERPRET(int32_t) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_f64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    ENOKI_REINTERPRET(uint32_t) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_f64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m{a1.coeff(0), a2.coeff(0)} { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(1)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return vaddq_f64(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return vsubq_f64(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return vmulq_f64(m, a.m); }
    ENOKI_INLINE Derived div_(Ref a) const { return vdivq_f64(m, a.m); }

#if defined(ENOKI_ARM_FMA)
    ENOKI_INLINE Derived fmadd_(Ref b, Ref c) const { return vfmaq_f64(c.m, m, b.m); }
    ENOKI_INLINE Derived fnmadd_(Ref b, Ref c) const { return vfmsq_f64(c.m, m, b.m); }
    ENOKI_INLINE Derived fmsub_(Ref b, Ref c) const { return vfmaq_f64(vnegq_f64(c.m), m, b.m); }
    ENOKI_INLINE Derived fnmsub_(Ref b, Ref c) const { return vfmsq_f64(vnegq_f64(c.m), m, b.m); }
#else
    ENOKI_INLINE Derived fmadd_(Ref b, Ref c) const { return vmlaq_f64(c.m, m, b.m); }
    ENOKI_INLINE Derived fnmadd_(Ref b, Ref c) const { return vmlsq_f64(c.m, m, b.m); }
    ENOKI_INLINE Derived fmsub_(Ref b, Ref c) const { return vmlaq_f64(vnegq_f64(c.m), m, b.m); }
    ENOKI_INLINE Derived fnmsub_(Ref b, Ref c) const { return vmlsq_f64(vnegq_f64(c.m), m, b.m); }
#endif

    template <typename T> ENOKI_INLINE Derived or_ (const T &a) const { return vreinterpretq_f64_s64(vorrq_s64(vreinterpretq_s64_f64(m), vreinterpretq_s64_f64(a.m))); }
    template <typename T> ENOKI_INLINE Derived and_(const T &a) const { return vreinterpretq_f64_s64(vandq_s64(vreinterpretq_s64_f64(m), vreinterpretq_s64_f64(a.m))); }
    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const { return vreinterpretq_f64_s64(vbicq_s64(vreinterpretq_s64_f64(m), vreinterpretq_s64_f64(a.m))); }
    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const { return vreinterpretq_f64_s64(veorq_s64(vreinterpretq_s64_f64(m), vreinterpretq_s64_f64(a.m))); }

    ENOKI_INLINE auto lt_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f64_u64(vcltq_f64(m, a.m))); }
    ENOKI_INLINE auto gt_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f64_u64(vcgtq_f64(m, a.m))); }
    ENOKI_INLINE auto le_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f64_u64(vcleq_f64(m, a.m))); }
    ENOKI_INLINE auto ge_ (Ref a) const { return mask_t<Derived>(vreinterpretq_f64_u64(vcgeq_f64(m, a.m))); }

    ENOKI_INLINE auto eq_ (Ref a) const {
        if constexpr (!IsMask_)
            return mask_t<Derived>(vreinterpretq_f64_u64(vceqq_f64(m, a.m)));
        else
            return mask_t<Derived>(vceqq_u64(vreinterpretq_f64_u64(m), vreinterpretq_f64_u64(a.m)));
    }

    ENOKI_INLINE auto neq_ (Ref a) const {
        if constexpr (!IsMask_)
            return mask_t<Derived>(vreinterpretq_f64_u64(vmvnq_u64(vceqq_f64(m, a.m))));
        else
            return mask_t<Derived>(vmvnq_u64(vceqq_u64(vreinterpretq_f64_u64(m), vreinterpretq_f64_u64(a.m))));
    }

    ENOKI_INLINE Derived abs_()      const { return vabsq_f64(m); }
    ENOKI_INLINE Derived neg_()      const { return vnegq_f64(m); }
    ENOKI_INLINE Derived not_()      const { return vreinterpretq_f64_s64(vmvnq_s64(vreinterpretq_s64_f64(m))); }

    ENOKI_INLINE Derived min_(Ref b) const { return vminq_f64(b.m, m); }
    ENOKI_INLINE Derived max_(Ref b) const { return vmaxq_f64(b.m, m); }

#if defined(ENOKI_ARM_64)
    ENOKI_INLINE Derived sqrt_()     const { return vsqrtq_f64(m);     }
    ENOKI_INLINE Derived round_()    const { return vrndnq_f64(m);     }
    ENOKI_INLINE Derived floor_()    const { return vrndmq_f64(m);     }
    ENOKI_INLINE Derived ceil_()     const { return vrndpq_f64(m);     }
#endif

    ENOKI_INLINE Derived rcp_() const {
        float64x2_t r = vrecpeq_f64(m);
        r = vmulq_f64(r, vrecpsq_f64(r, m));
        r = vmulq_f64(r, vrecpsq_f64(r, m));
        r = vmulq_f64(r, vrecpsq_f64(r, m));
        return r;
    }

    ENOKI_INLINE Derived rsqrt_() const {
        float64x2_t r = vrsqrteq_f64(m);
        r = vmulq_f64(r, vrsqrtsq_f64(vmulq_f64(r, r), m));
        r = vmulq_f64(r, vrsqrtsq_f64(vmulq_f64(r, r), m));
        r = vmulq_f64(r, vrsqrtsq_f64(vmulq_f64(r, r), m));
        return r;
    }

    template <typename Mask_>
    static ENOKI_INLINE Derived select_(const Mask_ &m, Ref t, Ref f) {
        return vbslq_f64(vreinterpretq_u64_f64(m.m), t.m, f.m);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hmax_() const { return vmaxvq_f64(m); }
    ENOKI_INLINE Value hmin_() const { return vminvq_f64(m); }
    ENOKI_INLINE Value hsum_() const { return vaddvq_f64(m); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        vst1q_f64((Value *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        vst1q_f64((Value *) ptr, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        assert((uintptr_t) ptr % 16 == 0);
        return vld1q_f64((const Value *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return vld1q_f64((const Value *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return vdupq_n_f64(0.0); }

    //! @}
    // -----------------------------------------------------------------------
};
#endif

/// Partial overload of StaticArrayImpl using ARM NEON intrinsics (32-bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<Value_, 4, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayBase<Value_, 4, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 4, uint32x4_t)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(vdupq_n_u32((uint32_t) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1, Value v2, Value v3)
        : m{(uint32_t) v0, (uint32_t) v1, (uint32_t) v2, (uint32_t) v3} { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(int32_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint32_t) : m(a.derived().m) { }
    ENOKI_CONVERT(float) : m(std::is_signed_v<Value> ?
          vreinterpretq_u32_s32(vcvtq_s32_f32(a.derived().m))
        : vcvtq_u32_f32(a.derived().m)) { }
#if defined(ENOKI_ARM_64)
    ENOKI_CONVERT(int64_t) : m(vmovn_high_u64(vmovn_u64(low(a).m), high(a).m)) { }
    ENOKI_CONVERT(uint64_t) : m(vmovn_high_u64(vmovn_u64(low(a).m), high(a).m)) { }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

#define ENOKI_REINTERPRET_BOOL(type, target)                   \
    ENOKI_REINTERPRET(type) {                                  \
        m = uint32x4_t {                                       \
            reinterpret_array<uint32_t>(a.derived().coeff(0)), \
            reinterpret_array<uint32_t>(a.derived().coeff(1)), \
            reinterpret_array<uint32_t>(a.derived().coeff(2)), \
            reinterpret_array<uint32_t>(a.derived().coeff(3))  \
        };                                                     \
    }

    ENOKI_REINTERPRET(int32_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint32_t) : m(a.derived().m) { }
#if defined(ENOKI_ARM_64)
    ENOKI_REINTERPRET(int64_t) : m(vcombine_u32(vmovn_u64(low(a).m), vmovn_u64(high(a).m))) { }
    ENOKI_REINTERPRET(uint64_t) : m(vcombine_u32(vmovn_u64(low(a).m), vmovn_u64(high(a).m))) { }
    ENOKI_REINTERPRET(double) : m(vcombine_u32(
        vmovn_u64(vreinterpretq_u64_f64(low(a).m)),
        vmovn_u64(vreinterpretq_u64_f64(high(a).m)))) { }
#else
    ENOKI_REINTERPRET_BOOL(int64_t, u32)
    ENOKI_REINTERPRET_BOOL(uint64_t, u32)
    ENOKI_REINTERPRET_BOOL(double, u32)
#endif
    ENOKI_REINTERPRET(float) : m(vreinterpretq_u32_f32(a.derived().m)) { }
    ENOKI_REINTERPRET_BOOL(bool, u32)

#undef ENOKI_REINTERPRET_BOOL

    //! @}
    // -----------------------------------------------------------------------


    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m{(uint32_t) a1.coeff(0), (uint32_t) a1.coeff(1), (uint32_t) a2.coeff(0), (uint32_t) a2.coeff(1)} { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0), coeff(1)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(2), coeff(3)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return vaddq_u32(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return vsubq_u32(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a) const { return vmulq_u32(m, a.m); }

    template <typename T> ENOKI_INLINE Derived or_ (const T &a) const { return vorrq_u32(m, a.m); }
    template <typename T> ENOKI_INLINE Derived and_(const T &a) const { return vandq_u32(m, a.m); }
    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const { return vbicq_u32(m, a.m); }
    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const { return veorq_u32(m, a.m); }

    ENOKI_INLINE auto lt_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcltq_s32(vreinterpretq_s32_u32(m), vreinterpretq_s32_u32(a.m)));
        else
            return mask_t<Derived>(vcltq_u32(m, a.m));
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcgtq_s32(vreinterpretq_s32_u32(m), vreinterpretq_s32_u32(a.m)));
        else
            return mask_t<Derived>(vcgtq_u32(m, a.m));
    }

    ENOKI_INLINE auto le_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcleq_s32(vreinterpretq_s32_u32(m), vreinterpretq_s32_u32(a.m)));
        else
            return mask_t<Derived>(vcleq_u32(m, a.m));
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcgeq_s32(vreinterpretq_s32_u32(m), vreinterpretq_s32_u32(a.m)));
        else
            return mask_t<Derived>(vcgeq_u32(m, a.m));
    }

    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>(vceqq_u32(m, a.m)); }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>(vmvnq_u32(vceqq_u32(m, a.m))); }

    ENOKI_INLINE Derived abs_() const {
        if (!std::is_signed<Value>())
            return m;
        return vreinterpretq_u32_s32(vabsq_s32(vreinterpretq_s32_u32(m)));
    }

    ENOKI_INLINE Derived neg_() const {
        return vreinterpretq_u32_s32(vnegq_s32(vreinterpretq_s32_u32(m)));
    }

    ENOKI_INLINE Derived not_()      const { return vmvnq_u32(m); }

    ENOKI_INLINE Derived max_(Ref b) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u32_s32(vmaxq_s32(vreinterpretq_s32_u32(b.m), vreinterpretq_s32_u32(m)));
        else
            return vmaxq_u32(b.m, m);
    }

    ENOKI_INLINE Derived min_(Ref b) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u32_s32(vminq_s32(vreinterpretq_s32_u32(b.m), vreinterpretq_s32_u32(m)));
        else
            return vminq_u32(b.m, m);
    }

    template <typename Mask_>
    static ENOKI_INLINE Derived select_(const Mask_ &m, Ref t, Ref f) {
        return vbslq_u32(m.m, t.m, f.m);
    }

    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        if constexpr (Imm == 0) {
            return derived();
        } else {
            if constexpr (std::is_signed_v<Value>)
                return vreinterpretq_u32_s32(
                    vshrq_n_s32(vreinterpretq_s32_u32(m), (int) Imm));
            else
                return vshrq_n_u32(m, (int) Imm);
        }
    }

    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        if constexpr (Imm == 0)
            return derived();
        else
            return vshlq_n_u32(m, (int) Imm);
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u32_s32(
                vshlq_s32(vreinterpretq_s32_u32(m), vdupq_n_s32(-(int) k)));
        else
            return vshlq_u32(m, vdupq_n_s32(-(int) k));
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return vshlq_u32(m, vdupq_n_s32((int) k));
    }

    ENOKI_INLINE Derived sr_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u32_s32(
                vshlq_s32(vreinterpretq_s32_u32(m),
                          vnegq_s32(vreinterpretq_s32_u32(a.m))));
        else
            return vshlq_u32(m, vnegq_s32(vreinterpretq_s32_u32(a.m)));
    }

    ENOKI_INLINE Derived sl_(Ref a) const {
        return vshlq_u32(m, vreinterpretq_s32_u32(a.m));
    }

#if defined(ENOKI_ARM_64)
    ENOKI_INLINE Derived mulhi_(Ref a) const {
    uint32x4_t ll, hh;
        if constexpr (std::is_signed_v<Value>) {
            int64x2_t l = vmull_s32(vreinterpret_s32_u32(vget_low_u32(m)),
                                    vreinterpret_s32_u32(vget_low_u32(a.m)));

            int64x2_t h = vmull_high_s32(vreinterpretq_s32_u32(m),
                                         vreinterpretq_s32_u32(a.m));

            ll = vreinterpretq_u32_s64(l);
            hh = vreinterpretq_u32_s64(h);
        } else {
            uint64x2_t l = vmull_u32(vget_low_u32(m),
                                     vget_low_u32(a.m));

            uint64x2_t h = vmull_high_u32(m, a.m);

            ll = vreinterpretq_u32_u64(l);
            hh = vreinterpretq_u32_u64(h);
        }
        return vuzp2q_u32(ll, hh);
    }
#endif

    ENOKI_INLINE Derived lzcnt_() const { return vclzq_u32(m); }
    ENOKI_INLINE Derived tzcnt_() const { return Value(32) - lzcnt(~derived() & (derived() - Value(1))); }
    ENOKI_INLINE Derived popcnt_() const { return vpaddlq_u16(vpaddlq_u8(vcntq_u8(vreinterpretq_u8_u32(m)))); }

    template <int I0, int I1, int I2, int I3>
    ENOKI_INLINE Derived shuffle_() const {
        /// Based on https://stackoverflow.com/a/32537433/1130282
        switch (I3 + I2*10 + I1*100 + I0*1000) {
            case 0123: return m;
            case 0000: return vdupq_lane_u32(vget_low_u32(m), 0);
            case 1111: return vdupq_lane_u32(vget_low_u32(m), 1);
            case 2222: return vdupq_lane_u32(vget_high_u32(m), 0);
            case 3333: return vdupq_lane_u32(vget_high_u32(m), 1);
            case 1032: return vrev64q_u32(m);
            case 0101: { uint32x2_t vt = vget_low_u32(m); return vcombine_u32(vt, vt); }
            case 2323: { uint32x2_t vt = vget_high_u32(m); return vcombine_u32(vt, vt); }
            case 1010: { uint32x2_t vt = vrev64_u32(vget_low_u32(m)); return vcombine_u32(vt, vt); }
            case 3232: { uint32x2_t vt = vrev64_u32(vget_high_u32(m)); return vcombine_u32(vt, vt); }
            case 0132: return vcombine_u32(vget_low_u32(m), vrev64_u32(vget_high_u32(m)));
            case 1023: return vcombine_u32(vrev64_u32(vget_low_u32(m)), vget_high_u32(m));
            case 2310: return vcombine_u32(vget_high_u32(m), vrev64_u32(vget_low_u32(m)));
            case 3201: return vcombine_u32(vrev64_u32(vget_high_u32(m)), vget_low_u32(m));
            case 3210: return vcombine_u32(vrev64_u32(vget_high_u32(m)), vrev64_u32(vget_low_u32(m)));
#if defined(ENOKI_ARM_64)
            case 0022: return vtrn1q_u32(m, m);
            case 1133: return vtrn2q_u32(m, m);
            case 0011: return vzip1q_u32(m, m);
            case 2233: return vzip2q_u32(m, m);
            case 0202: return vuzp1q_u32(m, m);
            case 1313: return vuzp2q_u32(m, m);
#endif
            case 1230: return vextq_u32(m, m, 1);
            case 2301: return vextq_u32(m, m, 2);
            case 3012: return vextq_u32(m, m, 3);

            default: {
                constexpr uint64_t prec0 = detail::arm_shuffle_helper_(I0) |
                                          (detail::arm_shuffle_helper_(I1) << 32);
                constexpr uint64_t prec1 = detail::arm_shuffle_helper_(I2) |
                                          (detail::arm_shuffle_helper_(I3) << 32);

                uint8x8x2_t tbl;
                tbl.val[0] = vreinterpret_u8_u32(vget_low_u32(m));
                tbl.val[1] = vreinterpret_u8_u32(vget_high_u32(m));

                uint8x8_t idx1 = vreinterpret_u8_u32(vcreate_u32(prec0));
                uint8x8_t idx2 = vreinterpret_u8_u32(vcreate_u32(prec1));

                uint32x2_t l = vreinterpret_u32_u8(vtbl2_u8(tbl, idx1));
                uint32x2_t h = vreinterpret_u32_u8(vtbl2_u8(tbl, idx2));

                return vcombine_u32(l, h);
            }
        }
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }


    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

#if defined(ENOKI_ARM_64)
    ENOKI_INLINE Value hmax_() const {
        if constexpr (std::is_signed_v<Value>)
            return Value(vmaxvq_s32(vreinterpretq_s32_u32(m)));
        else
            return Value(vmaxvq_u32(m));
    }

    ENOKI_INLINE Value hmin_() const {
        if constexpr (std::is_signed_v<Value>)
            return Value(vminvq_s32(vreinterpretq_s32_u32(m)));
        else
            return Value(vminvq_u32(m));
    }

    ENOKI_INLINE Value hsum_() const { return Value(vaddvq_u32(m)); }

    bool all_() const {
        if constexpr (Derived::Size == 4)
            return vmaxvq_s32(vreinterpretq_s32_u32(m)) < 0;
        else
            return Base::all_();
    }

    bool any_() const {
        if constexpr (Derived::Size == 4)
            return vminvq_s32(vreinterpretq_s32_u32(m)) < 0;
        else
            return Base::any_();
    }
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        assert((uintptr_t) ptr % 16 == 0);
        vst1q_u32((uint32_t *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        vst1q_u32((uint32_t *) ptr, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return vld1q_u32((const uint32_t *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return vld1q_u32((const uint32_t *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return vdupq_n_u32(0); }

    //! @}
    // -----------------------------------------------------------------------
};

#if defined(ENOKI_ARM_64)
/// Partial overload of StaticArrayImpl using ARM NEON intrinsics (64-bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<Value_, 2, IsMask_, Derived_, enable_if_int64_t<Value_>>
  : StaticArrayBase<Value_, 2, IsMask_, Derived_> {
    ENOKI_NATIVE_ARRAY(Value_, 2, uint64x2_t)

    // -----------------------------------------------------------------------
    //! @{ \name Value constructors
    // -----------------------------------------------------------------------

    ENOKI_INLINE StaticArrayImpl(Value value) : m(vdupq_n_u64((uint64_t) value)) { }
    ENOKI_INLINE StaticArrayImpl(Value v0, Value v1) : m{(uint64_t) v0, (uint64_t) v1} { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Type converting constructors
    // -----------------------------------------------------------------------

    ENOKI_CONVERT(int64_t) : m(a.derived().m) { }
    ENOKI_CONVERT(uint64_t) : m(a.derived().m) { }
    ENOKI_CONVERT(double) : m(std::is_signed_v<Value> ?
          vreinterpretq_u64_s64(vcvtq_s64_f64(a.derived().m))
        : vcvtq_u64_f64(a.derived().m)) { }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Reinterpreting constructors, mask converters
    // -----------------------------------------------------------------------

    ENOKI_REINTERPRET(int64_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(uint64_t) : m(a.derived().m) { }
    ENOKI_REINTERPRET(double) : m(vreinterpretq_u64_f64(a.derived().m)) { }
    ENOKI_REINTERPRET(bool) {
        m = uint64x2_t {
            reinterpret_array<uint64_t>(a.derived().coeff(0)),
            reinterpret_array<uint64_t>(a.derived().coeff(1))
        };
    }
    ENOKI_REINTERPRET(float) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_u64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    ENOKI_REINTERPRET(int32_t) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_u64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    ENOKI_REINTERPRET(uint32_t) {
        auto v0 = memcpy_cast<uint32_t>(a.derived().coeff(0)),
             v1 = memcpy_cast<uint32_t>(a.derived().coeff(1));
        m = vreinterpretq_u64_u32(uint32x4_t { v0, v0, v1, v1 });
    }

    //! @}
    // -----------------------------------------------------------------------


    // -----------------------------------------------------------------------
    //! @{ \name Converting from/to half size vectors
    // -----------------------------------------------------------------------

    StaticArrayImpl(const Array1 &a1, const Array2 &a2)
        : m{(uint64_t) a1.coeff(0), (uint64_t) a2.coeff(0)} { }

    ENOKI_INLINE Array1 low_()  const { return Array1(coeff(0)); }
    ENOKI_INLINE Array2 high_() const { return Array2(coeff(1)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Derived add_(Ref a) const { return vaddq_u64(m, a.m); }
    ENOKI_INLINE Derived sub_(Ref a) const { return vsubq_u64(m, a.m); }
    ENOKI_INLINE Derived mul_(Ref a_) const {
#if 1
        // Native ARM instructions + cross-domain penalities still
        // seem to be faster than the NEON approach below
        return Derived(
            coeff(0) * a_.coeff(0),
            coeff(1) * a_.coeff(1)
        );
#else
        // inp: [ah0, al0, ah1, al1], [bh0, bl0, bh1, bl1]
        uint32x4_t a = vreinterpretq_u32_u64(m),
                   b = vreinterpretq_u32_u64(a_.m);

        // uzp: [al0, al1, bl0, bl1], [bh0, bh1, ah0, ah1]
        uint32x4_t l = vuzp1q_u32(a, b);
        uint32x4_t h = vuzp2q_u32(b, a);

        uint64x2_t accum = vmull_u32(vget_low_u32(l), vget_low_u32(h));
        accum = vmlal_high_u32(accum, h, l);
        accum = vshlq_n_u64(accum, 32);
        accum = vmlal_u32(accum, vget_low_u32(l), vget_high_u32(l));

        return accum;
#endif
    }

    template <typename T> ENOKI_INLINE Derived or_ (const T &a) const { return vorrq_u64(m, a.m); }
    template <typename T> ENOKI_INLINE Derived and_(const T &a) const { return vandq_u64(m, a.m); }
    template <typename T> ENOKI_INLINE Derived andnot_(const T &a) const { return vbicq_u64(m, a.m); }
    template <typename T> ENOKI_INLINE Derived xor_(const T &a) const { return veorq_u64(m, a.m); }

    ENOKI_INLINE auto lt_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcltq_s64(vreinterpretq_s64_u64(m), vreinterpretq_s64_u64(a.m)));
        else
            return mask_t<Derived>(vcltq_u64(m, a.m));
    }

    ENOKI_INLINE auto gt_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcgtq_s64(vreinterpretq_s64_u64(m), vreinterpretq_s64_u64(a.m)));
        else
            return mask_t<Derived>(vcgtq_u64(m, a.m));
    }

    ENOKI_INLINE auto le_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcleq_s64(vreinterpretq_s64_u64(m), vreinterpretq_s64_u64(a.m)));
        else
            return mask_t<Derived>(vcleq_u64(m, a.m));
    }

    ENOKI_INLINE auto ge_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return mask_t<Derived>(vcgeq_s64(vreinterpretq_s64_u64(m), vreinterpretq_s64_u64(a.m)));
        else
            return mask_t<Derived>(vcgeq_u64(m, a.m));
    }

    ENOKI_INLINE auto eq_ (Ref a) const { return mask_t<Derived>(vceqq_u64(m, a.m)); }
    ENOKI_INLINE auto neq_(Ref a) const { return mask_t<Derived>(vmvnq_u64(vceqq_u64(m, a.m))); }

    ENOKI_INLINE Derived abs_() const {
        if (!std::is_signed<Value>())
            return m;
        return vreinterpretq_u64_s64(vabsq_s64(vreinterpretq_s64_u64(m)));
    }

    ENOKI_INLINE Derived neg_() const {
        return vreinterpretq_u64_s64(vnegq_s64(vreinterpretq_s64_u64(m)));
    }

    ENOKI_INLINE Derived not_()      const { return vmvnq_u64(m); }

    ENOKI_INLINE Derived min_(Ref b) const { return Derived(min(coeff(0), b.coeff(0)), min(coeff(1), b.coeff(1))); }
    ENOKI_INLINE Derived max_(Ref b) const { return Derived(max(coeff(0), b.coeff(0)), max(coeff(1), b.coeff(1))); }

    template <typename Mask_>
    static ENOKI_INLINE Derived select_(const Mask_ &m, Ref t, Ref f) {
        return vbslq_u64(m.m, t.m, f.m);
    }

    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        if constexpr (Imm == 0) {
            return derived();
        } else {
            if constexpr (std::is_signed_v<Value>)
                return vreinterpretq_u64_s64(
                    vshrq_n_s64(vreinterpretq_s64_u64(m), (int) Imm));
            else
                return vshrq_n_u64(m, (int) Imm);
        }
    }

    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        if constexpr (Imm == 0)
            return derived();
        else
            return vshlq_n_u64(m, (int) Imm);
    }

    ENOKI_INLINE Derived sr_(size_t k) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u64_s64(
                vshlq_s64(vreinterpretq_s64_u64(m), vdupq_n_s64(-(int) k)));
        else
            return vshlq_u64(m, vdupq_n_s64(-(int) k));
    }

    ENOKI_INLINE Derived sl_(size_t k) const {
        return vshlq_u64(m, vdupq_n_s64((int) k));
    }

    ENOKI_INLINE Derived sr_(Ref a) const {
        if constexpr (std::is_signed_v<Value>)
            return vreinterpretq_u64_s64(
                vshlq_s64(vreinterpretq_s64_u64(m),
                          vnegq_s64(vreinterpretq_s64_u64(a.m))));
        else
            return vshlq_u64(m, vnegq_s64(vreinterpretq_s64_u64(a.m)));
    }

    ENOKI_INLINE Derived sl_(Ref a) const {
        return vshlq_u64(m, vreinterpretq_s64_u64(a.m));
    }

    ENOKI_INLINE Derived popcnt_() const {
        return vpaddlq_u32(
            vpaddlq_u16(vpaddlq_u8(vcntq_u8(vreinterpretq_u8_u64(m)))));
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hsum_() const { return Value(vaddvq_u64(m)); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization, loading/writing data
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        vst1q_u64((uint64_t *) ENOKI_ASSUME_ALIGNED(ptr, 16), m);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        vst1q_u64((uint64_t *) ptr, m);
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return vld1q_u64((const uint64_t *) ENOKI_ASSUME_ALIGNED(ptr, 16));
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        return vld1q_u64((const uint64_t *) ptr);
    }

    static ENOKI_INLINE Derived zero_() { return vdupq_n_u64(0); }

    //! @}
    // -----------------------------------------------------------------------
};
#endif

/// Partial overload of StaticArrayImpl for the n=3 case (single precision)
template <bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<float, 3, IsMask_, Derived_>
  : StaticArrayImpl<float, 4, IsMask_, Derived_> {
    using Base = StaticArrayImpl<float, 4, IsMask_, Derived_>;

    ENOKI_DECLARE_3D_ARRAY(StaticArrayImpl)

    template <typename Derived2>
    ENOKI_INLINE StaticArrayImpl(
        const StaticArrayBase<half, 3, IsMask_, Derived2> &a) {
        float16x4_t value;
        memcpy(&value, a.data(), sizeof(uint16_t)*3);
        m = vcvt_f32_f16(value);
    }

    template <int I0, int I1, int I2>
    ENOKI_INLINE Derived shuffle_() const {
        return Derived(coeff(I0), coeff(I1), coeff(I2));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }


    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hmax_() const { return max(max(coeff(0), coeff(1)), coeff(2)); }
    ENOKI_INLINE Value hmin_() const { return min(min(coeff(0), coeff(1)), coeff(2)); }
    ENOKI_INLINE Value hsum_() const { return coeff(0) + coeff(1) + coeff(2); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Loading/writing data (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        memcpy(ptr, &m, sizeof(Value) * 3);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        store_(ptr);
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        Derived result;
        memcpy(&result.m, ptr, sizeof(Value) * 3);
        return result;
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return Base::load_unaligned_(ptr);
    }

    //! @}
    // -----------------------------------------------------------------------
};

/// Partial overload of StaticArrayImpl for the n=3 case (32 bit integers)
template <typename Value_, bool IsMask_, typename Derived_> struct ENOKI_MAY_ALIAS alignas(16)
    StaticArrayImpl<Value_, 3, IsMask_, Derived_, enable_if_int32_t<Value_>>
  : StaticArrayImpl<Value_, 4, IsMask_, Derived_> {
    using Base = StaticArrayImpl<Value_, 4, IsMask_, Derived_>;

    ENOKI_DECLARE_3D_ARRAY(StaticArrayImpl)

    template <int I0, int I1, int I2>
    ENOKI_INLINE Derived shuffle_() const {
        return Derived(coeff(I0), coeff(I1), coeff(I2));
    }

    template <typename Index>
    ENOKI_INLINE Derived shuffle_(const Index &index) const {
        return Base::shuffle_(index);
    }


    // -----------------------------------------------------------------------
    //! @{ \name Horizontal operations (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    ENOKI_INLINE Value hmax_() const { return max(max(coeff(0), coeff(1)), coeff(2)); }
    ENOKI_INLINE Value hmin_() const { return min(min(coeff(0), coeff(1)), coeff(2)); }
    ENOKI_INLINE Value hsum_() const { return coeff(0) + coeff(1) + coeff(2); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Loading/writing data (adapted for the n=3 case)
    // -----------------------------------------------------------------------

    using Base::load_;
    using Base::store_;
    using Base::load_unaligned_;
    using Base::store_unaligned_;

    ENOKI_INLINE void store_(void *ptr) const {
        memcpy(ptr, &m, sizeof(Value) * 3);
    }

    ENOKI_INLINE void store_unaligned_(void *ptr) const {
        store_(ptr);
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *ptr) {
        Derived result;
        memcpy(&result.m, ptr, sizeof(Value) * 3);
        return result;
    }

    static ENOKI_INLINE Derived load_(const void *ptr) {
        return Base::load_unaligned_(ptr);
    }

    //! @}
    // -----------------------------------------------------------------------
};

NAMESPACE_END(enoki)
