/*
    enoki/dynamic.h -- Dynamic heap-allocated array

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#pragma once

#include <enoki/array.h>

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#define ENOKI_DYNAMIC_H 1

NAMESPACE_BEGIN(enoki)

template <typename Packet_>
struct DynamicArrayReference : ArrayBase<value_t<Packet_>, DynamicArrayReference<Packet_>> {
    using Base = ArrayBase<value_t<Packet_>, DynamicArrayReference<Packet_>>;
    using Packet = Packet_;
    using ArrayType = DynamicArrayReference<array_t<Packet>>;
    using MaskType = DynamicArrayReference<mask_t<Packet>>;

    static constexpr size_t       PacketSize  = Packet::Size;
    static constexpr bool         IsMask      = Packet::IsMask;

    DynamicArrayReference(Packet *packets = nullptr) : m_packets(packets) { }

    ENOKI_INLINE Packet &packet(size_t i) {
        return ((Packet *) ENOKI_ASSUME_ALIGNED(m_packets, alignof(Packet)))[i];
    }

    ENOKI_INLINE const Packet &packet(size_t i) const {
        return ((const Packet *) ENOKI_ASSUME_ALIGNED(m_packets, alignof(Packet)))[i];
    }

    template <typename T>
    using ReplaceValue = DynamicArrayReference<replace_scalar_t<Packet, T>>;

private:
    Packet *m_packets;
};

template <typename Packet_, typename Derived_>
struct DynamicArrayImpl : ArrayBase<value_t<Packet_>, Derived_> {
    // -----------------------------------------------------------------------
    //! @{ \name Aliases and constants
    // -----------------------------------------------------------------------

    using Size                                = uint32_t;
    using Base                                = ArrayBase<value_t<Packet_>, Derived_>;
    using Packet                              = Packet_;
    using IndexPacket                         = uint_array_t<array_t<Packet_>, false>;
    using IndexScalar                         = scalar_t<IndexPacket>;
    using PacketHolder                        = std::unique_ptr<Packet[]>;

    static constexpr size_t       PacketSize  = Packet::Size;
    static constexpr bool         IsMask      = Packet::IsMask;

    using typename Base::Derived;
    using typename Base::Value;
    using typename Base::Scalar;
    using Base::derived;

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Constructors
    // -----------------------------------------------------------------------

    DynamicArrayImpl() = default;

    ENOKI_INLINE ~DynamicArrayImpl() {
        reset();
    }

    /// Initialize from a list of component values
    template <typename... Ts, enable_if_t<sizeof...(Ts) >= 2 &&
        std::conjunction_v<detail::is_constructible<Value, Ts>...>> = 0>
    ENOKI_INLINE DynamicArrayImpl(Ts... args) {
        Value storage[] = { (Value) args... };
        resize(sizeof...(Ts));
        memcpy(m_packets.get(), storage, sizeof(Value) * sizeof...(Ts));
    }

    DynamicArrayImpl(const DynamicArrayImpl &value) {
        operator=(value);
    }

    ENOKI_INLINE DynamicArrayImpl(DynamicArrayImpl &&value) {
        operator=(std::move(value));
    }

    template <typename Packet2, typename Derived2>
    DynamicArrayImpl(const DynamicArrayImpl<Packet2, Derived2> &value) {
        operator=(value);
    }

    template <typename Value2, typename Derived2>
    DynamicArrayImpl(const ArrayBase<Value2, Derived2> &value) {
        operator=(value);
    }

    template <typename Packet2, typename Derived2>
    DynamicArrayImpl(const DynamicArrayImpl<Packet2, Derived2> &other,
                     detail::reinterpret_flag) {
        static_assert(Packet2::Size == Packet::Size, "Packet sizes must match!");
        resize(other.size());
        for (size_t i = 0; i < other.packets(); ++i)
            packet(i) = reinterpret_array<Packet>(other.packet(i));
    }

#if defined(__GNUC__)
// Don't be so noisy about sign conversion in constructor
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

    template <typename T = Packet_, enable_if_mask_t<T> = 0>
    DynamicArrayImpl(bool value, detail::reinterpret_flag) {
        resize(1);
        packet(0) = Packet(value);
    }

    template <typename T, enable_if_t<is_scalar_v<T>> = 0>
    DynamicArrayImpl(const T &value) {
        using S = std::conditional_t<IsMask, bool, Scalar>;
        resize(1);
        packet(0) = Packet((S) value);
    }

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Assignment operators
    // -----------------------------------------------------------------------

    template <typename T, enable_if_t<is_scalar_v<T>> = 0>
    ENOKI_NOINLINE DynamicArrayImpl &operator=(const T &value) {
        resize(1);
        packet(0) = Packet(value);
        return derived();
    }

    ENOKI_NOINLINE DynamicArrayImpl &operator=(const DynamicArrayImpl &other) {
        resize(other.size());
        memcpy(m_packets.get(), other.m_packets.get(),
               packets() * sizeof(Packet));
        return derived();
    }

    template <typename Packet2, typename Derived2>
    ENOKI_NOINLINE DynamicArrayImpl &operator=(const DynamicArrayImpl<Packet2, Derived2> &other) {
        static_assert(Packet2::Size == Packet::Size, "Packet sizes must match!");
        resize(other.size());
        for (size_t i = 0; i < other.packets(); ++i)
            packet(i) = Packet(other.packet(i));
        return derived();
    }

    template <typename Value2, typename Derived2>
    ENOKI_NOINLINE DynamicArrayImpl &operator=(const ArrayBase<Value2, Derived2> &other) {
        resize(other.derived().size());
        for (size_t i = 0; i < other.derived().size(); ++i)
            coeff(i) = other.derived().coeff(i);
        return derived();
    }

    ENOKI_INLINE DynamicArrayImpl &operator=(DynamicArrayImpl &&other) {
        m_packets.swap(other.m_packets);
        std::swap(m_packets_allocated, other.m_packets_allocated);
        std::swap(m_size, other.m_size);
        return derived();
    }

    void reset() {
        if (is_mapped()) {
            m_packets.release();
        } else if (m_packets.get()) {
            ENOKI_TRACK_DEALLOC(m_packets.get(), packets_allocated() * sizeof(Packet));
            m_packets.reset();
        }

        m_size = m_packets_allocated = 0;
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Functions to access the array contents
    // -----------------------------------------------------------------------

    bool is_mapped() const { return (m_packets_allocated & 0x80000000u) != 0; }
    size_t size() const { return (size_t) m_size; }
    size_t packets() const { return ((size_t) m_size + PacketSize - 1) / PacketSize; }
    size_t packets_allocated() const { return (size_t) (m_packets_allocated & 0x7fffffffu); }
    size_t capacity() const { return packets_allocated() * Packet::Size; }

    bool empty() const { return m_size == 0; }

    size_t nbytes() const {
        return packets_allocated() * sizeof(Packet) + sizeof(Derived);
    }

    ENOKI_INLINE const Value *data() const {
        return (const Value *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet));
    }

    ENOKI_INLINE Value *data() {
        return (Value *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet));
    }

    ENOKI_INLINE const Packet *packet_ptr() const {
        return (const Packet *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet));
    }

    ENOKI_INLINE Packet *packet_ptr() {
        return (Packet *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet));
    }

    ENOKI_INLINE decltype(auto) coeff(size_t i) {
        return m_packets[i / PacketSize].coeff(i % PacketSize);
    }

    ENOKI_INLINE decltype(auto) coeff(size_t i) const {
        return m_packets[i / PacketSize].coeff(i % PacketSize);
    }

    ENOKI_INLINE Packet &packet(size_t i) {
        #if !defined(NDEBUG) && !defined(ENOKI_DISABLE_RANGE_CHECK)
            if (i >= packets())
                throw std::out_of_range(
                    "DynamicArray: out of range access (tried to access packet " +
                    std::to_string(i) + " in an array of size " +
                    std::to_string(packets()) + ")");
        #endif
        return ((Packet *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet)))[i];
    }

    ENOKI_INLINE const Packet &packet(size_t i) const {
        #if !defined(NDEBUG) && !defined(ENOKI_DISABLE_RANGE_CHECK)
            if (i >= packets())
                throw std::out_of_range(
                    "DynamicArray: out of range access (tried to access packet " +
                    std::to_string(i) + " in an array of size " +
                    std::to_string(packets()) + ")");
        #endif
        return ((const Packet *) ENOKI_ASSUME_ALIGNED(m_packets.get(), alignof(Packet)))[i];
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Vertical array operations
    // -----------------------------------------------------------------------

    #define ENOKI_FWD_UNARY_OPERATION(name, Return, op)                      \
        auto name##_() const {                                               \
            Return result;                                                   \
            result.resize(size());                                           \
            auto p1 = packet_ptr();                                          \
            auto pr = result.packet_ptr();                                   \
            for (size_t i = 0, n = result.packets();                         \
                 i < n; ++i, ++p1, ++pr) {                                   \
                Packet a = *p1;                                              \
                *pr = op;                                                    \
            }                                                                \
            return result;                                                   \
        }

    #define ENOKI_FWD_UNARY_OPERATION_IMM(name, Return, op)                  \
        template <size_t Imm> auto name##_() const {                         \
            Return result;                                                   \
            result.resize(size());                                           \
            auto p1 = packet_ptr();                                          \
            auto pr = result.packet_ptr();                                   \
            for (size_t i = 0, n = result.packets();                         \
                 i < n; ++i, ++p1, ++pr) {                                   \
                Packet a = *p1;                                              \
                *pr = op;                                                    \
            }                                                                \
            return result;                                                   \
        }

    #define ENOKI_FWD_BINARY_OPERATION(name, Return, op)                     \
        template <typename T>                                                \
        auto name##_(const T &d) const {                                     \
            Return result;                                                   \
            result.resize_like(*this, d);                                    \
            auto p1 = packet_ptr();                                          \
            auto p2 = d.packet_ptr();                                        \
            auto pr = result.packet_ptr();                                   \
            size_t s1 = size() == 1 ? 0 : 1,                                 \
                   s2 = d.size() == 1 ? 0 : 1;                               \
            for (size_t i = 0, n = result.packets(); i < n;                  \
                 ++i, ++pr, p1 += s1, p2 += s2) {                            \
                auto a1 = *p1;                                               \
                auto a2 = *p2;                                               \
                *pr = op;                                                    \
            }                                                                \
            return result;                                                   \
        }

    #define ENOKI_FWD_BINARY_OPERATION_SIZE(name, Return, op)                \
        auto name##_(size_t a2) const {                                      \
            Return result;                                                   \
            result.resize_like(*this);                                       \
            auto p1 = packet_ptr();                                          \
            auto pr = result.packet_ptr();                                   \
            for (size_t i = 0, n = result.packets(); i < n;                  \
                 ++i, ++pr, p1++) {                                          \
                auto a1 = *p1;                                               \
                *pr = op;                                                    \
            }                                                                \
            return result;                                                   \
        }

    #define ENOKI_FWD_TERNARY_OPERATION(name, Return, op)                    \
        template <typename T1, typename T2>                                  \
        auto name##_(const T1 &d1, const T2 &d2) const {                     \
            Return result;                                                   \
            result.resize_like(*this, d1, d2);                               \
            auto p1 = packet_ptr();                                          \
            auto p2 = d1.packet_ptr();                                       \
            auto p3 = d2.packet_ptr();                                       \
            auto pr = result.packet_ptr();                                   \
            size_t s1 = size() == 1 ? 0 : 1,                                 \
                   s2 = d1.size() == 1 ? 0 : 1,                              \
                   s3 = d2.size() == 1 ? 0 : 1;                              \
            for (size_t i = 0, n = result.packets(); i < n;                  \
                 ++i, ++pr, p1 += s1, p2 += s2, p3 += s3) {                  \
                auto a1 = *p1;                                               \
                auto a2 = *p2;                                               \
                auto a3 = *p3;                                               \
                *pr = op;                                                    \
            }                                                                \
            return result;                                                   \
        }

    #define ENOKI_FWD_MASKED_OPERATION(name, expr)                           \
        template <typename Mask>                                             \
        void m##name##_(const Derived &e, const Mask &m) {                   \
            resize_like(*this, e, m);                                        \
            auto pr = packet_ptr();                                          \
            auto p1 = e.packet_ptr();                                        \
            auto p2 = m.packet_ptr();                                        \
            size_t s1 = e.size() == 1 ? 0 : 1,                               \
                   s2 = m.size() == 1 ? 0 : 1;                               \
            for (size_t i = 0, n = packets(); i < n;                         \
                 ++i, ++pr, p1 += s1, p2 += s2)                              \
                (*pr).m##name##_(*p1, *p2);                                  \
        }

    ENOKI_FWD_BINARY_OPERATION(add, Derived, a1 + a2)
    ENOKI_FWD_BINARY_OPERATION(sub, Derived, a1 - a2)
    ENOKI_FWD_BINARY_OPERATION(mul, Derived, a1 * a2)
    ENOKI_FWD_BINARY_OPERATION(div, Derived, a1 / a2)
    ENOKI_FWD_BINARY_OPERATION(mod, Derived, a1 % a2)
    ENOKI_FWD_BINARY_OPERATION(sl,  Derived, a1 << a2)
    ENOKI_FWD_BINARY_OPERATION(sr,  Derived, a1 >> a2)
    ENOKI_FWD_BINARY_OPERATION(rol, Derived, rol(a1, a2))
    ENOKI_FWD_BINARY_OPERATION(ror, Derived, ror(a1, a2))
    ENOKI_FWD_BINARY_OPERATION(mulhi, Derived, mulhi(a1, a2))

    ENOKI_FWD_BINARY_OPERATION_SIZE(sl, Derived, a1 << a2)
    ENOKI_FWD_BINARY_OPERATION_SIZE(sr, Derived, a1 >> a2)

    ENOKI_FWD_UNARY_OPERATION_IMM(sl,  Derived, sl<Imm>(a))
    ENOKI_FWD_UNARY_OPERATION_IMM(sr,  Derived, sr<Imm>(a))
    ENOKI_FWD_UNARY_OPERATION_IMM(rol, Derived, rol<Imm>(a))
    ENOKI_FWD_UNARY_OPERATION_IMM(ror, Derived, ror<Imm>(a))

    ENOKI_FWD_UNARY_OPERATION(lzcnt, Derived, lzcnt(a))
    ENOKI_FWD_UNARY_OPERATION(tzcnt, Derived, tzcnt(a))
    ENOKI_FWD_UNARY_OPERATION(popcnt, Derived, popcnt(a))

    ENOKI_FWD_BINARY_OPERATION(or,     Derived, a1 | a2)
    ENOKI_FWD_BINARY_OPERATION(and,    Derived, a1 & a2)
    ENOKI_FWD_BINARY_OPERATION(andnot, Derived, andnot(a1, a2))
    ENOKI_FWD_BINARY_OPERATION(xor,    Derived, a1 ^ a2)

    ENOKI_FWD_UNARY_OPERATION(not, Derived, ~a);
    ENOKI_FWD_UNARY_OPERATION(neg, Derived, -a);

    ENOKI_FWD_BINARY_OPERATION(eq,  mask_t<Derived>, eq (a1, a2))
    ENOKI_FWD_BINARY_OPERATION(neq, mask_t<Derived>, neq(a1, a2))
    ENOKI_FWD_BINARY_OPERATION(gt,  mask_t<Derived>, a1 > a2)
    ENOKI_FWD_BINARY_OPERATION(ge,  mask_t<Derived>, a1 >= a2)
    ENOKI_FWD_BINARY_OPERATION(lt,  mask_t<Derived>, a1 < a2)
    ENOKI_FWD_BINARY_OPERATION(le,  mask_t<Derived>, a1 <= a2)

    ENOKI_FWD_TERNARY_OPERATION(fmadd,    Derived, fmadd(a1, a2, a3))
    ENOKI_FWD_TERNARY_OPERATION(fmsub,    Derived, fmsub(a1, a2, a3))
    ENOKI_FWD_TERNARY_OPERATION(fnmadd,   Derived, fnmadd(a1, a2, a3))
    ENOKI_FWD_TERNARY_OPERATION(fnmsub,   Derived, fnmsub(a1, a2, a3))
    ENOKI_FWD_TERNARY_OPERATION(fmsubadd, Derived, fmsubadd(a1, a2, a3))
    ENOKI_FWD_TERNARY_OPERATION(fmaddsub, Derived, fmaddsub(a1, a2, a3))

    ENOKI_FWD_BINARY_OPERATION(min, Derived, min(a1, a2))
    ENOKI_FWD_BINARY_OPERATION(max, Derived, max(a1, a2))

    ENOKI_FWD_UNARY_OPERATION(abs,   Derived, abs(a));
    ENOKI_FWD_UNARY_OPERATION(ceil,  Derived, ceil(a));
    ENOKI_FWD_UNARY_OPERATION(floor, Derived, floor(a));
    ENOKI_FWD_UNARY_OPERATION(sqrt,  Derived, sqrt(a));
    ENOKI_FWD_UNARY_OPERATION(round, Derived, round(a));
    ENOKI_FWD_UNARY_OPERATION(trunc, Derived, trunc(a));

    ENOKI_FWD_UNARY_OPERATION(rsqrt, Derived, rsqrt(a));
    ENOKI_FWD_UNARY_OPERATION(rcp,   Derived, rcp(a));

    ENOKI_FWD_MASKED_OPERATION(assign, b)
    ENOKI_FWD_MASKED_OPERATION(add, a + b)
    ENOKI_FWD_MASKED_OPERATION(sub, a - b)
    ENOKI_FWD_MASKED_OPERATION(mul, a * b)
    ENOKI_FWD_MASKED_OPERATION(div, a / b)
    ENOKI_FWD_MASKED_OPERATION(or, a | b)
    ENOKI_FWD_MASKED_OPERATION(and, a & b)
    ENOKI_FWD_MASKED_OPERATION(xor, a ^ b)

    #undef ENOKI_FWD_UNARY_OPERATION
    #undef ENOKI_FWD_UNARY_OPERATION_IMM
    #undef ENOKI_FWD_BINARY_OPERATION
    #undef ENOKI_FWD_TERNARY_OPERATION
    #undef ENOKI_FWD_MASKED_OPERATION

    template <typename Mask>
    static Derived select_(const Mask &mask, const Derived &t, const Derived &f) {
        if (ENOKI_UNLIKELY(f.empty())) {
            if (all(mask))
                return t;
            else
                throw std::runtime_error(
                    "DynamicArray::select(): array for false branch is empty, "
                    "and some entries were referenced.");
        }

        if (ENOKI_UNLIKELY(t.empty())) {
            if (none(mask))
                return f;
            else
                throw std::runtime_error(
                    "DynamicArray::select(): array for true branch is empty, "
                    "and some entries were referenced.");
        }

        Derived result;
        result.resize_like(mask, t, f);
        size_t i1 = 0, i1i = mask.size() == 1 ? 0 : 1,
               i2 = 0, i2i = t.size() == 1 ? 0 : 1,
               i3 = 0, i3i = f.size() == 1 ? 0 : 1;

        for (size_t i = 0; i < result.packets();
             ++i, i1 += i1i, i2 += i2i, i3 += i3i) {
            result.packet(i) = select(mask.packet(i1), t.packet(i2), f.packet(i3));
        }
        return result;
    }

    template <size_t Stride, typename Index, typename Mask>
    static Derived gather_(const void *mem, const Index &index, const Mask &mask) {
        Derived result;
        result.resize_like(index, mask);
        size_t i1 = 0, i1i = index.size() == 1 ? 0 : 1,
               i2 = 0, i2i = mask.size() == 1 ? 0 : 1,
               i = 0;
        if (!result.empty()) {
            for (; i < result.packets() - (PacketSize > 1 ? 1 : 0); ++i, i1 += i1i, i2 += i2i)
                result.packet(i) = gather<Packet, Stride>(mem, index.packet(i1), mask.packet(i2));
            if constexpr (PacketSize > 1) {
                auto mask2 = arange<IndexPacket>() <= IndexScalar((result.size() - 1) % PacketSize);
                result.packet(i) = gather<Packet, Stride>(mem, index.packet(i1), mask.packet(i2) & mask2);
                if (result.size() == 1)
                    result.packet(0) = result.coeff(0);
            }
        }
        return result;
    }

    template <size_t Stride, typename Index, typename Mask>
    void scatter_(void *mem, const Index &index, const Mask &mask) const {
        size_t i1 = 0, i1i = this->size() == 1 ? 0 : 1,
               i2 = 0, i2i = index.size() == 1 ? 0 : 1,
               i3 = 0, i3i = mask.size() == 1 ? 0 : 1,
               size = check_size(*this, index, mask),
               n_packets = (size + PacketSize - 1) / PacketSize,
               i = 0;

        if (n_packets > 0) {
            for (; i < n_packets - (PacketSize > 1 ? 1 : 0); ++i, i1 += i1i, i2 += i2i, i3 += i3i)
                scatter<Stride>(mem, packet(i1), index.packet(i2), mask.packet(i3));
            if constexpr (PacketSize > 1) {
                auto mask2 = arange<IndexPacket>() <= IndexScalar((size - 1) % PacketSize);
                scatter<Stride>(mem, packet(i1), index.packet(i2), mask.packet(i3) & mask2);
            }
        }
    }

    template <size_t Stride, typename Index, typename Mask>
    void scatter_add_(void *mem, const Index &index, const Mask &mask) const {
        size_t i1 = 0, i1i = this->size() == 1 ? 0 : 1,
               i2 = 0, i2i = index.size() == 1 ? 0 : 1,
               i3 = 0, i3i = mask.size() == 1 ? 0 : 1,
               size = check_size(*this, index, mask),
               n_packets = (size + PacketSize - 1) / PacketSize,
               i = 0;

        if (n_packets > 0) {
            for (; i < n_packets - (PacketSize > 1 ? 1 : 0); ++i, i1 += i1i, i2 += i2i, i3 += i3i)
                scatter_add<Stride>(mem, packet(i1), index.packet(i2), mask.packet(i3));
            if constexpr (PacketSize > 1) {
                auto mask2 = arange<IndexPacket>() <= IndexScalar((size - 1) % PacketSize);
                scatter_add<Stride>(mem, packet(i1), index.packet(i2), mask.packet(i3) & mask2);
            }
        }
    }

    template <size_t Stride, typename Index, typename Func, typename... Args, typename Mask>
    static ENOKI_INLINE void transform_(void *ptr, const Index &index, const Mask &mask,
                                        const Func &func, const Args &... args) {
        size_t size = check_size(index, mask, args...),
               n_packets = (size + PacketSize - 1) / PacketSize;

        if (n_packets > 0) {
            size_t i = 0;
            for (; i < n_packets - (PacketSize > 1 ? 1 : 0); ++i)
                transform<Packet, Stride>(
                    ptr,
                    enoki::packet(index, enoki::slices(index) <= 1 ? 0 : i),
                    func,
                    enoki::packet(args, enoki::slices(args) <= 1 ? 0 : i)...);

            if constexpr (PacketSize > 1) {
                auto mask2 = arange<IndexPacket>() <= IndexScalar((size - 1) % PacketSize);
                transform<Packet, Stride>(
                    ptr,
                    enoki::packet(index, enoki::slices(index) <= 1 ? 0 : i),
                    func,
                    enoki::packet(args, enoki::slices(args) <= 1 ? 0 : i) & mask2...);
            }
        }
    }

    template <typename Mask> Derived compress_(const Mask &mask) const {
        assert(mask.size() == size());
        size_t count = 0;
        Derived result;
        set_slices(result, size());
        Value *ptr = result.data();

        for (size_t i = 0; i < packets(); ++i)
            count += compress(ptr, packet(i), mask.packet(i));
        set_slices(result, count);
        return result;
    }

    template <typename T> T ceil2int_() const {
        T result;
        result.resize(size());
        auto p1 = packet_ptr();
        auto pr = result.packet_ptr();
        for (size_t i = 0, n = result.packets();
             i < n; ++i, ++p1, ++pr)
            *pr = ceil2int<typename T::Packet>(*p1);
        return result;
    }

    template <typename T> T floor2int_() const {
        T result;
        result.resize(size());
        auto p1 = packet_ptr();
        auto pr = result.packet_ptr();
        for (size_t i = 0, n = result.packets();
             i < n; ++i, ++p1, ++pr)
            *pr = floor2int<typename T::Packet>(*p1);
        return result;
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Horizontal array operations
    // -----------------------------------------------------------------------

    Derived reverse_() const {
        using CoeffValue = std::conditional_t<IsMask, bool, Value>;

        size_t n = size();

        Derived result;
        set_slices(result, n);

        for (size_t i = 0; i < n; ++i)
            result.coeff(i) = (CoeffValue) coeff(n - 1 - i);

        return result;
    }

    Derived psum_() const {
        Derived result;
        set_slices(result, size());

        if (!empty()) {
            // Difficult to vectorize this..
            result.coeff(0) = coeff(0);
            for (size_t i = 1; i < size(); ++i)
                result.coeff(i) = result.coeff(i - 1) + coeff(i);
        }

        return result;
    }

    Value hsum_() const {
        if (size() == 0) {
            return Value(Scalar(0));
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result = zero<Packet>();
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result += packet(i);

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] +=
                    packet(packets() - 1);
            }
            return hsum(result);
        }
    }

    Value hprod_() const {
        if (size() == 0) {
            return Value(Scalar(1));
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result = Scalar(1);
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result *= packet(i);

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] *=
                    packet(packets() - 1);
            }
            return hprod(result);
        }
    }

    Value hmin_() const {
        if (size() == 0) {
            return Value(std::numeric_limits<Scalar>::max());
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result = coeff(0);
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result = min(result, packet(i));

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] =
                    min(result, packet(packets() - 1));
            }
            return hmin(result);
        }
    }

    Value hmax_() const {
        if (size() == 0) {
            return Value(std::numeric_limits<Scalar>::min());
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result = coeff(0);
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result = max(result, packet(i));

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] =
                    max(result, packet(packets() - 1));
            }
            return hmax(result);
        }
    }

    bool any_() const {
        if (size() == 0) {
            return false;
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result(false);
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result |= packet(i);

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] |=
                    packet(packets() - 1);
            }
            return any(result);
        }
    }

    bool all_() const {
        if (size() == 0) {
            return true;
        } else if (size() == 1) {
            return coeff(0);
        } else {
            Packet result(true);
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result &= packet(i);

            if constexpr (PacketSize > 1) {
                result[arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize)] &=
                    packet(packets() - 1);
            }
            return all(result);
        }
    }

    size_t count_() const {
        size_t result = 0;
        if (!empty()) {
            for (size_t i = 0, count = packets() - (PacketSize > 1 ? 1 : 0); i < count; ++i)
                result += enoki::count(packet(i));

            if constexpr (PacketSize > 1) {
                auto mask = arange<IndexPacket>() <= IndexScalar((size() - 1) % PacketSize);
                result += enoki::count(packet(packets() - 1) & mask);
            }
        }
        return result;
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Initialization helper functions
    // -----------------------------------------------------------------------

    /**
     * \brief Resize the buffer to the desired size
     *
     * When the capacity is insufficient, the implementation destroys the
     * current contents and allocates a new (uninitialized) buffer
     *
     * When compiled in debug mode, newly allocated floating point arrays will
     * be initialized with NaNs.
     */
    ENOKI_NOINLINE void resize(size_t size) {
        if (size == (size_t) m_size)
            return;

        if (is_mapped())
            throw std::runtime_error("Can't resize a mapped dynamic array!");

        using CoeffValue = std::conditional_t<IsMask, bool, Value>;

        CoeffValue scalar = (m_size == 1) ? coeff(0) : zero<CoeffValue>();
        size_t n_packets = (size + PacketSize - 1) / PacketSize;

        if (n_packets > packets_allocated()) {
            if (!empty()) {
                ENOKI_TRACK_DEALLOC(m_packets.get(), packets_allocated() * sizeof(Packet));
            }
            m_packets = PacketHolder(new Packet[n_packets]);
            m_packets_allocated = (Size) n_packets;
            ENOKI_TRACK_ALLOC(m_packets.get(),
                              n_packets * sizeof(Packet));
        }

        if (m_size == 1) {
            /* Resizing a scalar array -- broadcast. */
            Packet p(scalar);
            for (size_t i = 0; i < n_packets; ++i)
                m_packets[i] = p;
        } else if (m_size == 0) {
            /* Potentially initialize array contents with NaNs */
            #if !defined(NDEBUG)
                for (size_t i = 0; i < n_packets; ++i)
                    new (&m_packets[i]) Packet();
            #endif
        }

        m_size = (Size) size;
        clean_trailing_();
    }

    // Clear the unused portion of a potential trailing partial packet
    void clean_trailing_() {
        IndexScalar remainder = (IndexScalar) (m_size % PacketSize);
        if (remainder > 0 && m_size != 1) {
            void *addr = m_packets.get() + packets_allocated() - 1;
            auto mask = arange<IndexPacket>() < IndexScalar(remainder);
            store(addr, load<Packet>(addr) & mask);
        }
    }

    static Derived map(void *ptr, size_t size, bool dealloc = false) {
        assert((uintptr_t) ptr % alignof(Packet) == 0);

        Derived r;
        r.m_packets = PacketHolder((Packet *) ptr);
        r.m_size = (Size) size;
        r.m_packets_allocated =
            (Size) ((size + PacketSize - 1) / PacketSize);

        if (!dealloc)
            r.m_packets_allocated |= 0x80000000u;

        return r;
    }

    static Derived copy(const void *ptr, size_t size) {
        Derived r;
        r.m_size = (Size) size;
        r.m_packets_allocated =
            (Size) ((size + PacketSize - 1) / PacketSize);
        r.m_packets = PacketHolder(new Packet[r.m_packets_allocated]);
        memcpy(r.m_packets.get(), ptr, size * sizeof(Value));
        return r;
    }

    Derived &managed() { return derived(); }
    Derived &eval() { return derived(); }
    Derived &managed() const { return derived(); }
    Derived &eval() const { return derived(); }

    template <typename... Args> void resize_like(const Args&... args) {
        resize(check_size(args...));
    }

private:

#if defined(__GNUC__)
// GCC 8.2: quench nonsensical warning in parameter pack expansion
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wparentheses" //  warning: suggest parentheses around ‘&&’ within ‘||’ [-Wparentheses]
#endif

    template <typename... Args> static size_t check_size(const Args&... args) {
        size_t max_size = std::max({ slices(args)... });
        if ((... || (slices(args) != max_size && slices(args) != 1))) {
            #if defined(NDEBUG)
                throw std::runtime_error(
                    "Incompatible sizes in dynamic array operation");
            #else
                std::string msg = "[";
                bool result[] = { ((msg += (std::to_string(slices(args)) + ", ")), false)... };
                (void) result;
                if (msg.size() > 2)
                    msg = msg.substr(0, msg.size() - 2);
                msg += "]";
                throw std::runtime_error(
                    "Incompatible sizes in dynamic array operation: " + msg);
            #endif
        }
        return max_size;
    }

#if defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

public:
    static Derived empty_(size_t size) {
        Derived result;
        result.resize(size);
        return result;
    }

    static Derived zero_(size_t size) {
        Derived result;
        result.resize(size);
        Packet value_p = zero<Packet>();
        for (size_t i = 0; i < result.packets(); ++i)
            result.packet(i) = value_p;
        return result;
    }

    static Derived full_(const Value &value, size_t size) {
        Derived result;
        result.resize(size);
        Packet value_p(value);
        for (size_t i = 0; i < result.packets(); ++i)
            result.packet(i) = value_p;
        return result;
    }

    /// Construct an evenly spaced integer sequence
    static Derived arange_(ssize_t start, ssize_t stop, ssize_t step) {
        Derived result;
        result.resize(size_t((stop - start + step - (step > 0 ? 1 : -1)) / step));
        Packet value_p = arange<Packet>(start, start + (ssize_t) Packet::Size * step, step),
               shift   = Value((ssize_t) PacketSize * step);
        for (size_t i = 0; i < result.packets(); ++i) {
            result.packet(i) = value_p;
            value_p += shift;
        }
        return result;
    }

    static Derived linspace_(Value min, Value max, size_t size) {
        Derived result;
        result.resize(size);

        Value step = (max - min) / Value(size - 1);

        Packet value_p = linspace<Packet>(min, min + step * (PacketSize - 1)),
               shift   = Value(step * PacketSize);

        for (size_t i = 0; i < result.packets(); ++i) {
            result.packet(i) = value_p;
            value_p += shift;
        }

        return result;
    }

    //! @}
    // -----------------------------------------------------------------------

    auto operator->() const {
        using BaseType = std::decay_t<std::remove_pointer_t<scalar_t<Derived_>>>;
        return call_support<BaseType, Derived_>(derived());
    }

    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        check_size(derived(), mask);
        for (size_t i = 0; i < mask.size(); ++i)
            if (mask.coeff(i))
                return coeff(i);
        return zero<Value>();
    }

    DynamicArrayReference<Packet> ref_wrap_() const {
        return m_packets.get();
    }
private:
    PacketHolder m_packets;
    Size m_size = 0;
    Size m_packets_allocated = 0;
};

template <typename Packet_>
struct DynamicArray : DynamicArrayImpl<Packet_, DynamicArray<Packet_>> {
    using Base = DynamicArrayImpl<Packet_, DynamicArray<Packet_>>;
    using Base::Base;
    using Base::operator=;

    using ArrayType = DynamicArray;
    using MaskType  = DynamicMask<mask_t<Packet_>>;

    template <typename T> using ReplaceValue =
        DynamicArray<typename Packet_::template ReplaceValue<T>>;

    DynamicArray(const DynamicArray &) = default;
    DynamicArray(DynamicArray &&) = default;
    DynamicArray &operator=(const DynamicArray &) = default;
    DynamicArray &operator=(DynamicArray &&) = default;
};

template <typename Packet_>
struct DynamicMask : DynamicArrayImpl<Packet_, DynamicMask<Packet_>> {
    using Base = DynamicArrayImpl<Packet_, DynamicMask<Packet_>>;

    using ArrayType = DynamicArray<array_t<Packet_>>;
    using MaskType  = DynamicMask;

    template <typename T> using ReplaceValue =
        DynamicMask<typename Packet_::template ReplaceValue<T>>;

    DynamicMask() = default;

    template <typename T> DynamicMask(T &&value)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }

    template <typename T> DynamicMask(T &&value, detail::reinterpret_flag)
        : Base(std::forward<T>(value), detail::reinterpret_flag()) { }
};

namespace detail {
    template <typename T> struct mutable_ref { using type = std::add_lvalue_reference_t<T>; };
    template <typename T> struct mutable_ref<const T &> { using type = T &; };
    template <typename T> using mutable_ref_t = typename mutable_ref<T>::type;

    /// Vectorized inner loop (void return value)
    template <typename Func, typename... Args, size_t... Index>
    ENOKI_INLINE void vectorize_inner_1(std::index_sequence<Index...>, Func &&f,
                                        size_t packet_count, Args &&... args) {
        ENOKI_NOUNROLL ENOKI_IVDEP for (size_t i = 0; i < packet_count; ++i)
            f(packet(args, i)...);
    }

    /// Vectorized inner loop (non-void return value)
    template <typename Func, typename Out, typename... Args, size_t... Index>
    ENOKI_INLINE void vectorize_inner_2(std::index_sequence<Index...>, Func &&f,
                                        size_t packet_count, Out &&out, Args &&... args) {
        ENOKI_NOUNROLL ENOKI_IVDEP for (size_t i = 0; i < packet_count; ++i)
            packet(out, i) = f(packet(args, i)...);
    }
}

template <bool Resize = false, typename Func, typename... Args>
auto vectorize(Func &&f, Args &&... args)
    -> make_dynamic_t<decltype(f(packet(args, 0)...))> /* LLVM bug #39326 */ {
#if defined(NDEBUG)
    constexpr bool Check = false;
#else
    constexpr bool Check = true;
#endif

    /** Determine the number of slices and packets of the input arrays,
        and broadcast scalar input arrays if requested */
    size_t packet_count = 0, slice_count = 0;

    bool unused1[] = { ((packet_count = !is_dynamic_v<Args> ? packet_count
        : (Resize ? std::max(packet_count, packets(args)) : packets(args))), false)... };

    bool unused2[] = { ((slice_count = !is_dynamic_v<Args> ? slice_count
        : (Resize ? std::max(slice_count, slices(args)) : slices(args))), false)... };

    (void) unused1; (void) unused2;

    if constexpr (Check || Resize) {
        size_t status[] = { (
            (!is_dynamic_v<Args> || array_size_v<Args> == 0) ||
            ((slice_count != 1 && slices(args) == 1 && Resize)
                 ? (set_slices((detail::mutable_ref_t<decltype(args)>) args, slice_count), true)
                 : (slices(args) == slice_count)))... };

        bool status_combined = true;
        for (bool s : status)
            status_combined &= s;

        if (!status_combined)
            throw std::runtime_error("vectorize(): vector arguments have incompatible lengths");
    }

    using Result = make_dynamic_t<decltype(f(packet(args, 0)...))>;
    if constexpr (std::is_void_v<Result>) {
        detail::vectorize_inner_1(std::make_index_sequence<sizeof...(Args)>(),
                                  f, packet_count, ref_wrap(args)...);
    } else {
        Result result;
        set_slices(result, slice_count);

        detail::vectorize_inner_2(std::make_index_sequence<sizeof...(Args)>(),
                                  f, packet_count, ref_wrap(result),
                                  ref_wrap(args)...);
        return result;
    }
}

template <typename Func, typename... Args>
auto vectorize_safe(Func &&f, Args &&... args)
    -> decltype(vectorize<true>(f, args...)) /* LLVM bug #39326 */ {
    return vectorize<true>(f, args...);
}

namespace detail {
    template <typename T>
    using reference_dynamic_t = std::conditional_t<
        is_dynamic_v<T>,
        std::add_lvalue_reference_t<T>,
        T
    >;

    /// Strip the class from a method type
    template <typename T> struct remove_class { };
    template <typename C, typename R, typename... A> struct remove_class<R (C::*)(A...)> { typedef R type(A...); };
    template <typename C, typename R, typename... A> struct remove_class<R (C::*)(A...) const> { typedef R type(A...); };
}

template <typename Func, typename Return, typename... Args>
auto vectorize_wrapper_detail(Func &&f_, Return (*)(Args...)) {
    return [f = std::forward<Func>(f_)](detail::reference_dynamic_t<enoki::make_dynamic_t<Args>>... args) {
        return vectorize_safe(f, args...);
    };
}

/// Vectorize a vanilla function pointer
template <typename Return, typename... Args>
auto vectorize_wrapper(Return (*f)(Args...)) {
    return vectorize_wrapper_detail(f, f);
}

/// Vectorize a lambda function method (possibly with internal state)
template <typename Func,
          typename FuncType = typename detail::remove_class<
              decltype(&std::remove_reference<Func>::type::operator())>::type>
auto vectorize_wrapper(Func &&f) {
    return vectorize_wrapper_detail(std::forward<Func>(f), (FuncType *) nullptr);
}

/// Vectorize a class method (non-const)
template <typename Return, typename Class, typename... Arg>
auto vectorize_wrapper(Return (Class::*f)(Arg...)) {
    return vectorize_wrapper_detail(
        [f](Class *c, Arg... args) -> Return { return (c->*f)(args...); },
        (Return(*)(Class *, Arg...)) nullptr);
}

/// Vectorize a class method (const)
template <typename Return, typename Class, typename... Arg>
auto vectorize_wrapper(Return (Class::*f)(Arg...) const) {
    return vectorize_wrapper_detail(
        [f](const Class *c, Arg... args) -> Return { return (c->*f)(args...); },
        (Return(*)(const Class *, Arg...)) nullptr);
}

#if defined(ENOKI_AUTODIFF_H) && !defined(ENOKI_BUILD)
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<DynamicArray<Packet<float>>>;
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<DynamicArray<Packet<float>>>;

    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT Tape<DynamicArray<Packet<double>>>;
    ENOKI_AUTODIFF_EXTERN template struct ENOKI_AUTODIFF_EXPORT DiffArray<DynamicArray<Packet<double>>>;
#endif

NAMESPACE_END(enoki)

#if defined(__GNUC__) && !defined(__clang__)
#  pragma GCC diagnostic pop
#endif
