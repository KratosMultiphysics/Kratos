#pragma once

#include "array_base.h"

NAMESPACE_BEGIN(enoki)

namespace detail {
    /// Compute binary OR of 'i' with right-shifted versions
    static constexpr size_t fill(size_t i) {
        return i != 0 ? i | fill(i >> 1) : 0;
    }

    /// Find the largest power of two smaller than 'i'
    static constexpr size_t lpow2(size_t i) {
        return i != 0 ? (fill(i-1) >> 1) + 1 : 0;
    }

    /// Compile-time integer logarithm
    static constexpr size_t clog2i(size_t value) {
        return (value > 1) ? 1 + clog2i(value >> 1) : 0;
    }
}

template <typename Value_, size_t Size_, bool IsMask_, typename Derived_>
struct StaticArrayBase : ArrayBase<Value_, Derived_> {
    using Base = ArrayBase<Value_, Derived_>;
    using typename Base::Derived;
    using typename Base::Value;
    using typename Base::Scalar;
    using Base::derived;

    // -----------------------------------------------------------------------
    //! @{ \name Basic declarations
    // -----------------------------------------------------------------------

    /// Number of array entries
    static constexpr size_t Size = Size_;

    /// Size of the low array part returned by low()
    static constexpr size_t Size1 = detail::lpow2(Size_);

    /// Size of the high array part returned by high()
    static constexpr size_t Size2 = Size_ - Size1;

    /// Size and ActualSize can be different, e.g. when representing 3D vectors using 4-wide registers
    static constexpr size_t ActualSize = Size;

    /// Is this a mask type?
    static constexpr bool IsMask = Base::IsMask || IsMask_;

    /// Does this array represent a fixed size vector?
    static constexpr bool IsVector = true;

    /// Type of the low array part returned by low()
    using Array1 = std::conditional_t<!IsMask_, Array<Value_, Size1>,
                                                Mask <Value_, Size1>>;

    /// Type of the high array part returned by high()
    using Array2 = std::conditional_t<!IsMask_, Array<Value_, Size2>,
                                                Mask <Value_, Size2>>;

    //! @}
    // -----------------------------------------------------------------------

    constexpr size_t size() const { return Derived::Size; }

    void resize(size_t size) {
        if (size != Derived::Size)
            throw std::length_error("Incompatible size for static array");
    }

    // -----------------------------------------------------------------------
    //! @{ \name Constructors
    // -----------------------------------------------------------------------

    StaticArrayBase() = default;
    StaticArrayBase(const StaticArrayBase &) = default;
    StaticArrayBase(StaticArrayBase &&) = default;
    StaticArrayBase &operator=(const StaticArrayBase &) = default;
    StaticArrayBase &operator=(StaticArrayBase &&) = default;

    /// Type cast fallback
    template <typename Value2, size_t Size2,
              typename Derived2, typename T = Derived,
              enable_if_t<Derived2::Size == T::Size> = 0>
    ENOKI_INLINE StaticArrayBase(
        const StaticArrayBase<Value2, Size2, IsMask_, Derived2> &a) {
        ENOKI_CHKSCALAR("Copy constructor (type cast)");
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) derived().coeff(i) = (const Value &) a.derived().coeff(i);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Fallback implementations of vertical operations
    // -----------------------------------------------------------------------

    /// Addition
    Derived add_(const Derived &a) const {
        ENOKI_CHKSCALAR("add");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) +
                                        (const Value &) a.coeff(i);
        return result;
    }

    /// Subtraction
    Derived sub_(const Derived &a) const {
        ENOKI_CHKSCALAR("sub");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) -
                                        (const Value &) a.coeff(i);
        return result;
    }

    /// Multiplication (low part)
    Derived mul_(const Derived &a) const {
        ENOKI_CHKSCALAR("mul");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) *
                                        (const Value &) a.coeff(i);
        return result;
    }

    /// Multiplication (high part)
    Derived mulhi_(const Derived &a) const {
        ENOKI_CHKSCALAR("mulhi");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = mulhi((const Value &) derived().coeff(i),
                                              (const Value &) a.coeff(i));
        return result;
    }

    /// Division
    Derived div_(const Derived &a) const {
        ENOKI_CHKSCALAR("div");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) /
                                        (const Value &) a.coeff(i);
        return result;
    }

    /// Modulo
    Derived mod_(const Derived &a) const {
        ENOKI_CHKSCALAR("mod");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) %
                                        (const Value &) a.coeff(i);
        return result;
    }

    /// Arithmetic NOT operation fallback
    ENOKI_INLINE Derived not_() const {
        if constexpr (!is_mask_v<Derived>) {
            const Scalar mask = memcpy_cast<Scalar>(int_array_t<Scalar>(-1));
            return derived() ^ mask;
        } else {
            return derived() ^ Derived(true);
        }
    }

    /// Arithmetic unary negation operation fallback
    ENOKI_INLINE Derived neg_() const {
        if constexpr (std::is_floating_point_v<Scalar>)
            return derived() ^ Scalar(-0.f);
        else
            return ~derived() + Scalar(1);
    }

    /// Arithmetic OR operation
    template <typename Array>
    ENOKI_INLINE Derived or_(const Array &d) const {
        Derived result;
        ENOKI_CHKSCALAR("or");
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                detail::or_((const Value &) derived().coeff(i), d.coeff(i));
        return result;
    }

    /// Arithmetic AND operation
    template <typename Array>
    ENOKI_INLINE Derived and_(const Array &d) const {
        ENOKI_CHKSCALAR("and");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                detail::and_((const Value &) derived().coeff(i), d.coeff(i));
        return result;
    }

    /// Arithmetic ANDNOT operation
    template <typename Array>
    ENOKI_INLINE Derived andnot_(const Array &d) const {
        ENOKI_CHKSCALAR("andnot");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                detail::andnot_((const Value &) derived().coeff(i), d.coeff(i));
        return result;
    }

    /// Arithmetic XOR operation
    template <typename Array>
    ENOKI_INLINE Derived xor_(const Array &d) const {
        ENOKI_CHKSCALAR("xor");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                detail::xor_((const Value &) derived().coeff(i), d.coeff(i));
        return result;
    }

    /// Left shift operator (uniform)
    ENOKI_INLINE Derived sl_(size_t value) const {
        ENOKI_CHKSCALAR("sl");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                (const Value &) derived().coeff(i) << value;
        return result;
    }

    /// Left shift operator (array)
    ENOKI_INLINE Derived sl_(const Derived &d) const {
        ENOKI_CHKSCALAR("sl");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) <<
                                        (const Value &) d.coeff(i);
        return result;
    }

    /// Left shift operator (immediate)
    template <size_t Imm> ENOKI_INLINE Derived sl_() const {
        ENOKI_CHKSCALAR("sl");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                sl<Imm>((const Value &) derived().coeff(i));
        return result;
    }

    /// Right shift operator (Uniform)
    ENOKI_INLINE Derived sr_(size_t value) const {
        ENOKI_CHKSCALAR("sr");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                (const Value &) derived().coeff(i) >> value;
        return result;
    }

    /// Right shift operator (Array)
    ENOKI_INLINE Derived sr_(const Derived &d) const {
        ENOKI_CHKSCALAR("sr");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = (const Value &) derived().coeff(i) >>
                                        (const Value &) d.coeff(i);
        return result;
    }

    /// Right shift operator (immediate)
    template <size_t Imm> ENOKI_INLINE Derived sr_() const {
        ENOKI_CHKSCALAR("sr");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                sr<Imm>((const Value &) derived().coeff(i));
        return result;
    }

    /// Equality comparison operation
    ENOKI_INLINE auto eq_(const Derived &d) const {
        ENOKI_CHKSCALAR("eq");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = eq((const Value &) derived().coeff(i),
                                 (const Value &) d.coeff(i));
        return result;
    }

    /// Inequality comparison operation
    ENOKI_INLINE auto neq_(const Derived &d) const {
        ENOKI_CHKSCALAR("neq");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = neq((const Value &) derived().coeff(i),
                                  (const Value &) d.coeff(i));
        return result;
    }

    /// Less than comparison operation
    ENOKI_INLINE auto lt_(const Derived &d) const {
        ENOKI_CHKSCALAR("lt");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) derived().coeff(i) <
                              (const Value &) d.coeff(i);
        return result;
    }

    /// Less than or equal comparison operation
    ENOKI_INLINE auto le_(const Derived &d) const {
        ENOKI_CHKSCALAR("le");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) derived().coeff(i) <=
                              (const Value &) d.coeff(i);
        return result;
    }

    /// Greater than comparison operation
    ENOKI_INLINE auto gt_(const Derived &d) const {
        ENOKI_CHKSCALAR("gt");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) derived().coeff(i) >
                              (const Value &) d.coeff(i);
        return result;
    }

    /// Greater than or equal comparison operation
    ENOKI_INLINE auto ge_(const Derived &d) const {
        ENOKI_CHKSCALAR("ge");
        mask_t<Derived> result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) derived().coeff(i) >=
                              (const Value &) d.coeff(i);
        return result;
    }

    /// Absolute value
    ENOKI_INLINE Derived abs_() const {
        ENOKI_CHKSCALAR("abs");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::abs((const Value &) derived().coeff(i));
        return result;
    }

    /// Square root
    ENOKI_INLINE Derived sqrt_() const {
        ENOKI_CHKSCALAR("sqrt");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::sqrt((const Value &) derived().coeff(i));
        return result;
    }

    /// Reciprocal fallback implementation
    ENOKI_INLINE Derived rcp_() const {
        return (Scalar) 1 / derived();
    }

    /// Reciprocal square root fallback implementation
    ENOKI_INLINE Derived rsqrt_() const {
        return (Scalar) 1 / sqrt(derived());
    }

    /// Round to smallest integral value not less than argument
    ENOKI_INLINE Derived ceil_() const {
        ENOKI_CHKSCALAR("ceil");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::ceil((const Value &) derived().coeff(i));
        return result;
    }

    /// Round to largest integral value not greater than argument
    ENOKI_INLINE Derived floor_() const {
        ENOKI_CHKSCALAR("floor");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::floor((const Value &) derived().coeff(i));
        return result;
    }

    /// Round to integral value
    ENOKI_INLINE Derived round_() const {
        ENOKI_CHKSCALAR("round");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::round((const Value &) derived().coeff(i));
        return result;
    }

    /// Round to zero
    ENOKI_INLINE Derived trunc_() const {
        ENOKI_CHKSCALAR("trunc");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                enoki::trunc((const Value &) derived().coeff(i));
        return result;
    }

    /// Element-wise maximum
    ENOKI_INLINE Derived max_(const Derived &d) const {
        ENOKI_CHKSCALAR("max");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = enoki::max((const Value &) derived().coeff(i),
                                                   (const Value &) d.coeff(i));
        return result;
    }

    /// Element-wise minimum
    ENOKI_INLINE Derived min_(const Derived &d) const {
        ENOKI_CHKSCALAR("min");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = enoki::min((const Value &) derived().coeff(i),
                                                   (const Value &) d.coeff(i));
        return result;
    }

    /// Fused multiply-add
    ENOKI_INLINE Derived fmadd_(const Derived &d1, const Derived &d2) const {
        if constexpr (array_depth_v<Value> > 0) {
            Derived result;
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) = fmadd((const Value &) derived().coeff(i),
                                                  (const Value &) d1.coeff(i),
                                                  (const Value &) d2.coeff(i));
            return result;
        } else {
            return derived() * d1 + d2;
        }
    }

    /// Fused negative multiply-add
    ENOKI_INLINE Derived fnmadd_(const Derived &d1, const Derived &d2) const {
        if constexpr (array_depth_v<Value> > 0) {
            Derived result;
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) = fnmadd((const Value &) derived().coeff(i),
                                                   (const Value &) d1.coeff(i),
                                                   (const Value &) d2.coeff(i));
            return result;
        } else {
            return -derived() * d1 + d2;
        }
    }

    /// Fused multiply-subtract
    ENOKI_INLINE Derived fmsub_(const Derived &d1, const Derived &d2) const {
        if constexpr (array_depth_v<Value> > 0) {
            Derived result;
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) = fmsub((const Value &) derived().coeff(i),
                                                  (const Value &) d1.coeff(i),
                                                  (const Value &) d2.coeff(i));
            return result;
        } else {
            return derived() * d1 - d2;
        }
    }

    /// Fused negative multiply-subtract
    ENOKI_INLINE Derived fnmsub_(const Derived &d1, const Derived &d2) const {
        if constexpr (array_depth_v<Value> > 0) {
            ENOKI_CHKSCALAR("fnmsub");
            Derived result;
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) = fnmsub((const Value &) derived().coeff(i),
                                                   (const Value &) d1.coeff(i),
                                                   (const Value &) d2.coeff(i));
            return result;
        } else {
            return -derived() * d1 - d2;
        }
    }

    /// Fused multiply-add/subtract fallback implementation
    ENOKI_INLINE Derived fmaddsub_(const Derived &b, const Derived &c) const {
        ENOKI_CHKSCALAR("fmaddsub");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i) {
            if (i % 2 == 0)
                (Value &) result.coeff(i) = fmsub((const Value &) derived().coeff(i),
                                                  (const Value &) b.coeff(i),
                                                  (const Value &) c.coeff(i));
            else
                (Value &) result.coeff(i) = fmadd((const Value &) derived().coeff(i),
                                                  (const Value &) b.coeff(i),
                                                  (const Value &) c.coeff(i));
        }
        return result;
    }

    /// Fused multiply-subtract/add fallback implementation
    ENOKI_INLINE Derived fmsubadd_(const Derived &b, const Derived &c) const {
        ENOKI_CHKSCALAR("fmsubadd");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i) {
            if (i % 2 == 0)
                (Value &) result.coeff(i) = fmadd((const Value &) derived().coeff(i),
                                                  (const Value &) b.coeff(i),
                                                  (const Value &) c.coeff(i));
            else
                (Value &) result.coeff(i) = fmsub((const Value &) derived().coeff(i),
                                                  (const Value &) b.coeff(i),
                                                  (const Value &) c.coeff(i));
        }
        return result;
    }

    /// Masked prefetch fallback
    template <bool Write, size_t Level, size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE void prefetch_(const void *mem, const Index &index, const Mask &mask) {
        ENOKI_CHKSCALAR("prefetch");
        for (size_t i = 0; i < Derived::Size; ++i)
            prefetch<Value, Write, Level, Stride>(mem, index.coeff(i), mask.coeff(i));
    }

    /// Masked gather fallback
    template <size_t Stride, typename Index, typename Mask>
    static ENOKI_INLINE Derived gather_(const void *mem, const Index &index, const Mask &mask) {
        ENOKI_CHKSCALAR("gather");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i) {
            if constexpr (!is_mask_v<Derived>) {
                (Value &) result.coeff(i) =
                    (const Value &) gather<Value, Stride>(mem, index.coeff(i), mask.coeff(i));
            } else {
                result.coeff(i) = gather<Value, Stride>(mem, index.coeff(i), mask.coeff(i));
            }
        }
        return result;
    }

    /// Masked scatter fallback
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_(void *mem, const Index &index, const Mask &mask) const {
        ENOKI_CHKSCALAR("scatter");
        for (size_t i = 0; i < Derived::Size; ++i)
            scatter<Stride>(mem, (const Value &) derived().coeff(i), index.coeff(i), mask.coeff(i));
    }

    /// Masked scatter_add-add fallback
    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_add_(void *mem, const Index &index, const Mask &mask) const {
        transform<Derived, Stride>(mem, index,
            [](auto &&a, auto &&b, auto &&) { a += b; },
            derived(),
            mask
        );
    }

    /// Ternary operator -- select between to values based on mask
    template <typename Mask>
    static ENOKI_INLINE auto select_(const Mask &m, const Derived &t, const Derived &f) {
        ENOKI_CHKSCALAR("select");
        Derived result;
        for (size_t i = 0; i < Size; ++i)
            (Value &) result.coeff(i) = select(m.coeff(i), (const Value &) t.coeff(i),
                                                           (const Value &) f.coeff(i));
        return result;
    }

    /// Shuffle operation fallback implementation
    template <size_t... Indices> ENOKI_INLINE Derived shuffle_() const {
        static_assert(sizeof...(Indices) == Size ||
                      sizeof...(Indices) == Derived::Size, "shuffle(): Invalid size!");
        ENOKI_CHKSCALAR("shuffle");
        Derived out;
        size_t idx = 0;
        bool result[] = { (out.coeff(idx++) = derived().coeff(Indices % Derived::Size), false)... };
        (void) idx; (void) result;
        return out;
    }

    template <typename Index> ENOKI_INLINE Derived shuffle_(const Index &index) const {
        ENOKI_CHKSCALAR("shuffle");
        Derived out;
        for (size_t i = 0; i < Derived::Size; ++i) {
            size_t idx = (size_t) index.coeff(i);
            out.coeff(i) = derived().coeff(idx % Derived::Size);
        }
        return out;
    }

    /// Rotate the entries of the array right
    template <size_t Imm> ENOKI_INLINE Derived ror_array_() const {
        return ror_array_<Imm>(std::make_index_sequence<Derived::Size>());
    }

    /// Rotate the entries of the array left
    template <size_t Imm>
    ENOKI_INLINE Derived rol_array_() const {
        return rol_array_<Imm>(std::make_index_sequence<Derived::Size>());
    }

    template <typename T> T floor2int_() const {
        if constexpr (array_depth_v<Value> > 0) {
            T result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) =
                    floor2int<value_t<T>>((const Value &) derived().coeff(i));
            return result;
        } else {
            return T(floor(derived()));
        }
    }

    template <typename T> T ceil2int_() const {
        if constexpr (array_depth_v<Value> > 0) {
            T result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) =
                    ceil2int<value_t<T>>((const Value &) derived().coeff(i));
            return result;
        } else {
            return T(enoki::ceil(derived()));
        }
    }

private:
    template <size_t Imm, size_t... Is>
    ENOKI_INLINE Derived ror_array_(std::index_sequence<Is...>) const {
        return shuffle<(Is + Derived::Size - Imm) % Derived::Size...>(derived());
    }

    template <size_t Imm, size_t... Is>
    ENOKI_INLINE Derived rol_array_(std::index_sequence<Is...>) const {
        return shuffle<(Is + Imm) % Derived::Size...>(derived());
    }

    template <typename T, size_t Offset, size_t... Is>
    ENOKI_INLINE T sub_array_(std::index_sequence<Is...>) const {
        return T((typename Derived::Value) derived().coeff(Offset + Is)...);
    }

public:
    /// Return the low array part (always a power of two)
    ENOKI_INLINE auto low_() const {
        return sub_array_<typename Derived::Array1, 0>(
            std::make_index_sequence<Derived::Size1>());
    }

    /// Return the high array part
    template <typename T = Derived, enable_if_t<T::Size2 != 0> = 0>
    ENOKI_INLINE auto high_() const {
        return sub_array_<typename Derived::Array2, Derived::Size1>(
            std::make_index_sequence<Derived::Size2>());
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Bit counting fallback implementations
    // -----------------------------------------------------------------------

    Derived popcnt_() const {
        using UInt = uint_array_t<Derived>;
        UInt w = reinterpret_array<UInt>(derived());
        using U = scalar_t<UInt>;

        if constexpr (sizeof(Scalar) <= 4) {
            w -= sr<1>(w) & U(0x55555555u);
            w = (w & U(0x33333333u)) + ((sr<2>(w)) & U(0x33333333u));
            w = (w + sr<4>(w)) & U(0x0F0F0F0Fu);
            w = sr<24>(w * U(0x01010101u));
        } else {
            w -= sr<1>(w) & U(0x5555555555555555ull);
            w = (w & U(0x3333333333333333ull)) + (sr<2>(w) & U(0x3333333333333333ull));
            w = (w + sr<4>(w)) & U(0x0F0F0F0F0F0F0F0Full);
            w = sr<56>(w * U(0x0101010101010101ull));
        }
        return Derived(w);
    }

    Derived lzcnt_() const {
        using UInt = uint_array_t<Derived>;
        UInt w = reinterpret_array<UInt>(derived());
        w |= sr<1>(w);
        w |= sr<2>(w);
        w |= sr<4>(w);
        w |= sr<8>(w);
        w |= sr<16>(w);
        if constexpr (sizeof(Scalar) > 4)
            w |= sr<32>(w);
        return popcnt(~w);
    }

    Derived tzcnt_() const {
        using UInt = uint_array_t<Derived>;
        UInt w = reinterpret_array<UInt>(derived());
        w |= sl<1>(w);
        w |= sl<2>(w);
        w |= sl<4>(w);
        w |= sl<8>(w);
        w |= sl<16>(w);
        if constexpr (sizeof(Scalar) > 4)
            w |= sl<32>(w);
        return popcnt(~w);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Fallback implementations of horizontal operations
    // -----------------------------------------------------------------------

    /// Reverse fallback
    ENOKI_INLINE Derived reverse_() const {
        ENOKI_CHKSCALAR("reverse");
        Derived result;
        for (size_t i = 0; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) derived().coeff(Derived::Size - 1 - i);
        return result;
    }

    /// Prefix sum fallback
    ENOKI_INLINE Derived psum_() const {
        ENOKI_CHKSCALAR("psum");
        Derived result;
        result.coeff(0) = (const Value &) derived().coeff(0);
        for (size_t i = 1; i < Derived::Size; ++i)
            result.coeff(i) = (const Value &) result.coeff(i - 1) +
                              (const Value &) derived().coeff(i);
        return result;
    }

    /// Prefix sum over innermost dimension
    ENOKI_INLINE auto psum_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(psum_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = psum_inner(derived().coeff(i));
            return result;
        } else {
            return psum(derived());
        }
    }

    /// Horizontal sum fallback
    ENOKI_INLINE Value hsum_() const {
        ENOKI_CHKSCALAR("hsum");
        Value result = (const Value &) derived().coeff(0);
        for (size_t i = 1; i < Derived::Size; ++i)
            result += (const Value &) derived().coeff(i);
        return result;
    }

    /// Horizontal sum over innermost dimension
    ENOKI_INLINE auto hsum_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(hsum_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = hsum_inner(derived().coeff(i));
            return result;
        } else {
            return hsum(derived());
        }
    }

    /// Horizontal product fallback
    ENOKI_INLINE Value hprod_() const {
        ENOKI_CHKSCALAR("hprod");
        Value result = (const Value &) derived().coeff(0);
        for (size_t i = 1; i < Derived::Size; ++i)
            result *= (const Value &) derived().coeff(i);
        return result;
    }

    /// Horizontal product over innermost dimension
    ENOKI_INLINE auto hprod_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(hprod_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = hprod_inner(derived().coeff(i));
            return result;
        } else {
            return hprod(derived());
        }
    }

    /// Horizontal maximum fallback
    ENOKI_INLINE Value hmax_() const {
        Value result = (const Value &) derived().coeff(0);
        ENOKI_CHKSCALAR("hmax");
        for (size_t i = 1; i < Derived::Size; ++i)
            result = max(result, (const Value &) derived().coeff(i));
        return result;
    }

    /// Horizontal maximum over innermost dimension
    ENOKI_INLINE auto hmax_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(hmax_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = hmax_inner(derived().coeff(i));
            return result;
        } else {
            return hmax(derived());
        }
    }

    /// Horizontal minimum fallback
    ENOKI_INLINE Value hmin_() const {
        Value result = (const Value &) derived().coeff(0);
        ENOKI_CHKSCALAR("hmin");
        for (size_t i = 1; i < Derived::Size; ++i)
            result = min(result, (const Value &) derived().coeff(i));
        return result;
    }

    /// Horizontal minimum over innermost dimension
    ENOKI_INLINE auto hmin_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(hmin_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = hmin_inner(derived().coeff(i));
            return result;
        } else {
            return hmin(derived());
        }
    }

    /// Horizontal mean over innermost dimension
    ENOKI_INLINE auto hmean_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(hmean_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = hmean_inner(derived().coeff(i));
            return result;
        } else {
            return hmean(derived());
        }
    }

    /// all() fallback implementation
    ENOKI_INLINE auto all_() const {
        ENOKI_CHKSCALAR("all");
        if constexpr (Derived::IsMask && std::is_scalar_v<Value_>) {
            bool result = derived().coeff(0);
            for (size_t i = 1; i < Derived::Size; ++i)
                result = result && derived().coeff(i);
            return result;
        } else {
            auto result = derived().coeff(0);
            for (size_t i = 1; i < Derived::Size; ++i)
                result &= derived().coeff(i);
            return result;
        }
    }

    /// all() over innermost dimension
    ENOKI_INLINE auto all_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(all_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = all_inner(derived().coeff(i));
            return result;
        } else {
            return all(derived());
        }
    }

    /// any() fallback implementation
    ENOKI_INLINE auto any_() const {
        ENOKI_CHKSCALAR("any");
        if constexpr (Derived::IsMask && std::is_scalar_v<Value_>) {
            bool result = derived().coeff(0);
            for (size_t i = 1; i < Derived::Size; ++i)
                result = result || derived().coeff(i);
            return result;
        } else {
            auto result = derived().coeff(0);
            for (size_t i = 1; i < Derived::Size; ++i)
                result |= derived().coeff(i);
            return result;
        }
    }

    /// any() over innermost dimension
    ENOKI_INLINE auto any_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(any_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = any_inner(derived().coeff(i));
            return result;
        } else {
            return any(derived());
        }
    }

    /// count() fallback implementation
    ENOKI_INLINE auto count_() const {
        ENOKI_CHKSCALAR("count");
        using Int = value_t<size_array_t<array_t<Derived>>>;
        const Int one(1);
        Int result(0);
        for (size_t i = 0; i < Derived::Size; ++i)
            masked(result, derived().coeff(i)) += one;
        return result;
    }

    /// count() over innermost dimension
    ENOKI_INLINE auto count_inner_() const {
        if constexpr (is_array_v<Value>) {
            using Value = decltype(count_inner(derived().coeff(0)));
            using Result = typename Derived::template ReplaceValue<Value>;
            Result result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = count_inner(derived().coeff(i));
            return result;
        } else {
            return count(derived());
        }
    }

    /// Dot product fallback implementation
    ENOKI_INLINE Value dot_(const Derived &a) const {
        ENOKI_CHKSCALAR("dot");
        if constexpr (is_array_v<Value>) {
            Value result = (const Value &) derived().coeff(0) *
                           (const Value &) a.coeff(0);
            for (size_t i = 1; i < Size; ++i)
                result = fmadd((const Value &) derived().coeff(i),
                               (const Value &) a.coeff(i), result);
            return result;
        } else {
            return hsum(derived() * a);
        }
    }

    /// Extract fallback implementation
    template <typename Mask>
    ENOKI_INLINE Value extract_(const Mask &mask) const {
        ENOKI_CHKSCALAR("extract");
        for (size_t i = 0; i < Derived::Size; ++i)
            if (mask.coeff(i))
                return (const Value &) derived().coeff(i);
        return zero<Value>();
    }

    template <typename Mask>
    ENOKI_INLINE size_t compress_(Scalar *&mem, const Mask &mask) const {
        ENOKI_CHKSCALAR("compress");
        size_t result = 0;
        for (size_t i = 0; i < Derived::Size; ++i)
            result += compress(mem, (const Value &) derived().coeff(i), mask.coeff(i));
        return result;
    }

    /// Combined gather-modify-scatter operation without conflicts (fallback implementation)
    template <size_t Stride, typename Index, typename Func, typename Mask,
              typename... Args>
    static ENOKI_INLINE void transform_(void *mem, const Index &index,
                                        const Mask &, const Func &func,
                                        const Args &... args) {
        ENOKI_CHKSCALAR("transform");
        for (size_t i = 0; i < Derived::Size; ++i)
            transform<Value, Stride>(
                mem, index.coeff(i), func, args.coeff(i)...);
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Miscellaneous
    // -----------------------------------------------------------------------

    /// Return the size in bytes
    ENOKI_INLINE size_t nbytes() const {
        if constexpr (is_dynamic_v<Derived>) {
            size_t result = 0;
            for (size_t i = 0; i < Derived::Size; ++i)
                result += derived().coeff(i).nbytes();
            return result;
        } else {
            return sizeof(Derived);
        }
    }

    static ENOKI_INLINE Derived load_(const void *mem) {
        Derived result;
        if constexpr (is_scalar_v<Value>) {
            memcpy(result.data(), mem, sizeof(const Value &) * Derived::Size);
        } else {
            ENOKI_CHKSCALAR("load");
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) = load<Value>(static_cast<const Value *>(mem) + i);
        }
        return result;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_(const void *mem, const Mask &mask) {
        Derived result;
        ENOKI_CHKSCALAR("load");
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) =
                load<Value>(static_cast<const Value *>(mem) + i, mask.coeff(i));
        return result;
    }

    static ENOKI_INLINE Derived load_unaligned_(const void *mem) {
        Derived result;
        if constexpr (is_scalar_v<Value>) {
            memcpy(result.data(), mem, sizeof(const Value &) * Derived::Size);
        } else {
            ENOKI_CHKSCALAR("load_unaligned");
            for (size_t i = 0; i < Derived::Size; ++i)
                (Value &) result.coeff(i) =
                    load_unaligned<Value>(static_cast<const Value *>(mem) + i);
        }
        return result;
    }

    template <typename Mask>
    static ENOKI_INLINE Derived load_unaligned_(const void *mem, const Mask &mask) {
        Derived result;
        ENOKI_CHKSCALAR("load_unaligned");
        for (size_t i = 0; i < Derived::Size; ++i)
            (Value &) result.coeff(i) = load_unaligned<Value>(
                static_cast<const Value *>(mem) + i, mask.coeff(i));
        return result;
    }

    void store_(void *mem) const {
        if constexpr (is_scalar_v<Value>) {
            memcpy(mem, derived().data(), sizeof(const Value &) * derived().size());
        } else {
            ENOKI_CHKSCALAR("store");
            for (size_t i = 0; i < derived().size(); ++i)
                store<Value>(static_cast<Value *>(mem) + i, derived().coeff(i));
        }
    }

    template <typename Mask>
    void store_(void *mem, const Mask &mask) const {
        ENOKI_CHKSCALAR("store");
        for (size_t i = 0; i < derived().size(); ++i)
            store<Value>(static_cast<Value *>(mem) + i, derived().coeff(i),
                         mask.coeff(i));
    }

    void store_unaligned_(void *mem) const {
        if constexpr (is_scalar_v<Value>) {
            memcpy(mem, derived().data(), sizeof(const Value &) * derived().size());
        } else {
            ENOKI_CHKSCALAR("store_unaligned");
            for (size_t i = 0; i < derived().size(); ++i)
                store_unaligned<Value>(static_cast<Value *>(mem) + i,
                                       derived().coeff(i));
        }
    }

    template <typename Mask>
    void store_unaligned_(void *mem, const Mask &mask) const {
        ENOKI_CHKSCALAR("store_unaligned");
        for (size_t i = 0; i < derived().size(); ++i)
            store_unaligned<Value>(static_cast<Value *>(mem) + i,
                                   derived().coeff(i), mask.coeff(i));
    }

    static ENOKI_INLINE Derived zero_() { return Derived(zero<Value>()); }

    template <typename T> static Derived full_(const T &value, size_t size) {
        ENOKI_MARK_USED(size);

        if constexpr (array_depth_v<T> > array_depth_v<Value> ||
                      (array_depth_v<T> == array_depth_v<Value> &&
                       (is_dynamic_array_v<Value> || is_scalar_v<Value>))) {
            return Derived(value);
        } else {
            Derived result;
            for (size_t i = 0; i < Derived::Size; ++i)
                result.coeff(i) = Value::full_(value, size);
            return result;
        }
    }

    /// Construct an evenly spaced integer sequence
    static ENOKI_INLINE Derived arange_(ssize_t start, ssize_t stop, ssize_t step) {
        (void) stop;
        return linspace_(std::make_index_sequence<Derived::Size>(),
                         start, step);
    }

    /// Construct an array that linearly interpolates from min..max
    static ENOKI_INLINE Derived linspace_(Scalar min, Scalar max) {
        if constexpr (Derived::Size == 0) {
            return Derived();
        } else if constexpr (Derived::Size == 1) {
            return Derived(min);
        } else {
            return linspace_(std::make_index_sequence<Derived::Size>(), min,
                (max - min) / (Scalar) (Derived::Size - 1));
        }
    }

    /// Return an unitialized array
    static ENOKI_INLINE Derived empty_() { Derived result; return result; }

private:
    template <typename T, size_t... Is>
    static ENOKI_INLINE auto linspace_(std::index_sequence<Is...>, T offset, T step) {
        ENOKI_MARK_USED(step);
        if constexpr (sizeof...(Is) == 1)
            return Derived((Scalar) offset);
        else
            return Derived(((Scalar) ((T) Is * step + offset))...);
    }

public:

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Component access
    // -----------------------------------------------------------------------

    ENOKI_INLINE decltype(auto) x() const {
        static_assert(Derived::ActualSize >= 1, "StaticArrayBase::x(): requires Size >= 1");
        return derived().coeff(0);
    }

    ENOKI_INLINE decltype(auto) x() {
        static_assert(Derived::ActualSize >= 1, "StaticArrayBase::x(): requires Size >= 1");
        return derived().coeff(0);
    }

    ENOKI_INLINE decltype(auto) y() const {
        static_assert(Derived::ActualSize >= 2, "StaticArrayBase::y(): requires Size >= 2");
        return derived().coeff(1);
    }

    ENOKI_INLINE decltype(auto) y() {
        static_assert(Derived::ActualSize >= 2, "StaticArrayBase::y(): requires Size >= 2");
        return derived().coeff(1);
    }

    ENOKI_INLINE decltype(auto) z() const {
        static_assert(Derived::ActualSize >= 3, "StaticArrayBase::z(): requires Size >= 3");
        return derived().coeff(2);
    }

    ENOKI_INLINE decltype(auto) z() {
        static_assert(Derived::ActualSize >= 3, "StaticArrayBase::z(): requires Size >= 3");
        return derived().coeff(2);
    }

    ENOKI_INLINE decltype(auto) w() const {
        static_assert(Derived::ActualSize >= 4, "StaticArrayBase::w(): requires Size >= 4");
        return derived().coeff(3);
    }

    ENOKI_INLINE decltype(auto) w() {
        static_assert(Derived::ActualSize >= 4, "StaticArrayBase::w(): requires Size >= 4");
        return derived().coeff(3);
    }

    ENOKI_INLINE decltype(auto) data() { return &derived().coeff(0); }
    ENOKI_INLINE decltype(auto) data() const { return &derived().coeff(0); }

    ENOKI_INLINE Derived& managed() {
        if constexpr (is_cuda_array_v<Value_>) {
            for (size_t i = 0; i < Derived::Size; ++i)
                derived().coeff(i).managed();
        }
        return derived();
    }

    ENOKI_INLINE const Derived& managed() const {
        if constexpr (is_cuda_array_v<Value_>) {
            for (size_t i = 0; i < Derived::Size; ++i)
                derived().coeff(i).managed();
        }
        return derived();
    }

    ENOKI_INLINE Derived& eval() {
        if constexpr (is_cuda_array_v<Value_>) {
            for (size_t i = 0; i < Derived::Size; ++i)
                derived().coeff(i).eval();
        }
        return derived();
    }

    ENOKI_INLINE const Derived& eval() const {
        if constexpr (is_cuda_array_v<Value_>) {
            for (size_t i = 0; i < Derived::Size; ++i)
                derived().coeff(i).eval();
        }
        return derived();
    }

    //! @}
    // -----------------------------------------------------------------------
};

NAMESPACE_END(enoki)
