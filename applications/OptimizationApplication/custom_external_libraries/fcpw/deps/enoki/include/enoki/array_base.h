/*
    enoki/array_base.h -- Base class of all Enoki arrays

    Enoki is a C++ template library that enables transparent vectorization
    of numerical kernels using SIMD instruction sets available on current
    processor architectures.

    Copyright (c) 2019 Wenzel Jakob <wenzel.jakob@epfl.ch>

    All rights reserved. Use of this source code is governed by a BSD-style
    license that can be found in the LICENSE file.
*/

#include <enoki/array_router.h>
#include <enoki/array_masked.h>
#include <enoki/array_struct.h>

NAMESPACE_BEGIN(enoki)

template <typename Value_, typename Derived_> struct ArrayBase {
    // -----------------------------------------------------------------------
    //! @{ \name Curiously Recurring Template design pattern
    // -----------------------------------------------------------------------

    /// Alias to the derived type
    using Derived = Derived_;

    /// Cast to derived type
    ENOKI_INLINE Derived &derived()             { return (Derived &) *this; }

    /// Cast to derived type (const version)
    ENOKI_INLINE const Derived &derived() const { return (Derived &) *this; }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Basic declarations
    // -----------------------------------------------------------------------

    /// Actual type underlying the derived array
    using Value = Value_;

    /// Scalar data type all the way at the lowest level
    using Scalar = scalar_t<Value_>;

    /// Specifies how deeply nested this array is
    static constexpr size_t Depth = 1 + array_depth_v<Value>;

    /// Is this a mask type?
    static constexpr bool IsMask = is_mask_v<Value_>;

    /// Is this a dynamically allocated array (no by default)
    static constexpr bool IsDynamic = is_dynamic_v<Value_>;

    /// Does this array compute derivatives using automatic differentation?
    static constexpr bool IsDiff = is_diff_array_v<Value_>;

    /// Does this array reside on the GPU? (via CUDA)
    static constexpr bool IsCUDA = is_cuda_array_v<Value_>;

    /// Does this array map operations onto native vector instructions?
    static constexpr bool IsNative = false;

    /// Is this an AVX512-style 'k' mask register?
    static constexpr bool IsKMask = false;

    /// Is the storage representation of this array implemented recursively?
    static constexpr bool IsRecursive = false;

    /// Always prefer broadcasting to the outer dimensions of a N-D array
    static constexpr bool BroadcastPreferOuter = true;

    /// Does this array represent a fixed size vector?
    static constexpr bool IsVector = false;

    /// Does this array represent a complex number?
    static constexpr bool IsComplex = false;

    /// Does this array represent a quaternion?
    static constexpr bool IsQuaternion = false;

    /// Does this array represent a matrix?
    static constexpr bool IsMatrix = false;

    /// Does this array represent the result of a 'masked(...)' epxpression?
    static constexpr bool IsMaskedArray = false;

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Iterators
    // -----------------------------------------------------------------------

    ENOKI_INLINE auto begin() const { return derived().data(); }
    ENOKI_INLINE auto begin()       { return derived().data(); }
    ENOKI_INLINE auto end()   const { return derived().data() + derived().size(); }
    ENOKI_INLINE auto end()         { return derived().data() + derived().size(); }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Element access
    // -----------------------------------------------------------------------

    /// Array indexing operator with bounds checks in debug mode
    ENOKI_INLINE decltype(auto) operator[](size_t i) {
        #if !defined(NDEBUG) && !defined(ENOKI_DISABLE_RANGE_CHECK)
            if (i >= derived().size())
                throw std::out_of_range(
                    "ArrayBase: out of range access (tried to access index " +
                    std::to_string(i) + " in an array of size " +
                    std::to_string(derived().size()) + ")");
        #endif
        return derived().coeff(i);
    }

    /// Array indexing operator with bounds checks in debug mode, const version
    ENOKI_INLINE decltype(auto) operator[](size_t i) const {
        #if !defined(NDEBUG) && !defined(ENOKI_DISABLE_RANGE_CHECK)
            if (i >= derived().size())
                throw std::out_of_range(
                    "ArrayBase: out of range access (tried to access index " +
                    std::to_string(i) + " in an array of size " +
                    std::to_string(derived().size()) + ")");
        #endif
        return derived().coeff(i);
    }

    template <typename Mask, enable_if_mask_t<Mask> = 0>
    ENOKI_INLINE auto operator[](const Mask &m) {
        return detail::MaskedArray<Derived>{ derived(), (const mask_t<Derived> &) m };
    }

    //! @}
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    //! @{ \name Fallback implementations for masked operations
    // -----------------------------------------------------------------------

    #define ENOKI_MASKED_OPERATOR_FALLBACK(name, expr)                      \
        template <typename T, typename Mask>                                \
        ENOKI_INLINE void m##name##_(const T &e, const Mask &m) {           \
            derived() = select(m, expr, derived());                         \
        }

    ENOKI_MASKED_OPERATOR_FALLBACK(assign, e)
    ENOKI_MASKED_OPERATOR_FALLBACK(add, derived() + e)
    ENOKI_MASKED_OPERATOR_FALLBACK(sub, derived() - e)
    ENOKI_MASKED_OPERATOR_FALLBACK(mul, derived() * e)
    ENOKI_MASKED_OPERATOR_FALLBACK(div, derived() / e)
    ENOKI_MASKED_OPERATOR_FALLBACK(or, derived() | e)
    ENOKI_MASKED_OPERATOR_FALLBACK(and, derived() & e)
    ENOKI_MASKED_OPERATOR_FALLBACK(xor, derived() ^ e)

    #undef ENOKI_MASKED_OPERATOR_FALLBACK

    //! @}
    // -----------------------------------------------------------------------

    /// Dot product fallback implementation
    ENOKI_INLINE auto dot_(const Derived &a) const { return hsum(derived() * a); }

    /// Horizontal mean fallback implementation
    ENOKI_INLINE auto hmean_() const {
        return hsum(derived()) * (1.f / derived().size());
    }

    template <size_t Stride, typename Index, typename Mask>
    ENOKI_INLINE void scatter_add_(void *mem, const Index &index,
                                   const Mask &mask) const {
        transform<Derived, Stride>(
            mem, index, [](auto &a, auto &b, auto &) { a += b; },
            derived(), mask);
    }
};

namespace detail {
    template <typename T>
    ENOKI_INLINE bool convert_mask(T value) {
        if constexpr (std::is_same_v<T, bool>)
            return value;
        else
            return memcpy_cast<typename type_chooser<sizeof(T)>::UInt>(value) != 0;
    }

    template <typename Stream, typename Array, size_t N, typename... Indices>
    void print(Stream &os, const Array &a, bool abbrev,
               const std::array<size_t, N> &size, Indices... indices) {
        ENOKI_MARK_USED(size);
        ENOKI_MARK_USED(abbrev);
        if constexpr (sizeof...(Indices) == N) {
            os << a.derived().coeff(indices...);
        } else {
            constexpr size_t k = N - sizeof...(Indices) - 1;
            os << "[";
            for (size_t i = 0; i < size[k]; ++i) {
                if constexpr (is_dynamic_array_v<Array>) {
                    if (size[k] > 20 && i == 5 && abbrev) {
                        if (k > 0) {
                            os << ".. " << size[k] - 10 << " skipped ..,\n";
                            for (size_t j = 0; j <= sizeof...(Indices); ++j)
                                os << " ";
                        } else {
                            os << ".. " << size[k] - 10 << " skipped .., ";
                        }
                        i = size[k] - 6;
                        continue;
                    }
                }
                print(os, a, abbrev, size, i, indices...);
                if (i + 1 < size[k]) {
                    if constexpr (k == 0) {
                        os << ", ";
                    } else {
                        os << ",\n";
                        for (size_t j = 0; j <= sizeof...(Indices); ++j)
                            os << " ";
                    }
                }
            }
            os << "]";
        }
    }
}

template <typename Value, typename Derived>
ENOKI_NOINLINE std::ostream &operator<<(std::ostream &os, const ArrayBase<Value, Derived> &a) {
    if (ragged(a))
        os << "[ragged array]";
    else
        detail::print(os, a, true, shape(a));
    return os;
}


NAMESPACE_END(enoki)
