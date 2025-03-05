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

#include <enoki/array_generic.h>
#include <enoki/array_idiv.h>

NAMESPACE_BEGIN(enoki)

/// Analagous to meshgrid() in NumPy or MATLAB; for dynamic arrays
template <typename T, enable_if_dynamic_array_t<T> = 0>
Array<T, 2> meshgrid(const T &x, const T &y) {
    if constexpr (is_cuda_array_v<T> || is_diff_array_v<T>) {
        x.eval(); y.eval();

        if (x.size() == 1) {
            T x2(x);
            set_slices(x2, slices(y));
            return Array<T, 2>(
                std::move(x2),
                y
            );
        }

        uint32_t n = (uint32_t) x.size() * (uint32_t) y.size();
        divisor<uint32_t> div((uint32_t) x.size());

        using UInt32 = uint32_array_t<T>;
        UInt32 index = arange<UInt32>(n),
               yi    = div(index),
               xi    = index - yi * (uint32_t) x.size();

        return Array<T, 2>(
            gather<T>(x, xi),
            gather<T>(y, yi)
        );
    } else {
        T X, Y;
        set_slices(X, x.size() * y.size());
        set_slices(Y, x.size() * y.size());

        size_t pos = 0;

        if (x.size() % T::PacketSize == 0) {
            /* Fast path */

            for (size_t i = 0; i < y.size(); ++i) {
                for (size_t j = 0; j < packets(x); ++j) {
                    packet(X, pos) = packet(x, j);
                    packet(Y, pos) = y.coeff(i);
                    pos++;
                }
            }
        } else {
            for (size_t i = 0; i < y.size(); ++i) {
                for (size_t j = 0; j < x.size(); ++j) {
                    X.coeff(pos) = x.coeff(j);
                    Y.coeff(pos) = y.coeff(i);
                    pos++;
                }
            }
        }

        return Array<T, 2>(std::move(X), std::move(Y));
    }
}

/// Vectorized N-dimensional 'range' iterable with automatic mask computation
template <typename Value> struct range {
    static constexpr size_t Dimension = array_depth_v<Value> == 2 ?
        array_size_v<Value> : 1;
    static constexpr size_t PacketSize = array_depth_v<Value> == 2 ?
        array_size_v<value_t<Value>> : array_size_v<Value>;

    using Scalar = scalar_t<Value>;
    using Packet = Array<Scalar, PacketSize>;
    using Size   = Array<Scalar, Dimension>;

    struct iterator {
        iterator(size_t index) : index(index) { }
        iterator(size_t index, Size size)
            : index(index), index_p(arange<Packet>()), size(size) {
            for (size_t i = 0; i < Dimension - 1; ++i)
                div[i] = size[i];
        }

        bool operator==(const iterator &it) const { return it.index == index; }
        bool operator!=(const iterator &it) const { return it.index != index; }

        iterator &operator++() {
            index += 1;
            index_p += Scalar(Packet::Size);
            return *this;
        }

        std::pair<Value, mask_t<Packet>> operator*() const {
            if constexpr (array_depth_v<Value> == 1) {
                return { index_p, index_p < size[0] };
            } else {
                Value value;
                value[0] = index_p;
                ENOKI_UNROLL for (size_t i = 0; i < Dimension - 1; ++i)
                    value[i + 1] = div[i](value[i]);
                Packet offset = zero<Packet>();
                ENOKI_UNROLL for (size_t i = Dimension - 2; ; --i) {
                    offset = size[i] * (value[i + 1] + offset);
                    value[i] -= offset;
                    if (i == 0)
                        break;
                }

                return { value, value[Dimension - 1] < size[Dimension - 1] };
            }
        }

    private:
        size_t index;
        Packet index_p;
        Size size;
        divisor<Scalar> div[Dimension > 1 ? (Dimension - 1) : 1];
    };

    template <typename... Args>
    range(Args&&... args) : size(args...) { }

    iterator begin() {
        return iterator(0, size);
    }

    iterator end() {
        return iterator((hprod(size) + Packet::Size - 1) / Packet::Size);
    }

private:
    Size size;
};

template <typename Predicate,
          typename Args  = typename function_traits<Predicate>::Args,
          typename Index = std::decay_t<std::tuple_element_t<0, Args>>>
Index binary_search(scalar_t<Index> start_,
                    scalar_t<Index> end_,
                    const Predicate &pred) {
    Index start(start_), end(end_);

    scalar_t<Index> iterations = (start_ < end_) ?
        (log2i(end_ - start_) + 1) : 0;

    for (size_t i = 0; i < iterations; ++i) {
        Index middle = sr<1>(start + end);

        mask_t<Index> cond = pred(middle);

        masked(start,  cond) = min(middle + 1, end);
        masked(end,   !cond) = middle;
    }

    return start;
}

// -----------------------------------------------------------------------
//! @{ \name Stack memory allocation
// -----------------------------------------------------------------------

/**
 * \brief Wrapper around alloca(), which returns aligned (and, optionally,
 * zero-initialized) memory
 */
#define ENOKI_ALIGNED_ALLOCA(Array, Count, Clear)                             \
    enoki::detail::alloca_helper<Array, Clear>((uint8_t *) alloca(            \
        sizeof(Array) * (Count) + enoki::max_packet_size - 4),                \
        sizeof(Array) * (Count))

namespace detail {
    template <typename Array, bool Clear>
    ENOKI_INLINE Array *alloca_helper(uint8_t *ptr, size_t size) {
        (uintptr_t &) ptr +=
            ((max_packet_size - (uintptr_t) ptr) % max_packet_size);
        if constexpr (Clear)
            memset(ptr, 0, size);
        return (Array *) ptr;
    }
}

//! @}
// -----------------------------------------------------------------------

NAMESPACE_END(enoki)
