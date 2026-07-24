//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
//
//
//  License:         MIT License
//                   meshio++ default license: LICENSE
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#pragma once

/**
 * @file value_io.hpp
 * @brief Shared helpers to read scalar values out of an `NDArray` regardless
 * of its runtime dtype, used mainly by the ASCII format writers.
 *
 * ASCII writers need a single numeric value at a time (to format as text)
 * without caring whether the underlying `NDArray` holds `float`, `double`,
 * or any integer width — `read_double`/`read_int` do that dispatch once per
 * call. For hot loops where the same dispatch would otherwise happen inside
 * every iteration, `dispatch_dtype` hoists the `switch` on `DType` *outside*
 * the loop: it instantiates a caller-supplied templated lambda once per
 * concrete C++ type and lets the loop body run with a statically-typed
 * pointer, avoiding a per-element branch.
 */

// System includes
#include <cstddef>
#include <cstdint>

// Project includes
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Whether a dtype is one of the two floating-point kinds.
 * @param dt The dtype to test.
 * @return `true` for `Float32`/`Float64`, `false` for any integer dtype.
 */
inline bool is_float_dtype(DType dt) {
    return dt == DType::Float32 || dt == DType::Float64;
}

/**
 * @brief Reads element `i` of `a` as a `double`, regardless of `a`'s dtype.
 *
 * Dispatches on `a.dtype()` and `static_cast`s the underlying element
 * (narrowing for large 64-bit integers is possible but consistent with how
 * this codebase already treats "read as double" for display/ASCII purposes).
 * @param a Source array.
 * @param i Flat (linear) element index into `a`'s buffer.
 * @return `a`'s `i`-th element converted to `double`.
 */
inline double read_double(const NDArray& rA, std::size_t i) {
    switch (rA.Dtype()) {
        case DType::Float32:
            return static_cast<double>(rA.As<float>()[i]);
        case DType::Float64:
            return rA.As<double>()[i];
        case DType::Int8:
            return static_cast<double>(rA.As<std::int8_t>()[i]);
        case DType::Int16:
            return static_cast<double>(rA.As<std::int16_t>()[i]);
        case DType::Int32:
            return static_cast<double>(rA.As<std::int32_t>()[i]);
        case DType::Int64:
            return static_cast<double>(rA.As<std::int64_t>()[i]);
        case DType::UInt8:
            return static_cast<double>(rA.As<std::uint8_t>()[i]);
        case DType::UInt16:
            return static_cast<double>(rA.As<std::uint16_t>()[i]);
        case DType::UInt32:
            return static_cast<double>(rA.As<std::uint32_t>()[i]);
        case DType::UInt64:
            return static_cast<double>(rA.As<std::uint64_t>()[i]);
    }
    return 0.0;
}

/**
 * @brief Reads element `i` of `a` as an `int64_t`, regardless of `a`'s dtype.
 *
 * For any integer dtype this is a plain widening/narrowing cast; for a
 * floating-point dtype it falls back to `read_double` and truncates toward
 * zero via the `static_cast<int64_t>`.
 * @param a Source array.
 * @param i Flat (linear) element index into `a`'s buffer.
 * @return `a`'s `i`-th element converted to `int64_t`.
 */
inline std::int64_t read_int(const NDArray& rA, std::size_t i) {
    switch (rA.Dtype()) {
        case DType::Int8:
            return rA.As<std::int8_t>()[i];
        case DType::Int16:
            return rA.As<std::int16_t>()[i];
        case DType::Int32:
            return rA.As<std::int32_t>()[i];
        case DType::Int64:
            return rA.As<std::int64_t>()[i];
        case DType::UInt8:
            return rA.As<std::uint8_t>()[i];
        case DType::UInt16:
            return rA.As<std::uint16_t>()[i];
        case DType::UInt32:
            return rA.As<std::uint32_t>()[i];
        case DType::UInt64:
            return static_cast<std::int64_t>(rA.As<std::uint64_t>()[i]);
        default:
            return static_cast<std::int64_t>(read_double(rA, i));
    }
}

/**
 * @brief Number of rows (first-dimension extent) of `a`.
 * @param a Array to query.
 * @return `a.shape()[0]`, or 0 if `a` has no shape.
 */
inline std::size_t rows(const NDArray& rA) {
    return rA.Shape().empty() ? 0 : rA.Shape()[0];
}

/**
 * @brief Number of columns (second-dimension extent) of `a`, treating a
 * 1-D (or shapeless) array as having exactly one column.
 * @param a Array to query.
 * @return `a.shape()[1]` if `a` has at least 2 dimensions, else 1.
 */
inline std::size_t cols(const NDArray& rA) {
    return rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
}

/**
 * @brief Hoists a per-element `DType` switch out of a hot loop.
 *
 * Invokes the C++20 templated lambda `f.template operator()<T>()` with `T`
 * bound to the concrete C++ scalar type corresponding to `dt`, so the caller
 * writes the loop body once, generically, and gets a statically-typed
 * pointer (`a.as<T>()`) inside — the `switch` on `dt` happens exactly once,
 * not once per element:
 * @code
 * detail::dispatch_dtype(a.dtype(), [&]<class T>() {
 *     const T* src = a.as<T>();
 *     // ... plain, typed loop over src ...
 * });
 * @endcode
 *
 * @tparam F A callable with a templated `operator()<T>()` (a C++20 generic
 *           lambda with an explicit template parameter).
 * @param dt Runtime dtype selecting which instantiation of `f` to invoke.
 * @param f The generic callable to instantiate and invoke.
 * @return Whatever `f.template operator()<T>()` returns (perfectly forwarded
 *         via `decltype(auto)`).
 */
template <class F>
decltype(auto) dispatch_dtype(DType dt, F&& f) {
    switch (dt) {
        case DType::Float32:
            return f.template operator()<float>();
        case DType::Float64:
            return f.template operator()<double>();
        case DType::Int8:
            return f.template operator()<std::int8_t>();
        case DType::Int16:
            return f.template operator()<std::int16_t>();
        case DType::Int32:
            return f.template operator()<std::int32_t>();
        case DType::Int64:
            return f.template operator()<std::int64_t>();
        case DType::UInt8:
            return f.template operator()<std::uint8_t>();
        case DType::UInt16:
            return f.template operator()<std::uint16_t>();
        case DType::UInt32:
            return f.template operator()<std::uint32_t>();
        case DType::UInt64:
            return f.template operator()<std::uint64_t>();
    }
    return f.template operator()<double>();  // unreachable
}

}  // namespace detail
}  // namespace meshioplusplus
