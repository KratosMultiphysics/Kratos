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
 * @file vtk_xml.hpp
 * @brief Shared VTK-XML `<DataArray>` helpers used by the VTU
 * (UnstructuredGrid) and VTP (PolyData) readers/writers.
 *
 * Both formats are the same XML container — `<VTKFile>` with per-array
 * `<DataArray>` elements in ascii or base64 "binary" encoding (raw or
 * zlib-compressed, framed by `detail/vtu_binary.hpp`) — differing only in
 * the grid element in between. This header holds the container-level pieces:
 * the DType <-> VTK type-name mapping (`vtu_type_str`/`dtype_from_vtu`),
 * ASCII float/array emission (`vtu_ascii_double`/`vtu_ascii_ndarray`),
 * DataArray text parsing for both encodings (`vtu_parse_ascii` /
 * `vtu_parse_binary`), and the `vtu_to_int64` widening helper.
 *
 * Deliberately pugixml-free (the single-header amalgamation only bundles
 * pugixml in its implementation section): the readers each keep a thin local
 * wrapper that pulls the `format`/`type`/`NumberOfComponents` attributes off
 * a `<DataArray>` node and dispatches to `vtu_parse_ascii`/`vtu_parse_binary`.
 */

// System includes
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/detail/value_io.hpp"
#include "meshioplusplus/detail/vtu_binary.hpp"
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Map an `NDArray` dtype to its VTK-XML type-name string.
 * @param dt The dtype.
 * @return The VTK type name (e.g. `"Float64"`, `"Int32"`).
 */
inline const char* vtu_type_str(DType dt) {
    switch (dt) {
        case DType::Float32:
            return "Float32";
        case DType::Float64:
            return "Float64";
        case DType::Int8:
            return "Int8";
        case DType::Int16:
            return "Int16";
        case DType::Int32:
            return "Int32";
        case DType::Int64:
            return "Int64";
        case DType::UInt8:
            return "UInt8";
        case DType::UInt16:
            return "UInt16";
        case DType::UInt32:
            return "UInt32";
        case DType::UInt64:
            return "UInt64";
    }
    return "Float64";
}

/**
 * @brief Map a VTK-XML type-name string to the `NDArray` dtype.
 * @param rS The VTK type name.
 * @return The matching dtype.
 * @throws ReadError on an unknown type name.
 */
inline DType dtype_from_vtu(const std::string& rS) {
    if (rS == "Float32")
        return DType::Float32;
    if (rS == "Float64")
        return DType::Float64;
    if (rS == "Int8")
        return DType::Int8;
    if (rS == "Int16")
        return DType::Int16;
    if (rS == "Int32")
        return DType::Int32;
    if (rS == "Int64")
        return DType::Int64;
    if (rS == "UInt8")
        return DType::UInt8;
    if (rS == "UInt16")
        return DType::UInt16;
    if (rS == "UInt32")
        return DType::UInt32;
    if (rS == "UInt64")
        return DType::UInt64;
    throw ReadError("Illegal VTU data type '" + rS + "'");
}

/**
 * @brief Emit one float in VTK's `%.11e` ASCII format followed by a newline.
 * @param rOs Output stream.
 * @param v The value.
 */
inline void vtu_ascii_double(std::ostream& rOs, double v) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%.11e", v);
    rOs << buf << '\n';
}

/**
 * @brief Emit a whole `NDArray` in ASCII, one value per line (floats via
 * `vtu_ascii_double`, integers as plain decimals).
 * @param rOs Output stream.
 * @param rA The array.
 */
inline void vtu_ascii_ndarray(std::ostream& rOs, const NDArray& rA) {
    const bool flt = is_float_dtype(rA.Dtype());
    const std::size_t n = rA.Size();
    for (std::size_t i = 0; i < n; ++i) {
        if (flt)
            vtu_ascii_double(rOs, read_double(rA, i));
        else
            rOs << read_int(rA, i) << '\n';
    }
}

/**
 * @brief Store one parsed value into a dtype-erased array slot.
 * @param rA Destination array.
 * @param i Flat index.
 * @param d The value when `rA` is a float dtype.
 * @param v The value when `rA` is an integer dtype.
 */
inline void vtu_store(NDArray& rA, std::size_t i, double d, std::int64_t v) {
    switch (rA.Dtype()) {
        case DType::Float32:
            rA.As<float>()[i] = static_cast<float>(d);
            break;
        case DType::Float64:
            rA.As<double>()[i] = d;
            break;
        case DType::Int8:
            rA.As<std::int8_t>()[i] = static_cast<std::int8_t>(v);
            break;
        case DType::Int16:
            rA.As<std::int16_t>()[i] = static_cast<std::int16_t>(v);
            break;
        case DType::Int32:
            rA.As<std::int32_t>()[i] = static_cast<std::int32_t>(v);
            break;
        case DType::Int64:
            rA.As<std::int64_t>()[i] = v;
            break;
        case DType::UInt8:
            rA.As<std::uint8_t>()[i] = static_cast<std::uint8_t>(v);
            break;
        case DType::UInt16:
            rA.As<std::uint16_t>()[i] = static_cast<std::uint16_t>(v);
            break;
        case DType::UInt32:
            rA.As<std::uint32_t>()[i] = static_cast<std::uint32_t>(v);
            break;
        case DType::UInt64:
            rA.As<std::uint64_t>()[i] = static_cast<std::uint64_t>(v);
            break;
    }
}

/**
 * @brief Parse whitespace-separated ASCII DataArray text into a flat array.
 * @param pText The element text (may be null).
 * @param dt Target dtype (drives float vs integer parsing).
 * @return A 1-D owning array of every parsed value.
 */
inline NDArray vtu_parse_ascii(const char* pText, DType dt) {
    const bool isflt = is_float_dtype(dt);
    std::vector<double> dv;
    std::vector<std::int64_t> iv;
    const char* p = pText ? pText : "";
    while (*p) {
        while (*p && std::isspace(static_cast<unsigned char>(*p)))
            ++p;
        if (!*p)
            break;
        char* endp = nullptr;
        if (isflt) {
            double x = std::strtod(p, &endp);
            if (endp == p)
                break;
            dv.push_back(x);
        } else {
            long long x = std::strtoll(p, &endp, 10);
            if (endp == p)
                break;
            iv.push_back(static_cast<std::int64_t>(x));
        }
        p = endp;
    }
    std::size_t n = isflt ? dv.size() : iv.size();
    NDArray a(dt, {n});
    for (std::size_t i = 0; i < n; ++i)
        vtu_store(a, i, isflt ? dv[i] : 0.0, isflt ? 0 : iv[i]);
    return a;
}

/**
 * @brief Trim leading/trailing whitespace from a C string.
 * @param pS The string (may be null).
 * @return The trimmed copy.
 */
inline std::string vtu_strip(const char* pS) {
    std::string t = pS ? pS : "";
    std::size_t b = 0, e = t.size();
    while (b < e && std::isspace(static_cast<unsigned char>(t[b])))
        ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(t[e - 1])))
        --e;
    return t.substr(b, e - b);
}

/**
 * @brief Decode a base64 "binary" DataArray payload into a flat array.
 * @param rText The stripped base64 text.
 * @param dt Target dtype.
 * @param compression 0 = none, 1 = zlib.
 * @param hsz Header integer size in bytes (4 for UInt32, 8 for UInt64).
 * @return A 1-D owning array over the decoded bytes.
 */
inline NDArray vtu_parse_binary(const std::string& rText, DType dt, int compression,
                                std::size_t hsz) {
    std::vector<unsigned char> bytes;
    if (compression == 0)
        bytes = vtu_decode_uncompressed(rText.c_str(), rText.size(), hsz);
    else
        bytes = vtu_decode_zlib(rText.c_str(), rText.size(), hsz);
    std::size_t isz = dtype_size(dt);
    std::size_t n = isz ? bytes.size() / isz : 0;
    NDArray a(dt, {n});
    if (n)
        std::memcpy(a.Data(), bytes.data(), n * isz);
    return a;
}

/**
 * @brief Widen a dtype-erased integer array to a `std::int64_t` vector.
 * @param rA The array.
 * @return The widened values.
 */
inline std::vector<std::int64_t> vtu_to_int64(const NDArray& rA) {
    std::vector<std::int64_t> v(rA.Size());
    for (std::size_t i = 0; i < rA.Size(); ++i)
        v[i] = read_int(rA, i);
    return v;
}

}  // namespace detail
}  // namespace meshioplusplus
