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
 * @file byteswap.hpp
 * @brief Endianness conversion primitives used by every binary format reader
 * and writer that has to swap between file and host byte order.
 *
 * Every conversion goes through this header's `bswap16`/`32`/`64` (each a
 * single compiler intrinsic on GCC/Clang/MSVC, with a portable
 * shift-and-mask fallback) or the generic `bswap_copy`/`bswap_inplace`
 * helpers — never a hand-written per-byte reversal loop. Callers that need
 * to byte-swap many elements should combine these with `parallel_for_bw`
 * (parallel.hpp), since byte-swapping is a memory-bandwidth-bound operation.
 */

// System includes
#include <cstdint>
#include <cstring>

#ifdef _MSC_VER
#include <stdlib.h>
#endif

namespace meshioplusplus {
namespace detail {

/**
 * @brief Reverses the byte order of a 16-bit value.
 * @param v Value in one byte order.
 * @return `v` with its two bytes swapped.
 */
inline std::uint16_t bswap16(std::uint16_t v) {
#if defined(_MSC_VER)
    return _byteswap_ushort(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap16(v);
#else
    return static_cast<std::uint16_t>((v << 8) | (v >> 8));
#endif
}

/**
 * @brief Reverses the byte order of a 32-bit value.
 * @param v Value in one byte order.
 * @return `v` with its four bytes reversed.
 */
inline std::uint32_t bswap32(std::uint32_t v) {
#if defined(_MSC_VER)
    return _byteswap_ulong(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap32(v);
#else
    return ((v & 0x000000FFu) << 24) | ((v & 0x0000FF00u) << 8) | ((v & 0x00FF0000u) >> 8) |
           ((v & 0xFF000000u) >> 24);
#endif
}

/**
 * @brief Reverses the byte order of a 64-bit value.
 *
 * On the portable fallback path, implemented as two 32-bit swaps of the
 * high/low halves rather than a per-byte loop.
 * @param v Value in one byte order.
 * @return `v` with its eight bytes reversed.
 */
inline std::uint64_t bswap64(std::uint64_t v) {
#if defined(_MSC_VER)
    return _byteswap_uint64(v);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(v);
#else
    return (static_cast<std::uint64_t>(bswap32(static_cast<std::uint32_t>(v))) << 32) |
           bswap32(static_cast<std::uint32_t>(v >> 32));
#endif
}

/**
 * @brief Reverses `n` bytes (`n` in `{1,2,4,8}`) from `src` into `dst`, using
 * the matching intrinsic (`bswap16`/`32`/`64`) rather than a per-byte loop.
 *
 * `dst == src` is fine (the swap goes through a stack temporary), but `dst`
 * and `src` must not otherwise *partially* overlap.
 * @param dst Destination buffer, at least `n` bytes.
 * @param src Source buffer, at least `n` bytes.
 * @param n Element width in bytes: 1, 2, 4, or 8 (1 — or any other value —
 *          degrades to a plain copy, since a single byte has no order to
 *          reverse).
 */
inline void bswap_copy(char* pDst, const char* pSrc, int n) {
    switch (n) {
        case 8: {
            std::uint64_t v;
            std::memcpy(&v, pSrc, 8);
            v = bswap64(v);
            std::memcpy(pDst, &v, 8);
            break;
        }
        case 4: {
            std::uint32_t v;
            std::memcpy(&v, pSrc, 4);
            v = bswap32(v);
            std::memcpy(pDst, &v, 4);
            break;
        }
        case 2: {
            std::uint16_t v;
            std::memcpy(&v, pSrc, 2);
            v = bswap16(v);
            std::memcpy(pDst, &v, 2);
            break;
        }
        default:  // n == 1 (or unexpected): plain copy
            if (pDst != pSrc)
                std::memcpy(pDst, pSrc, static_cast<std::size_t>(n));
            break;
    }
}

/**
 * @brief In-place variant of `bswap_copy`: reverses `n` bytes at `p`.
 * @param pP Buffer to reverse in place, at least `n` bytes.
 * @param n Element width in bytes: 1, 2, 4, or 8.
 */
inline void bswap_inplace(char* pP, int n) {
    bswap_copy(pP, pP, n);
}

}  // namespace detail
}  // namespace meshioplusplus
