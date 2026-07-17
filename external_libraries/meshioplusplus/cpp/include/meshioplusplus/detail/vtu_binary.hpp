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
 * @file vtu_binary.hpp
 * @brief Base64 and VTU "binary" `DataArray` codecs (raw and zlib-compressed),
 * shared helpers behind the VTU (VTK XML) reader/writer's binary I/O.
 *
 * VTU's binary encoding wraps raw little-endian bytes (optionally
 * zlib-deflated in fixed-size blocks) as base64 text inside the XML. This
 * header provides both halves: plain base64 encode/decode
 * (`b64encode`/`b64decode`), and the VTU-specific framing on top of it —
 * `vtu_decode_uncompressed`/`vtu_encode_binary(zlib_compress=false)` for the
 * uncompressed scheme (a little-endian byte-count header followed by raw
 * data) and `vtu_decode_zlib`/`vtu_encode_binary(zlib_compress=true)` for the
 * compressed block scheme (num_blocks / max_block_size / last_block_size
 * header, then each block's compressed size, then the concatenated deflated
 * blocks). zlib support is conditionally compiled on
 * `MESHIOPLUSPLUS_HAS_ZLIB`; without it, the zlib-specific functions throw
 * rather than compiling out entirely, since they're still callable — the
 * absence is discovered at runtime and routes the caller to the Python
 * fallback. Both directions parallelize per independent unit of work
 * (base64 3-byte groups; zlib blocks) via `parallel_for`, since base64/zlib
 * are genuinely compute-bound (unlike the memory-bandwidth-bound gather/
 * byteswap work elsewhere, which uses `parallel_for_bw` instead).
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// External includes
#ifdef MESHIOPLUSPLUS_HAS_ZLIB
#include <zlib.h>
#endif

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/parallel.hpp"

namespace meshioplusplus {
namespace detail {

/** @brief The standard base64 alphabet (RFC 4648), indexed by 6-bit value. */
inline const char* b64_table() {
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
}

/**
 * @brief Base64-encodes `len` bytes of `data`.
 *
 * Every full 3-byte input group maps to exactly 4 output characters at a
 * fixed, independently-computable offset, so the output string is
 * pre-sized once and each group is encoded directly into its slot via
 * `parallel_for` — no synchronization or intermediate buffering needed. Any
 * trailing 1- or 2-byte group is handled afterward, sequentially, with the
 * standard `'='` padding.
 * @param pData Bytes to encode.
 * @param len Number of bytes in `pData`.
 * @return The base64-encoded text, `'='`-padded to a multiple of 4 characters.
 */
inline std::string b64encode(const unsigned char* pData, std::size_t len) {
    const char* tbl = b64_table();
    // Every 3-byte group maps to 4 output chars at a deterministic offset:
    // pre-size the output and write by index -> parallel over groups.
    const std::size_t ngroups = len / 3;  // full groups
    std::string out(((len + 2) / 3) * 4, '\0');
    parallel_for(ngroups, [&](std::size_t g) {
        const std::size_t i = g * 3;
        unsigned n =
            (unsigned(pData[i]) << 16) | (unsigned(pData[i + 1]) << 8) | unsigned(pData[i + 2]);
        char* o = out.data() + g * 4;
        o[0] = tbl[(n >> 18) & 63];
        o[1] = tbl[(n >> 12) & 63];
        o[2] = tbl[(n >> 6) & 63];
        o[3] = tbl[n & 63];
    });
    const std::size_t i = ngroups * 3;
    if (i < len) {  // trailing 1- or 2-byte group with '=' padding
        const bool two = (i + 1 < len);
        unsigned n = unsigned(pData[i]) << 16;
        if (two)
            n |= unsigned(pData[i + 1]) << 8;
        char* o = out.data() + ngroups * 4;
        o[0] = tbl[(n >> 18) & 63];
        o[1] = tbl[(n >> 12) & 63];
        o[2] = two ? tbl[(n >> 6) & 63] : '=';
        o[3] = '=';
    }
    return out;
}

/**
 * @brief Base64-decodes `len` characters of `s`.
 *
 * Builds (and caches, in a function-local `static`) an inverse lookup table
 * from ASCII byte to 6-bit value on first call. Silently skips `'='`
 * padding and whitespace (`\n \r space \t`), and silently ignores any other
 * character outside the base64 alphabet, rather than treating either as an
 * error — VTU-embedded base64 can be split across lines.
 * @param pS Base64 text to decode (need not be NUL-terminated; length is explicit).
 * @param len Number of characters in `pS` to consider.
 * @return The decoded raw bytes.
 */
inline std::vector<unsigned char> b64decode(const char* pS, std::size_t len) {
    static int8_t inv[256];
    static bool init = false;
    if (!init) {
        for (int i = 0; i < 256; ++i)
            inv[i] = -1;
        const char* tbl = b64_table();
        for (int i = 0; i < 64; ++i)
            inv[(unsigned char)tbl[i]] = static_cast<int8_t>(i);
        init = true;
    }
    std::vector<unsigned char> out;
    out.reserve(len / 4 * 3);
    int buf = 0, bits = 0;
    for (std::size_t i = 0; i < len; ++i) {
        char ch = pS[i];
        if (ch == '=' || ch == '\n' || ch == '\r' || ch == ' ' || ch == '\t')
            continue;
        int v = inv[(unsigned char)ch];
        if (v < 0)
            continue;
        buf = (buf << 6) | v;
        bits += 6;
        if (bits >= 8) {
            bits -= 8;
            out.push_back(static_cast<unsigned char>((buf >> bits) & 0xFF));
        }
    }
    return out;
}

#ifdef MESHIOPLUSPLUS_HAS_ZLIB
/**
 * @brief Compresses one block with zlib's default `compress()` (a single
 * deflate call, no streaming).
 * @param pSrc Bytes to compress.
 * @param n Number of bytes in `pSrc`.
 * @return The compressed bytes (sized to zlib's actual output, not the bound).
 * @throws WriteError if zlib does not return `Z_OK`.
 */
inline std::vector<unsigned char> zlib_compress_block(const unsigned char* pSrc, std::size_t n) {
    uLongf bound = compressBound(static_cast<uLong>(n));
    std::vector<unsigned char> out(bound);
    uLongf destLen = bound;
    int r = compress(out.data(), &destLen, pSrc, static_cast<uLong>(n));
    if (r != Z_OK)
        throw WriteError("zlib compression failed");
    out.resize(destLen);
    return out;
}

/**
 * @brief Decompresses one zlib-compressed block whose decompressed size is
 * already known.
 * @param pSrc Compressed bytes.
 * @param n Number of compressed bytes in `pSrc`.
 * @param expected Exact expected decompressed size (from the VTU block header).
 * @return The decompressed bytes.
 * @throws ReadError if zlib does not return `Z_OK`.
 */
inline std::vector<unsigned char> zlib_decompress(const unsigned char* pSrc, std::size_t n,
                                                  std::size_t expected) {
    std::vector<unsigned char> out(expected);
    uLongf destLen = static_cast<uLongf>(expected);
    int r = uncompress(out.data(), &destLen, pSrc, static_cast<uLong>(n));
    if (r != Z_OK)
        throw ReadError("zlib decompression failed");
    out.resize(destLen);
    return out;
}
#endif  // MESHIOPLUSPLUS_HAS_ZLIB

/**
 * @brief Reads a little-endian unsigned integer of `isz` bytes from `p`.
 * @param pP Buffer to read from, at least `isz` bytes.
 * @param isz Width in bytes of the integer to read (typically 4 or 8, the
 *            VTU header_type item size).
 * @return The decoded value, widened to `uint64_t`.
 */
inline std::uint64_t read_uint_le(const unsigned char* pP, std::size_t isz) {
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < isz; ++i)
        v |= static_cast<std::uint64_t>(pP[i]) << (8 * i);
    return v;
}

/**
 * @brief Decodes an uncompressed VTU "binary" `DataArray`: base64 text of a
 * little-endian byte-count header followed by the raw payload.
 * @param pText Base64-encoded DataArray text.
 * @param len Length of `pText` in characters.
 * @param hsz `header_type` item size in bytes (4 for `UInt32`, 8 for `UInt64`).
 * @return The decoded raw payload bytes (header stripped).
 * @throws ReadError if the decoded data is shorter than the header, or
 *         shorter than the header declares.
 */
inline std::vector<unsigned char> vtu_decode_uncompressed(const char* pText, std::size_t len,
                                                          std::size_t hsz) {
    std::vector<unsigned char> all = b64decode(pText, len);
    if (all.size() < hsz)
        throw ReadError("VTU binary data too short");
    std::uint64_t total = read_uint_le(all.data(), hsz);
    if (all.size() < hsz + total)
        throw ReadError("VTU binary data truncated");
    return std::vector<unsigned char>(all.begin() + hsz, all.begin() + hsz + total);
}

/**
 * @brief Decodes a zlib-compressed VTU "binary" `DataArray` (the VTK block
 * compression scheme).
 *
 * The format, all base64-encoded: a header of `num_blocks`, `max_block`,
 * `last_block` (each `hsz` bytes), then `num_blocks` compressed-size
 * entries, then the concatenated deflated blocks themselves (each block
 * `max_block` bytes decompressed, except the last which is `last_block`).
 * Decoded in three passes: decode just enough base64 to learn
 * `num_blocks`, decode the rest of the header to get each block's
 * compressed size, then base64-decode the block data. Input offsets are a
 * cheap sequential prefix sum of the per-block compressed sizes; output
 * offsets are `k * max_block` by construction, so with both known up front
 * the per-block `inflate` calls are independent and run under
 * `parallel_for` with `grain=1` (each ~32 KiB block is a full unit of
 * inflate work, so per-block dispatch is exactly right — this is
 * compute-bound work, unlike the memory-gather use of `parallel_for_bw`
 * elsewhere).
 *
 * @param pText Base64-encoded DataArray text.
 * @param len Length of `pText` in characters.
 * @param hsz `header_type` item size in bytes (4 for `UInt32`, 8 for `UInt64`).
 * @return The decoded, decompressed raw payload bytes (all blocks concatenated).
 * @throws ReadError if built without `MESHIOPLUSPLUS_HAS_ZLIB`, or if the
 *         header/data is truncated, or if any block fails to decompress.
 */
inline std::vector<unsigned char> vtu_decode_zlib(const char* pText, std::size_t len,
                                                  std::size_t hsz) {
#ifndef MESHIOPLUSPLUS_HAS_ZLIB
    (void)pText;
    (void)len;
    (void)hsz;
    throw ReadError("VTU zlib decompression requires a zlib-enabled build");
#else
    std::size_t first_chars = ((hsz + 2) / 3) * 4;
    if (len < first_chars)
        throw ReadError("VTU zlib header too short");
    std::vector<unsigned char> hb = b64decode(pText, first_chars);
    std::uint64_t num_blocks = read_uint_le(hb.data(), hsz);

    std::size_t num_header_bytes = hsz * (3 + static_cast<std::size_t>(num_blocks));
    std::size_t num_header_chars = ((num_header_bytes + 2) / 3) * 4;
    if (len < num_header_chars)
        throw ReadError("VTU zlib header truncated");
    std::vector<unsigned char> header = b64decode(pText, num_header_chars);

    std::uint64_t max_block = read_uint_le(header.data() + hsz, hsz);
    std::uint64_t last_block = read_uint_le(header.data() + 2 * hsz, hsz);
    std::vector<std::uint64_t> comp_sizes(num_blocks);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        comp_sizes[k] = read_uint_le(header.data() + (3 + k) * hsz, hsz);

    std::vector<unsigned char> blockdata =
        b64decode(pText + num_header_chars, len - num_header_chars);

    // Input offsets are a (cheap, sequential) prefix sum of comp_sizes; the
    // output offset of block k is k*max_block per the VTU block scheme -> the
    // per-block inflate runs in parallel into a pre-sized buffer.
    std::vector<std::size_t> in_off(static_cast<std::size_t>(num_blocks) + 1, 0);
    for (std::uint64_t k = 0; k < num_blocks; ++k)
        in_off[static_cast<std::size_t>(k) + 1] =
            in_off[static_cast<std::size_t>(k)] + static_cast<std::size_t>(comp_sizes[k]);

    const std::size_t total = num_blocks ? static_cast<std::size_t>(num_blocks - 1) *
                                                   static_cast<std::size_t>(max_block) +
                                               static_cast<std::size_t>(last_block)
                                         : 0;
    std::vector<unsigned char> out(total);
    parallel_for(
        static_cast<std::size_t>(num_blocks),
        [&](std::size_t k) {
            std::size_t expected = (k + 1 == num_blocks) ? static_cast<std::size_t>(last_block)
                                                         : static_cast<std::size_t>(max_block);
            auto dec = zlib_decompress(blockdata.data() + in_off[k],
                                       static_cast<std::size_t>(comp_sizes[k]), expected);
            std::memcpy(out.data() + k * static_cast<std::size_t>(max_block), dec.data(),
                        std::min(dec.size(), expected));
        },
        /*grain=*/1);  // each block is 32 KB of inflate work
    return out;
#endif  // MESHIOPLUSPLUS_HAS_ZLIB
}

/**
 * @brief Encodes raw little-endian bytes as a VTU "binary" `DataArray` text,
 * either uncompressed or zlib-compressed (block scheme).
 *
 * Uncompressed (`zlib_compress == false`): a 4-byte little-endian length
 * header followed by the raw bytes, base64-encoded as one unit.
 *
 * Compressed (`zlib_compress == true`): splits `data` into fixed 32 KiB
 * blocks, deflates each independently under `parallel_for` with `grain=1`
 * (each block is a full, sizeable unit of compute — one whole deflate call
 * — so per-block dispatch is ideal; this is compute-bound, unlike the
 * memory-gather work that uses `parallel_for_bw`), then emits the
 * `num_blocks`/`max_block`/`last_block_size`/per-block-compressed-size
 * header followed by the concatenated compressed blocks, all base64-encoded.
 *
 * @param pData Raw bytes to encode (already in the file's target byte order).
 * @param nbytes Number of bytes in `pData`.
 * @param zlib_compress Whether to zlib-compress (block scheme) or emit raw.
 * @return The base64-encoded VTU `DataArray` text.
 * @throws WriteError if `zlib_compress` is requested but the build lacks
 *         `MESHIOPLUSPLUS_HAS_ZLIB`.
 */
inline std::string vtu_encode_binary(const unsigned char* pData, std::size_t nbytes,
                                     bool zlib_compress) {
    if (!zlib_compress) {
        std::vector<unsigned char> buf(4 + nbytes);
        std::uint32_t header = static_cast<std::uint32_t>(nbytes);
        std::memcpy(buf.data(), &header, 4);
        if (nbytes)
            std::memcpy(buf.data() + 4, pData, nbytes);
        return b64encode(buf.data(), buf.size());
    }

#ifndef MESHIOPLUSPLUS_HAS_ZLIB
    throw WriteError("VTU zlib compression requires a zlib-enabled build");
#else
    const std::uint32_t max_block = 32768;
    std::uint32_t num_blocks = static_cast<std::uint32_t>((nbytes + max_block - 1) / max_block);
    std::uint32_t last_block_size =
        num_blocks ? static_cast<std::uint32_t>(nbytes - std::size_t(num_blocks - 1) * max_block)
                   : max_block;

    // Blocks are independent -> compress in parallel into pre-sized slots.
    std::vector<std::vector<unsigned char> > blocks(num_blocks);
    parallel_for(
        num_blocks,
        [&](std::size_t b) {
            std::size_t off = b * max_block;
            std::size_t len = std::min<std::size_t>(max_block, nbytes - off);
            blocks[b] = zlib_compress_block(pData + off, len);
        },
        /*grain=*/1);  // each block is 32 KB of deflate work

    std::vector<std::uint32_t> header;
    header.reserve(3 + num_blocks);
    header.push_back(num_blocks);
    header.push_back(max_block);
    header.push_back(last_block_size);
    for (const auto& b : blocks)
        header.push_back(static_cast<std::uint32_t>(b.size()));

    std::string out = b64encode(reinterpret_cast<const unsigned char*>(header.data()),
                                header.size() * sizeof(std::uint32_t));
    std::vector<unsigned char> concat;
    for (const auto& b : blocks)
        concat.insert(concat.end(), b.begin(), b.end());
    out += b64encode(concat.data(), concat.size());
    return out;
#endif  // MESHIOPLUSPLUS_HAS_ZLIB
}

}  // namespace detail
}  // namespace meshioplusplus
