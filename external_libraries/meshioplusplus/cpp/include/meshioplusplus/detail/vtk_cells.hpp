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
 * @file vtk_cells.hpp
 * @brief Shared reconstruction of meshio cell blocks from the VTK/VTU
 * connectivity + end-offsets + types representation.
 *
 * Both the VTU reader and the VTK 5.1 legacy reader store cells in the same
 * layout — a flat `connectivity` array of node indices, an `offsets` array
 * giving each cell's end position within it, and a `types` array giving each
 * cell's VTK type id — so this header's `detail::reconstruct_cells` (ported
 * from `vtk_cells_from_data` in `_vtk_common.py`) is the single place that
 * turns that layout back into meshio's per-type cell-block representation
 * (appended straight onto the output `Mesh`),
 * grouping consecutive same-type runs and further splitting runs of
 * variable-node-count types (polygon, VTK_LAGRANGE_*) by per-cell size.
 * It leans heavily on `parallel_for_bw`/`parallel_copy_i64` (memory-gather
 * and memory-fault-bound work) since reconstructing connectivity is pure
 * data movement, not compute.
 */

// System includes
#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/parallel.hpp"
#include "meshioplusplus/types.hpp"
#include "meshioplusplus/vtk_common.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief Copies `n` `int64_t` elements from `src` to `dst`, splitting the
 * copy into large contiguous chunks run across `parallel_for_bw`'s
 * bandwidth-capped threads.
 *
 * `dst` is assumed to be a fresh allocation, so most of the wall-clock cost
 * is first-touch page faults (the OS zeroing/mapping pages on first write)
 * rather than the memcpy itself — servicing those faults concurrently across
 * a few threads beats one thread doing a single serial `memcpy`. Falls back
 * to a single sequential `memcpy` when `n` doesn't even fill one 4 MiB chunk
 * (`nchunks <= 1`). Uses `grain=1` so every chunk (already coarse at 512Ki
 * elements) dispatches individually rather than being batched further by the
 * default grain.
 * @param pDst Destination buffer, at least `n` elements, ideally freshly
 *            allocated (unfaulted) memory.
 * @param pSrc Source buffer, at least `n` elements.
 * @param n Number of `int64_t` elements to copy.
 */
inline void parallel_copy_i64(std::int64_t* pDst, const std::int64_t* pSrc, std::size_t n) {
    constexpr std::size_t kChunk = 1u << 19;  // 512Ki elements (4 MiB) per task
    const std::size_t nchunks = (n + kChunk - 1) / kChunk;
    if (nchunks <= 1) {
        std::memcpy(pDst, pSrc, n * sizeof(std::int64_t));
        return;
    }
    // grain=1: each chunk is already coarse (4 MiB), so dispatch per chunk —
    // otherwise the default grain (2048) would run these few chunks serially.
    parallel_for_bw(
        nchunks,
        [&](std::size_t c) {
            const std::size_t off = c * kChunk;
            const std::size_t len = std::min(kChunk, n - off);
            std::memcpy(pDst + off, pSrc + off, len * sizeof(std::int64_t));
        },
        1);
}

/**
 * @brief Extracts rows `[r0, r1)` of a 2-D (or column-vector) `NDArray` into
 * a new, freshly-allocated `NDArray`.
 *
 * The output buffer is allocated via `NDArray::Uninit` (skipping the
 * zero-fill) since the single `memcpy` below fully overwrites it.
 * @param rA Source array; row size is `rA.Shape()[1]` if 2-D, else 1.
 * @param r0 First row to include (inclusive).
 * @param r1 One past the last row to include (exclusive).
 * @return A new owning `NDArray` with `r1 - r0` rows, same dtype/row-width as `rA`.
 */
inline NDArray slice_rows(const NDArray& rA, std::size_t r0, std::size_t r1) {
    std::size_t nc = rA.Shape().size() >= 2 ? rA.Shape()[1] : 1;
    std::size_t isz = dtype_size(rA.Dtype());
    std::size_t rowbytes = nc * isz;
    std::vector<std::size_t> shape = rA.Shape();
    if (shape.empty())
        shape = {0};
    shape[0] = r1 - r0;
    NDArray out = NDArray::Uninit(rA.Dtype(), shape);  // fully overwritten below
    if (r1 > r0)
        std::memcpy(out.Data(), rA.Data() + r0 * rowbytes, (r1 - r0) * rowbytes);
    return out;
}

/**
 * @brief Reconstructs meshio cell blocks (and the matching per-block
 * `cell_data`) from the VTK/VTU flat connectivity + end-offsets + types
 * representation, appending them to @p rMesh.
 *
 * Ported from `vtk_cells_from_data` in `_vtk_common.py`; shared by the VTU
 * reader and the VTK 5.1 legacy reader, which store cells identically.
 * Walks `types` and groups consecutive cells of the same VTK type into a
 * run; a run of a *fixed*-node-count type becomes one rectangular cell
 * block (data gathered per-row via `vtk_to_meshio_order`, or
 * block-copied via `parallel_copy_i64` when the run is contiguous in
 * `conn` with no reordering needed); a run of a *variable*-node-count type
 * (`is_special_cell`, e.g. polygon or VTK_LAGRANGE_*) is further split into
 * sub-runs of a single common node count each, since a rectangular cell
 * block still requires a `(num_cells, n)` layout — each such sub-run
 * is emitted as its own separate block sharing the same meshio type
 * name. Matching slices of every array in `cell_data_raw` are appended to
 * @p rMesh's cell data in lockstep with the appended blocks, via
 * `slice_rows`.
 *
 * @param pConn Flat node-index connectivity buffer. Passed as a raw
 *             `int64_t*` (rather than an `NDArray`) so callers can hand in
 *             an `NDArray`'s buffer directly — VTK 5.1 connectivity is
 *             already `vtktypeint64` — without an intermediate
 *             to-int64 copy.
 * @param rOffsets End offsets, one per cell: `rOffsets[i]` is the index in
 *                `pConn` just past cell `i`'s last node (so cell `i`'s nodes
 *                are `pConn[rOffsets[i-1] .. rOffsets[i])`, with `rOffsets[-1]`
 *                treated as 0).
 * @param rTypes VTK cell type id for each cell, same length as `rOffsets`.
 * @param rCellDataRaw Per-name cell-data arrays covering the whole mesh
 *                      (all cells concatenated), to be re-sliced per output
 *                      block.
 * @param rMesh Mesh appended to: one rectangular cell block per contiguous
 *              same-type (and, for special types, same-size) run, plus — in
 *              lockstep — for each name in `rCellDataRaw` one sliced
 *              `NDArray` per new block.
 * @throws ReadError if a cell's VTK type id is 42 (polyhedron — unsupported
 *         by the C++ reader) or is otherwise not in `vtk_to_meshio_type()`,
 *         or if a resolved meshio type has no entry in `num_nodes_per_cell()`.
 */
inline void reconstruct_cells(const std::int64_t* pConn, const std::vector<std::int64_t>& rOffsets,
                              const std::vector<std::int64_t>& rTypes,
                              const std::unordered_map<std::string, NDArray>& rCellDataRaw,
                              Mesh& rMesh) {
    const auto& vmap = vtk_to_meshio_type();
    const std::size_t ncells = rTypes.size();

    auto add_cd = [&](std::size_t start, std::size_t end) {
        for (const auto& kv : rCellDataRaw)
            rMesh.AppendCellData(kv.first, slice_rows(kv.second, start, end));
    };

    std::size_t start = 0;
    while (start < ncells) {
        std::size_t end = start + 1;
        while (end < ncells && rTypes[end] == rTypes[start])
            ++end;

        int vtk_type = static_cast<int>(rTypes[start]);
        if (vtk_type == 42)
            throw ReadError("polyhedron cells are not supported by the C++ reader");
        auto it = vmap.find(vtk_type);
        if (it == vmap.end())
            throw ReadError("VTK cell type " + std::to_string(vtk_type) +
                            " not supported by the C++ reader");
        const std::string& meshio_type = it->second;

        if (is_special_cell(meshio_type)) {
            std::int64_t first_node = (start == 0) ? 0 : rOffsets[start - 1];
            std::vector<std::int64_t> start_cn;
            start_cn.reserve(end - start + 1);
            start_cn.push_back(first_node);
            for (std::size_t i = start; i < end; ++i)
                start_cn.push_back(rOffsets[i]);
            std::vector<std::int64_t> sizes(end - start);
            for (std::size_t i = 0; i < sizes.size(); ++i)
                sizes[i] = start_cn[i + 1] - start_cn[i];

            std::size_t i = 0;
            while (i < sizes.size()) {
                std::size_t j = i;
                while (j < sizes.size() && sizes[j] == sizes[i])
                    ++j;
                std::int64_t sz = sizes[i];
                std::size_t m = j - i;
                NDArray data = NDArray::Uninit(DType::Int64, {m, static_cast<std::size_t>(sz)});
                std::int64_t* out = data.As<std::int64_t>();
                const std::size_t ii = i;
                // Contiguous uniform-size sub-run -> block memcpy.
                const std::int64_t sub_first = start_cn[ii];
                bool sub_regular = true;
                for (std::size_t r = 0; sub_regular && r < m; ++r)
                    if (rOffsets[start + ii + r] !=
                        sub_first + static_cast<std::int64_t>(r + 1) * sz)
                        sub_regular = false;
                if (sub_regular) {
                    parallel_copy_i64(out, pConn + sub_first, m * static_cast<std::size_t>(sz));
                } else {
                    parallel_for_bw(m, [&](std::size_t r) {
                        std::int64_t endoff = rOffsets[start + ii + r];
                        std::int64_t base = endoff - sz;
                        for (std::int64_t c = 0; c < sz; ++c)
                            out[r * sz + c] = pConn[base + c];
                    });
                }
                rMesh.AddCellBlock(meshio_type, std::move(data));
                add_cd(start + i, start + j);
                i = j;
            }
        } else {
            auto nit = num_nodes_per_cell().find(meshio_type);
            if (nit == num_nodes_per_cell().end())
                throw ReadError("Unknown node count for cell type " + meshio_type);
            int n = nit->second;
            std::vector<int> order = vtk_to_meshio_order(vtk_type);
            std::size_t m = end - start;
            NDArray data = NDArray::Uninit(DType::Int64, {m, static_cast<std::size_t>(n)});
            std::int64_t* out = data.As<std::int64_t>();
            const int* ord = order.empty() ? nullptr : order.data();
            const std::size_t ss = start;
            // Regular run (offsets advance by exactly n per cell) with identity
            // node order -> the run's connectivity is one contiguous slice:
            // block memcpy instead of a per-row gather.
            const std::int64_t first = (ss == 0) ? 0 : rOffsets[ss - 1];
            bool regular = true;
            for (std::size_t r = 0; regular && r < m; ++r)
                if (rOffsets[ss + r] !=
                    first + static_cast<std::int64_t>((r + 1) * static_cast<std::size_t>(n)))
                    regular = false;
            if (!ord && regular) {
                // Contiguous slice -> parallel block copy (fault-bound).
                parallel_copy_i64(out, pConn + first, m * static_cast<std::size_t>(n));
            } else {
                parallel_for_bw(m, [&](std::size_t r) {
                    std::int64_t endoff = rOffsets[ss + r];
                    std::int64_t base = endoff - n;
                    for (int j = 0; j < n; ++j) {
                        int col = ord ? ord[j] : j;
                        out[r * n + j] = pConn[base + col];
                    }
                });
            }
            rMesh.AddCellBlock(meshio_type, std::move(data));
            add_cd(start, end);
        }
        start = end;
    }
}

}  // namespace detail
}  // namespace meshioplusplus
