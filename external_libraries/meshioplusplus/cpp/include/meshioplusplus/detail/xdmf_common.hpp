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
 * @file xdmf_common.hpp
 * @brief XDMF cell-type-name maps and cell-data raw<->blocks conversion,
 * shared between the XDMF format implementation and HMF (which reuses
 * XDMF's topology names and raw cell-data layout).
 *
 * Ported from `src/meshio/xdmf/common.py` and the `raw_from_cell_data` /
 * `cell_data_from_raw` helpers in `src/meshio/_common.py`. XDMF (and HMF)
 * store per-cell-type-block data as one array *per cell type name string* in
 * the XML/XDMF Topology, and store cell_data for a mixed mesh as one
 * concatenated raw array per data name (all cell blocks laid end-to-end)
 * rather than one array per block — `concat_cell_data`/`split_raw_cell_data`
 * are what let this header's callers go between meshio's per-block
 * `cell_data` representation and that concatenated-raw representation.
 */

// System includes
#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

// Project includes
#include "meshioplusplus/exceptions.hpp"
#include "meshioplusplus/mesh.hpp"
#include "meshioplusplus/ndarray.hpp"

namespace meshioplusplus {
namespace xdmfcommon {

/**
 * @brief Maps a meshio cell-type name to its XDMF topology type name.
 * @param t meshio cell-type name (e.g. `"triangle"`, `"tetra10"`).
 * @return The corresponding XDMF `TopologyType` string (e.g. `"Triangle"`).
 * @throws WriteError if `t` has no XDMF equivalent.
 */
inline const char* meshio_to_xdmf(const std::string& rT) {
    static const std::unordered_map<std::string, const char*> m = {
        {"vertex", "Polyvertex"},
        {"line", "Polyline"},
        {"line3", "Edge_3"},
        {"quad", "Quadrilateral"},
        {"quad8", "Quadrilateral_8"},
        {"quad9", "Quadrilateral_9"},
        {"pyramid", "Pyramid"},
        {"pyramid13", "Pyramid_13"},
        {"tetra", "Tetrahedron"},
        {"triangle", "Triangle"},
        {"triangle6", "Triangle_6"},
        {"tetra10", "Tetrahedron_10"},
        {"wedge", "Wedge"},
        {"wedge15", "Wedge_15"},
        {"wedge18", "Wedge_18"},
        {"hexahedron", "Hexahedron"},
        {"hexahedron20", "Hexahedron_20"},
        {"hexahedron24", "Hexahedron_24"},
        {"hexahedron27", "Hexahedron_27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw WriteError("XDMF: unsupported cell type " + rT);
    return it->second;
}

/**
 * @brief Maps an XDMF topology type name to a meshio cell-type name.
 *
 * Accepts both the canonical XDMF spelling and common abbreviations some
 * writers emit (e.g. both `"Hexahedron_20"` and `"Hex_20"` map to
 * `"hexahedron20"`).
 * @param t XDMF `TopologyType` string as found in the file.
 * @return The corresponding meshio cell-type name.
 * @throws ReadError if `t` is not a recognized topology type.
 */
inline std::string xdmf_to_meshio(const std::string& rT) {
    static const std::unordered_map<std::string, std::string> m = {
        {"Polyvertex", "vertex"},
        {"Polyline", "line"},
        {"Edge_3", "line3"},
        {"Quadrilateral", "quad"},
        {"Quadrilateral_8", "quad8"},
        {"Quad_8", "quad8"},
        {"Quadrilateral_9", "quad9"},
        {"Quad_9", "quad9"},
        {"Pyramid", "pyramid"},
        {"Pyramid_13", "pyramid13"},
        {"Tetrahedron", "tetra"},
        {"Triangle", "triangle"},
        {"Triangle_6", "triangle6"},
        {"Tri_6", "triangle6"},
        {"Tetrahedron_10", "tetra10"},
        {"Tet_10", "tetra10"},
        {"Wedge", "wedge"},
        {"Wedge_15", "wedge15"},
        {"Wedge_18", "wedge18"},
        {"Hexahedron", "hexahedron"},
        {"Hexahedron_20", "hexahedron20"},
        {"Hex_20", "hexahedron20"},
        {"Hexahedron_24", "hexahedron24"},
        {"Hex_24", "hexahedron24"},
        {"Hexahedron_27", "hexahedron27"},
        {"Hex_27", "hexahedron27"}};
    auto it = m.find(rT);
    if (it == m.end())
        throw ReadError("XDMF: unsupported topology type " + rT);
    return it->second;
}

/**
 * @brief Concatenates one cell-data name's per-block arrays along axis 0
 * into a single raw array, matching Python's `raw_from_cell_data`.
 *
 * Used when writing XDMF/HMF cell data for a mixed-cell-type mesh: XDMF
 * stores cell data as one flat array per data name (all blocks' rows
 * back-to-back) rather than one array per block.
 * @param rMesh The mesh whose cell data to concatenate.
 * @param rName The cell-data name; must have at least one block, all blocks
 *              sharing dtype/trailing shape.
 * @return A new array with the same trailing shape as the first block and
 *         first dimension equal to the sum of each block's row count.
 */
inline NDArray concat_cell_data(const Mesh& rMesh, const std::string& rName) {
    const std::size_t nblocks = rMesh.CellDataNumBlocks(rName);
    std::size_t total_rows = 0;
    std::vector<std::size_t> shape = rMesh.CellData(rName, 0).Shape();
    for (std::size_t b = 0; b < nblocks; ++b) {
        const auto& bshape = rMesh.CellData(rName, b).Shape();
        total_rows += bshape.empty() ? 0 : bshape[0];
    }
    shape[0] = total_rows;
    NDArray out(rMesh.CellData(rName, 0).Dtype(), shape);
    std::size_t off = 0;
    for (std::size_t b = 0; b < nblocks; ++b) {
        const NDArray& blk = rMesh.CellData(rName, b);
        std::memcpy(out.Data() + off, blk.Data(), blk.Nbytes());
        off += blk.Nbytes();
    }
    return out;
}

/**
 * @brief Splits a raw, whole-mesh cell-data array (as read from XDMF/HMF)
 * back into one `NDArray` per cell block, matching Python's
 * `cell_data_from_raw`.
 *
 * Inverse of `concat_cell_data`.
 * @param raw The concatenated array covering every cell block's rows,
 *            in cell-block order.
 * @param sizes Row count of each cell block, in the same order the blocks
 *              appear in `raw`; must sum to `raw`'s row count.
 * @return One `NDArray` per entry in `sizes`, each holding that block's slice.
 */
inline std::vector<NDArray> split_raw_cell_data(const NDArray& rRaw,
                                                const std::vector<std::size_t>& rSizes) {
    std::size_t ncols = rRaw.Ndim() >= 2 ? rRaw.Shape()[1] : 1;
    std::size_t off = 0;
    std::vector<NDArray> blocks;
    for (std::size_t bs : rSizes) {
        std::vector<std::size_t> bshape = rRaw.Shape();
        if (!bshape.empty())
            bshape[0] = bs;
        NDArray b(rRaw.Dtype(), bshape);
        std::size_t elems = bs * ncols;
        std::memcpy(b.Data(), rRaw.Data() + off * ncols * dtype_size(rRaw.Dtype()),
                    elems * dtype_size(rRaw.Dtype()));
        off += bs;
        blocks.push_back(std::move(b));
    }
    return blocks;
}

}  // namespace xdmfcommon
}  // namespace meshioplusplus
