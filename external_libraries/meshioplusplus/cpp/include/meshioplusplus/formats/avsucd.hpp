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
 * @file avsucd.hpp
 * @brief AVS-UCD (.avs) ASCII C++ reader/writer.
 *
 * AVS Unstructured Cell Data: `#`-comment lines, a header line of 5 integers
 * (`num_nodes num_cells num_node_data num_cell_data 0`), then `num_nodes`
 * rows of `id x y z` (id is an **arbitrary integer**, not necessarily
 * sequential or 1-based — both read and write maintain explicit id<->index
 * maps), `num_cells` rows of `id material_id avsucd_type_name node_id0
 * node_id1 ...` (node count per row is simply "everything after the first 3
 * fields", no fixed per-type table on read), and, when the corresponding
 * header counts are nonzero, node-data / cell-data sections (a component-
 * count header line, `", real"`-suffixed label lines, then per-entity data
 * rows resolved through the same id map).
 *
 * Node types map through the `pt`/`line`/`tri`/`quad`/`tet`/`pyr`/`prism`/
 * `hex` <-> `vertex`/`line`/`triangle`/`quad`/`tetra`/`pyramid`/`wedge`/
 * `hexahedron` table with fixed node-order permutations for `tetra`
 * (`[0,1,3,2]`), `pyramid` (`[4,0,1,2,3]`), `wedge` (`[3,4,5,0,1,2]`), and
 * `hexahedron` (`[4,5,6,7,0,1,2,3]`) on write; the read-side inverse is the
 * same table for the involutions (tetra/wedge/hexahedron) but a distinct
 * `[1,2,3,4,0]` for pyramid. See doc/formats/avsucd.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as an AVS-UCD .avs file.
 *
 * Renumbers all node and cell ids sequentially from 1, regardless of any
 * original ids. `cell_data["avsucd:material"]` (the first integer-typed
 * cell_data array found, if any; others are silently dropped in the C++
 * writer) is written as each cell row's material id; other cell_data/
 * point_data arrays become additional labeled data sections. 2D points are
 * promoted to 3D. Floats use `%.17g` for points and `%.14e` for data.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if a cell block's type has no AVS-UCD type-name mapping
 * @note reads/writes `cell_data["avsucd:material"]`; other point_data/
 *       cell_data names pass through as-is (post-strip(), spaces replaced
 *       with underscores — not reversible)
 */
void write_avsucd(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an AVS-UCD .avs file.
 *
 * Builds an id->index map while reading nodes and cells so that arbitrary,
 * sparse, or non-contiguous file ids resolve correctly; applies the AVS-UCD
 * -> meshio++ node-order permutation per cell type; splits any multi-block
 * cell_data array back into per-block pieces using cumulative block-length
 * offsets (assumes blocks are contiguous in read order).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `cell_data["avsucd:material"]` set from each
 *         cell row's material id
 * @throws ReadError if the file can't be opened or a cell row names an
 *         unknown AVS-UCD type
 */
Mesh read_avsucd(const std::string& rPath);

}  // namespace meshioplusplus
