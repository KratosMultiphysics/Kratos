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
 * @file tetgen.hpp
 * @brief TetGen (.node/.ele) C++ reader/writer — a shared-stem file pair.
 *
 * TetGen stores a mesh as two sibling files sharing a stem: `<stem>.node`
 * (header `npoints dim nattrs nbmarkers`, `dim` must be 3, rows
 * `idx x y z attr1..attrN marker1..markerM`) and `<stem>.ele` (header
 * `ntets 4 nattrs`, rows `idx n0 n1 n2 n3 attr1..attrK`). Either path
 * (`.node` or `.ele`) selects the pair. The `.node` file's node index base
 * (0 or 1) is auto-detected from its first row and all indices must then be
 * **exactly consecutive** from that base (ReadError on any gap); the `.ele`
 * connectivity is shifted by that same detected base so files using either
 * numbering read correctly. `tetra` is the only representable cell type —
 * TetGen only ever describes tetrahedra.
 *
 * On write, attribute/marker keys are partitioned into at most one "ref" key
 * (the first key containing the substring `:ref`, or else the first key
 * present) plus the remaining plain attributes; the C++ writer special-cases
 * exact-integer ref values to print as plain integers (falling back to
 * `%.16e` otherwise), which can format float-valued refs slightly
 * differently than the Python writer's plain `str()` formatting. Each
 * `tetra` cell block written to `.ele` restarts its element-id counter at 0
 * — a mesh with multiple `tetra` blocks would produce duplicate element ids
 * (a genuine round-trip risk, though TetGen conventionally emits exactly one
 * block). The format cannot be read from or written to an in-memory buffer.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a TetGen `<stem>.node` / `<stem>.ele` file pair.
 *
 * `path` may be either sibling path; the stem is derived and both files are
 * written. Point attribute/marker columns come from point_data (first
 * `:ref`-containing key floats to the front as the boundary-marker column,
 * `%.16e` or integer formatting per value); cell attribute/ref columns come
 * from cell_data the same way. Only `tetra` cells are written; each tetra
 * block gets its own `.ele` numbering restarting at 0.
 *
 * @param rPath filesystem path to either the `.node` or `.ele` sibling
 * @param rMesh the mesh to write (must contain only `tetra` cells)
 * @throws WriteError if either output file cannot be opened, or the mesh
 *         contains non-tetra cells
 * @note point_data keys produced: `"tetgen:attr{k}"`, `"tetgen:ref"`,
 *       `"tetgen:ref2"`, ...; cell_data keys: `"tetgen:ref"`, `"tetgen:ref2"`, ...
 */
void write_tetgen(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a TetGen `.node`/`.ele` file pair.
 *
 * `path` may name either sibling; the other is derived from the shared
 * stem. Detects the node index base from the `.node` file's first row and
 * requires strictly consecutive indices; `.ele` connectivity is rebased by
 * the same amount.
 *
 * @param rPath filesystem path to either the `.node` or `.ele` sibling
 * @return the read Mesh (a single `tetra` CellBlock)
 * @throws ReadError if the sibling file is missing, `dim != 3`, or node
 *         indices are non-consecutive from the detected base
 * @note point_data keys produced: `"tetgen:attr{k}"` (node attribute
 *       columns) and `"tetgen:ref"`/`"tetgen:ref2"`/... (boundary marker
 *       columns); cell_data key: `"tetgen:ref"`/... (region attribute
 *       columns, one array per column since TetGen has one cell block).
 */
Mesh read_tetgen(const std::string& rPath);

}  // namespace meshioplusplus
