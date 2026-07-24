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
 * @file ansysinp.hpp
 * @brief Ansys MAPDL "coded database" (.cdb / .inp) C++ reader/writer.
 *
 * An autonomous format distinct from the unrelated Fluent `.msh` format in
 * ansys.hpp (both are named "ansys" in meshio++). Mirrors
 * `src/meshioplusplus/ansysInp/_ansysInp.py`. Parses whitespace/keyword-
 * delimited MAPDL command blocks directly: `ET`/`ETBLOCK` (element-type
 * declarations), `NBLOCK` (fixed-width node rows, field widths parsed from
 * the format-spec line such as `(3i9,6e20.13)` rather than hardcoded),
 * `EBLOCK` (element rows `(mat, type, real, secnum, esys, birth, death,
 * solkey, nodes_per_elem, ..., elem_id, node_ids...)`, with a continuation
 * line when there are more than 8 node ids), and `CMBLOCK` (named
 * components: `NODE` -> point set, `ELEM*` -> cell set, with negative
 * values expanding a range `-k` after base `b` into `range(b+1, k+1)`).
 *
 * Ansys element type ids group into 4 families (`solid`, `shell`, `plane`,
 * `line`) and combine with the actual node count read to resolve a meshio++
 * type (e.g. (solid, 10) -> `tetra10`, (shell/plane, 8) -> `quad8`); see
 * doc/formats/ansysinp.md for the full family/(family,nodes) tables and the
 * fixed meshio++ -> Ansys-type-id reverse map used on write (one id per
 * meshio++ type, e.g. `tetra10->187`, regardless of the id the file was
 * originally read with).
 *
 * `CMBLOCK` point/cell sets are custom attributes on the Python `Mesh`, not
 * carried by the Mesh conversion layer, so they travel out-of-band through
 * the @ref AnsysInfo side-channel struct (the same pattern as `MedInfo`) that
 * the binding layer `setattr`s onto the Python Mesh as `point_sets`/
 * `cell_sets`.
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel carrying CMBLOCK-derived point/cell sets across the
 *        Mesh conversion boundary (the `point_sets`/`cell_sets` Python Mesh
 *        attributes are not part of the C++ Mesh/NDArray conversion layer).
 */
struct AnsysInfo {
    /** Component name -> node indices (0-based), from `CMBLOCK ...,NODE`. */
    std::map<std::string, std::vector<std::int64_t>> mPointSets;
    /**
     * Component name -> per-cell-block lists of local cell indices
     * (0-based), one inner list per mesh cell block in block order (the
     * order blocks were first encountered while reading `EBLOCK`), from
     * `CMBLOCK ...,ELEM`.
     */
    std::map<std::string, std::vector<std::vector<std::int64_t>>> mCellSets;
};

/**
 * @brief Read an Ansys MAPDL coded-database (.cdb/.inp) file.
 *
 * Parses `ET`/`ETBLOCK` element-type declarations, `NBLOCK` node rows,
 * `EBLOCK` element rows (resolving each row's meshio++ type from the
 * element's family + node count), and `CMBLOCK` named components. A line
 * matching the exclusion list (known keywords, `KEYWORD,` syntax, `!`/`/`
 * comments) stops a block's row-reading loop early.
 *
 * @param rPath filesystem path to read
 * @param[out] rInfo receives `CMBLOCK` point/cell sets (0-based indices),
 *        keyed by component name
 * @return the read Mesh (no point_data/cell_data/field_data — only
 *         geometry, connectivity, and named sets are represented)
 * @throws ReadError if the file can't be opened, no `NBLOCK`/`EBLOCK`/
 *         `CMBLOCK` is found at all, a `CMBLOCK` negative range value
 *         appears before any base value, or an `EBLOCK` row's (family,
 *         node-count) pair has no meshio++ type mapping
 */
Mesh read_ansysinp(const std::string& rPath, AnsysInfo& rInfo);

/**
 * @brief Write `mesh` (plus `info`'s named sets) as an Ansys MAPDL
 *        coded-database file.
 *
 * Always emits exactly one `NBLOCK`/`EBLOCK` pair with fixed `i9`/`e20.13`
 * field widths (not preserving an original file's exact layout or element
 * type ids) — a read-write round trip is semantically but not byte-
 * identical. 2D input meshes are padded to 3D with a zero z-column (MAPDL
 * has no native 2D coordinate concept). Each meshio++ cell type is written
 * with the fixed reverse element-type-id map from the file-level doc
 * comment.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param rInfo point/cell sets to emit as `CMBLOCK` components
 * @throws WriteError if a cell block's meshio++ type has no entry in the
 *         reverse element-type map ("Unhandled meshio type")
 * @note point_sets/cell_sets travel via `info`, not via `mesh` — the Python
 *       binding setattrs them onto/from the Mesh object separately
 */
void write_ansysinp(const std::string& rPath, const Mesh& rMesh, const AnsysInfo& rInfo);

}  // namespace meshioplusplus
