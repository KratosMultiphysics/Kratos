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
 * @file gmsh.hpp
 * @brief Gmsh mesh format (.msh, versions 2.2 and 4.1) C++ reader/writer.
 *
 * `$MeshFormat` (`version filetype datasize`; `filetype` 0=ascii, 1=binary,
 * with a 4-byte endianness-detection integer `1` for binary) is read first
 * and picks the reader: `"2"`/`"2.2"` -> the 2.2 path, `"4"`/`"4.1"` -> the
 * 4.1 path. The C++ reader (@ref read_gmsh) handles **versions 2.2 and 4.1
 * only** — version 4.0 (which needs `$Entities`) and `$Periodic` records
 * always throw and defer to the Python reader (see
 * doc/formats/gmsh.md#quirks-limitations).
 *
 * **Version 2.2**: `$PhysicalNames`; `$Nodes` (ascii `id x y z` rows, or
 * binary `(int32 id, 3xdouble)`); `$Elements` (ascii `id type ntags
 * tag1..tagN node1..nodeK` per line, binary per-block `elem_type num_elems
 * num_tags` header then flat int32 rows). The first two element tags become
 * `cell_data["gmsh:physical"]`/`cell_data["gmsh:geometrical"]`.
 *
 * **Version 4.1** restructures node/element blocks: `$Nodes` header
 * `numEntityBlocks numNodes minNodeTag maxNodeTag`, per block `entityDim
 * entityTag parametric numNodesInBlock` followed by a node-tag list then a
 * matching coordinate list (tags may be sparse/out of order, requiring a
 * tag->index remap); `$Elements` likewise groups per entity block, with rows
 * of `elementTag node1..nodeK`. `point_data["gmsh:dim_tags"]` (an `(N,2)`
 * `(entity_dim, entity_tag)` array) and `cell_sets["gmsh:bounding_entities"]`
 * are v4.1-only concepts.
 *
 * Five element types need a node-order permutation between Gmsh and
 * meshio++ (`tetra10`, `hexahedron20`, `hexahedron27`, `wedge15`,
 * `pyramid13` — see doc/formats/gmsh.md for the exact permutation arrays);
 * everything else uses natural order. The C++ type table covers a curated
 * subset up through roughly `hexahedron125`/`tetra286` — not the full
 * ~110-entry Python table — so a file referencing a higher-order type
 * outside that subset falls back to Python transparently. `field_data` maps
 * from `$PhysicalNames` as `[phys_num, phys_dim]`; `mesh.gmsh_periodic` (a
 * mesh-level attribute, not a data-dict key) is only ever populated by the
 * Python reader.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` to `path` as a Gmsh 2.2 .msh file (ascii or binary).
 *
 * Emits `$MeshFormat` (version "2.2"), `$PhysicalNames` (from
 * `field_data`), `$Nodes`, and `$Elements` with `gmsh:physical`/
 * `gmsh:geometrical` as the first two element tags. Applies the gmsh <->
 * meshio++ node-order permutation for `tetra10`/`hexahedron20`/
 * `hexahedron27`/`wedge15`/`pyramid13`.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write node/element bodies as binary (`true`, with the
 *        endianness-detection integer) or ASCII (`false`)
 * @throws WriteError if a cell block's type has no Gmsh type-code mapping
 * @note reads/writes `cell_data["gmsh:physical"]`/`cell_data["gmsh:geometrical"]`
 *       and `field_data` (as `$PhysicalNames`)
 * @note the shim only attempts this C++ path when `float_fmt == ".16e"` and
 *       `mesh.gmsh_periodic` is unset
 */
void write_gmsh22(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Write `mesh` to `path` as a Gmsh 4.1 .msh file (ascii or binary).
 *
 * Intended for meshes without entity information (no
 * `point_data["gmsh:dim_tags"]`); `$Entities` is not emitted, so more than
 * one cell type cannot be written this way (Gmsh 4.1 requires `$Entities`
 * to disambiguate cell-to-entity assignment for mixed meshes).
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write node/element bodies as binary (`true`) or ASCII
 *        (`false`)
 * @throws WriteError if `mesh` has more than one cell type (since
 *         `$Entities` is never emitted here) or a cell block's type has no
 *         Gmsh type-code mapping
 * @note the shim only attempts this C++ path when `float_fmt == ".16e"`, no
 *       `gmsh_periodic`, and no `gmsh:dim_tags` in `point_data` — any v4.1
 *       write carrying `gmsh:dim_tags` or periodic data always goes through
 *       Python instead
 */
void write_gmsh41(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read a Gmsh .msh file (versions 2.2 and 4.1 only).
 *
 * Dispatches on the `$MeshFormat` version string; parses `$PhysicalNames`,
 * `$Nodes`, `$Elements` (applying the gmsh <-> meshio++ node-order
 * permutation where needed), and, for 4.1, `$Entities`/per-entity node and
 * element blocks with node-tag->index remapping.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `cell_data["gmsh:physical"]`/
 *         `cell_data["gmsh:geometrical"]` from the first two element tags,
 *         `point_data["gmsh:dim_tags"]` and `cell_sets["gmsh:bounding_entities"]`
 *         (v4.1 only), and `field_data` from `$PhysicalNames`
 * @throws ReadError for anything not handled by the C++ path — version not
 *         2.2/4.1 (e.g. 4.0, which needs `$Entities`), `$Periodic` records,
 *         a Gmsh element type outside the curated type-code subset, or
 *         parametric nodes — so the Python reader can take over
 * @note the C++ reader never populates `mesh.gmsh_periodic`; only the
 *       Python fallback does, for files containing `$Periodic`
 */
Mesh read_gmsh(const std::string& rPath);

}  // namespace meshioplusplus
