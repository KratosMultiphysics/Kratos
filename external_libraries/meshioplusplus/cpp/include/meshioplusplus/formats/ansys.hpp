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
 * @file ansys.hpp
 * @brief Ansys/Fluent mesh (.msh) C++ reader/writer.
 *
 * Not to be confused with the unrelated Ansys MAPDL "coded database" format
 * handled by ansysinp.hpp. This is the Fluent `.msh` format: fully
 * parenthesis-nested "Scheme-like" sections `(<index> ...)`, where the index
 * may be a bare decimal (ASCII payload) or prefixed `20`/`30` for a binary
 * payload (`20xx` = float32 nodes / int32 cells, `30xx` = float64 / int64).
 * All connectivity and zone-header integers in both ASCII and binary bodies
 * are **hexadecimal** — the format's defining quirk. Section `10` gives node
 * blocks (`zone-id first last type ND`), section `12` gives cell blocks
 * (`zone-id first last zone-type element-type`; `zone-type == 0` is a dead
 * zone producing no cells; `element-type == 0` is a "mixed" zone that is
 * structurally skipped, Fluent's own heterogeneous-cell encoding being
 * unresolved here), and section `13` gives boundary faces. All zones are
 * folded into one flat cell list, then every connectivity array has the
 * first point-zone's `first` index subtracted so numbering normalizes to 0.
 * `point_data`/`cell_data`/`field_data` are always empty for this format —
 * it carries geometry and zone/boundary structure only.
 *
 * See doc/formats/ansys.md for the full section grammar and the
 * element-type/face-type code tables.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a Fluent .msh file.
 *
 * Emits, in order: a `(1 "...")` header, `(2 DIM)`, a `(10 (0 1 N 0))` node-
 * count declaration, a `(12 (0 1 N 0))` cell-count declaration, one node
 * block (`10` ascii or `3010` binary), then one cell block per meshio++ cell
 * type using the fixed reverse map `triangle:1, tetra:2, quad:3,
 * hexahedron:4, pyramid:5, wedge:6` (ascii section `12`, binary `2012`
 * int32 or `3012` int64). No face (`13`) sections, and no `mixed`/polyhedral
 * cell support, are ever emitted.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary write node/cell bodies as binary (`true`, `20xx`/`30xx`
 *        prefixed sections) or ASCII (`false`)
 * @throws WriteError if `mesh` is not 2D or 3D, or if a cell block's type
 *         has no entry in the meshio++ -> Ansys type-code map ("illegal
 *         cell type")
 */
void write_ansys(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read a Fluent .msh file.
 *
 * Parses `(0 ...)`/`(1 ...)`/`(2 ...)` header/comment/dimension sections
 * (bracket-skipped), `(<pfx>10 ...)` node sections (ascii one point per
 * line, binary a raw float32/float64 block), and `(<pfx>12 ...)` cell
 * sections (dead zones -> no cells; `mixed` zones structurally skipped, body
 * not decoded). All hexadecimal header/body integers are converted; the
 * result is one flat `Mesh.mCells` list with the first point-zone's `first`
 * index subtracted from every connectivity array.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (point_data/cell_data/field_data always empty)
 * @throws ReadError if the file can't be opened, a section header is
 *         malformed or truncated, a cell zone's `element-type` isn't one of
 *         the known volume codes (0/1/.../6), or a face (`13`) section
 *         carries a data body — **any** real face section always defers the
 *         whole file to the Python fallback, so files with real boundary
 *         face zones (a common real-world case) are never handled here;
 *         binary "mixed" faces additionally raise unconditionally even in
 *         the Python path
 * @note point_data/cell_data/field_data are never produced — this format
 *       carries no per-node/per-cell field values, only geometry and zone
 *       structure
 */
Mesh read_ansys(const std::string& rPath);

}  // namespace meshioplusplus
