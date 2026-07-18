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
 * @file ensight.hpp
 * @brief EnSight Gold (.case/.geo) C++ reader/writer — geometry only.
 *
 * EnSight Gold stores a dataset as a small `.case` index file plus a
 * geometry file (conventionally `.geo`). Only the mesh-geometry subset is
 * handled: the `.case` `FORMAT`/`GEOMETRY` sections (the file must declare
 * `type: ensight gold`; `VARIABLE`/`TIME` sections are ignored) and the Gold
 * geometry file in both ASCII and C-binary form (the leading `"C Binary"`
 * 80-char record selects binary; `"Fortran Binary"` is rejected). Binary
 * files use 32-bit ints/floats in the writing machine's byte order; the
 * reader auto-detects a foreign byte order from the plausibility of the
 * part-number/node-count records and byte-swaps accordingly.
 *
 * Element keywords `point`, `bar2/3`, `tria3/6`, `quad4/8`, `tetra4/10`,
 * `pyramid5/13`, `penta6/15`, `hexa8/20` map to the corresponding meshio
 * cell types; ragged `nsided`/`nfaced` sections are read into polygon /
 * polyhedron blocks (grouped by node count, the openfoam convention). Node
 * ordering matches meshio for every type except `penta15`, which differs
 * from meshio's `wedge15` by the involution
 * `{0,2,1,3,5,4, 8,7,6, 11,10,9, 12,14,13}` (the same map VTK's EnSight
 * readers apply). Per the Gold specification connectivity is **positional**
 * (1-based index into the part's coordinate list); `node id given/ignore`
 * id arrays are present in the file but skipped.
 *
 * Multi-part files are concatenated into one point array; every per-part
 * element section becomes its own cell block, and when the file has two or
 * more parts the owning part number is recorded as the integer cell_data
 * field `"ensight:part"`. The writer emits a single part (`node id assign`,
 * `element id assign`) and drops point/cell/field data (mesh-only scope).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as an EnSight Gold `.case` + `.geo` sibling pair.
 *
 * `rPath` may name either sibling (`.case` or `.geo`); the other is derived
 * from the shared stem and both files are written. The geometry is a single
 * part (id 1, description "Mesh") with `node id assign` / `element id
 * assign` numbering; coordinates are always emitted with three components
 * (z padded with 0 for 2D meshes). ASCII floats use `%12.5e` (EnSight caps
 * usable precision at 6 significant digits); binary emits 32-bit
 * ints/floats in host byte order with 80-char string records.
 *
 * @param rPath filesystem path to either the `.case` or `.geo` sibling
 * @param rMesh the mesh to write
 * @param binary true for C-binary geometry, false for ASCII
 * @throws WriteError if a file cannot be opened, the mesh contains a cell
 *         type without an EnSight keyword, a ragged (polygon/polyhedron)
 *         block is present (not written in v1), points have more than three
 *         components, or a binary mesh exceeds 32-bit counts
 * @note point/cell/field data are not written (geometry-only scope).
 */
void write_ensight(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read an EnSight Gold `.case` file (or a Gold geometry file directly).
 *
 * A `.case` path is parsed for its `GEOMETRY`/`model:` entry (resolved
 * relative to the case file's directory; transient wildcard names are
 * rejected); any other path is treated as a Gold geometry file. ASCII and
 * C-binary geometries are both handled, including foreign-endian binaries.
 *
 * @param rPath filesystem path to a `.case` file or a Gold `.geo` file
 * @return the read Mesh (parts concatenated; one cell block per per-part
 *         element section; `"ensight:part"` cell_data when the file has
 *         two or more parts)
 * @throws ReadError on malformed input, a non-Gold case file, Fortran-binary
 *         geometry, an unknown element keyword, or out-of-range connectivity
 */
Mesh read_ensight(const std::string& rPath);

}  // namespace meshioplusplus
