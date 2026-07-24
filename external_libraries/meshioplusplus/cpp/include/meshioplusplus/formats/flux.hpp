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
 * @file flux.hpp
 * @brief Altair FLUX mesh (.pf3) C++ reader/writer.
 *
 * ASCII with French keyword headers (as handled by FEconv). Header lines
 * (`dim`, `nel`, `nnod`) are located by substring search for their French
 * label (e.g. `"NOMBRE DE DIMENSIONS"`) rather than fixed line position, so
 * header ordering is tolerant. The element block (after "DESCRIPTEUR DE
 * TOPOLOGIE") holds `nel` records as a continuous token stream: a 12-integer
 * header (field 3 = region reference -> `cell_data["pf3:ref"]`; field 6 =
 * `desc3`, the type code selecting the meshio++ type; field 7 = node count)
 * followed by 1-based connectivity. The coordinate block (after
 * "COORDONNEES DES NOEUDS") holds `nnod` rows of `node_index x1 ... x_dim`
 * (the leading index is discarded; rows assumed already in file order).
 * Unlike UNV/gmsh/mphtxt, **no node-order permutation** is applied — ids
 * pass through in file order directly. Hybrid (multi-type) meshes are
 * supported. See doc/formats/flux.md for the full `desc3` <-> meshio++ type
 * table and the meshio++ -> `(desc1, desc2, desc3)` reverse table.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a FLUX .pf3 file.
 *
 * Emits the dimension/element-count/node-count header, then per-element
 * 12-integer records via the fixed meshio++ -> `(desc1, desc2, desc3)`
 * table, `cell_data["pf3:ref"]` as each element's region reference (field
 * 3; defaults if absent), followed by 1-based connectivity, then the
 * "COORDONNEES DES NOEUDS" coordinate block. Several header fields are
 * always-placeholder: region counts are hardcoded `1 0 0 0 0 0`, and both
 * "max nodes per element" and "max integration points" fields are hardcoded
 * to `20` regardless of actual mesh content.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write (hybrid/multi-type meshes are supported)
 * @throws WriteError if a cell block's type has no entry in the meshio++ ->
 *         `desc3` table
 * @note reads/writes `cell_data["pf3:ref"]`
 */
void write_flux(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a FLUX .pf3 file.
 *
 * Locates the `dim`/`nel`/`nnod` header fields by French-label substring
 * search, then parses `nel` element records (12-int header + 1-based
 * connectivity, `desc3` selecting the meshio++ type) and `nnod` coordinate
 * rows, with no node-order permutation applied.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `cell_data["pf3:ref"]` set from each
 *         element's region-reference field
 * @throws ReadError if the file can't be opened, the element/coordinate
 *         section markers are missing, an element header is truncated, or
 *         an element's `desc3` type code is unrecognized
 * @note region *names* (which FLUX may store separately) are never read —
 *       only the numeric per-element reference in `cell_data["pf3:ref"]`
 */
Mesh read_flux(const std::string& rPath);

}  // namespace meshioplusplus
