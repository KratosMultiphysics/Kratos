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
 * @file netgen.hpp
 * @brief Netgen neutral mesh (.vol) C++ reader/writer — common path only.
 *
 * A `.vol` file starts with a literal `mesh3d` line, then keyword blocks
 * in arbitrary order (`dimension`, `geomtype`, `points`, `pointelements`,
 * `edgesegments`/`edgesegmentsgi`, `surfaceelements*`, `volumeelements`),
 * each `<keyword>\n<count>\n<count data rows>`; blank/`#`-comment lines are
 * skipped anywhere. The element type per row is inferred from a
 * fixed-column node count (variable per row: e.g. `surfaceelements`
 * triangle=3/quad=4/triangle6=6/quad8=8; `volumeelements`
 * tetra=4/pyramid=5/wedge=6/hexahedron=8/tetra10=10/pyramid13=13/
 * wedge15=15/hexahedron20=20). A single per-element region/material index
 * is stored across every block -> `cell_data["netgen:index"]`.
 *
 * Indices are **1-based** on disk (`+1`/`-1` applied). The full
 * Netgen->meshio++ node permutation table (`meshio++[i] = netgen[table[i]]`,
 * meshio++->Netgen uses the exact per-entry inverse) covers `triangle6`
 * `[0,1,2,5,3,4]`, `quad8` `[0,1,2,3,4,7,5,6]`, `tetra` `[0,2,1,3]`,
 * `tetra10` `[0,2,1,3,5,7,4,6,9,8]`, `pyramid` `[0,3,2,1,4]`, `pyramid13`
 * `[0,3,2,1,4,7,6,8,5,9,12,11,10]`, `wedge` `[0,2,1,3,5,4]`, `wedge15`
 * `[0,2,1,3,5,4,7,8,6,13,14,12,9,11,10]`, `hexahedron` `[0,3,2,1,4,7,6,5]`,
 * `hexahedron20`
 * `[0,3,2,1,4,7,6,5,10,9,11,8,16,19,18,17,14,13,15,12]`
 * (`line`/`triangle`/`quad`/`vertex` use natural order).
 *
 * **Deferred to Python** (the reader throws when it meets any of these
 * tokens, and the writer is gated off by the shim when the mesh carries
 * the corresponding data): the `identifications`/`identificationtypes`
 * periodic node-pair tables (stored in `mesh.info`, which has no C++-core
 * representation), `materials`/`bcnames`/`cd2names`/`cd3names` codimension
 * name tables (-> non-empty `field_data`), the two-physical-line
 * `edgesegmentsgi2` variant, `face_colours`/`singular_*` sections, and the
 * gzip `.vol.gz` container (the C++ reader/writer explicitly refuse the
 * `.gz` suffix; Python handles it via `gzip.open`).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a Netgen neutral mesh (.vol) file, ascii,
 *        common-path only.
 *
 * Emits `mesh3d`, `dimension`, `points`, and per-dimension element blocks
 * (`pointelements`/edge/`surfaceelements`/`volumeelements` as applicable)
 * with 1-based connectivity after applying the meshio++->Netgen node
 * permutation. The single per-cell region/material marker is taken from
 * `cell_data["netgen:index"]` if present, else the first integer-dtype
 * cell_data array found (Netgen has no way to store the array's name).
 * Refuses (via the shim) meshes carrying `mesh.info` entries or non-empty
 * `field_data`, and never handles the `.vol.gz` suffix.
 *
 * @param rPath filesystem path to the .vol file to create/overwrite
 * @param rMesh the mesh to write
 * @param rFloatFmt coordinate format string (e.g. `".16e"`)
 * @throws WriteError on an unsupported cell type, mixed content this path
 *         doesn't implement, or a `.gz` path
 * @note reads `cell_data["netgen:index"]` if present
 */
void write_netgen(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt);

/**
 * @brief Read a Netgen neutral mesh (.vol) file into a Mesh, ascii,
 *        common-path only.
 *
 * Parses `dimension`, `geomtype` (unexpected values only warn),
 * `points`, and the point/edge/surface/volume element blocks, inferring
 * each row's cell type from its node count at a fixed column position and
 * converting 1-based indices to 0-based via the inverse Netgen->meshio++
 * node permutation. The single per-element region marker becomes
 * `cell_data["netgen:index"]`.
 *
 * @param rPath filesystem path to the .vol file to read
 * @return the read Mesh, with `cell_data["netgen:index"]` populated
 * @throws ReadError on `identifications`/`identificationtypes`,
 *         `materials`/`bcnames`/`cd2names`/`cd3names`, the two-line
 *         `edgesegmentsgi2` variant, `face_colours`/`singular_*` sections,
 *         a `.gz` path, or a malformed file — all of which route to the
 *         Python fallback
 */
Mesh read_netgen(const std::string& rPath);

}  // namespace meshioplusplus
