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
 * @file flac3d.hpp
 * @brief Itasca FLAC3D grid (.f3grid) C++ reader/writer — common path
 *        (ASCII + binary), excluding cell groups.
 *
 * Format is auto-detected by checking the first 8 bytes for a null byte
 * (binary if found, else ASCII). Binary (little-endian): an 8-byte header
 * pair of undocumented meaning on read but reproduced verbatim on write
 * (`1375135718, 3`), then `uint32` node count and per-node `(point_id:
 * uint32, x,y,z: float64x3)`, then for `zone` and `face` in that order:
 * `uint32` cell count and per-cell `(cell_id: uint32, num_verts: uint32,
 * node_ids: uint32 x num_verts)` — `num_verts == 7` is a degenerate "B7"
 * hexahedron-as-7-node encoding, handled by duplicating the last node to
 * make 8. ASCII mirrors the same structure with `G`/`Z`/`F` record lines.
 * Cells are grouped into blocks of **consecutive same-typed cells in file
 * order** (not merged globally by type).
 *
 * The format's central quirk is the **right-handed zone reorder**: FLAC3D
 * requires each zone's first four corner nodes to form a right-handed
 * system, so on write the C++ core computes the scalar triple product of
 * the first three edge vectors (from the primary meshio++->FLAC3D node
 * order) and picks the primary order if positive, else a pre-tabulated
 * "flipped" alternate order (only `tetra`/`pyramid`/`wedge`/`hexahedron`
 * have a flipped variant; `triangle`/`quad` do not need one). This
 * determinant check only happens on write — the read-side reorder is a
 * fixed, unconditional permutation, assuming a well-formed file already
 * stores correctly-handed zones. `ZGROUP`/`FGROUP` cell-group sections are
 * always deferred to the Python fallback (see @ref read_flac3d). See
 * doc/formats/flac3d.md for the full node-order tables.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a FLAC3D .f3grid file (ASCII or binary), gridpoints
 *        and zone/face cells only.
 *
 * Applies the meshio++ -> FLAC3D node-order table per cell type, choosing
 * between the primary and flipped order per zone based on the sign of the
 * first-three-edge-vectors scalar triple product (the right-handed
 * reorder). Emits `* ZONES` before `* FACES` in the ASCII layout.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param rFloatFmt coordinate format string (ASCII only; ignored for
 *        binary)
 * @param binary write the binary FLAC3D layout (`true`) or ASCII (`false`)
 * @throws WriteError if a file cannot be opened for writing
 * @note the shim only attempts this C++ path when `mesh.cell_sets` is
 *       empty — `ZGROUP`/`FGROUP` are always written by the Python fallback,
 *       which also hardcodes group slots (`SLOT 1` ASCII / `"Default"`
 *       binary) rather than preserving an original slot name
 */
void write_flac3d(const std::string& rPath, const Mesh& rMesh, const std::string& rFloatFmt,
                  bool binary);

/**
 * @brief Read a FLAC3D .f3grid file (ASCII or binary), gridpoints and
 *        zone/face cells only.
 *
 * Detects ASCII vs. binary from the first 8 bytes, then parses gridpoints
 * and, in faces-then-zones internal order (a structural asymmetry relative
 * to the writer's ZONES-then-FACES section order — harmless since the
 * reader doesn't depend on write order), zone/face cell blocks, applying the
 * fixed FLAC3D -> meshio++ node-order permutation per type and expanding
 * degenerate 7-node "B7" hexahedra to 8 nodes.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `cell_data["cell_ids"]` set to each cell's
 *         original FLAC3D global id (split per block, faces numbered before
 *         zones)
 * @throws ReadError if the file can't be opened, the file ends
 *         unexpectedly, a cell's node count doesn't match any known FLAC3D
 *         type, or the file contains a `ZGROUP`/`FGROUP` section (ASCII) or
 *         binary group section — always deferring the whole file to the
 *         Python fallback, since `cell_sets` (built from those groups) is
 *         not carried by the Mesh conversion layer
 * @note point_data/field_data are never produced; `cell_data["cell_ids"]` is
 *       the only key this reader sets
 */
Mesh read_flac3d(const std::string& rPath);

}  // namespace meshioplusplus
