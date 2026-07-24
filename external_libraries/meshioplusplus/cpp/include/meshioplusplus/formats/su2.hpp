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
 * @file su2.hpp
 * @brief SU2 (.su2) ascii mesh C++ reader/writer.
 *
 * Line-oriented `KEY= value` records (`%` starts a comment; blank/malformed
 * lines are skipped with a warning). `NDIME= 2|3` sets the dimension;
 * `NPOIN= n` is followed by `n` coordinate rows (auto-detecting and
 * stripping an optional trailing global-index column some SU2 files add);
 * `NELEM= n` (volume cells, tagged `su2:tag = 0`) and `MARKER_ELEMS= n`
 * (boundary cells under a preceding `MARKER_TAG=`, tagged with that marker's
 * id) rows are parsed as one integer block and binned by VTK-style numeric
 * type code into cell blocks (mixed types in one section split into
 * separate blocks). Boundary blocks of the same cell type from separate
 * `MARKER_ELEMS` sections are merged into one block per type after the full
 * scan. `NMARK= n` is only soft-checked against the actual marker count
 * (mismatch just warns).
 *
 * Cell type codes: 3=line(2), 5=triangle(3), 9=quad(4), 10=tetra(4),
 * 12=hexahedron(8), 13=wedge(6), 14=pyramid(5) — `(nodes)` per cell.
 *
 * @note cell_data key produced/consumed: `"su2:tag"` (volume cells always 0;
 *       boundary cells get their marker's tag id, auto-incremented from 1
 *       for non-numeric string tags).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh to an SU2 file.
 *
 * Emits `NDIME=` from `points.shape[1]`, `NPOIN=` + coordinates, volume
 * cells (triangle/quad in 2D, tetra/hexahedron/wedge/pyramid in 3D) under
 * `NELEM=` each prefixed by its numeric type code, and boundary markers
 * grouped by the first integer-typed `cell_data` array found (the
 * "first-int-array" convention shared with several other formats), emitting
 * one `MARKER_TAG=`/`MARKER_ELEMS=` pair per distinct tag value. Unsupported
 * cell types warn and are skipped. Only one integer cell_data array can
 * drive the markers; additional candidates are dropped with a warning.
 *
 * @param path filesystem path to write
 * @param mesh the mesh to write
 * @throws WriteError on an unopenable output path or an unwritable geometry
 * @note reads `cell_data["su2:tag"]` to build boundary markers.
 */
void write_su2(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an SU2 mesh file.
 *
 * Parses `NDIME`/`NPOIN`/`NELEM`/`NMARK`/`MARKER_TAG`/`MARKER_ELEMS` records
 * per the grammar above, reconstructing volume and boundary cell blocks with
 * `su2:tag` cell_data.
 *
 * @param path filesystem path to read
 * @return the read Mesh
 * @throws ReadError on a malformed record (e.g. an element row that doesn't
 *         match a known VTK-style type code)
 * @note cell_data key produced: `"su2:tag"`.
 */
Mesh read_su2(const std::string& rPath);

}  // namespace meshioplusplus
