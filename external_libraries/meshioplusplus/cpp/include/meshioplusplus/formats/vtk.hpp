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
 * @file vtk.hpp
 * @brief Legacy VTK (.vtk) `UNSTRUCTURED_GRID` C++ reader/writer, versions
 *        4.2 and 5.1, ascii and binary.
 *
 * Only `DATASET UNSTRUCTURED_GRID` is handled by the C++ core; any other
 * dataset type (`STRUCTURED_POINTS`, `STRUCTURED_GRID`, `RECTILINEAR_GRID`)
 * always falls back to Python, which converts those into unstructured
 * line/quad/hex cells in Fortran (column-major) order. Binary numeric data
 * is **always big-endian on disk** regardless of host platform — an
 * explicit VTK-wiki convention, not a meshio++ choice — so binary I/O
 * byte-swaps through `detail/byteswap.hpp` intrinsics whenever the host is
 * little-endian, and always builds one pre-sized buffer for a single
 * `os.write`/bulk read rather than per-element stream operations.
 *
 * The two on-disk `CELLS` layouts differ completely between versions:
 * - **4.2**: interleaved — `CELLS <num_cells> <total_ints>` then, per cell,
 *   `[n, p0, ..., p_{n-1}]`, followed by a separate `CELL_TYPES <num_cells>`
 *   section. Reconstructed with a list-based per-block append.
 * - **5.1** (no official published spec; reverse-engineered from real files
 *   and a ParaView forum thread): `CELLS <num_offsets> <num_conn_items>` /
 *   `OFFSETS <dtype>` / offsets array (first entry 0, last equals
 *   `len(connectivity)`) / `CONNECTIVITY <dtype>` / flat connectivity array,
 *   `<dtype>` a literal token like `vtktypeint64`. Reconstructed via the
 *   shared offset-diff/vectorized helper in `detail/vtk_cells.hpp` (also
 *   used by the VTU reader) — zero-copy-friendly block reconstruction when
 *   the connectivity is contiguous and node order is identity.
 *
 * Cell types are shared with VTU (see vtu.hpp / doc/formats/vtk.md); only
 * `wedge` needs a node-order permutation relative to VTK (`[0,2,1,3,5,4]`,
 * self-inverse) — every other type uses natural order.
 *
 * `_cpp_ok(mesh)` (Python-side gate, not in this header) skips the C++ write
 * path for meshes with polyhedron cells or any 2-component vector data,
 * because the Python writer pads 2-component vectors to 3 components
 * **in place** on the input mesh, a mutation the C++ writer deliberately
 * does not replicate.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a VTK legacy `UNSTRUCTURED_GRID` file.
 *
 * Emits the version header line (`# vtk DataFile Version 4.2` or `5.1`),
 * `ASCII`/`BINARY`, `DATASET UNSTRUCTURED_GRID`, the `POINTS` block, the
 * `CELLS`/`CELL_TYPES` (4.2) or `CELLS`/`OFFSETS`/`CONNECTIVITY` (5.1)
 * blocks, then `POINT_DATA`/`CELL_DATA` `SCALARS`/`VECTORS`/`TENSORS`/
 * `FIELD` sections. Binary output is always big-endian regardless of host.
 * `wedge` cells are permuted `[0,2,1,3,5,4]` to VTK's node order; every other
 * type is written in meshio++'s natural order.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary true for big-endian binary numeric data, false for ascii
 * @param v51 true selects the version 5.1 `OFFSETS`+`CONNECTIVITY` `CELLS`
 *        layout, false selects the legacy 4.2 interleaved `[n,p0,...]`
 *        `CELLS`+`CELL_TYPES` layout
 * @throws WriteError on a field name containing spaces (VTK doesn't support
 *         them), a polyhedron cell block (unsupported by the C++ writer), an
 *         unknown cell type, or an unopenable output path
 * @note point_data/cell_data map generically to `SCALARS`/`VECTORS`/
 *       `TENSORS`/`FIELD` blocks; no reserved key names.
 */
void write_vtk(const std::string& rPath, const Mesh& rMesh, bool binary, bool v51);

/**
 * @brief Read a VTK legacy file.
 *
 * The version string on line 1 (`# vtk DataFile Version <ver>`) selects
 * between the 4.2 and 5.1 sub-reader (only the literal value `"5.1"`
 * triggers the 5.1 path; anything else, including genuinely older version
 * strings, goes through the 4.2 path). `COLOR_SCALARS` sections are read and
 * discarded (only to advance the file cursor correctly). `LOOKUP_TABLE`
 * entries after a `SCALARS` line are consumed but discarded.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if `DATASET` is anything other than
 *         `UNSTRUCTURED_GRID` (structured points/grid, rectilinear grid all
 *         fall back to Python), on a truncated binary section, an ascii
 *         parse failure, an unknown VTK data-type token, or an unrecognized
 *         section keyword
 * @note point_data/cell_data map generically from `SCALARS`/`VECTORS`/
 *       `TENSORS`/`FIELD` blocks; `point_sets`/`cell_sets` round-trip as
 *       extra data arrays (5.1 files only), same convention as VTU.
 */
Mesh read_vtk(const std::string& rPath);

}  // namespace meshioplusplus
