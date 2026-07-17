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
 * @file ply.hpp
 * @brief PLY (Polygon File Format / Stanford triangle format) C++
 *        reader/writer, ascii and binary (little/big endian).
 *
 * A PLY file is a header (`format ascii 1.0` / `format binary_little_endian
 * 1.0` / `format binary_big_endian 1.0`, then `element vertex <N>` +
 * `property <type> <name>` lines, optionally `element face <M>` + a
 * `property list <count_type> <index_type> vertex_indices`) followed by the
 * vertex rows and face rows. Unlike VTU/VTK, endianness is read straight off
 * the `format` line rather than from a separate byte-order attribute.
 *
 * Vertex properties beyond `x`/`y`/`z` (e.g. `nx,ny,nz` normals,
 * `confidence`, `intensity`) become `point_data`. Faces are grouped into
 * cell blocks by vertex count (1=vertex, 2=line, 3=triangle, 4=quad, else
 * polygon); in **binary** mode, since a face row's length is only known by
 * reading its own leading list-count, the reader first walks the buffer
 * computing per-row byte offsets and then groups **consecutive
 * constant-length runs** into separate cell blocks (position in the file
 * matters here, unlike most other formats where same-typed cells always
 * merge into one block). The C++ reader rejects any face `property` beyond
 * the single index list, and any list-typed *vertex* property, forcing the
 * Python fallback (which supports arbitrary extra face properties as
 * `cell_data`). On write, 64-bit integer cell data is silently downcast to
 * int32 (PLY has no 64-bit integer property type); only
 * vertex/line/triangle/quad/polygon cell types are writable.
 *
 * See doc/formats/ply.md for the full grammar and quirks (e.g. the
 * deliberate `uchar`-as-signed-1-byte quirk in the list-count parsing path).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh to a PLY file, ascii or binary.
 *
 * Writes `element vertex`/`element face` sections; point_data columns beyond
 * x/y/z become extra vertex properties, cell blocks are grouped by type into
 * one `property list` face element. Requires all cell blocks to share one
 * dtype; 64-bit integer cell_data is downcast to int32 with a warning (PLY
 * has no 64-bit int property type). Only vertex/line/triangle/quad/polygon
 * cell types are written; other types are skipped with a warning.
 * Multi-dimensional point_data is silently filtered.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary true for `binary_little_endian`/`binary_big_endian`
 *        (host-endian) output, false for `format ascii 1.0`
 * @throws WriteError on an unopenable output path
 */
void write_ply(const std::string& rPath, const Mesh& rMesh, bool binary);

/**
 * @brief Read a PLY file (ascii or binary, either endianness).
 *
 * Parses the header, decodes vertex rows into points plus point_data (any
 * vertex property beyond x/y/z), and decodes face rows into cell blocks
 * grouped by vertex count and, in binary mode, by contiguous constant-length
 * run. `obj_info` header lines (a MeshLab convention) are skipped without
 * being parsed.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if a face carries list-typed vertex properties or extra
 *         (non-index) face properties beyond `vertex_indices` — the shim
 *         then falls back to the Python reader, which does support those.
 * @note point_data keys are the raw PLY property names (no `ply:` prefix).
 */
Mesh read_ply(const std::string& rPath);

}  // namespace meshioplusplus
