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
 * @file stl.hpp
 * @brief STL (stereolithography) C++ reader/writer, ascii and binary.
 *
 * STL has no shared-vertex table: every triangle facet repeats its own
 * three vertex coordinates. ASCII layout is
 * `solid [name] / facet normal nx ny nz / outer loop / vertex x y z (x3) /
 * endloop / endfacet ... endsolid`. Binary layout is an 80-byte free-form
 * header, a little-endian `uint32` triangle count, then that many fixed
 * 50-byte records (`float32[3]` normal, `float32[3][3]` vertices, `int16`
 * attribute byte count, conventionally 0 and not enforced).
 *
 * Binary-vs-ascii detection: files under 80 bytes are ascii; otherwise the
 * header + triangle count are read and the expected size
 * `84 + num_triangles*50` is compared against the actual file size — this
 * deliberately avoids the naive "starts with the literal word `solid`"
 * check, since binary STL headers sometimes also start with that word.
 *
 * On read, all raw (possibly duplicate) triangle vertices are uniquified in
 * **first-occurrence order** into a shared point table, so point indices are
 * not preserved across a round-trip (only geometry is). ASCII coordinates
 * parse as float64; binary coordinates parse as float32, matching the
 * on-disk precision. Only `triangle` cells are supported; a mesh with other
 * cell types written to STL warns and drops everything else. An empty STL
 * (0 triangles) yields a Mesh with no cell blocks at all.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh's triangle facets to an STL file, ascii or binary.
 *
 * When `skin` is true (the default) and the mesh contains supported 3D
 * volume cells (tetra/hexahedron/wedge/pyramid and their common
 * higher-order variants), the boundary skin is extracted first
 * (`extract_skin`, linearized to corner nodes) with quads triangulated as
 * `(0,1,2)`/`(0,2,3)`, and only that skin is written; any pre-existing
 * surface blocks are dropped with a warning (a volume mesh's surface blocks
 * are usually its boundary — writing both would duplicate every facet).
 * With `skin=false` (or no volume cells) the legacy behavior applies:
 * non-triangle cell blocks are dropped and only existing `triangle` blocks
 * are written. If `cell_data["facet_normals"]` is present it is written
 * verbatim; otherwise normals are computed per-facet from the cross product
 * of two edge vectors.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param binary true for the 80-byte-header + 50-byte-record binary layout,
 *        false for the `solid`/`facet`/`endfacet` ascii layout
 * @param skin extract and write the boundary skin when volume cells are
 *        present (default true); false keeps the legacy triangles-only path
 * @throws WriteError on an unopenable output path
 * @note cell_data key produced/consumed: `"facet_normals"`.
 */
void write_stl(const std::string& rPath, const Mesh& rMesh, bool binary, bool skin = true);

/**
 * @brief Read an STL file, auto-detecting ascii vs. binary.
 *
 * Uses the file-size heuristic described above (never the "starts with
 * `solid`" check) to pick ascii or binary parsing, then de-duplicates raw
 * facet vertices in first-occurrence order to build the point table and a
 * single `triangle` cell block. ASCII parsing uses a fast custom line reader
 * that only inspects the last 3 whitespace-separated tokens per line
 * (discarding any leading keyword such as `vertex`).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (a single `triangle` CellBlock, or none if the file
 *         has zero triangles); points are float64 for ascii input, float32
 *         for binary input
 * @throws ReadError on a malformed/truncated file
 */
Mesh read_stl(const std::string& rPath);

}  // namespace meshioplusplus
