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
 * @file obj_off.hpp
 * @brief Wavefront OBJ (.obj) and Geomview OFF (.off) ASCII C++
 *        readers/writers — two distinct surface formats sharing one
 *        header.
 *
 * **OBJ**: a line-oriented format with `v x y z` (points), `vn`/`vt`
 * (vertex normals / texture coordinates, stored raw with arbitrary column
 * count), `s` (smooth-shading toggle, ignored), `f i1[/t1[/n1]] ...` (faces
 * — only the leading vertex index of each `i/t/n` token is used; a run of
 * faces stays in one cell block until the per-face vertex count changes or
 * a new `g` line appears), and `g <name>` (starts a new group, incrementing
 * a running group-id counter starting at -1; only the numeric id is kept,
 * not the name string). Faces are grouped by vertex count into `triangle`
 * (3), `quad` (4), or `polygon` (else). Blank trailing groups are dropped.
 *
 * **OFF**: a minimal format — a literal `"OFF"` first line, a
 * `<nverts> <nfaces> <nedges>` header (edge count parsed but discarded),
 * `nverts` coordinate rows, then `nfaces` rows each `3 i j k` (a leading
 * vertex count that **must** be 3 — any other value is a hard `ReadError`,
 * since only triangular faces are supported). OFF carries no point_data,
 * cell_data, or field_data at all.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a Geomview OFF (.off) file.
 *
 * Emits the `"OFF"` header line, `<nverts> <nfaces> 0` (edge count always
 * 0), vertex coordinate rows, then one `3 i j k` row per triangle. Only
 * `triangle` cells are representable.
 *
 * @param rPath filesystem path to the .off file to create/overwrite
 * @param rMesh the mesh to write (triangle cells only)
 * @throws WriteError on any non-triangle cell type
 */
void write_off(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Geomview OFF (.off) file into a Mesh.
 *
 * Validates the `"OFF"` first line, reads the vertex/face/edge counts
 * (edge count discarded), then `nverts` coordinate rows and `nfaces`
 * triangle rows (each row's leading count must be exactly 3).
 *
 * @param rPath filesystem path to the .off file to read
 * @return the read Mesh (points + a single `triangle` cell block only —
 *         no point_data/cell_data/field_data)
 * @throws ReadError if the first line isn't `"OFF"`, or any face row's
 *         leading vertex count isn't 3 ("Can only read triangular faces")
 */
Mesh read_off(const std::string& rPath);

/**
 * @brief Write a Mesh to a Wavefront OBJ (.obj) file.
 *
 * Emits `v x y z` rows, `vn`/`vt` rows from `point_data["obj:vn"]`/
 * `["obj:vt"]` if present, and `f` face rows grouped by cell block
 * (triangle/quad/polygon), 1-based indices. Group (`g`) lines are emitted
 * per distinct value found in `cell_data["obj:group_ids"]`, if present.
 *
 * @param rPath filesystem path to the .obj file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on an unsupported cell type
 * @note reads `point_data["obj:vn"]`, `point_data["obj:vt"]`,
 *       `cell_data["obj:group_ids"]` if present
 */
void write_obj(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Wavefront OBJ (.obj) file into a Mesh.
 *
 * Parses `v` (points), `vn`/`vt` (stored raw, arbitrary column count),
 * and `f` (faces; only the leading vertex index of each `i[/t[/n]]` token
 * is kept — the `/vt`/`/vn` index references are discarded entirely).
 * Faces are grouped by vertex count into `triangle`/`quad`/`polygon` cell
 * blocks; a run of same-count faces breaks into a new block whenever the
 * count changes **or** a `g` line appears, even if the count didn't
 * change. `g <name>` increments a running group-id counter (starting at
 * -1 for faces before any `g` line); only the id survives, not the name.
 * Empty trailing groups are dropped after the full file is scanned.
 *
 * @param rPath filesystem path to the .obj file to read
 * @return the read Mesh, with `point_data["obj:vn"]` (if any `vn` lines
 *         were seen), `point_data["obj:vt"]` (if any `vt` lines were
 *         seen), and `cell_data["obj:group_ids"]` (one int array per cell
 *         block, the originating group id, `-1` if before the first `g`)
 * @throws ReadError on a malformed file
 */
Mesh read_obj(const std::string& rPath);

}  // namespace meshioplusplus
