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
 * @file triangle.hpp
 * @brief Triangle (.node/.ele/.poly) C++ reader/writer — Shewchuk's 2D
 * mesh generator, the planar analogue of TetGen.
 *
 * A `.node`/`.ele` path selects the shared-stem sibling pair: `<stem>.node`
 * (header `npoints 2 nattrs nbmarkers`, rows `idx x y attrs.. markers..`;
 * `dim` must be 2) plus, when present, `<stem>.ele` (header
 * `ntriangles 3|6 nattrs`, giving `triangle` or `triangle6` cells; a lone
 * `.node` file reads as a point cloud). A `.poly` path reads the PSLG file:
 * a vertex section (inline, or the sibling `.node` when the vertex count is
 * 0) plus a segment section that becomes a `line` cell block; holes and
 * regional attributes are skipped with a warning. The node index base (0 or
 * 1) is auto-detected and indices must be exactly consecutive, mirroring
 * the tetgen reader.
 *
 * Data naming mirrors tetgen: vertex attribute columns become
 * `point_data["triangle:attr<k>"]`, boundary-marker columns
 * `"triangle:ref"`/`"triangle:ref2"`/..., and element attribute / segment
 * marker columns become the matching `cell_data` keys.
 *
 * Note that `.node`/`.ele` default to the tetgen format in the extension
 * registry; the Python dispatcher falls through to this format when tetgen
 * rejects a 2D file, while the flat bindings need an explicit
 * `format="triangle"` (only `.poly` defaults here).
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as Triangle files.
 *
 * A `.node`/`.ele` path writes the sibling pair — points (2D only) with
 * attribute/marker columns from point_data, and every `triangle`/
 * `triangle6` block (which must all share one type) with attribute columns
 * from cell_data; other cell types are skipped. A `.poly` path writes a
 * PSLG: inline vertices, the `line` blocks as segments (markers from the
 * first `:ref` cell_data key), and zero holes.
 *
 * @param rPath filesystem path ending in `.node`, `.ele`, or `.poly`
 * @param rMesh the mesh to write (points must be 2D)
 * @throws WriteError if a file cannot be opened, points are not 2D, or
 *         `triangle` and `triangle6` blocks are mixed
 */
void write_triangle(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read Triangle files (`.node`/`.ele` pair or `.poly`).
 *
 * @param rPath filesystem path ending in `.node`, `.ele`, or `.poly`
 * @return the read Mesh (2-column points; `triangle`/`triangle6` cells from
 *         `.ele`, `line` cells from `.poly` segments)
 * @throws ReadError on malformed input, `dim != 2`, non-consecutive node
 *         indices, or an unsupported nodes-per-triangle count
 */
Mesh read_triangle(const std::string& rPath);

}  // namespace meshioplusplus
