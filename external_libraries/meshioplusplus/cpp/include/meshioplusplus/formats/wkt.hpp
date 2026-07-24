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
 * @file wkt.hpp
 * @brief WKT (Well-Known Text) Triangulated Irregular Network C++
 *        reader/writer.
 *
 * A WKT TIN is a single `TIN (((x y z, x y z, x y z, x y z)), ...)`
 * expression: each triangle is one closed 4-point ring (`((p0, p1, p2,
 * p0))`) — the 4th point repeats the 1st to close the ring. The C++ reader
 * parses this by tracking **parenthesis depth** rather than matching
 * literal substrings (the point list sits at depth 3: `TIN`->1, the
 * triangle polygon->2, its linestring->3), which makes it naturally
 * tolerant of arbitrary whitespace/newlines between and inside the nested
 * parentheses. Points are de-duplicated by **exact** floating-point value
 * (no epsilon tolerance) in first-occurrence order; the repeated closing
 * point of each ring is dropped once the 3 unique corner indices are
 * recovered. A ring whose last point doesn't equal its first is a parse
 * error. `triangle` is the only cell type WKT can produce, and no
 * point_data/cell_data/field_data are read or written.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh's triangles as a WKT `TIN (...)` expression.
 *
 * Emits one closed 4-point ring per `triangle` cell (re-appending each
 * triangle's first point to close the ring). Only `triangle` cells are
 * representable; no point_data/cell_data is emitted (WKT carries none).
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write (only `triangle` cells contribute)
 * @throws WriteError on an unopenable output path
 */
void write_wkt(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a WKT TIN file into a single-`triangle`-block Mesh.
 *
 * Parses the `TIN (((...)), ...)` expression by tracking parenthesis depth,
 * de-duplicating points by exact value in first-occurrence order and
 * dropping each ring's repeated closing point.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (points plus a single `triangle` CellBlock; no
 *         point_data/cell_data/field_data)
 * @throws ReadError if a ring's last point does not equal its first (not a
 *         closed linestring), or the file doesn't parse as `TIN (...)`
 */
Mesh read_wkt(const std::string& rPath);

}  // namespace meshioplusplus
