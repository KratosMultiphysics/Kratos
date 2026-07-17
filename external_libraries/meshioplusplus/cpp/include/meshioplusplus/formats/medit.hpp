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
 * @file medit.hpp
 * @brief Medit / GMF (.mesh) ASCII C++ reader/writer.
 *
 * Medit (INRIA "libMeshb" GMF) is a keyword-section format. This header
 * covers only the **ASCII** `.mesh` variant: whitespace/`#`-comment
 * tokenized, `MeshVersionFormatted`/`Dimension` header, then `Vertices`
 * (`x1 x2 [x3] ref` rows, 1-based) and element keyword sections
 * (`Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/`Prisms`/`Pyramids`/
 * `Hexahedra`/`Hexaedra`), each a count followed by that many
 * `<node ids> ref` rows. The trailing per-row integer becomes
 * `point_data["medit:ref"]`/`cell_data["medit:ref"]`; Medit stores at most
 * one such integer column, so only the first int-dtype point/cell data
 * array is kept on write (extras are dropped with a warning).
 *
 * The binary `.meshb` GMF variant (position-indexed records, endianness
 * flip via a leading magic code, version-dependent int/float widths) is
 * **not implemented here** — dispatch on the `"b"` filename suffix always
 * routes to the Python fallback for that variant.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to an ASCII Medit (.mesh) file.
 *
 * Emits `MeshVersionFormatted 2` (float64 coordinates) and `Dimension`,
 * then one keyword section per cell type present (`Vertices` is always
 * written; `Edges`/`Triangles`/`Quadrilaterals`/`Tetrahedra`/`Prisms`/
 * `Pyramids`/`Hexahedra` as applicable), each row 1-based node ids
 * followed by a trailing reference integer. At most one int-dtype
 * point_data array and one int-dtype cell_data array are used to populate
 * the `ref` columns (Medit's single-reference-column limitation); if
 * `point_data`/`cell_data` contain more than one integer candidate, the
 * first is used and the rest are dropped with a warning. Ends with `End`.
 *
 * @param path filesystem path to the .mesh file to create/overwrite
 * @param mesh the mesh to write
 * @throws WriteError on an unsupported cell type
 * @note reads `point_data`/`cell_data` key `"medit:ref"` if present
 *       (preferred over other int-dtype arrays for the ref column)
 */
void write_medit_ascii(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an ASCII Medit (.mesh) file into a Mesh.
 *
 * Parses `MeshVersionFormatted` (0/1 -> float32 coords, 2 -> float64) and
 * `Dimension`, then every recognized element keyword section, converting
 * 1-based node ids to 0-based. Sections such as `Corners`, `Normals`,
 * `NormalAtVertices`, `SubDomainFromMesh`, `VertexOnGeometricVertex`/
 * `Edge`, `EdgeOnGeometricEdge`, `Identifier`, `Geometry`,
 * `RequiredVertices`, `TangentAtVertices`, `Tangents`, `Ridges` are
 * recognized only enough to be token-skipped, and are otherwise discarded.
 *
 * @param path filesystem path to the .mesh file to read
 * @return the read Mesh, with `point_data["medit:ref"]` and
 *         `cell_data["medit:ref"]` (one array per cell block) populated
 *         from each row's trailing reference integer
 * @throws ReadError on a malformed file or unrecognized required section
 */
Mesh read_medit_ascii(const std::string& rPath);

}  // namespace meshioplusplus
