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
 * @file freefem.hpp
 * @brief FreeFem++ mesh (.msh) C++ reader/writer (as handled by FEconv).
 *
 * ASCII: a 3-integer header `nver n_el1 n_el2` (vertex count, then the two
 * element-block counts), then `nver` rows of `x y [z] ref`, `n_el1` volume-
 * element rows, and `n_el2` boundary-element rows, each row ending in an
 * integer region/boundary label. All connectivity is 1-based. The spatial
 * dimension is **inferred** from the first vertex row's token count minus
 * one (must resolve to 2 or 3) — there is no explicit dimension field. In
 * 2D, volume elements are `triangle` (3 nodes) and boundary elements are
 * `line` (2 nodes); in 3D, volume elements are `tetra` (4 nodes) and
 * boundary elements are `triangle` (3 nodes) — the only cell types this
 * format supports. Blank lines are ignored. The `.msh` extension is shared
 * with `ansys` and `gmsh`; on extension-based auto-detection `freefem` is
 * tried last, so pass `file_format="freefem"` explicitly. See
 * doc/formats/freefem.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a FreeFem++ .msh file.
 *
 * Emits the two dimension-appropriate cell types only (triangle+line for
 * 2D, tetra+triangle for 3D), each vertex/element row ending in its
 * `point_data`/`cell_data["freefem:ref"]` label (defaulting to zero when
 * absent for cell_data). Points use `%.16e` formatting (vs. full Python
 * `repr()` precision in the Python writer — same effective precision,
 * different string form).
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if `rMesh` is not 2D or 3D, or if it contains a cell
 *         type other than the two appropriate for its dimension (forcing
 *         the Python fallback, which performs a warn-and-skip instead of
 *         hard-failing)
 * @note reads/writes `point_data["freefem:ref"]` and
 *       `cell_data["freefem:ref"]`
 */
void write_freefem(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a FreeFem++ .msh file.
 *
 * Reads the 3-integer header, infers the spatial dimension from the first
 * vertex row's token count, then parses `n_el1` volume-element rows and
 * `n_el2` boundary-element rows (1-based connectivity, dimension-dependent
 * types as described in the file-level doc comment).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `point_data["freefem:ref"]` (per-vertex label)
 *         and `cell_data["freefem:ref"]` (per-element label, one array per
 *         cell block)
 * @throws ReadError if the file can't be opened, the header isn't 3
 *         integers, the inferred vertex dimension isn't 2 or 3, or a
 *         vertex/element section is truncated
 */
Mesh read_freefem(const std::string& rPath);

}  // namespace meshioplusplus
