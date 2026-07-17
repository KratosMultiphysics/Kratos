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
 * @file tecplot.hpp
 * @brief Tecplot ASCII finite-element (.dat/.tec) C++ reader/writer,
 *        single-zone only.
 *
 * A Tecplot FE file is a `VARIABLES = "X" "Y" "Z" ...` list followed by one
 * or more `ZONE T="..." N=<nodes> E=<elements> F=FEPOINT|FEBLOCK
 * ET=TRIANGLE|FEQUADRILATERAL|FETETRAHEDRON|FEBRICK [VARLOCATION=(...)]`
 * blocks. meshio++ only reads/writes a **single** FE zone: on read, only the
 * first zone is parsed and any subsequent zones are silently ignored (not
 * merged or errored on). `VARLOCATION=([a-b]=CELLCENTERED)` (1-based,
 * inclusive ranges) marks cell-centered variables, otherwise cell-centered-
 * ness is inferred from `NV=`. `FEBLOCK` packing reads one variable's full
 * array before the next; `FEPOINT` reads one full-variable-tuple row per
 * node. `X`/`x` (and optional `Y`/`Z`) become point coordinates; everything
 * else becomes point_data or cell_data, keyed by the raw variable name (no
 * `tecplot:` prefix).
 *
 * Zone type -> meshio++ type: LINESEG/FELINESEG->line,
 * TRIANGLE/FETRIANGLE->triangle, QUADRILATERAL/FEQUADRILATERAL->quad,
 * TETRAHEDRON/FETETRAHEDRON->tetra, BRICK/FEBRICK->hexahedron. On write,
 * pyramid/wedge/hexahedron all degrade to an 8-node FEBRICK zone, padding
 * with duplicated corner nodes (pyramid: `[0,1,2,3,4,4,4,4]`; wedge:
 * `[0,1,4,3,2,2,5,5]`).
 *
 * The multi-cell-type write path (Python degrades everything into a single
 * FEQUADRILATERAL/FEBRICK zone via "order_2" padding tables) exists **only
 * in the Python writer**: the C++ writer throws WriteError as soon as more
 * than one distinct cell type is present, forcing the Python fallback.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a mesh as a single Tecplot FE zone.
 *
 * Emits `VARIABLES`/`ZONE` headers for the one supported cell type present
 * (line/triangle/quad/tetra/hexahedron, with pyramid/wedge padded into
 * FEBRICK), then FEBLOCK-packed coordinate/point_data/cell_data columns
 * (data wrapped at 20 values per line) and 1-based connectivity.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if the mesh contains **more than one** distinct cell
 *         type (the Python fallback handles that case by degrading
 *         everything into one FEQUADRILATERAL/FEBRICK zone)
 */
void write_tecplot(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a Tecplot ASCII file's first FE zone.
 *
 * Parses the `VARIABLES` list and the first `ZONE` header (tolerating
 * multi-line continuation and a quoted `T="..."` title), then its FEBLOCK or
 * FEPOINT data body and 1-based connectivity. Any zones after the first are
 * silently ignored.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @throws ReadError if `X`/`x` is missing, the zone header uses an
 *         unsupported `F=`/`ZONETYPE=` combination, or the header/data
 *         otherwise doesn't parse (e.g. an adversarial zone title that is
 *         literally the string `"VARLOCATION"`) — the shim then falls back
 *         to the more tolerant Python reader.
 * @note point_data/cell_data keys are the raw Tecplot variable names (no
 *       prefix); `X`/`Y`/`Z` are reserved for coordinates.
 */
Mesh read_tecplot(const std::string& rPath);

}  // namespace meshioplusplus
