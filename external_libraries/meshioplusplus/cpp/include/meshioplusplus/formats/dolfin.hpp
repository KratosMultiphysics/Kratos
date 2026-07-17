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
 * @file dolfin.hpp
 * @brief Legacy DOLFIN/FEniCS XML (.xml) C++ reader/writer.
 *
 * A DOLFIN XML file holds exactly one mesh
 * (`<dolfin><mesh celltype="triangle|tetrahedron" dim="2|3"><vertices/>
 * <cells/></mesh></dolfin>`), `triangle`/`tetra` only, with meshio++'s own
 * node order (no permutation needed). Vertices and cells are placed by their
 * `index` attribute, not document order. Each `cell_data` array lives in a
 * **separate sibling file** `<stem>_<name>.xml`
 * (`<dolfin><mesh_function type="int|uint|float" dim="D" size="N"><entity
 * index="..." value="..."/></mesh_function></dolfin>`), matched by scanning
 * the mesh file's directory for the regex `"{stem}_([^.]+)\.xml"`; the `dim`
 * attribute there is **not** the topological dimension — it is a z-flatness
 * check (`2` if the mesh is 2D or all point z-coordinates are ~0, else `3`).
 * Implemented via the vendored pugixml plus `std::filesystem` for the
 * directory scan; no Python fallback is needed for this format. See
 * doc/formats/dolfin.md.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a DOLFIN XML mesh file (plus one sibling
 *        `<stem>_<name>.xml` per cell_data array).
 *
 * If the mesh has both `triangle` and `tetra` cells, `tetra` is preferred
 * and every other cell type is discarded (DOLFIN XML stores exactly one
 * cell type per mesh) — this call always emits the legacy-format warning.
 *
 * @param path filesystem path to write (sibling cell_data files are placed
 *        next to it, named from its stem)
 * @param mesh the mesh to write
 * @throws WriteError if, after preferring tetra, the mesh has neither
 *         triangle nor tetra cells, or if `mesh`'s dimension is not 2 or 3,
 *         or if a file cannot be opened for writing
 * @note writes one `<stem>_<key>.xml` file per `cell_data` key; no
 *       point_data or field_data is ever written
 */
void write_dolfin(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a DOLFIN XML mesh file (plus any sibling cell_data files).
 *
 * Parses the `<mesh>` element (triangle or tetrahedron only) by `index`
 * attribute via pugixml, then scans the file's directory for sibling
 * `<stem>_<name>.xml` files and reads each as one `cell_data[name]` array.
 *
 * @param path filesystem path to read
 * @return the read Mesh (cell_data from sibling files; no point_data, no
 *         field_data)
 * @throws ReadError if the main file can't be parsed, is missing `<dolfin>`/
 *         `<mesh>`, names an unsupported cell type, or if a sibling
 *         cell-data file contains more than one `<mesh_function>`
 */
Mesh read_dolfin(const std::string& rPath);

}  // namespace meshioplusplus
