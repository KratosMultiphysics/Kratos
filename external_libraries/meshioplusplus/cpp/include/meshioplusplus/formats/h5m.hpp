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
 * @file h5m.hpp
 * @brief MOAB H5M (.h5m) HDF5-backed C++ reader/writer.
 *
 * MOAB stores its mesh under a `tstt` root group in an HDF5 file: node
 * coordinates at `tstt/nodes/coordinates` (1-based `start_id`), one
 * `tstt/elements/<H5M_TYPE>/connectivity` dataset per cell block (e.g.
 * `Tet4`, `Tri3`, `Edge2`, `Hex8`, `Prism6`, `Pyramid5`, `Quad4`), a global
 * tag registry under `tstt/tags/<name>`, and per-node tag data under
 * `tstt/nodes/tags/<name>`. 2D point-data arrays are stored as `(n,)`
 * datasets of `k`-tuples via an HDF5 ARRAY/compound datatype (created with
 * `H5Tarray_create2`), not as `(n,k)` datasets. Every element/node index is
 * 1-based on disk; the reader/writer apply the `+1`/`-1` shift.
 *
 * The writer only supports three cell types on the way out
 * (`line`->`Edge2`, `triangle`->`Tri3`, `tetra`->`Tet4`); any other type is
 * silently skipped. Element/cell tags and MOAB "sets" are not read at all.
 * There is **no cell_data support end-to-end**: the reference Python
 * writer's cell-data path has a pre-existing bug (it misattributes the last
 * `elements` sub-group from a prior loop to every cell type), which the C++
 * writer deliberately does not replicate — it simply never writes cell
 * data, and the shim only attempts the C++ write path when
 * `mesh.cell_data` is empty (see doc/formats/h5m.md quirks).
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a MOAB H5M (.h5m) file.
 *
 * Points are written under `tstt/nodes/coordinates` (1-based `start_id`
 * tracked in a running global-id counter shared with every element block).
 * Only `line`, `triangle`, and `tetra` cell blocks are emitted (as `Edge2`,
 * `Tri3`, `Tet4` respectively); any other cell type present in the mesh is
 * silently skipped (no warning, unlike the Python fallback). Arbitrary
 * `point_data` keys are written as tag datasets under `nodes/tags/<name>`
 * plus a registry entry under `tstt/tags/<name>`. `cell_data` is never
 * written by this function.
 *
 * @param rPath filesystem path to the .h5m file to create/overwrite
 * @param rMesh the mesh to write
 * @param add_global_ids if true, write a conventional `GLOBAL_ID` node tag
 *        (values `1..n`) when the mesh doesn't already carry one
 * @param gzip_level HDF5 gzip compression level (0 = none) applied to the
 *        written datasets
 * @throws WriteError on an unsupported layout
 */
void write_h5m(const std::string& rPath, const Mesh& rMesh, bool add_global_ids, int gzip_level);

/**
 * @brief Read a MOAB H5M (.h5m) file into a Mesh.
 *
 * Reads `tstt/nodes/coordinates` and every `tstt/elements/<H5M_TYPE>/`
 * connectivity block, mapping H5M type names to meshio++ types (`Edge2`->
 * `line`, `Tri3`->`triangle`, `Tet4`->`tetra`, `Prism6`->`wedge`,
 * `Pyramid5`->`pyramid`, `Quad4`->`quad`, `Hex8`->`hexahedron`).
 * Connectivity is 1-based on disk and shifted to 0-based. Per-node tag
 * datasets under `nodes/tags/<name>` become `point_data`. Element/cell tags
 * and the `sets` group are ignored entirely (MOAB supports them; this
 * reader does not read them).
 *
 * @param rPath filesystem path to the .h5m file to read
 * @return the read Mesh (points, cells, point_data only — no cell_data)
 * @throws ReadError on a malformed/unsupported HDF5 layout
 */
Mesh read_h5m(const std::string& rPath);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
