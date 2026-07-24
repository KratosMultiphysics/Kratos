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
 * @file exodus.hpp
 * @brief Exodus II (.e/.exo/.ex2) C++ reader/writer, stored in netCDF using
 *        its classic variable/dimension conventions.
 *
 * Key variables: `coord(num_dim, num_nodes)` (transposed relative to
 * meshio++'s `(n, dim)` layout) or separate `coordx`/`coordy`/`coordz`
 * (both accepted on read); `eb_prop1(num_el_blk)` arbitrary distinct block
 * ids; `connect{k}(num_el_in_blk{k}, num_nod_per_el{k})` per element block
 * with a text `elem_type` attribute and 1-based node indices;
 * `name_nod_var`/`vals_nod_var{k}` point-data (first timestep only);
 * `name_elem_var`/`vals_elem_var{idx}[eb{block}]` cell data, concatenated
 * across blocks then re-split by target cell-block size. Compiled in only
 * when `MESHIOPLUSPLUS_HAS_NETCDF` is defined; otherwise the Python
 * `netCDF4` fallback handles this format. See doc/formats/exodus.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_NETCDF

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as an Exodus II (netCDF classic) file.
 *
 * Writes global attrs (`title`, `version=5.1f`, `api_version=5.1f`,
 * `floating_point_word_size=8`), a dummy single `0.0` `time_whole` step, one
 * `connect{k}` variable per cell block (element type mapped through the
 * canonical meshio++ -> Exodus reverse table, e.g. `hexahedron -> HEX8`,
 * `tetra -> TETRA`, `tetra4 -> TET4` as a distinct entry from plain
 * `tetra`), and point_data/cell_data as `vals_nod_var`/`vals_elem_var`
 * variables.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @throws WriteError if a cell block's type has no entry in the meshio++ ->
 *         Exodus type table, or if the connectivity dtype is unsupported
 * @note the shim only attempts this C++ path when `mesh.point_sets` is
 *       empty — the C++ writer has no support for Exodus node sets at all
 */
void write_exodus(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read an Exodus II (netCDF classic) file.
 *
 * Reads coordinates (either `coord` or `coordx`/`coordy`/`coordz`), one cell
 * block per `connect{k}` variable (via its `elem_type` attribute and the
 * Exodus -> meshio++ type table), and point_data/cell_data — with the
 * point-data name recombination quirk `categorize()` reproduces on purpose
 * from the reference implementation: names ending `X`/`Y`/`Z` (or `_R`/`_Z`)
 * are stacked into a 3- (or 2-) component vector when a sibling exists, but
 * the "sibling found" check uses Python truthiness on the found variable
 * index, so index `0` is treated the same as "not found" — a latent
 * reference-implementation edge case deliberately preserved rather than
 * fixed, so the two implementations agree. Only the first timestep is ever
 * read (a warning is emitted if more exist, matching a known ParaView writer
 * limitation).
 *
 * @param rPath filesystem path to read
 * @return the read Mesh, with `mesh.point_sets` from node sets (1-based in
 *         file) and `mesh.info` from `info_records`/`qa_records`
 * @throws ReadError if a variable has an unsupported netCDF type, point-data
 *         names are inconsistent, a `connect{k}` names an unknown Exodus
 *         element type, the connectivity dtype is unsupported, or the file
 *         contains `info_records`/`qa_records`/`ns_names`/`node_ns*` — any
 *         of the latter always defers the whole file to the Python fallback
 *         since node sets/info strings aren't carried by the conversion
 *         layer
 * @note point_data keys ending X/Y/Z or _R/_Z may be recombined into vector
 *       arrays; cell_data is split per cell block by node count
 */
Mesh read_exodus(const std::string& rPath);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_NETCDF
