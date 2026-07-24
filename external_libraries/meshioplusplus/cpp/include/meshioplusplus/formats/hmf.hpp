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
 * @file hmf.hpp
 * @brief HMF (.hmf) — meshio++'s experimental HDF5 mesh container.
 *
 * HMF is not a third-party format: it is meshio++'s own HDF5-backed
 * container, reusing the XDMF topology-name vocabulary
 * (`meshio_to_xdmf_type`/`xdmf_to_meshio_type`, see xdmf.hpp) for its
 * `TopologyType` attribute. Layout: `domain/grid/Geometry` (points, with a
 * `GeometryType` attribute "X"/"XY"/"XYZ" — asserted but otherwise unused
 * after validation), one `domain/grid/Topology{k}` dataset per cell block,
 * and `domain/grid/NodeAttributes/<name>` / `CellAttributes/<name>` groups
 * for point_data/cell_data keyed by name verbatim. Only one `domain`/`grid`
 * pair is supported per file. **The format may change at any time** — the
 * writer always emits a warning to that effect.
 *
 * If two `Topology{k}` datasets resolve to the same meshio++ type, the
 * reader deliberately replicates the Python "later entry wins" semantics
 * (accumulation is keyed by meshio++ type name) rather than merging or
 * erroring. Unlike the reference `h5py` reader (which has a known
 * correctness issue here), the C++ reader correctly round-trips
 * **multi-block** cell data sharing the same `CellAttributes` name.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to meshio++'s HMF (.hmf) HDF5 container.
 *
 * Emits file attributes `type="hmf"`, `version="0.1-alpha"`, then
 * `domain/grid/Geometry` (points) with a `GeometryType` matching the point
 * dimensionality, one `Topology{k}` dataset per cell block (named via the
 * XDMF type-name table), and `NodeAttributes`/`CellAttributes` subgroups
 * for every point_data/cell_data key (any key name is preserved verbatim;
 * multiple cell blocks contributing to the same cell_data name are
 * concatenated). Always logs a warning that the format may change.
 *
 * @param path filesystem path to the .hmf file to create/overwrite
 * @param mesh the mesh to write
 * @param gzip_level HDF5 gzip compression level (0 = none) applied to the
 *        written datasets
 * @throws WriteError on an unsupported cell type or mesh layout
 */
void write_hmf(const std::string& rPath, const Mesh& rMesh, int gzip_level);

/**
 * @brief Read a meshio++ HMF (.hmf) HDF5 container into a Mesh.
 *
 * Reads `domain/grid/Geometry` as points and every `Topology{k}` dataset as
 * a cell block, resolving its meshio++ type from the `TopologyType`
 * attribute via the shared XDMF type table. If two `Topology{k}` datasets
 * map to the same meshio++ type, the later one (by dataset index) replaces
 * the earlier one rather than merging — this exact "later entry wins"
 * semantics is deliberately kept for parity with the reference Python
 * reader. `NodeAttributes`/`CellAttributes` datasets become point_data/
 * cell_data keyed by their stored name; multi-block cell_data under one
 * name round-trips correctly here even though the reference `h5py` reader
 * has a known bug for that case.
 *
 * @param path filesystem path to the .hmf file to read
 * @return the read Mesh
 * @throws ReadError if `GeometryType` is not one of "X"/"XY"/"XYZ", or on a
 *         malformed/unsupported HDF5 layout
 */
Mesh read_hmf(const std::string& rPath);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
