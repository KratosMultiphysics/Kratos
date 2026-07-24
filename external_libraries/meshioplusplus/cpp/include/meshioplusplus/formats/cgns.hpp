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
 * @file cgns.hpp
 * @brief CGNS (.cgns) C++ reader/writer — a minimal tetrahedra-only subset
 *        stored in HDF5, not the full CGNS/SIDS specification.
 *
 * On-disk layout: `Base/Zone1/GridCoordinates/{CoordinateX,Y,Z}/" data"` and
 * `Base/Zone1/GridElements/{ElementRange," "ElementConnectivity}/" data"`
 * (the leading-space dataset name `" data"` in every leaf group is this
 * implementation's own ad hoc convention, not part of the real CGNS/HDF5
 * spec, but shared identically between the Python and C++ writers).
 * `ElementRange` is `[1, n_cells]` (1-based inclusive) and
 * `ElementConnectivity` is flat 1-based tetra connectivity; `+-1` is applied
 * on read/write while preserving the connectivity array's original integer
 * dtype. This is the least complete format meshio++ supports: `tetra` is the
 * only cell type accepted or emitted, and no point_data/cell_data/field_data
 * is read or written at all. Compiled in only when
 * `MESHIOPLUSPLUS_HAS_HDF5` is defined; otherwise the Python `h5py` fallback
 * handles this format with identical on-disk behavior. See
 * doc/formats/cgns.md.
 */

#ifdef MESHIOPLUSPLUS_HAS_HDF5

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write `mesh` as a minimal CGNS/HDF5 file (tetrahedra only).
 *
 * Emits `Base/Zone1/GridCoordinates` (CoordinateX/Y/Z) and
 * `Base/Zone1/GridElements` (ElementRange = `[1,n]`, ElementConnectivity =
 * flat 1-based node ids), converting 0-based to 1-based indices while
 * preserving the connectivity array's integer dtype.
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write — only its `"tetra"` cell block (if any) is
 *        emitted; any other cell type present is silently ignored, not
 *        warned
 * @param gzip_level HDF5 gzip compression level applied to every dataset
 *        (CoordinateX/Y/Z, ElementRange, ElementConnectivity); write-only —
 *        HDF5 decompresses transparently on read regardless of the level
 *        used to write
 * @throws WriteError if the connectivity array's dtype is unsupported
 */
void write_cgns(const std::string& rPath, const Mesh& rMesh, int gzip_level);

/**
 * @brief Read a CGNS/HDF5 file written by @ref write_cgns (or a compatible
 *        file following the same minimal layout).
 *
 * Reads `Base/Zone1/GridCoordinates` and `GridElements`, converting 1-based
 * connectivity to 0-based.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh (points + one `"tetra"` cell block; no point_data/
 *         cell_data/field_data)
 * @throws ReadError if `"Base"` or `"Base/Zone1"` is missing ("Malformed
 *         CGNS?"), `ElementRange`/`ElementConnectivity` are malformed, the
 *         connectivity doesn't reshape to exactly 4 columns per cell ("Can
 *         only read tetrahedra."), or the connectivity dtype is unsupported
 */
Mesh read_cgns(const std::string& rPath);

}  // namespace meshioplusplus

#endif  // MESHIOPLUSPLUS_HAS_HDF5
