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
 * @file mphtxt.hpp
 * @brief COMSOL text mesh (.mphtxt) C++ reader/writer.
 *
 * `.mphtxt` (as handled by FEconv) is a flat ASCII token stream (comments
 * from `#` to end of line): a version pair, tag-name and type-name tables
 * (discarded), then one or more "object" records — **only the first mesh
 * object in the file is parsed**; the rest are ignored entirely. That
 * object holds `sdim` (spatial dimension), `n_points`, `lowest` (the file's
 * actual lowest node index, not assumed to be 1), node coordinates, and a
 * sequence of element-type blocks (hybrid/multi-type meshes are
 * supported), each: a COMSOL type name (`"tet"`, `"tri2"`, …), node/element
 * counts, connectivity shifted by `-lowest` to 0-based, discarded
 * parameter tokens, a per-element **geometric entity index** ->
 * `cell_data["mphtxt:geom"]`, and discarded up/down topology-link pairs.
 *
 * COMSOL <-> meshio++ type map: `vtx`->`vertex`, `edg`/`edg2`->`line`/
 * `line3`, `tri`/`tri2`->`triangle`/`triangle6`, `quad`/`quad2`->`quad`/
 * `quad9`, `tet`/`tet2`->`tetra`/`tetra10`, `prism`/`prism2`->`wedge`/
 * `wedge18`, `pyr`->`pyramid`, `hex`/`hex2`->`hexahedron`/`hexahedron27`.
 * Node-order permutation is applied for `quad` (`[0,1,3,2]`, self-inverse)
 * and `hexahedron` (`[0,1,3,2,4,5,7,6]`, self-inverse); every other type
 * uses natural order.
 */

// System includes
#include <string>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Write a Mesh to a COMSOL text mesh (.mphtxt) file.
 *
 * Emits the version/tag/type-name header tables, then a single mesh
 * object with 1-based-equivalent (`lowest = 1`) connectivity, applying the
 * `quad`/`hexahedron` node-order permutation on the way out. Per-element
 * geometric entity indices come from `cell_data["mphtxt:geom"]` (one array
 * per block); parameter and up/down-link sections are always written
 * empty/zero.
 *
 * @param rPath filesystem path to the .mphtxt file to create/overwrite
 * @param rMesh the mesh to write
 * @throws WriteError on a cell type with no COMSOL equivalent (the C++
 *         writer raises here, forcing the Python fallback, whereas the
 *         Python writer merely warns and skips the type)
 * @note reads `cell_data["mphtxt:geom"]` if present
 */
void write_mphtxt(const std::string& rPath, const Mesh& rMesh);

/**
 * @brief Read a COMSOL text mesh (.mphtxt) file into a Mesh.
 *
 * Parses only the **first** mesh object in the file (subsequent objects
 * are ignored). Reads node coordinates, then each element-type block,
 * converting connectivity from the file's `lowest`-based indexing to
 * 0-based and applying the inverse `quad`/`hexahedron` node-order
 * permutation. Element parameter values and up/down topology-link pairs
 * are discarded.
 *
 * @param rPath filesystem path to the .mphtxt file to read
 * @return the read Mesh, with `cell_data["mphtxt:geom"]` populated (one
 *         array per cell block) from each element's geometric entity index
 * @throws ReadError on a malformed file or an unrecognized COMSOL type name
 */
Mesh read_mphtxt(const std::string& rPath);

}  // namespace meshioplusplus
