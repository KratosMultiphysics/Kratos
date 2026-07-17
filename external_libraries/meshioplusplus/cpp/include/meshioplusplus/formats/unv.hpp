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
 * @file unv.hpp
 * @brief I-DEAS Universal (.unv) C++ reader/writer — datasets 2411 (nodes)
 *        and 2412 (elements) only.
 *
 * A UNV file is a sequence of datasets, each delimited by a line containing
 * only `-1`, followed by a numeric dataset-id line and the dataset body.
 * Dataset **2411**: two-line node records (`label CS1 CS2 color` then a
 * coordinate line, Fortran `D`/`d` exponents normalized before parsing);
 * node labels are arbitrary integers, so a `label -> 0-based index` map is
 * built while reading. Dataset **2412**: a 6-integer record (`label fedesc
 * pid ... ... num_nodes`) selects the meshio++ type from the FE-descriptor
 * id (11/21->line, 22/24->line3, 41/81/91->triangle, 42/82/92->triangle6,
 * 44/84/94/122->quad, 45/85/95->quad8, 111->tetra, 118->tetra10,
 * 112->wedge, 115->hexahedron, 116->hexahedron20), followed by an extra
 * discarded 3-integer orientation record for beam descriptors (11/21/22/24)
 * — beam orientation is a genuinely lossy round-trip, always rewritten as
 * `0 0 0` on write — then the node-label records themselves.
 *
 * Parabolic (second-order) types use the Salome/Code-Aster mid-node
 * "sandwich" ordering (corner, mid-node, corner, mid-node, ...), converted
 * to meshio++'s "all corners then all edge nodes" convention via a fixed
 * permutation table per type (line3 `[0,2,1]`, triangle6
 * `[0,3,1,4,2,5]`, quad8 `[0,4,1,5,2,6,3,7]`, tetra10
 * `[0,4,1,5,2,6,7,8,9,3]`, hexahedron20 20-entry table) — applied directly
 * on read and inverted on write.
 *
 * Field/results datasets (2414 and legacy 55, 56, 57) are read
 * and written by the C++ core: data at nodes (location 1) -> `point_data`,
 * data on elements (location 2) -> `cell_data`; the field name becomes the
 * data key (de-duplicated on collision), and the component count (1/3/6/9)
 * is the array's inner dimension. On write, the default emits dataset 2414;
 * with `code_aster=true` it emits dataset 55 for `point_data` and 57 for
 * `cell_data` (the Code-Aster convention). Complex data and the
 * nodes-on-elements location (3) are skipped with a warning.
 *
 * Permanent-group datasets (2467, 2477, 2452, 2435, 2432, 2430 ->
 * point_sets/cell_sets) are decoded by the `UnvInfo` overloads of read_unv /
 * write_unv (a side-channel, since point_sets/cell_sets are not part of the
 * Mesh/NDArray conversion layer); the group-less overloads ignore them.
 */

// System includes
#include <cstdint>
#include <map>
#include <string>
#include <vector>

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Side-channel carrying permanent-group-derived point/cell sets across
 *        the Mesh conversion boundary (the `point_sets`/`cell_sets` Python
 *        Mesh attributes are not part of the C++ Mesh/NDArray layer).
 *
 * Mirrors `AnsysInfo`: node groups (UNV entity type 8) become `mPointSets`
 * (0-based node indices); element groups (entity type 7) become `mCellSets`
 * (per-cell-block lists of 0-based local cell indices, one inner list per
 * mesh cell block in block order).
 */
struct UnvInfo {
    std::map<std::string, std::vector<std::int64_t>> mPointSets;
    std::map<std::string, std::vector<std::vector<std::int64_t>>> mCellSets;
};

/**
 * @brief Write a mesh as a UNV file (datasets 2411 + 2412 only).
 *
 * Emits node records (dataset 2411, labels = 1-based row index) and element
 * records (dataset 2412), choosing one canonical FE descriptor per meshio++
 * type (line->21, line3->24, triangle->91, triangle6->92, quad->94,
 * quad8->95, tetra->111, tetra10->118, wedge->112, hexahedron->115,
 * hexahedron20->116), applying the inverse sandwich permutation for
 * parabolic types, and always writing a placeholder `0 0 0` beam
 * orientation record for line/line3 elements.
 *
 * Also emits field datasets from `point_data` (dataset 2414 location 1, or
 * dataset 55 in Code-Aster mode) and `cell_data` (dataset 2414 location 2, or
 * dataset 57 in Code-Aster mode); the reserved key `unv:pid` is excluded (it
 * is the per-element property id carried by dataset 2412, not a field).
 *
 * @param rPath filesystem path to write
 * @param rMesh the mesh to write
 * @param code_aster emit legacy datasets 55/57 for fields instead of 2414
 * @param node_dataset node dataset id to emit — `2411` (default) or `781`
 * @throws WriteError if the mesh carries `point_sets`/`cell_sets` (no
 *         dataset-2467 writer in C++ — the shim falls back to Python)
 * @note unsupported cell types are warned about and skipped (matching the
 *       Python writer); reads `cell_data["unv:pid"]` for the per-element
 *       property id (defaults to `1` if absent).
 */
void write_unv(const std::string& rPath, const Mesh& rMesh, bool code_aster = false,
               int node_dataset = 2411);

/**
 * @brief Write a mesh plus permanent groups (dataset 2467) as a UNV file.
 *
 * Same as the group-less overload, additionally emitting `rInfo`'s point sets
 * (node groups, entity type 8) and cell sets (element groups, entity type 7)
 * as dataset-2467 records after the field datasets.
 *
 * @param rInfo point/cell sets to emit as dataset-2467 groups
 */
void write_unv(const std::string& rPath, const Mesh& rMesh, const UnvInfo& rInfo,
               bool code_aster = false, int node_dataset = 2411);

/**
 * @brief Read a UNV file's node (2411) and element (2412) datasets.
 *
 * Splits the file into datasets on `-1` delimiter lines, builds a node
 * label->index map from dataset 2411, then decodes dataset 2412 element
 * records into typed cell blocks using the FE-descriptor table and the
 * sandwich-order permutation for parabolic types.
 *
 * Field datasets (2414/55/56/57) are decoded into `point_data`/`cell_data`.
 *
 * This group-less overload discards any permanent groups; use the `UnvInfo`
 * overload to receive them.
 *
 * @param rPath filesystem path to read
 * @return the read Mesh
 * @note cell_data key produced: `"unv:pid"` (element property id, dataset-
 *       2412 record-1 field 2); field datasets add point_data/cell_data keyed
 *       by field name.
 */
Mesh read_unv(const std::string& rPath);

/**
 * @brief Read a UNV file, additionally decoding permanent-group datasets
 *        (2467/2477/2452/2435/2432/2430) into `rInfo`.
 *
 * @param[out] rInfo receives node groups as `mPointSets` and element groups
 *        as `mCellSets` (0-based indices).
 */
Mesh read_unv(const std::string& rPath, UnvInfo& rInfo);

}  // namespace meshioplusplus
