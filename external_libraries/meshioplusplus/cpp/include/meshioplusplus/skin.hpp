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
 * @file skin.hpp
 * @brief Boundary-skin extraction: derive the surface mesh of a 3D volume
 * mesh (the algorithm of Kratos Multiphysics' `SkinDetectionProcess`).
 *
 * A face of a volume cell is a *boundary* face when its sorted corner-node
 * key appears exactly once across all volume cells — faces shared by two
 * cells are interior and cancel out. The per-cell-type face topology comes
 * from `detail/cell_faces.hpp`; supported volume types are tetra,
 * hexahedron, wedge, pyramid and their common higher-order variants
 * (tetra10, hexahedron20/27, wedge15, pyramid13/14), whose boundary faces
 * come out as `triangle6`/`quad8`/`quad9`.
 *
 * Used by the STL/PLY writers (skin-by-default for volume meshes), by the
 * SVG/TikZ 3D projection path, and exposed to Python as
 * `meshioplusplus.extract_skin`.
 */

// Project includes
#include "meshioplusplus/mesh.hpp"

namespace meshioplusplus {

/**
 * @brief Extract the boundary skin of a volume mesh as a new surface mesh.
 *
 * Walks every supported, non-ragged 3D cell block, keeps the faces whose
 * sorted corner-node key occurs exactly once (Kratos `SkinDetectionProcess`
 * hashing), and builds a new mesh containing only those faces. Output cell
 * blocks appear in the fixed canonical order `triangle`, `triangle6`,
 * `quad`, `quad8`, `quad9` (only non-empty blocks are emitted), with the
 * faces of each type in volume-block enumeration order. Points are
 * *compacted*: only nodes referenced by a boundary face are kept, in
 * ascending original-index order, and `point_data` rows are subset the same
 * way. `cell_data`, `field_data`, sets, and `info` are dropped (a boundary
 * face has no canonical 1:1 source cell).
 *
 * Unsupported 3D blocks (polyhedron/ragged blocks, Lagrange and other
 * very-high-order types) are skipped with a warning; existing 0/1/2-D
 * blocks are ignored (the skin is derived from volume cells only).
 *
 * @param rMesh the volume mesh to extract the skin from
 * @param linearize when true, only the corner nodes of each boundary face
 *        are emitted (`triangle`/`quad` output even for higher-order volume
 *        cells) and unused mid-nodes are compacted away — what the STL/PLY/
 *        SVG/TikZ writers need
 * @return a new surface mesh holding the boundary faces
 * @throws std::invalid_argument if `rMesh` contains no supported volume
 *         cell block at all
 */
Mesh extract_skin(const Mesh& rMesh, bool linearize = false);

/**
 * @brief Whether a mesh has at least one cell block the skin extractor
 * supports (a non-ragged 3D block with a known face table).
 * @param rMesh the mesh to test
 * @return true when `extract_skin(rMesh)` would have input to work on
 */
bool has_skinnable_cells(const Mesh& rMesh);

}  // namespace meshioplusplus
