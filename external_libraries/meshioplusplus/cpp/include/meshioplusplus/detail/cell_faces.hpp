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
 * @file cell_faces.hpp
 * @brief Per-cell-type boundary-face topology tables (local node indices of
 * each face of a 3D volume cell), used by the skin extractor (`skin.hpp`).
 *
 * Every face row lists **corner nodes first** (wound so the face normal
 * points *outward*, away from the element interior), then the mid-edge nodes
 * (mid of corner k → corner k+1), then the face-center node — i.e. each row
 * is itself a valid meshio/VTK `triangle6`/`quad8`/`quad9` node ordering.
 *
 * Node numbering is meshio's (= VTK's). Conventions baked in, matching the
 * rest of this repo (see `openfoam.cpp`'s `build_*` orientation checks and
 * the `tests/helpers.py` fixtures):
 *  - tetra: base `(0,1,2)` normal points toward apex 3 → outward base is
 *    `(0,2,1)`.
 *  - hexahedron: base `(0,1,2,3)` normal points toward the top `(4,5,6,7)`.
 *  - wedge: base `(0,1,2)` normal points toward the top `(3,4,5)` (the
 *    gmsh-like layout this codebase uses; do NOT trust vtkWedge's own
 *    doc-comment/face array, whose orientation is famously inconsistent).
 *  - pyramid: base `(0,1,2,3)` normal points toward apex 4.
 *  - wedge15 mid-node numbering is pure VTK_QUADRATIC_WEDGE — EnSight's
 *    penta15 involution (see CLAUDE.md) must NOT be applied here.
 *
 * The outward winding of every row is enforced by a gtest invariant
 * (`cpp/tests/test_skin.cpp`): on the reference element, the Newell normal
 * of each face's corner ring must point away from the cell centroid.
 *
 * KEEP IN SYNC: `src/meshioplusplus/_skin.py` carries the Python twin of
 * these tables for the pure-Python fallback — any change here must be
 * mirrored there.
 */

// System includes
#include <array>
#include <cstdint>
#include <vector>

// Project includes
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {
namespace detail {

/**
 * @brief One boundary face of a volume cell: its meshio face type and the
 * local node indices into the cell's connectivity row.
 */
struct CellFaceDef {
    /// Face cell type: `Triangle`/`Quad`/`Triangle6`/`Quad8`/`Quad9`.
    CellType mFaceType;
    /// Number of corner nodes (3 or 4) — the leading entries of `mNodes`.
    std::uint8_t mNumCorners;
    /// Total node count of the face (3, 4, 6, 8, or 9).
    std::uint8_t mNumNodes;
    /// Local node indices, corners first (outward winding), then mid-edge
    /// nodes, then the face center; only the first `mNumNodes` are valid.
    std::array<std::uint8_t, 9> mNodes;
};

/**
 * @brief The boundary faces of a volume cell type, or an empty list for
 * types the skin extractor does not support.
 * @param VolumeType The volume cell type to query.
 * @return Reference to the process-wide face table (empty if unsupported).
 */
inline const std::vector<CellFaceDef>& cell_faces(CellType VolumeType) {
    using CT = CellType;
    static const std::vector<CellFaceDef> empty = {};
    static const std::vector<CellFaceDef> tetra = {
        {CT::Triangle, 3, 3, {0, 1, 3}},
        {CT::Triangle, 3, 3, {1, 2, 3}},
        {CT::Triangle, 3, 3, {2, 0, 3}},
        {CT::Triangle, 3, 3, {0, 2, 1}},
    };
    static const std::vector<CellFaceDef> tetra10 = {
        {CT::Triangle6, 3, 6, {0, 1, 3, 4, 8, 7}},
        {CT::Triangle6, 3, 6, {1, 2, 3, 5, 9, 8}},
        {CT::Triangle6, 3, 6, {2, 0, 3, 6, 7, 9}},
        {CT::Triangle6, 3, 6, {0, 2, 1, 6, 5, 4}},
    };
    static const std::vector<CellFaceDef> hexahedron = {
        {CT::Quad, 4, 4, {0, 4, 7, 3}}, {CT::Quad, 4, 4, {1, 2, 6, 5}},
        {CT::Quad, 4, 4, {0, 1, 5, 4}}, {CT::Quad, 4, 4, {3, 7, 6, 2}},
        {CT::Quad, 4, 4, {0, 3, 2, 1}}, {CT::Quad, 4, 4, {4, 5, 6, 7}},
    };
    static const std::vector<CellFaceDef> hexahedron20 = {
        {CT::Quad8, 4, 8, {0, 4, 7, 3, 16, 15, 19, 11}},
        {CT::Quad8, 4, 8, {1, 2, 6, 5, 9, 18, 13, 17}},
        {CT::Quad8, 4, 8, {0, 1, 5, 4, 8, 17, 12, 16}},
        {CT::Quad8, 4, 8, {3, 7, 6, 2, 19, 14, 18, 10}},
        {CT::Quad8, 4, 8, {0, 3, 2, 1, 11, 10, 9, 8}},
        {CT::Quad8, 4, 8, {4, 5, 6, 7, 12, 13, 14, 15}},
    };
    // VTK face-center numbering: 20=(0,1,5,4), 21=(1,2,6,5), 22=(2,3,7,6),
    // 23=(3,0,4,7), 24=bottom (0,1,2,3), 25=top (4,5,6,7); 26 = body center.
    static const std::vector<CellFaceDef> hexahedron27 = {
        {CT::Quad9, 4, 9, {0, 4, 7, 3, 16, 15, 19, 11, 23}},
        {CT::Quad9, 4, 9, {1, 2, 6, 5, 9, 18, 13, 17, 21}},
        {CT::Quad9, 4, 9, {0, 1, 5, 4, 8, 17, 12, 16, 20}},
        {CT::Quad9, 4, 9, {3, 7, 6, 2, 19, 14, 18, 10, 22}},
        {CT::Quad9, 4, 9, {0, 3, 2, 1, 11, 10, 9, 8, 24}},
        {CT::Quad9, 4, 9, {4, 5, 6, 7, 12, 13, 14, 15, 25}},
    };
    static const std::vector<CellFaceDef> wedge = {
        {CT::Triangle, 3, 3, {0, 2, 1}}, {CT::Triangle, 3, 3, {3, 4, 5}},
        {CT::Quad, 4, 4, {0, 1, 4, 3}},  {CT::Quad, 4, 4, {1, 2, 5, 4}},
        {CT::Quad, 4, 4, {2, 0, 3, 5}},
    };
    static const std::vector<CellFaceDef> wedge15 = {
        {CT::Triangle6, 3, 6, {0, 2, 1, 8, 7, 6}},
        {CT::Triangle6, 3, 6, {3, 4, 5, 9, 10, 11}},
        {CT::Quad8, 4, 8, {0, 1, 4, 3, 6, 13, 9, 12}},
        {CT::Quad8, 4, 8, {1, 2, 5, 4, 7, 14, 10, 13}},
        {CT::Quad8, 4, 8, {2, 0, 3, 5, 8, 12, 11, 14}},
    };
    static const std::vector<CellFaceDef> pyramid = {
        {CT::Quad, 4, 4, {0, 3, 2, 1}},  {CT::Triangle, 3, 3, {0, 1, 4}},
        {CT::Triangle, 3, 3, {1, 2, 4}}, {CT::Triangle, 3, 3, {2, 3, 4}},
        {CT::Triangle, 3, 3, {3, 0, 4}},
    };
    static const std::vector<CellFaceDef> pyramid13 = {
        {CT::Quad8, 4, 8, {0, 3, 2, 1, 8, 7, 6, 5}}, {CT::Triangle6, 3, 6, {0, 1, 4, 5, 10, 9}},
        {CT::Triangle6, 3, 6, {1, 2, 4, 6, 11, 10}}, {CT::Triangle6, 3, 6, {2, 3, 4, 7, 12, 11}},
        {CT::Triangle6, 3, 6, {3, 0, 4, 8, 9, 12}},
    };
    // pyramid14 = pyramid13 plus node 13 at the base-face center.
    static const std::vector<CellFaceDef> pyramid14 = {
        {CT::Quad9, 4, 9, {0, 3, 2, 1, 8, 7, 6, 5, 13}},
        {CT::Triangle6, 3, 6, {0, 1, 4, 5, 10, 9}},
        {CT::Triangle6, 3, 6, {1, 2, 4, 6, 11, 10}},
        {CT::Triangle6, 3, 6, {2, 3, 4, 7, 12, 11}},
        {CT::Triangle6, 3, 6, {3, 0, 4, 8, 9, 12}},
    };
    switch (VolumeType) {
        case CT::Tetra:
            return tetra;
        case CT::Tetra10:
            return tetra10;
        case CT::Hexahedron:
            return hexahedron;
        case CT::Hexahedron20:
            return hexahedron20;
        case CT::Hexahedron27:
            return hexahedron27;
        case CT::Wedge:
            return wedge;
        case CT::Wedge15:
            return wedge15;
        case CT::Pyramid:
            return pyramid;
        case CT::Pyramid13:
            return pyramid13;
        case CT::Pyramid14:
            return pyramid14;
        default:
            return empty;
    }
}

/**
 * @brief Whether the skin extractor supports a volume cell type.
 * @param Type The cell type to test.
 * @return `true` when `cell_faces(Type)` is non-empty.
 */
inline bool skin_supported(CellType Type) {
    return !cell_faces(Type).empty();
}

}  // namespace detail
}  // namespace meshioplusplus
