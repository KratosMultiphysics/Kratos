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
 * @file vtk_common.hpp
 * @brief VTK cell-type metadata shared by the VTU and VTK legacy format
 * implementations, ported from `src/meshio/_vtk_common.py`.
 *
 * Holds the meshio-type <-> VTK-cell-type-id maps (`meshio_to_vtk_type`/
 * `vtk_to_meshio_type`), the one node-order quirk that differs between the
 * two conventions (`meshio_to_vtk_order`/`vtk_to_meshio_order`, for the
 * linear wedge), and `is_special_cell`, which flags the cell types whose
 * per-cell node count is not fixed (`"polygon"` and the VTK_LAGRANGE_*
 * family) and therefore need the offsets-based reconstruction in
 * `detail/vtk_cells.hpp` rather than a plain fixed-width connectivity slice.
 */

// System includes
#include <string>
#include <unordered_map>
#include <vector>

namespace meshioplusplus {

/**
 * @brief Maps a meshio cell-type name to its VTK cell type id.
 *
 * Inverse of `vtk_to_meshio_type()`. Lazily constructed once (function-local
 * `static`) and returned by `const&`. Covers linear and quadratic standard
 * VTK cells plus the Lagrange (68-74) and Bezier (75-81) high-order
 * families.
 * @return Reference to the process-wide singleton lookup table.
 */
inline const std::unordered_map<std::string, int>& meshio_to_vtk_type() {
    static const std::unordered_map<std::string, int> m = {
        {"empty", 0},
        {"vertex", 1},
        {"line", 3},
        {"triangle", 5},
        {"polygon", 7},
        {"pixel", 8},
        {"quad", 9},
        {"tetra", 10},
        {"hexahedron", 12},
        {"wedge", 13},
        {"pyramid", 14},
        {"penta_prism", 15},
        {"hexa_prism", 16},
        {"line3", 21},
        {"triangle6", 22},
        {"quad8", 23},
        {"tetra10", 24},
        {"hexahedron20", 25},
        {"wedge15", 26},
        {"pyramid13", 27},
        {"quad9", 28},
        {"hexahedron27", 29},
        {"quad6", 30},
        {"wedge12", 31},
        {"wedge18", 32},
        {"hexahedron24", 33},
        {"triangle7", 34},
        {"line4", 35},
        {"polyhedron", 42},
        {"VTK_LAGRANGE_CURVE", 68},
        {"VTK_LAGRANGE_TRIANGLE", 69},
        {"VTK_LAGRANGE_QUADRILATERAL", 70},
        {"VTK_LAGRANGE_TETRAHEDRON", 71},
        {"VTK_LAGRANGE_HEXAHEDRON", 72},
        {"VTK_LAGRANGE_WEDGE", 73},
        {"VTK_LAGRANGE_PYRAMID", 74},
        {"VTK_BEZIER_CURVE", 75},
        {"VTK_BEZIER_TRIANGLE", 76},
        {"VTK_BEZIER_QUADRILATERAL", 77},
        {"VTK_BEZIER_TETRAHEDRON", 78},
        {"VTK_BEZIER_HEXAHEDRON", 79},
        {"VTK_BEZIER_WEDGE", 80},
        {"VTK_BEZIER_PYRAMID", 81},
    };
    return m;
}

/**
 * @brief Node-index permutation applied when writing a meshio cell block's
 * connectivity out in VTK order.
 *
 * Only the linear `"wedge"` differs between the two conventions (meshio/gmsh
 * prism ordering vs. `vtkWedge`'s); every other supported type has identical
 * ordering, signaled by returning an empty vector (callers should treat
 * empty as "no permutation needed", not as an error).
 * @param meshio_type The meshio cell-type name.
 * @return `result[j]` = the meshio-order index to place at VTK-order
 *         position `j`; empty if the ordering is already identical.
 */
inline std::vector<int> meshio_to_vtk_order(const std::string& rMeshioType) {
    if (rMeshioType == "wedge")
        return {0, 2, 1, 3, 5, 4};
    return {};
}

/**
 * @brief Maps a VTK cell type id to a meshio cell-type name.
 *
 * Covers only the subset meshio itself can represent (matches
 * `vtk_to_meshio_type` in `_vtk_common.py`); ids meshio has no equivalent
 * for are simply absent from the map, and callers (e.g.
 * `detail::reconstruct_cells`) must treat a failed lookup as an unsupported
 * cell type. Lazily constructed once (function-local `static`) and returned
 * by `const&`.
 * @return Reference to the process-wide singleton lookup table.
 */
inline const std::unordered_map<int, std::string>& vtk_to_meshio_type() {
    static const std::unordered_map<int, std::string> m = {
        {0, "empty"},
        {1, "vertex"},
        {3, "line"},
        {5, "triangle"},
        {7, "polygon"},
        {8, "pixel"},
        {9, "quad"},
        {10, "tetra"},
        {12, "hexahedron"},
        {13, "wedge"},
        {14, "pyramid"},
        {15, "penta_prism"},
        {16, "hexa_prism"},
        {21, "line3"},
        {22, "triangle6"},
        {23, "quad8"},
        {24, "tetra10"},
        {25, "hexahedron20"},
        {26, "wedge15"},
        {27, "pyramid13"},
        {28, "quad9"},
        {29, "hexahedron27"},
        {30, "quad6"},
        {31, "wedge12"},
        {32, "wedge18"},
        {33, "hexahedron24"},
        {34, "triangle7"},
        {35, "line4"},
        {42, "polyhedron"},
        {68, "VTK_LAGRANGE_CURVE"},
        {69, "VTK_LAGRANGE_TRIANGLE"},
        {70, "VTK_LAGRANGE_QUADRILATERAL"},
        {71, "VTK_LAGRANGE_TETRAHEDRON"},
        {72, "VTK_LAGRANGE_HEXAHEDRON"},
        {73, "VTK_LAGRANGE_WEDGE"},
        {74, "VTK_LAGRANGE_PYRAMID"},
        {75, "VTK_BEZIER_CURVE"},
        {76, "VTK_BEZIER_TRIANGLE"},
        {77, "VTK_BEZIER_QUADRILATERAL"},
        {78, "VTK_BEZIER_TETRAHEDRON"},
        {79, "VTK_BEZIER_HEXAHEDRON"},
        {80, "VTK_BEZIER_WEDGE"},
        {81, "VTK_BEZIER_PYRAMID"},
    };
    return m;
}

/**
 * @brief Inverse of `meshio_to_vtk_order`, applied when reading VTK
 * connectivity back into meshio order.
 *
 * Only the linear wedge (VTK type id 13) differs; its permutation
 * `[0,2,1,3,5,4]` happens to be its own inverse, so the same literal serves
 * both directions. Empty means no permutation needed.
 * @param vtk_type The VTK cell type id being read.
 * @return `result[j]` = the VTK-order index to place at meshio-order
 *         position `j`; empty if the ordering is already identical.
 */
inline std::vector<int> vtk_to_meshio_order(int vtk_type) {
    if (vtk_type == 13)
        return {0, 2, 1, 3, 5, 4};
    return {};
}

/**
 * @brief Whether a meshio cell type has a variable node count per cell in a
 * VTK/VTU connectivity+offsets representation.
 *
 * True for `"polygon"` and every `VTK_LAGRANGE_*` type. These cannot be
 * described by a single fixed nodes-per-cell count, so
 * `detail::reconstruct_cells` (vtk_cells.hpp) reconstructs them from the
 * end-offsets array (grouping same-size runs) instead of slicing a uniform
 * `(num_cells, n)` block.
 * @param meshio_type The meshio cell-type name to test.
 * @return `true` if `meshio_type` needs offsets-based reconstruction.
 */
inline bool is_special_cell(const std::string& rMeshioType) {
    return rMeshioType == "polygon" || rMeshioType.rfind("VTK_LAGRANGE_", 0) == 0;
}

}  // namespace meshioplusplus
