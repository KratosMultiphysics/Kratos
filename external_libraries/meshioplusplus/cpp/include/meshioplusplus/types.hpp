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
 * @file types.hpp
 * @brief Cell-type metadata tables, ported 1:1 from the Python reference so
 * C++ and Python agree on a single definition.
 *
 * `num_nodes_per_cell()` is ported from `src/meshio/_common.py` and
 * `topological_dimension()` from `src/meshio/_mesh.py`. Both are keyed by
 * meshio's own cell-type name strings (e.g. `"triangle"`, `"tetra10"`,
 * `"hexahedron20"`) rather than any per-format native name — each format
 * module maps its own names to/from these before consulting these tables.
 * See <https://github.com/nschloe/meshio/wiki/Node-ordering-in-cells> for
 * the node-ordering convention these types assume.
 */

// System includes
#include <string>
#include <unordered_map>

namespace meshioplusplus {

/**
 * @brief Table mapping a meshio cell-type name to its fixed node count.
 *
 * Lazily constructed once (function-local `static`) and returned by
 * `const&`; only rectangular, fixed-node-count cell types appear here — cell
 * types whose node count varies per cell (`"polygon"`, the VTK_LAGRANGE_*
 * family) are represented via `CellBlock`'s ragged storage instead and are
 * intentionally absent. Covers meshio's linear through high-order elements
 * (e.g. `"line"` through `"line11"`, `"tetra"` through `"tetra286"`).
 * @return Reference to the process-wide singleton lookup table.
 */
inline const std::unordered_map<std::string, int>& num_nodes_per_cell() {
    static const std::unordered_map<std::string, int> m = {
        {"vertex", 1},
        {"line", 2},
        {"triangle", 3},
        {"quad", 4},
        {"quad8", 8},
        {"tetra", 4},
        {"hexahedron", 8},
        {"hexahedron20", 20},
        {"hexahedron24", 24},
        {"wedge", 6},
        {"pyramid", 5},
        //
        {"line3", 3},
        {"triangle6", 6},
        {"quad9", 9},
        {"tetra10", 10},
        {"hexahedron27", 27},
        {"wedge15", 15},
        {"wedge18", 18},
        {"pyramid13", 13},
        {"pyramid14", 14},
        //
        {"line4", 4},
        {"triangle10", 10},
        {"quad16", 16},
        {"tetra20", 20},
        {"wedge40", 40},
        {"hexahedron64", 64},
        //
        {"line5", 5},
        {"triangle15", 15},
        {"quad25", 25},
        {"tetra35", 35},
        {"wedge75", 75},
        {"hexahedron125", 125},
        //
        {"line6", 6},
        {"triangle21", 21},
        {"quad36", 36},
        {"tetra56", 56},
        {"wedge126", 126},
        {"hexahedron216", 216},
        //
        {"line7", 7},
        {"triangle28", 28},
        {"quad49", 49},
        {"tetra84", 84},
        {"wedge196", 196},
        {"hexahedron343", 343},
        //
        {"line8", 8},
        {"triangle36", 36},
        {"quad64", 64},
        {"tetra120", 120},
        {"wedge288", 288},
        {"hexahedron512", 512},
        //
        {"line9", 9},
        {"triangle45", 45},
        {"quad81", 81},
        {"tetra165", 165},
        {"wedge405", 405},
        {"hexahedron729", 729},
        //
        {"line10", 10},
        {"triangle55", 55},
        {"quad100", 100},
        {"tetra220", 220},
        {"wedge550", 550},
        {"hexahedron1000", 1000},
        {"hexahedron1331", 1331},
        //
        {"line11", 11},
        {"triangle66", 66},
        {"quad121", 121},
        {"tetra286", 286},
    };
    return m;
}

/**
 * @brief Table mapping a meshio cell-type name to its topological dimension
 * (0 = vertex, 1 = line/curve, 2 = surface, 3 = volume).
 *
 * Lazily constructed once (function-local `static`) and returned by
 * `const&`. Includes the standard meshio types plus the VTK Lagrange
 * high-order family (`"VTK_LAGRANGE_CURVE"`, `..._TRIANGLE`,
 * `..._QUADRILATERAL`, `..._TETRAHEDRON`, `..._HEXAHEDRON`, `..._WEDGE`,
 * `..._PYRAMID`), which carry a variable node count per cell (see
 * `vtk_common.hpp`'s `is_special_cell`) but still have a fixed dimension.
 * @return Reference to the process-wide singleton lookup table.
 */
inline const std::unordered_map<std::string, int>& topological_dimension() {
    static const std::unordered_map<std::string, int> m = {
        {"line", 1},
        {"polygon", 2},
        {"triangle", 2},
        {"quad", 2},
        {"tetra", 3},
        {"hexahedron", 3},
        {"wedge", 3},
        {"pyramid", 3},
        {"line3", 1},
        {"triangle6", 2},
        {"quad9", 2},
        {"tetra10", 3},
        {"hexahedron27", 3},
        {"wedge18", 3},
        {"pyramid14", 3},
        {"vertex", 0},
        {"quad8", 2},
        {"hexahedron20", 3},
        {"triangle10", 2},
        {"triangle15", 2},
        {"triangle21", 2},
        {"line4", 1},
        {"line5", 1},
        {"line6", 1},
        {"tetra20", 3},
        {"tetra35", 3},
        {"tetra56", 3},
        {"quad16", 2},
        {"quad25", 2},
        {"quad36", 2},
        {"triangle28", 2},
        {"triangle36", 2},
        {"triangle45", 2},
        {"triangle55", 2},
        {"triangle66", 2},
        {"quad49", 2},
        {"quad64", 2},
        {"quad81", 2},
        {"quad100", 2},
        {"quad121", 2},
        {"line7", 1},
        {"line8", 1},
        {"line9", 1},
        {"line10", 1},
        {"line11", 1},
        {"tetra84", 3},
        {"tetra120", 3},
        {"tetra165", 3},
        {"tetra220", 3},
        {"tetra286", 3},
        {"wedge40", 3},
        {"wedge75", 3},
        {"hexahedron64", 3},
        {"hexahedron125", 3},
        {"hexahedron216", 3},
        {"hexahedron343", 3},
        {"hexahedron512", 3},
        {"hexahedron729", 3},
        {"hexahedron1000", 3},
        {"wedge126", 3},
        {"wedge196", 3},
        {"wedge288", 3},
        {"wedge405", 3},
        {"wedge550", 3},
        {"VTK_LAGRANGE_CURVE", 1},
        {"VTK_LAGRANGE_TRIANGLE", 2},
        {"VTK_LAGRANGE_QUADRILATERAL", 2},
        {"VTK_LAGRANGE_TETRAHEDRON", 3},
        {"VTK_LAGRANGE_HEXAHEDRON", 3},
        {"VTK_LAGRANGE_WEDGE", 3},
        {"VTK_LAGRANGE_PYRAMID", 3},
    };
    return m;
}

}  // namespace meshioplusplus
