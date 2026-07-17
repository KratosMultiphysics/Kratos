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
 * @file cell_type.hpp
 * @brief `CellType`: a compact enum for meshio cell-type names, with
 * name/node-count/dimension lookup tables.
 *
 * The format layer identifies cell types by meshio's name strings
 * (`"triangle"`, `"tetra10"`, ...; see `types.hpp`). The NATIVE and KRATOS
 * mesh backends store cell types as this enum instead — an integer compare
 * beats a string compare in per-block hot paths, and the KRATOS backend's
 * geometry-name tables (`backends/kratos_names.hpp`) key off it. The enum
 * covers every fixed-node-count type in `num_nodes_per_cell()` plus the
 * variable-node-count families (`Polygon`, `Polyhedron`, the VTK Lagrange
 * types); anything else maps to `CellType::Custom` and keeps its name
 * out-of-band (see `NativeCellBlock::mTypeName`).
 *
 * The single source of truth is the `MESHIOPLUSPLUS_CELL_TYPES` X-macro
 * below: `(EnumName, "meshio name", nodes-per-cell or -1 if variable,
 * topological dimension)`. The tables in `types.hpp` remain the reference
 * the entries were transcribed from.
 */

// System includes
#include <cstdint>
#include <string>
#include <unordered_map>

namespace meshioplusplus {

// X(EnumName, MeshioName, NumNodes /* -1 = variable */, TopologicalDim)
#define MESHIOPLUSPLUS_CELL_TYPES(X)                                 \
    X(Vertex, "vertex", 1, 0)                                        \
    X(Line, "line", 2, 1)                                            \
    X(Line3, "line3", 3, 1)                                          \
    X(Line4, "line4", 4, 1)                                          \
    X(Line5, "line5", 5, 1)                                          \
    X(Line6, "line6", 6, 1)                                          \
    X(Line7, "line7", 7, 1)                                          \
    X(Line8, "line8", 8, 1)                                          \
    X(Line9, "line9", 9, 1)                                          \
    X(Line10, "line10", 10, 1)                                       \
    X(Line11, "line11", 11, 1)                                       \
    X(Triangle, "triangle", 3, 2)                                    \
    X(Triangle6, "triangle6", 6, 2)                                  \
    X(Triangle10, "triangle10", 10, 2)                               \
    X(Triangle15, "triangle15", 15, 2)                               \
    X(Triangle21, "triangle21", 21, 2)                               \
    X(Triangle28, "triangle28", 28, 2)                               \
    X(Triangle36, "triangle36", 36, 2)                               \
    X(Triangle45, "triangle45", 45, 2)                               \
    X(Triangle55, "triangle55", 55, 2)                               \
    X(Triangle66, "triangle66", 66, 2)                               \
    X(Quad, "quad", 4, 2)                                            \
    X(Quad8, "quad8", 8, 2)                                          \
    X(Quad9, "quad9", 9, 2)                                          \
    X(Quad16, "quad16", 16, 2)                                       \
    X(Quad25, "quad25", 25, 2)                                       \
    X(Quad36, "quad36", 36, 2)                                       \
    X(Quad49, "quad49", 49, 2)                                       \
    X(Quad64, "quad64", 64, 2)                                       \
    X(Quad81, "quad81", 81, 2)                                       \
    X(Quad100, "quad100", 100, 2)                                    \
    X(Quad121, "quad121", 121, 2)                                    \
    X(Tetra, "tetra", 4, 3)                                          \
    X(Tetra10, "tetra10", 10, 3)                                     \
    X(Tetra20, "tetra20", 20, 3)                                     \
    X(Tetra35, "tetra35", 35, 3)                                     \
    X(Tetra56, "tetra56", 56, 3)                                     \
    X(Tetra84, "tetra84", 84, 3)                                     \
    X(Tetra120, "tetra120", 120, 3)                                  \
    X(Tetra165, "tetra165", 165, 3)                                  \
    X(Tetra220, "tetra220", 220, 3)                                  \
    X(Tetra286, "tetra286", 286, 3)                                  \
    X(Hexahedron, "hexahedron", 8, 3)                                \
    X(Hexahedron20, "hexahedron20", 20, 3)                           \
    X(Hexahedron24, "hexahedron24", 24, 3)                           \
    X(Hexahedron27, "hexahedron27", 27, 3)                           \
    X(Hexahedron64, "hexahedron64", 64, 3)                           \
    X(Hexahedron125, "hexahedron125", 125, 3)                        \
    X(Hexahedron216, "hexahedron216", 216, 3)                        \
    X(Hexahedron343, "hexahedron343", 343, 3)                        \
    X(Hexahedron512, "hexahedron512", 512, 3)                        \
    X(Hexahedron729, "hexahedron729", 729, 3)                        \
    X(Hexahedron1000, "hexahedron1000", 1000, 3)                     \
    X(Hexahedron1331, "hexahedron1331", 1331, 3)                     \
    X(Wedge, "wedge", 6, 3)                                          \
    X(Wedge15, "wedge15", 15, 3)                                     \
    X(Wedge18, "wedge18", 18, 3)                                     \
    X(Wedge40, "wedge40", 40, 3)                                     \
    X(Wedge75, "wedge75", 75, 3)                                     \
    X(Wedge126, "wedge126", 126, 3)                                  \
    X(Wedge196, "wedge196", 196, 3)                                  \
    X(Wedge288, "wedge288", 288, 3)                                  \
    X(Wedge405, "wedge405", 405, 3)                                  \
    X(Wedge550, "wedge550", 550, 3)                                  \
    X(Pyramid, "pyramid", 5, 3)                                      \
    X(Pyramid13, "pyramid13", 13, 3)                                 \
    X(Pyramid14, "pyramid14", 14, 3)                                 \
    X(Polygon, "polygon", -1, 2)                                     \
    X(Polyhedron, "polyhedron", -1, 3)                               \
    X(VtkLagrangeCurve, "VTK_LAGRANGE_CURVE", -1, 1)                 \
    X(VtkLagrangeTriangle, "VTK_LAGRANGE_TRIANGLE", -1, 2)           \
    X(VtkLagrangeQuadrilateral, "VTK_LAGRANGE_QUADRILATERAL", -1, 2) \
    X(VtkLagrangeTetrahedron, "VTK_LAGRANGE_TETRAHEDRON", -1, 3)     \
    X(VtkLagrangeHexahedron, "VTK_LAGRANGE_HEXAHEDRON", -1, 3)       \
    X(VtkLagrangeWedge, "VTK_LAGRANGE_WEDGE", -1, 3)                 \
    X(VtkLagrangePyramid, "VTK_LAGRANGE_PYRAMID", -1, 3)

/**
 * @brief Compact identifier for a meshio cell type.
 *
 * `Custom` is the catch-all for names not in the table (parameterized types
 * like `"polyhedron12"` keep their exact spelling out-of-band alongside the
 * enum value).
 */
enum class CellType : std::uint16_t {
#define MESHIOPLUSPLUS_CELL_TYPE_ENUM(Name, Str, N, Dim) Name,
    MESHIOPLUSPLUS_CELL_TYPES(MESHIOPLUSPLUS_CELL_TYPE_ENUM)
#undef MESHIOPLUSPLUS_CELL_TYPE_ENUM
        Custom,
};

/**
 * @brief The meshio name for a `CellType` (e.g. `CellType::Tetra10` →
 * `"tetra10"`).
 * @param type The cell type to convert; `Custom` yields `""` (the caller is
 *             expected to carry the real name out-of-band).
 * @return Reference to the process-wide name string.
 */
inline const std::string& cell_type_name(CellType type) {
    static const std::string names[] = {
#define MESHIOPLUSPLUS_CELL_TYPE_NAME(Name, Str, N, Dim) Str,
        MESHIOPLUSPLUS_CELL_TYPES(MESHIOPLUSPLUS_CELL_TYPE_NAME)
#undef MESHIOPLUSPLUS_CELL_TYPE_NAME
            "",  // Custom
    };
    return names[static_cast<std::size_t>(type)];
}

/**
 * @brief The `CellType` for a meshio cell-type name; `CellType::Custom` for
 * anything not in the table.
 * @param rName The meshio cell-type name (e.g. `"triangle"`).
 * @return The matching enum value, or `Custom`.
 */
inline CellType cell_type_from_name(const std::string& rName) {
    static const std::unordered_map<std::string, CellType> m = {
#define MESHIOPLUSPLUS_CELL_TYPE_LOOKUP(Name, Str, N, Dim) {Str, CellType::Name},
        MESHIOPLUSPLUS_CELL_TYPES(MESHIOPLUSPLUS_CELL_TYPE_LOOKUP)
#undef MESHIOPLUSPLUS_CELL_TYPE_LOOKUP
    };
    auto it = m.find(rName);
    return it == m.end() ? CellType::Custom : it->second;
}

/**
 * @brief Fixed nodes-per-cell of a `CellType`, or -1 for variable-node-count
 * types (`Polygon`, `Polyhedron`, the VTK Lagrange family) and `Custom`.
 * @param type The cell type to query.
 * @return The node count, or -1.
 */
inline int cell_type_num_nodes(CellType type) {
    static const int counts[] = {
#define MESHIOPLUSPLUS_CELL_TYPE_NODES(Name, Str, N, Dim) N,
        MESHIOPLUSPLUS_CELL_TYPES(MESHIOPLUSPLUS_CELL_TYPE_NODES)
#undef MESHIOPLUSPLUS_CELL_TYPE_NODES
            - 1,  // Custom
    };
    return counts[static_cast<std::size_t>(type)];
}

/**
 * @brief Topological dimension (0 = vertex, 1 = curve, 2 = surface,
 * 3 = volume) of a `CellType`, or -1 for `Custom`.
 * @param type The cell type to query.
 * @return The dimension, or -1.
 */
inline int cell_type_dimension(CellType type) {
    static const int dims[] = {
#define MESHIOPLUSPLUS_CELL_TYPE_DIM(Name, Str, N, Dim) Dim,
        MESHIOPLUSPLUS_CELL_TYPES(MESHIOPLUSPLUS_CELL_TYPE_DIM)
#undef MESHIOPLUSPLUS_CELL_TYPE_DIM
            - 1,  // Custom
    };
    return dims[static_cast<std::size_t>(type)];
}

}  // namespace meshioplusplus
