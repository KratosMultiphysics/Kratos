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
 * @file kratos_names.hpp
 * @brief Kratos Multiphysics entity/geometry name tables mapped to
 * `CellType`, ported from this repo's own MIT `src/meshioplusplus/mdpa/_mdpa.py`
 * (`_kratos_elements_to_meshio_type`, `_kratos_conditions_to_meshio_type`,
 * `_kratos_geometries_to_meshio_type`, and the default
 * `_meshio_to_kratos_element/condition_type` pick tables).
 *
 * Used by the `ModelPart` backend (`model_part.hpp`, `kratos_mesh.hpp`) to
 * resolve Kratos entity names on creation and to pick default Kratos names
 * when converting meshio cell blocks into Elements/Conditions, and by the
 * templated bridge (`kratos_bridge.hpp`) to name entities it creates in a
 * real `Kratos::ModelPart`. Backend-independent, header-only.
 */

// System includes
#include <string>
#include <unordered_map>

// Project includes
#include "meshioplusplus/cell_type.hpp"

namespace meshioplusplus {

/**
 * @brief Kratos name -> `CellType` lookup covering element names
 * (`Element3D4N`, `SurfaceElement3D3N`, ...), condition names
 * (`SurfaceCondition3D3N`, ...), geometry names (`Tetrahedra3D4`,
 * `Triangle2D3`, ...), and plain meshio cell-type names (`"tetra"`).
 * @param rName The name to resolve.
 * @return The matching `CellType`, or `CellType::Custom` if unknown.
 */
inline CellType cell_type_from_kratos_name(const std::string& rName) {
    static const std::unordered_map<std::string, CellType> m = {
        // Elements (ported from _kratos_elements_to_meshio_type).
        {"Element2D1N", CellType::Vertex},
        {"Element2D2N", CellType::Line},
        {"Element2D3N", CellType::Triangle},
        {"Element2D6N", CellType::Triangle6},
        {"Element2D4N", CellType::Quad},
        {"Element2D8N", CellType::Quad8},
        {"Element2D9N", CellType::Quad9},
        {"Element3D1N", CellType::Vertex},
        {"Element3D2N", CellType::Line},
        {"Element3D3N", CellType::Triangle},
        {"Element3D4N", CellType::Tetra},
        {"Element3D5N", CellType::Pyramid},
        {"Element3D6N", CellType::Wedge},
        {"Element3D8N", CellType::Hexahedron},
        {"Element3D10N", CellType::Tetra10},
        {"Element3D15N", CellType::Wedge15},
        {"Element3D20N", CellType::Hexahedron20},
        {"Element3D27N", CellType::Hexahedron27},
        {"PointElement2D1N", CellType::Vertex},
        {"PointElement3D1N", CellType::Vertex},
        {"LineElement2D2N", CellType::Line},
        {"LineElement2D3N", CellType::Line3},
        {"LineElement3D2N", CellType::Line},
        {"LineElement3D3N", CellType::Line3},
        {"SurfaceElement3D3N", CellType::Triangle},
        {"SurfaceElement3D6N", CellType::Triangle6},
        {"SurfaceElement3D4N", CellType::Quad},
        {"SurfaceElement3D8N", CellType::Quad8},
        {"SurfaceElement3D9N", CellType::Quad9},
        // Conditions (ported from _kratos_conditions_to_meshio_type).
        {"PointCondition2D1N", CellType::Vertex},
        {"PointCondition3D1N", CellType::Vertex},
        {"LineCondition2D2N", CellType::Line},
        {"LineCondition2D3N", CellType::Line3},
        {"LineCondition3D2N", CellType::Line},
        {"LineCondition3D3N", CellType::Line3},
        {"SurfaceCondition3D3N", CellType::Triangle},
        {"SurfaceCondition3D6N", CellType::Triangle6},
        {"SurfaceCondition3D4N", CellType::Quad},
        {"SurfaceCondition3D8N", CellType::Quad8},
        {"SurfaceCondition3D9N", CellType::Quad9},
        {"PrismCondition2D4N", CellType::Quad},
        {"PrismCondition3D6N", CellType::Wedge},
        // Geometries (ported from _kratos_geometries_to_meshio_type).
        {"Point2D", CellType::Vertex},
        {"Point3D", CellType::Vertex},
        {"Line2D2", CellType::Line},
        {"Line3D2", CellType::Line},
        {"Line2D3", CellType::Line3},
        {"Line3D3", CellType::Line3},
        {"Triangle2D3", CellType::Triangle},
        {"Triangle3D3", CellType::Triangle},
        {"Triangle2D6", CellType::Triangle6},
        {"Triangle3D6", CellType::Triangle6},
        {"Quadrilateral2D4", CellType::Quad},
        {"Quadrilateral3D4", CellType::Quad},
        {"Quadrilateral2D8", CellType::Quad8},
        {"Quadrilateral3D8", CellType::Quad8},
        {"Quadrilateral2D9", CellType::Quad9},
        {"Quadrilateral3D9", CellType::Quad9},
        {"Tetrahedra3D4", CellType::Tetra},
        {"Tetrahedra3D10", CellType::Tetra10},
        {"Prism3D6", CellType::Wedge},
        {"Prism3D15", CellType::Wedge15},
        {"Pyramid3D5", CellType::Pyramid},
        {"Pyramid3D13", CellType::Pyramid13},
        {"Hexahedra3D8", CellType::Hexahedron},
        {"Hexahedra3D20", CellType::Hexahedron20},
        {"Hexahedra3D27", CellType::Hexahedron27},
    };
    auto it = m.find(rName);
    if (it != m.end())
        return it->second;
    return cell_type_from_name(rName);  // plain meshio names; Custom if unknown
}

/**
 * @brief Default Kratos *element* name for a cell type (ported from
 * `_meshio_to_kratos_element_type`), falling back to the Kratos geometry
 * name, then the meshio name itself for types with no Kratos equivalent.
 * @param type The cell type.
 * @return The Kratos name (resolvable back via `cell_type_from_kratos_name`).
 */
inline const std::string& kratos_element_name(CellType type) {
    static const std::unordered_map<CellType, std::string> m = {
        {CellType::Vertex, "Element3D1N"},
        {CellType::Line, "Element3D2N"},
        {CellType::Triangle, "Element3D3N"},
        {CellType::Tetra, "Element3D4N"},
        {CellType::Pyramid, "Element3D5N"},
        {CellType::Wedge, "Element3D6N"},
        {CellType::Hexahedron, "Element3D8N"},
        {CellType::Line3, "LineElement3D3N"},
        {CellType::Triangle6, "Element2D6N"},
        {CellType::Quad, "Element2D4N"},
        {CellType::Quad8, "Element2D8N"},
        {CellType::Quad9, "Element2D9N"},
        {CellType::Tetra10, "Element3D10N"},
        {CellType::Hexahedron20, "Element3D20N"},
        {CellType::Hexahedron27, "Element3D27N"},
        // No Element* default in Kratos conventions -> geometry names.
        {CellType::Wedge15, "Prism3D15"},
        {CellType::Pyramid13, "Pyramid3D13"},
    };
    auto it = m.find(type);
    if (it != m.end())
        return it->second;
    return cell_type_name(type);  // meshio name; ResolveEntityType understands it
}

/**
 * @brief Default Kratos *condition* name for a cell type (ported from
 * `_meshio_to_kratos_condition_type`), with the same fallbacks as
 * `kratos_element_name`.
 * @param type The cell type.
 * @return The Kratos name (resolvable back via `cell_type_from_kratos_name`).
 */
inline const std::string& kratos_condition_name(CellType type) {
    static const std::unordered_map<CellType, std::string> m = {
        {CellType::Vertex, "PointCondition3D1N"},      {CellType::Line, "LineCondition3D2N"},
        {CellType::Line3, "LineCondition3D3N"},        {CellType::Triangle, "SurfaceCondition3D3N"},
        {CellType::Triangle6, "SurfaceCondition3D6N"}, {CellType::Quad, "SurfaceCondition3D4N"},
        {CellType::Quad8, "SurfaceCondition3D8N"},     {CellType::Quad9, "SurfaceCondition3D9N"},
        {CellType::Wedge, "PrismCondition3D6N"},
    };
    auto it = m.find(type);
    if (it != m.end())
        return it->second;
    return kratos_element_name(type);  // geometry/meshio-name fallback chain
}

}  // namespace meshioplusplus
