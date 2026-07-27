//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include "meshioplusplus/cell_type.hpp"

// Project includes
#include "geometries/geometry_data.h"
#include "custom_utilities/meshioplusplus_conversion_utilities.h"
#include "meshioplusplus_fast_suite.h"

namespace Kratos::Testing {

namespace {
namespace mio = meshioplusplus;
using KGT = GeometryData::KratosGeometryType;
} // namespace

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusConversionUtilitiesCellTypeMapping, KratosMeshioPlusPlusFastSuite)
{
    // Every (Kratos geometry type -> meshio++ cell type) pair the conversion supports,
    // exercised directly rather than only through whichever two or three element types the
    // other tests happen to build a ModelPart with.
    const std::vector<std::pair<KGT, mio::CellType>> expected_mapping = {
        {KGT::Kratos_Point2D, mio::CellType::Vertex},
        {KGT::Kratos_Point3D, mio::CellType::Vertex},
        {KGT::Kratos_Line2D2, mio::CellType::Line},
        {KGT::Kratos_Line3D2, mio::CellType::Line},
        {KGT::Kratos_Line2D3, mio::CellType::Line3},
        {KGT::Kratos_Line3D3, mio::CellType::Line3},
        {KGT::Kratos_Line2D4, mio::CellType::Line4},
        {KGT::Kratos_Line2D5, mio::CellType::Line5},
        {KGT::Kratos_Triangle2D3, mio::CellType::Triangle},
        {KGT::Kratos_Triangle3D3, mio::CellType::Triangle},
        {KGT::Kratos_Triangle2D6, mio::CellType::Triangle6},
        {KGT::Kratos_Triangle3D6, mio::CellType::Triangle6},
        {KGT::Kratos_Triangle2D10, mio::CellType::Triangle10},
        {KGT::Kratos_Triangle2D15, mio::CellType::Triangle15},
        {KGT::Kratos_Quadrilateral2D4, mio::CellType::Quad},
        {KGT::Kratos_Quadrilateral3D4, mio::CellType::Quad},
        {KGT::Kratos_Quadrilateral2D8, mio::CellType::Quad8},
        {KGT::Kratos_Quadrilateral3D8, mio::CellType::Quad8},
        {KGT::Kratos_Quadrilateral2D9, mio::CellType::Quad9},
        {KGT::Kratos_Quadrilateral3D9, mio::CellType::Quad9},
        {KGT::Kratos_Tetrahedra3D4, mio::CellType::Tetra},
        {KGT::Kratos_Tetrahedra3D10, mio::CellType::Tetra10},
        {KGT::Kratos_Prism3D6, mio::CellType::Wedge},
        {KGT::Kratos_Prism3D15, mio::CellType::Wedge15},
        {KGT::Kratos_Pyramid3D5, mio::CellType::Pyramid},
        {KGT::Kratos_Pyramid3D13, mio::CellType::Pyramid13},
        {KGT::Kratos_Hexahedra3D8, mio::CellType::Hexahedron},
        {KGT::Kratos_Hexahedra3D20, mio::CellType::Hexahedron20},
        {KGT::Kratos_Hexahedra3D27, mio::CellType::Hexahedron27},
    };

    for (const auto& r_case : expected_mapping) {
        KRATOS_EXPECT_EQ(static_cast<int>(Internals::MeshioCellTypeFromKratosGeometry(r_case.first)),
                         static_cast<int>(r_case.second));
    }
}

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusConversionUtilitiesCellTypeMappingThrowsOnUnsupported, KratosMeshioPlusPlusFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        Internals::MeshioCellTypeFromKratosGeometry(GeometryData::KratosGeometryType::Kratos_Sphere3D1),
        "is not supported by MeshioPlusPlusIO");
}

} // namespace Kratos::Testing
