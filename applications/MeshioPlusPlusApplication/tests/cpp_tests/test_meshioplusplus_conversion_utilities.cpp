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
#include <cstdint>

// External includes
#include "meshioplusplus/cell_type.hpp"
#include "meshioplusplus/ndarray.hpp"

// Project includes
#include "containers/model.h"
#include "geometries/geometry_data.h"
#include "includes/variables.h"
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

/***********************************************************************************/
/***********************************************************************************/

KRATOS_TEST_CASE_IN_SUITE(MeshioPlusPlusConversionUtilitiesMeshToModelPartAcceptsInt64Arrays, KratosMeshioPlusPlusFastSuite)
{
    // ApplyDataArrayToEntities used to reject every non-Float64 array outright, dropping
    // meshio++'s own Int64 bookkeeping arrays (refine:level, refine:parent_cell, ...)
    // regardless of name. No *real* operation currently produces an Int64 array under a name
    // that also happens to be a registered Kratos Variable (those arrays are all
    // colon-namespaced, e.g. "refine:level", which no Kratos Variable ever is), so this is
    // exercised directly against the mesh rather than through a real operation: staging a
    // synthetic Int64 point_data array under "DOMAIN_SIZE" (a registered Variable<int>, also
    // used this way in test_meshioplusplus_io.cpp).
    Model model;
    auto& r_source = model.CreateModelPart("source");
    r_source.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_source.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_source.CreateNewNode(3, 0.0, 1.0, 0.0);

    // Mirrors Internals::ModelPartToMeshWithData's own sequence: fill the topology, stage
    // data directly on the model part view, then InvalidateBlocks() once at the end.
    meshioplusplus::Mesh mesh;
    Internals::FillMeshioModelPart(r_source, mesh.GetModelPart(), true, true, false);
    meshioplusplus::NDArray domain_size_array =
        meshioplusplus::NDArray::Uninit(meshioplusplus::DType::Int64, {3});
    std::int64_t* p_values = domain_size_array.As<std::int64_t>();
    p_values[0] = 2;
    p_values[1] = 3;
    p_values[2] = 3;
    mesh.GetModelPart().SetNodalData("DOMAIN_SIZE", domain_size_array);
    mesh.InvalidateBlocks();

    auto& r_destination = model.CreateModelPart("destination");
    Internals::MeshToModelPart(mesh, r_destination);

    KRATOS_EXPECT_EQ(r_destination.GetNode(1).GetValue(DOMAIN_SIZE), 2);
    KRATOS_EXPECT_EQ(r_destination.GetNode(2).GetValue(DOMAIN_SIZE), 3);
    KRATOS_EXPECT_EQ(r_destination.GetNode(3).GetValue(DOMAIN_SIZE), 3);
}

} // namespace Kratos::Testing
