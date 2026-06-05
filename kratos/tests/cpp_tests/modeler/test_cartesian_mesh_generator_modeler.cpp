//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/element.h"
#include "includes/kratos_components.h"
#include "modeler/cartesian_mesh_generator_modeler.h"

namespace Kratos::Testing
{

namespace {

/// Build a square 2x2 bounding-box model part with 8 corner nodes only (no elements),
/// so the node-based BB fallback is exercised.
ModelPart& CreateBoxSourceModelPart(Model& rModel, const std::string& rName,
    double xMin, double yMin, double zMin,
    double xMax, double yMax, double zMax)
{
    auto& r_mp = rModel.CreateModelPart(rName);
    r_mp.CreateNewNode(1, xMin, yMin, zMin);
    r_mp.CreateNewNode(2, xMax, yMin, zMin);
    r_mp.CreateNewNode(3, xMax, yMax, zMin);
    r_mp.CreateNewNode(4, xMin, yMax, zMin);
    r_mp.CreateNewNode(5, xMin, yMin, zMax);
    r_mp.CreateNewNode(6, xMax, yMin, zMax);
    r_mp.CreateNewNode(7, xMax, yMax, zMax);
    r_mp.CreateNewNode(8, xMin, yMax, zMax);
    return r_mp;
}

} // anonymous namespace

/***********************************************************************************/
/***********************************************************************************/

/// Test that GenerateMesh (direct API) fills a 3D model part with the correct
/// node and element counts for a unit-cube source with a given element size.
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModeler3DNodeCount, KratosCoreFastSuite)
{
    Model model;
    // 1x1x1 box, element size 0.5 → 2 segments per direction → 3^3=27 nodes, 2^3*6=48 tets
    auto& r_source = CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    auto& r_output = model.CreateModelPart("output");
    r_output.CreateNewProperties(0);

    CartesianMeshGeneratorModeler modeler(r_source, 0.5);
    modeler.GenerateMesh(r_output, KratosComponents<Element>::Get("Element3D4N"));

    // 3 nodes per direction → 27 total
    KRATOS_EXPECT_EQ(r_output.NumberOfNodes(), 27u);
    // 2 cells per direction → 8 hex cells × 6 tets = 48 tetrahedra
    KRATOS_EXPECT_EQ(r_output.NumberOfElements(), 48u);
}

/***********************************************************************************/
/***********************************************************************************/

/// Test that every generated node lies inside (or on the boundary of) the source BB.
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModeler3DNodesInsideBB, KratosCoreFastSuite)
{
    Model model;
    auto& r_source = CreateBoxSourceModelPart(model, "source", -1.0, -1.0, -1.0, 1.0, 1.0, 1.0);
    auto& r_output = model.CreateModelPart("output");
    r_output.CreateNewProperties(0);

    CartesianMeshGeneratorModeler modeler(r_source, 0.5);
    modeler.GenerateMesh(r_output, KratosComponents<Element>::Get("Element3D4N"));

    for (const auto& r_node : r_output.Nodes()) {
        KRATOS_EXPECT_GE(r_node.X(), -1.0 - 1e-12);
        KRATOS_EXPECT_LE(r_node.X(),  1.0 + 1e-12);
        KRATOS_EXPECT_GE(r_node.Y(), -1.0 - 1e-12);
        KRATOS_EXPECT_LE(r_node.Y(),  1.0 + 1e-12);
        KRATOS_EXPECT_GE(r_node.Z(), -1.0 - 1e-12);
        KRATOS_EXPECT_LE(r_node.Z(),  1.0 + 1e-12);
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// Test that each element has exactly 4 nodes (tetrahedron) and that all node IDs
/// reference nodes that are present in the output model part.
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModeler3DElementConnectivity, KratosCoreFastSuite)
{
    Model model;
    auto& r_source = CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 1.0, 1.0, 1.0);
    auto& r_output = model.CreateModelPart("output");
    r_output.CreateNewProperties(0);

    CartesianMeshGeneratorModeler modeler(r_source, 1.0);
    modeler.GenerateMesh(r_output, KratosComponents<Element>::Get("Element3D4N"));

    // 1 cell per direction, 6 tets
    KRATOS_EXPECT_EQ(r_output.NumberOfElements(), 6u);

    for (const auto& r_elem : r_output.Elements()) {
        KRATOS_EXPECT_EQ(r_elem.GetGeometry().size(), 4u);
        for (unsigned int n = 0; n < 4; ++n) {
            const auto node_id = r_elem.GetGeometry()[n].Id();
            KRATOS_EXPECT_TRUE(r_output.HasNode(node_id));
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

/// Test the factory / registry constructor (Model + Parameters) via SetupModelPart().
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModelerSetupModelPart, KratosCoreFastSuite)
{
    Model model;
    CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 2.0, 2.0, 2.0);

    Parameters params(R"({
        "input_model_part_name"  : "source",
        "output_model_part_name" : "output",
        "element_name"           : "Element3D4N",
        "element_size"           : 1.0
    })");

    CartesianMeshGeneratorModeler modeler(model, params);
    modeler.SetupModelPart();

    const auto& r_output = model.GetModelPart("output");
    // 2x2x2 box, element size 1 → 3^3=27 nodes, 2^3*6=48 tets
    KRATOS_EXPECT_EQ(r_output.NumberOfNodes(), 27u);
    KRATOS_EXPECT_EQ(r_output.NumberOfElements(), 48u);
}

/***********************************************************************************/
/***********************************************************************************/

/// Test that the bounding box helper correctly handles a node-only model part
/// (no elements), which is the typical output of StlIO (conditions only).
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModelerBoundingBoxNodeOnly, KratosCoreFastSuite)
{
    Model model;
    auto& r_source = CreateBoxSourceModelPart(model, "source", -2.0, -3.0, -4.0, 2.0, 3.0, 4.0);
    auto& r_output = model.CreateModelPart("output");
    r_output.CreateNewProperties(0);

    CartesianMeshGeneratorModeler modeler(r_source, 2.0);
    Point min_pt, max_pt;
    modeler.CalculateBoundingBox(r_source, min_pt, max_pt);

    KRATOS_EXPECT_NEAR(min_pt.X(), -2.0, 1e-12);
    KRATOS_EXPECT_NEAR(min_pt.Y(), -3.0, 1e-12);
    KRATOS_EXPECT_NEAR(min_pt.Z(), -4.0, 1e-12);
    KRATOS_EXPECT_NEAR(max_pt.X(),  2.0, 1e-12);
    KRATOS_EXPECT_NEAR(max_pt.Y(),  3.0, 1e-12);
    KRATOS_EXPECT_NEAR(max_pt.Z(),  4.0, 1e-12);
}

/***********************************************************************************/
/***********************************************************************************/

/// Test CalculateDivisionNumbers: exact multiples produce the right segment count.
KRATOS_TEST_CASE_IN_SUITE(CartesianMeshGeneratorModelerDivisionNumbers, KratosCoreFastSuite)
{
    Model model;
    // 3x2x1 box, element size 0.5 → 6, 4, 2 segments
    auto& r_source = CreateBoxSourceModelPart(model, "source", 0.0, 0.0, 0.0, 3.0, 2.0, 1.0);
    auto& r_output = model.CreateModelPart("output");
    r_output.CreateNewProperties(0);

    CartesianMeshGeneratorModeler modeler(r_source, 0.5);
    modeler.GenerateMesh(r_output, KratosComponents<Element>::Get("Element3D4N"));

    // (6+1)*(4+1)*(2+1) = 7*5*3 = 105 nodes
    KRATOS_EXPECT_EQ(r_output.NumberOfNodes(), 105u);
    // 6*4*2*6 = 288 tetrahedra
    KRATOS_EXPECT_EQ(r_output.NumberOfElements(), 288u);
}

} // namespace Kratos::Testing
