// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "containers/model.h"
#include "custom_processes/apply_constant_interpolate_line_pressure_process.h"
#include "includes/model_part.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "testing/testing.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include "includes/smart_pointers.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

ModelPart& CreateTestModelPart(Model& rModel)
{
    auto& result = rModel.CreateModelPart("TestPart");
    result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_properties = result.CreateNewProperties(0);

    // Create nodes directly in the model part
    result.CreateNewNode(1, 0.0, 0.0, 0.0);
    result.CreateNewNode(2, 1.0, 0.0, 0.0);
    result.CreateNewNode(3, 1.0, 1.0, 0.0);
    result.CreateNewNode(4, 0.0, 1.0, 0.0);

    // Create elements using the node IDs
    std::vector<ModelPart::IndexType> elem1_nodes = {1, 2, 3};
    std::vector<ModelPart::IndexType> elem2_nodes = {1, 3, 4};
    result.CreateNewElement("Element2D3N", 1, elem1_nodes, p_properties);
    result.CreateNewElement("Element2D3N", 2, elem2_nodes, p_properties);

    return result;
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_Construction, KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    // Act & Assert
    EXPECT_NO_THROW(ApplyConstantInterpolateLinePressureProcess process(r_model_part, params));
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_ThrowsOnInvalidDirections, KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 1,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        ApplyConstantInterpolateLinePressureProcess process(r_model_part, params),
        "Gravity direction cannot be the same as Out-of-Plane directions");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_ExecuteInitializeSolutionStep,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    // Set initial pressure values for boundary nodes
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = -10.0 * static_cast<double>(r_node.Id());
    }

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 100.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    for (auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), -10.0 * r_node.Id());
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_SeepageBranch, KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    // Set up so that CalculatePressure returns a value less than the cut-off for node 1,
    // and a value greater than or equal to the cut-off for node 2.
    // We'll do this by setting the initial WATER_PRESSURE and using a high cut-off.

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    }

    // Set a high cut-off so all nodes will be below the cut-off
    const auto params_below = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": true,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 100.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process_below(r_model_part, params_below);

    // Act and Assert
    process_below.ExecuteInitializeSolutionStep();

    // All nodes should be set and fixed (since pressure < cut-off)
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
        // The value is set by CalculatePressure, which in this test setup is 0.0
        KRATOS_EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    }

    // Now set a low cut-off so all nodes will be above the cut-off
    const auto params_above = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": true,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": -100.0,
        "table": 1
    })");

    // Reset nodes to unfixed
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Free(WATER_PRESSURE);
    }

    ApplyConstantInterpolateLinePressureProcess process_above(r_model_part, params_above);
    process_above.ExecuteInitializeSolutionStep();

    // All nodes should be free (since pressure >= cut-off)
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.IsFixed(WATER_PRESSURE))
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_Info, KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": false,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 0.0,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act & Assert
    KRATOS_EXPECT_EQ(process.Info(), "ApplyConstantInterpolateLinePressureProcess");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_DoesNotFreeWhenIsFixedIsNotProvided,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
        r_node.Fix(WATER_PRESSURE);
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = -10.0 * static_cast<double>(r_node.Id());
    }

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 1.0e9,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE))
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_FreesWhenIsFixedIsExplicitlyFalse,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateTestModelPart(model);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
        r_node.Fix(WATER_PRESSURE);
        r_node.FastGetSolutionStepValue(WATER_PRESSURE) = -10.0 * static_cast<double>(r_node.Id());
    }

    const auto params = Parameters(R"({
        "model_part_name": "TestPart",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": false,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 1.0e9,
        "table": 1
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    for (const auto& r_node : r_model_part.Nodes()) {
        KRATOS_EXPECT_FALSE(r_node.IsFixed(WATER_PRESSURE))
    }
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_ExecuteInitializeSolutionStepRunsOnlyOnce,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: boundary nodes on top (y=10) and bottom (y=0), plus an interior node (y=5).
    // CalculatePressure interpolates the interior node to 50.0, which is intentionally
    // different from the sentinel value set after the first call.  This allows reliable
    // detection of any accidental re-execution of the process.
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Top boundary (y = 10)
    r_model_part.CreateNewNode(1, 0.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    r_model_part.CreateNewNode(2, 10.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;

    // Bottom boundary (y = 0)
    r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    r_model_part.CreateNewNode(4, 10.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;

    // Interior node — not part of any element, so not a boundary node.
    // CalculatePressure will interpolate it to 50.0.
    auto interior_node = r_model_part.CreateNewNode(5, 5.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{3, 4}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "is_fixed": true,
        "is_seepage": false,
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "pressure_tension_cut_off": 1.0e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act: first call sets the interior node to the interpolated value (-50.0).
    process.ExecuteInitializeSolutionStep();
    constexpr auto tolerance                   = Defaults::absolute_tolerance * 100.0;
    constexpr auto expected_interpolated_value = -50;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE),
                       expected_interpolated_value, tolerance);

    // Reset the interior node to a sentinel that clearly differs from -50.0.
    constexpr auto sentinel = -99999.0;
    interior_node->Free(WATER_PRESSURE);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = sentinel;

    // Act: second call must be a no-op (guarded by mIsInitialized).
    process.ExecuteInitializeSolutionStep();

    // Assert: interior node retains the sentinel value, not the interpolated -50.0.
    KRATOS_EXPECT_FALSE(interior_node->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), sentinel, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_InterpolatesInteriorFromTopAndBottomBoundaries,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");

    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // TOP boundary (y = 10)
    auto top_node_1 = r_model_part.CreateNewNode(1, 0.0, 10.0, 0.0);
    auto top_node_2 = r_model_part.CreateNewNode(2, 5.0, 10.0, 0.0);
    auto top_node_3 = r_model_part.CreateNewNode(3, 10.0, 10.0, 0.0);

    top_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    top_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = -50.0;
    top_node_3->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;

    // BOTTOM boundary (y = 0)
    auto bottom_node_1 = r_model_part.CreateNewNode(4, 0.0, 0.0, 0.0);
    auto bottom_node_2 = r_model_part.CreateNewNode(5, 5.0, 0.0, 0.0);
    auto bottom_node_3 = r_model_part.CreateNewNode(6, 10.0, 0.0, 0.0);

    bottom_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    bottom_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = -150.0;
    bottom_node_3->FastGetSolutionStepValue(WATER_PRESSURE) = -200.0;

    auto interior_node = r_model_part.CreateNewNode(7, 7.5, 5.0, 0.0);

    // Minimal elements to detect boundary nodes
    Properties::Pointer p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{2, 3}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 3, std::vector<ModelPart::IndexType>{4, 5}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 4, std::vector<ModelPart::IndexType>{5, 6}, p_props);

    const auto params = Parameters(R"(
    {
        "model_part_name" : "Main",
        "variable_name"   : "WATER_PRESSURE",
        "gravity_direction" : 1,
        "out_of_plane_direction" : 2,
        "is_fixed" : false,
        "pressure_tension_cut_off" : 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_top = -75.0;  // from top boundary
    constexpr auto expected_bot = -175.0; // from bottom boundary
    constexpr auto expected_value = (expected_top - expected_bot) / (10.0 - 0.0) * (5.0 - 0.0) + expected_bot;

    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value,
                       Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_NoBoundaryNodes, KratosGeoMechanicsFastSuite)
{
    // Arrange
    Model      model;
    ModelPart& mp = model.CreateModelPart("Main");
    mp.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Create nodes but NO elements → no boundary detection
    mp.CreateNewNode(1, 0.0, 0.0, 0.0);
    mp.CreateNewNode(2, 1.0, 0.0, 0.0);

    const auto params = Parameters(R"(
    {
        "model_part_name" : "Main",
        "variable_name"   : "WATER_PRESSURE"
    })");

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(ApplyConstantInterpolateLinePressureProcess(mp, params),
                                      "No boundary node is found");
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_ExtrapolatesWhenNodeIsRightOfAllBoundaryNodes,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: boundary nodes only at x=[0,5]; interior node at x=10 lies beyond the right edge.
    // CalculateBoundaryPressure will take the "only-left" branch and call
    // InterpolateBoundaryPressureWithOneContainer, which extrapolates the pressure.
    //
    // Hand-calculated expected pressure for interior node (x=10, y=5):
    //   Top extrapolation  (from x=0,P=0 and x=5,P=50):  slope=10 → P_top = 10*(10-5)+50 = 100
    //   Coord_top = y(node2) = 10
    //   Bottom extrapolation (from x=0,P=100 and x=5,P=150): slope=10 → P_bot = 10*(10-5)+150 = 200
    //   Coord_bot = 0
    //   Final: slope_p = (100-200)/(10-0) = -10;  P = -10*(5-0)+200 = 150
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    r_model_part.CreateNewNode(1, 0.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    r_model_part.CreateNewNode(2, 5.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -50.0;
    r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    r_model_part.CreateNewNode(4, 5.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -150.0;

    auto interior_node = r_model_part.CreateNewNode(5, 10.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{3, 4}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_value = -150.0;
    constexpr auto tolerance      = Defaults::absolute_tolerance * 100.0;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_ExtrapolatesWhenNodeIsLeftOfAllBoundaryNodes,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: boundary nodes only at x=[5,10]; interior node at x=0 lies beyond the left edge.
    // CalculateBoundaryPressure will take the "only-right" branch and call
    // InterpolateBoundaryPressureWithOneContainer, which extrapolates the pressure.
    //
    // Hand-calculated expected pressure for interior node (x=0, y=5):
    //   Top extrapolation  (from x=5,P=50 and x=10,P=100): slope=10 → P_top = 10*(0-5)+50 = 0
    //   Coord_top = 10
    //   Bottom extrapolation (from x=5,P=150 and x=10,P=200): slope=10 → P_bot = 10*(0-5)+150 = 100
    //   Coord_bot = 0
    //   Final: slope_p = (0-100)/(10-0) = -10;  P = -10*(5-0)+100 = 50
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    r_model_part.CreateNewNode(1, 5.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -50.0;
    r_model_part.CreateNewNode(2, 10.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    r_model_part.CreateNewNode(3, 5.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -150.0;
    r_model_part.CreateNewNode(4, 10.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -200.0;

    auto interior_node = r_model_part.CreateNewNode(5, 0.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{3, 4}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_value = -50.0;
    constexpr auto tolerance      = Defaults::absolute_tolerance * 100.0;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_NodeAtSameXAsBoundaryNode, KratosGeoMechanicsFastSuite)
{
    // Arrange: interior node shares the same x-coordinate as the leftmost boundary node.
    // FindClosestNodeOnBoundaryNodes returns the same node for both the "left" and "right"
    // candidate sets, so InterpolateBoundaryPressure takes its else-branch (horizontal
    // distance == 0) and uses the boundary node's pressure directly.
    //
    // Hand-calculated expected pressure for interior node (x=0, y=5):
    //   Top "left" and "right" both resolve to node1(x=0,P=0).
    //     Else-branch → P_top=0, Coord_top=10.
    //   Bottom "left" and "right" both resolve to node3(x=0,P=100).
    //     Else-branch → P_bot=100, Coord_bot=0.
    //   Final: slope_p = (0-100)/(10-0) = -10;  P = -10*(5-0)+100 = 50
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    r_model_part.CreateNewNode(1, 0.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
    r_model_part.CreateNewNode(2, 5.0, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -50.0;
    r_model_part.CreateNewNode(3, 0.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    r_model_part.CreateNewNode(4, 5.0, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -150.0;

    // Interior node at x=0: same horizontal coordinate as node1 and node3.
    auto interior_node = r_model_part.CreateNewNode(5, 0.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{3, 4}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_value = -50.0;
    constexpr auto tolerance      = Defaults::absolute_tolerance * 100.0;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_OneContainerVerticalFallback_WhenHorizontalDifferenceIsTiny,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: force one-container interpolation by placing the interior node to the right of all
    // boundary nodes with negative tiny x-values. This avoids cancellation in the closest-node
    // search while keeping the boundary-node horizontal delta <= TINY so
    // InterpolateBoundaryPressureWithOneContainer takes its vertical fallback branch:
    //   rPressure   = pressureLeft;
    //   rCoordinate = CoordinatesLeft[mGravityDirection];
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // GeoMechanics TINY is 1e-60, so use a smaller spacing to force the fallback branch.
    constexpr auto dx = 1.0e-80;

    // Top boundary (y=10): both nodes have equal pressure, making the expected value
    // independent of FoundNodes ordering inside the fallback branch.
    r_model_part.CreateNewNode(1, -2.0 * dx, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -120.0;
    r_model_part.CreateNewNode(2, -dx, 10.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -120.0;

    // Bottom boundary (y=0): same rationale as top boundary.
    r_model_part.CreateNewNode(3, -2.0 * dx, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -220.0;
    r_model_part.CreateNewNode(4, -dx, 0.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -220.0;

    auto interior_node = r_model_part.CreateNewNode(5, 0.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);
    r_model_part.CreateNewElement("Element2D2N", 2, std::vector<ModelPart::IndexType>{3, 4}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert: top value is 120 at y=10 and bottom value is 220 at y=0.
    // Linear interpolation at y=5 gives 170.
    constexpr auto expected_value = -170.0;
    constexpr auto tolerance      = Defaults::absolute_tolerance * 100.0;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressure_FlatBoundary_ReturnsBoundaryPressure,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: all boundary nodes are on a single horizontal line (y=5).
    // Both FindTopBoundaryNodes and FindBottomBoundaryNodes return the same set, so
    // CoordinateTop == CoordinateBottom.  CalculatePressure takes its else-branch
    // (|CoordinateTop - CoordinateBottom| <= TINY) and returns PressureBottom directly.
    //
    // Hand-calculated expected pressure for interior node (x=5, y=5):
    //   Top and bottom both interpolate between node1(x=0,P=100) and node2(x=10,P=200):
    //     slope=10, P=10*(5-0)+100=150, Coord=5.
    //   Else-branch: P_bottom=150 is returned directly.
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    r_model_part.CreateNewNode(1, 0.0, 5.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    r_model_part.CreateNewNode(2, 10.0, 5.0, 0.0)->FastGetSolutionStepValue(WATER_PRESSURE) = -200.0;

    // Interior node on the same horizontal level as the boundary.
    auto interior_node = r_model_part.CreateNewNode(3, 5.0, 5.0, 0.0);
    interior_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    constexpr auto expected_value = -150.0;
    constexpr auto tolerance      = Defaults::absolute_tolerance * 100.0;
    KRATOS_EXPECT_NEAR(interior_node->FastGetSolutionStepValue(WATER_PRESSURE), expected_value, tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_FillListOfBoundaryNodesFast_DynamicAllocation,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: Create a node shared by 11 elements to trigger dynamic push_back in FillListOfBoundaryNodesFast.
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    // Create 12 nodes in a line
    for (std::size_t i = 1; i <= 12; ++i) {
        r_model_part.CreateNewNode(i, static_cast<double>(i), 0.0, 0.0);
    }
    auto p_props = r_model_part.CreateNewProperties(0);
    // Create 11 elements, all sharing node 1
    for (std::size_t i = 2; i <= 12; ++i) {
        r_model_part.CreateNewElement(
            "Element2D2N", i - 1,
            std::vector<ModelPart::IndexType>{static_cast<ModelPart::IndexType>(1),
                                              static_cast<ModelPart::IndexType>(i)},
            p_props);
    }

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    // Act & Assert: Should not throw, and should cover the dynamic allocation branch
    EXPECT_NO_THROW(ApplyConstantInterpolateLinePressureProcess process(r_model_part, params));
}

KRATOS_TEST_CASE_IN_SUITE(ApplyConstantInterpolateLinePressureProcess_CalculateBoundaryPressure_ErrorBranch,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange: Build valid boundary nodes, plus one isolated node above them to make top boundary search empty.
    Model      model;
    ModelPart& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);

    auto p_node_1 = r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    auto p_node_2 = r_model_part.CreateNewNode(2, 0.0, 1.0, 0.0);
    auto p_node_3 = r_model_part.CreateNewNode(3, 0.0, 2.0, 0.0); // isolated node

    p_node_1->FastGetSolutionStepValue(WATER_PRESSURE) = -100.0;
    p_node_2->FastGetSolutionStepValue(WATER_PRESSURE) = -50.0;
    p_node_3->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;

    auto p_props = r_model_part.CreateNewProperties(0);
    r_model_part.CreateNewElement("Element2D2N", 1, std::vector<ModelPart::IndexType>{1, 2}, p_props);

    const auto params = Parameters(R"({
        "model_part_name": "Main",
        "variable_name": "WATER_PRESSURE",
        "gravity_direction": 1,
        "out_of_plane_direction": 2,
        "is_fixed": false,
        "pressure_tension_cut_off": 1e9
    })");

    ApplyConstantInterpolateLinePressureProcess process(r_model_part, params);

    // Act & Assert: The isolated node triggers the CalculateBoundaryPressure error path.
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.ExecuteInitializeSolutionStep(),
                                      "There is not enough points around interpolation, node Id");
}

} // namespace