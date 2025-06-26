// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_processes/fix_water_pressures_above_phreatic_line.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_FixesAllWaterPressuresToZeroWhenAllNodesAbovePhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, -1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [-2.0],
                "move_mesh": false
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();

    for (auto& r_node : r_model_part.Nodes()) {
        EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
        EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_DoesNothingWhenAllNodesBelowPhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, -1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [1.0],
                "move_mesh": false
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 1.0);
    EXPECT_FALSE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_OnlyFixesNodesAbovePhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, -1.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [-0.5],
                "move_mesh": false
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_FreesNodesWhenTheyGetBelowPhreaticLine_WithMoveMesh,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, -2.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [-0.5],
                "move_mesh": true
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));

    r_model_part.GetNode(1).Coordinates()[1] = -1.0; // move node 1 below phreatic line
    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_FALSE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_FreesNodesWhenTheyGetBelowPhreaticLine_WithoutMoveMesh,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 0.0, -2.0, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [-0.5],
                "move_mesh": false
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));

    // Move node 1 below phreatic line using total displacements (without moving the mesh, since move_mesh is false)
    r_model_part.GetNode(1).FastGetSolutionStepValue(TOTAL_DISPLACEMENT)[1] = -1.0;
    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_FALSE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_InterpolatesMultiLineCorrectly,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(TOTAL_DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.4, -0.4, 0.0);
    r_model_part.CreateNewNode(2, 0.6, -0.6, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    auto                                      test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [-1.0, 0.0],
                "move_mesh": false
            }  )"};
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    process.ExecuteInitializeSolutionStep();
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

} // namespace Kratos::Testing