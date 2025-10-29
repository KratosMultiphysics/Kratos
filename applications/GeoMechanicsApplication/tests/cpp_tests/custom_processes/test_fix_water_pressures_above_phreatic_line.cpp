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

namespace
{
using namespace Kratos;

Parameters CreateParametersWithConstantPhreaticLineAtHeight(double height)
{
    return Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0],
                "y_coordinates": [)" +
                      std::to_string(height) + R"(]
            }  )"};
}

ModelPart& CreateModelPartWithTwoNodesAtHeights(Model& rModel, double y_node_1, double y_node_2)
{
    auto& r_result = rModel.CreateModelPart("foo");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.CreateNewNode(1, 0.0, y_node_1, 0.0);
    r_result.CreateNewNode(2, 0.0, y_node_2, 0.0);

    for (auto& r_node : r_result.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    return r_result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_ThrowsUponConstructionWhenXOrYCoordinatesAreEmpty,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -2.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);

    const auto test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [],
                "y_coordinates": []
            }  )"};

    // Act & Assert
    const std::string expected_error_message = "The x_coordinates and y_coordinates "
                                               "of the phreatic line must be non-empty vectors.";
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        FixWaterPressuresAbovePhreaticLineProcess(r_model_part, test_parameters), expected_error_message);
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_ThrowsUponConstructionWhenXYCoordinatesAreNotEqualLength,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -2.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);

    const auto test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [1.0],
                "y_coordinates": [2.0, 3.0]
            }  )"};

    // Act & Assert
    const std::string expected_error_message = "The lengths of the x_coordinates and y_coordinates "
                                               "of the phreatic line must be equal.";
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        FixWaterPressuresAbovePhreaticLineProcess(r_model_part, test_parameters), expected_error_message);
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_FixesAllWaterPressuresToZeroWhenAllNodesAboveOrOnPhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -2.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);
    const auto     test_parameters = CreateParametersWithConstantPhreaticLineAtHeight(y_node_2);
    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    for (auto& r_node : r_model_part.Nodes()) {
        EXPECT_DOUBLE_EQ(r_node.FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
        EXPECT_TRUE(r_node.IsFixed(WATER_PRESSURE));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_DoesNothingWhenAllNodesBelowPhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -1.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);
    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    const auto test_parameters = CreateParametersWithConstantPhreaticLineAtHeight(1.0);
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 1.0);
    EXPECT_FALSE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_OnlyFixesWaterPressureDoFsAbovePhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -1.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    const auto test_parameters = CreateParametersWithConstantPhreaticLineAtHeight(-0.5);
    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_FreesWaterPressureDoFsWhenTheyAreBelowPhreaticLine,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    constexpr auto y_node_1     = 0.0;
    constexpr auto y_node_2     = -2.0;
    auto           model        = Model{};
    auto&          r_model_part = CreateModelPartWithTwoNodesAtHeights(model, y_node_1, y_node_2);

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    const auto test_parameters = CreateParametersWithConstantPhreaticLineAtHeight(-0.5);

    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));

    // Act
    // Move node 1 below phreatic line
    r_model_part.GetNode(1).Y() = -1.0;
    process.ExecuteInitializeSolutionStep();

    // Assert
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_FALSE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(TestFixWaterPressureAbovePhreaticLine_InterpolatesMultiLineCorrectly,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("foo");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.CreateNewNode(1, 0.4, -0.4, 0.0);
    r_model_part.CreateNewNode(2, 0.6, -0.6, 0.0);

    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }

    r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;
    r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;

    const auto test_parameters = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [-1.0, 0.0]
            }  )"};

    FixWaterPressuresAbovePhreaticLineProcess process(r_model_part, test_parameters);

    // Act
    process.ExecuteInitializeSolutionStep();

    // Assert
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(1).FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    EXPECT_TRUE(r_model_part.GetNode(1).IsFixed(WATER_PRESSURE));
    EXPECT_DOUBLE_EQ(r_model_part.GetNode(2).FastGetSolutionStepValue(WATER_PRESSURE), 2.0);
    EXPECT_FALSE(r_model_part.GetNode(2).IsFixed(WATER_PRESSURE));
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoFixWaterPressureAbovePhreaticLine, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto       model              = Model{};
    auto&      r_empty_model_part = model.CreateModelPart("foo");
    const auto test_parameters    = Parameters{R"(
            {
                "model_part_name": "foo",
                "x_coordinates": [0.0, 1.0],
                "y_coordinates": [-1.0, 0.0]
            }  )"};
    const FixWaterPressuresAbovePhreaticLineProcess process(r_empty_model_part, test_parameters);

    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "FixWaterPressuresAbovePhreaticLineProcess");
}
} // namespace Kratos::Testing