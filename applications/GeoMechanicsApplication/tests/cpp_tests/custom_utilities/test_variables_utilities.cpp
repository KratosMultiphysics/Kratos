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
#include "custom_utilities/variables_utilities.hpp"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TestVariablesUtilitiesGetsCorrectComponents, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto& component = VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "X");

    KRATOS_EXPECT_EQ(component, ACCELERATION_X);
}

KRATOS_TEST_CASE_IN_SUITE(TestVariablesUtilitiesThrowsWhenComponentDoesNotExist,
                          KratosGeoMechanicsFastSuiteWithoutKernel){KRATOS_EXPECT_EXCEPTION_IS_THROWN(
    VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "?"),
    "Error: The component \"ACCELERATION_?\" is not registered!")}

KRATOS_TEST_CASE_IN_SUITE(GeometryUtilities_ReturnsGetNodalValues, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.CreateNewNode(11, 0.0, 0.0, 0.0));
    nodes.push_back(r_model_part.CreateNewNode(30, 0.0, 0.0, 0.0));
    const Geometry geometry(1, nodes);

    for (auto& r_node : geometry) {
        r_node.FastGetSolutionStepValue(TEMPERATURE) = 4.5;
    }

    // Act
    const auto temperatures = VariablesUtilities::GetNodalValuesOf<2>(TEMPERATURE, geometry);

    // Assert
    auto expected_temperatures = std::vector({4.5, 4.5});
    KRATOS_EXPECT_VECTOR_EQ(temperatures, expected_temperatures);

    // Act
    const auto more_temperatures = VariablesUtilities::GetNodalValues(r_model_part.Nodes(), TEMPERATURE);

    // Assert
    KRATOS_EXPECT_VECTOR_EQ(temperatures, expected_temperatures);
}

} // namespace Kratos::Testing
