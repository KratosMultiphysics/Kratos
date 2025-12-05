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
#include "containers/model.h"
#include "custom_utilities/variables_utilities.hpp"
#include "includes/expect.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestVariablesUtilitiesGetsCorrectComponents)
{
    const auto& component = VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "X");

    EXPECT_EQ(component, ACCELERATION_X);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       TestVariablesUtilitiesThrowsWhenComponentDoesNotExist){KRATOS_EXPECT_EXCEPTION_IS_THROWN(
    VariablesUtilities::GetComponentFromVectorVariable(ACCELERATION.Name(), "?"),
    "Error: The component \"ACCELERATION_?\" is not registered!")}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, GeometryUtilities_ReturnsGetNodalValuesOf)
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
    KRATOS_EXPECT_VECTOR_EQ(temperatures, std::vector({4.5, 4.5}));
}

} // namespace Kratos::Testing
