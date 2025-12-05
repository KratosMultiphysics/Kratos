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
//                   Gennady Markelov
//

#include "custom_elements/integration_coefficients_calculator.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, IntegrationCoefficientsCalculatorWithoutModifier_ReturnsCorrectValue)
{
    // Set
    const auto integration_coefficient_calculator = IntegrationCoefficientsCalculator{};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    const auto                                       detJs = Vector{ScalarVector{1, 2.0}};

    // Act
    const auto calculated_coefficients = integration_coefficient_calculator.Run<>(integration_points, detJs);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    EXPECT_EQ(calculated_coefficients.size(), 1);
    EXPECT_NEAR(calculated_coefficients[0], 1.0, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, IntegrationCoefficientsCalculatorWithoutModifier_CloneReturnsNullptr)
{
    // Set
    const auto integration_coefficient_calculator = IntegrationCoefficientsCalculator{};

    // Act
    const auto clone_modifier = integration_coefficient_calculator.CloneModifier();

    // Assert
    EXPECT_EQ(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
