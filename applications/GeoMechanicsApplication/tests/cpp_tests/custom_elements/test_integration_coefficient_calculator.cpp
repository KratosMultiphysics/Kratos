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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(IntegrationCoefficientsCalculator_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto integration_coefficient_calculator = IntegrationCoefficientsCalculator{};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    auto                                             detJs = Vector{ScalarVector{1, 2.0}};

    // Act
    const auto calculated_coefficients = integration_coefficient_calculator.Run<>(integration_points, detJs);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(IntegrationCoefficientsCalculator_CloneReturnsNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto integration_coefficient_calculator = IntegrationCoefficientsCalculator{};

    // Act
    const auto clone_modifier = integration_coefficient_calculator.CloneModifier();

    // Assert
    KRATOS_EXPECT_EQ(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
