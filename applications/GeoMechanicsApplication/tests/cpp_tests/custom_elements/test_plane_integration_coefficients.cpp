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
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(PlaneIntegrationCoefficients_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto calculate_integration_coefficients = CalculateIntegrationCoefficients0{};
    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    const Geometry<Node>::IntegrationPointType       integration_point(0.5, 0.3, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    Vector                                           detJs(1);
    detJs <<= 2.0;

    // Act
    const auto calculated_coefficients = calculate_integration_coefficients.Run<>(integration_points, detJs);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(PlaneIntegrationCoefficients_ClobeReturnsNotNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto plane_integration_coefficients = CalculateIntegrationCoefficients0{};

    // Act
    const auto clone_modifier = plane_integration_coefficients.CloneModifier();

    // Assert
    KRATOS_EXPECT_EQ(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
