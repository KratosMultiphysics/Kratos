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

#include "custom_elements/pw_line_integration_coefficients.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(PwLineIntegrationCoefficients_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const std::unique_ptr<IntegrationCoefficientsCalculator> p_pw_line_integration_coefficients =
        std::make_unique<PwLineIntegrationCoefficients>();
    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    const Geometry<Node>::IntegrationPointType       integration_point(0.5, 0.3, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    Vector                                           detJs(1);
    detJs <<= 2.0;
    const auto cross_area = 0.5;
    // Act
    const auto calculated_coefficients = p_pw_line_integration_coefficients->CalculateIntegrationCoefficients(
        integration_points, detJs, cross_area);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 0.5 (cross area) = 0.5
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 0.5, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(PwLineIntegrationCoefficients_ClobeReturnsNotNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const std::unique_ptr<IntegrationCoefficientsCalculator> p_pw_line_integration_coefficients =
        std::make_unique<PwLineIntegrationCoefficients>();

    // Act
    const auto clone = p_pw_line_integration_coefficients->Clone();

    // Assert
    KRATOS_EXPECT_NE(clone, nullptr);
}
} // namespace Kratos::Testing
