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
#include "custom_geometries/line_interface_geometry.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficients_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto interface_integration_coefficients = IntegrationCoefficientsCalculator{};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    Vector                                           detJs(1);
    detJs <<= 2.0;

    // Act
    const auto calculated_coefficients = interface_integration_coefficients.Run<>(integration_points, detJs);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficients_CloneReturnsNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto interface_integration_coefficients = IntegrationCoefficientsCalculator{};

    // Act
    const auto clone_modifier = interface_integration_coefficients.CloneModifier();

    // Assert
    KRATOS_EXPECT_EQ(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
