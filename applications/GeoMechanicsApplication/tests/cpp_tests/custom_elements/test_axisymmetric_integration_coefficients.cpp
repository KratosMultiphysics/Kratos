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
//                   Marjan Fathian
//                   Gennady Markelov
//

#include "custom_elements/integration_coefficient_modifier_for_axisymmetric_element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricIntegrationCoefficientsCalculator_ReturnsCorrectValue,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto axisymmetric_integration_coefficient_calculator = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForAxisymmetricElement>()};
    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    const Geometry<Node>::IntegrationPointType       integration_point(0.5, 0.3, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    Vector                                           detJs(1);
    detJs <<= 2.0;
    const auto p_element = ElementSetupUtilities::Create2D3NElement();

    // Act
    const auto calculated_coefficients = axisymmetric_integration_coefficient_calculator.Run<>(
        integration_points, detJs, p_element.get());

    // Assert
    // The expected number is calculated as follows:
    // 2.0 * pi * 0.8 (radius) * 2.0 (detJ) * 0.5 (weight) = 5.02655
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 5.02655, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricIntegrationCoefficientsCalculator_CloneReturnsNotNullptr,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto axisymmetric_integration_coefficient_calculator = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForAxisymmetricElement>()};

    // Act
    const auto clone_modifier = axisymmetric_integration_coefficient_calculator.CloneModifier();

    // Assert
    KRATOS_EXPECT_NE(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
