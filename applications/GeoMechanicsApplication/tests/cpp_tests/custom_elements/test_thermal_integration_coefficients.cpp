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

#include "custom_elements/thermal_integration_coefficients.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <boost/numeric/ublas/assignment.hpp>
#include <structural_mechanics_application_variables.h>
#include <tests/cpp_tests/test_utilities/element_setup_utilities.h>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThermalIntegrationCoefficients_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto calculator = CalculateIntegrationCoefficients0{
        std::make_unique<IntegrationCoefficientModifierForThermalElement>()};
    // The shape function values for this integration point are 0.2, 0.5 and 0.3 for nodes 1, 2 and 3 respectively
    const Geometry<Node>::IntegrationPointType       integration_point(0.5, 0.3, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    Vector                                           detJs(1);
    detJs <<= 2.0;
    auto p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CROSS_AREA, 0.5);
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    const auto p_geometry   = std::make_shared<Line2D2<Node>>(nodes);
    const auto line_element = Element{1, p_geometry, p_properties};

    // Act
    auto calculated_coefficients = calculator.Run<>(integration_points, detJs, &line_element);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 0.5 (cross area) = 0.5
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 0.5, 1e-5);

    const auto plane_element = ElementSetupUtilities::Create2D3NElement();
    calculated_coefficients  = calculator.Run<>(integration_points, detJs, plane_element.get());

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(ThermalIntegrationCoefficients_ClobeReturnsNotNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const std::unique_ptr<IntegrationCoefficientsCalculator> p_pw_line_integration_coefficients =
        std::make_unique<ThermalIntegrationCoefficients>();

    // Act
    const auto clone = p_pw_line_integration_coefficients->Clone();

    // Assert
    KRATOS_EXPECT_NE(clone, nullptr);
}

} // namespace Kratos::Testing
