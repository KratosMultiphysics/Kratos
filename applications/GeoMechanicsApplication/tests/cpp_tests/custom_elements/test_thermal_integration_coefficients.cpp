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
#include "structural_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"

#include <boost/numeric/ublas/assignment.hpp>

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ThermalIntegrationCoefficients_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto thermal_integration_coefficients = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForThermalElement>()};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
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

    // Act and Assert
    auto calculated_coefficients =
        thermal_integration_coefficients.Run<>(integration_points, detJs, &line_element);

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 0.5 (cross area) = 0.5
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 0.5, 1e-5);

    nodes.push_back(make_intrusive<Node>(2, 1.0, 1.0, 0.0));
    const auto plane_element = ElementSetupUtilities::Create2D3NElement(nodes, p_properties);
    calculated_coefficients =
        thermal_integration_coefficients.Run<>(integration_points, detJs, plane_element.get());

    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) = 1.0
    // cross area is not taken into account
    KRATOS_EXPECT_NEAR(plane_element.get()->GetProperties()[CROSS_AREA], 0.5, 1e-5);
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 1.0, 1e-5);
}

KRATOS_TEST_CASE_IN_SUITE(ThermalIntegrationCoefficients_CloneReturnsNotNullptr, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto thermal_integration_coefficients = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForThermalElement>()};

    // Act
    const auto clone_modifier = thermal_integration_coefficients.CloneModifier();

    // Assert
    KRATOS_EXPECT_NE(clone_modifier, nullptr);
}

} // namespace Kratos::Testing
