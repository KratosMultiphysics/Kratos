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

#include "custom_elements/integration_coefficient_modifier_for_pw_line_element.h"
#include "structural_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(PwLineIntegrationCoefficientsCalculator_ReturnsCorrectValue, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto pw_line_integration_coefficient_calculator = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>()};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    const auto                                       detJs        = Vector{ScalarVector{1, 2.0}};
    auto                                             p_properties = std::make_shared<Properties>();
    p_properties->SetValue(CROSS_AREA, 0.5);
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element    = Element{1, p_geometry, p_properties};

    // Act
    const auto calculated_coefficients =
        pw_line_integration_coefficient_calculator.Run<>(integration_points, detJs, &element);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 0.5 (cross area) = 0.5
    KRATOS_EXPECT_EQ(calculated_coefficients.size(), 1);
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 0.5, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(PwLineIntegrationCoefficientsCalculator_CloneReturnsNotNullptr,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Set
    const auto pw_line_integration_coefficient_calculator = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>()};

    // Act
    const auto clone_modifier = pw_line_integration_coefficient_calculator.CloneModifier();

    // Assert
    KRATOS_EXPECT_NE(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
