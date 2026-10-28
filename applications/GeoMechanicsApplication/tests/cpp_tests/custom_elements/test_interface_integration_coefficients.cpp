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

#include "custom_elements/integration_coefficient_modifier_for_interface_element.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficientModifier_WithPerimeterProperty_ReturnsScaledCoefficient,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto modifier = IntegrationCoefficientModifierForInterfaceElement{};
    const Geometry<Node>::IntegrationPointType integration_point(0.0, 0.0, 0.0, 0.5);
    const auto                                 integration_coefficient = 2.0;
    auto                                       p_properties = std::make_shared<Properties>();
    p_properties->SetValue(GEO_CONTACT_SURFACE_PERIMETER, 5.0);
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element    = Element{1, p_geometry, p_properties};

    // Act
    const auto modified_coefficient = modifier(integration_coefficient, integration_point, element);

    // Assert
    // Expected: 2.0 * 5.0 = 10.0
    KRATOS_EXPECT_NEAR(modified_coefficient, 10.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficientModifier_WithoutPerimeterProperty_ReturnsUnchangedCoefficient,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto modifier = IntegrationCoefficientModifierForInterfaceElement{};
    const Geometry<Node>::IntegrationPointType integration_point(0.0, 0.0, 0.0, 0.5);
    const double                               integration_coefficient = 2.0;
    auto                                       p_properties = std::make_shared<Properties>();
    PointerVector<Node>                        nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element    = Element{1, p_geometry, p_properties};

    // Act
    const auto modified_coefficient = modifier(integration_coefficient, integration_point, element);

    // Assert
    KRATOS_EXPECT_NEAR(modified_coefficient, 2.0, Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficientModifier_CloneReturnsNotNullptr,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto modifier = IntegrationCoefficientModifierForInterfaceElement{};

    // Act
    const auto cloned_modifier = modifier.Clone();

    // Assert
    KRATOS_EXPECT_NE(cloned_modifier, nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceIntegrationCoefficientCalculator_WithModifier_ReturnsScaledCoefficients,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto interface_integration_coefficient_calculator = IntegrationCoefficientsCalculator{
        std::make_unique<IntegrationCoefficientModifierForInterfaceElement>()};
    const Geometry<Node>::IntegrationPointType       integration_point(0.0, 0.0, 0.0, 0.5);
    const Geometry<Node>::IntegrationPointsArrayType integration_points{integration_point};
    const auto                                       detJs        = Vector{ScalarVector{1, 2.0}};
    auto                                             p_properties = std::make_shared<Properties>();
    p_properties->SetValue(GEO_CONTACT_SURFACE_PERIMETER, 3.0);
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(0, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(1, 1.0, 0.0, 0.0));
    const auto p_geometry = std::make_shared<Line2D2<Node>>(nodes);
    const auto element    = Element{1, p_geometry, p_properties};

    // Act
    const auto calculated_coefficients =
        interface_integration_coefficient_calculator.Run<>(integration_points, detJs, &element);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 3.0 (perimeter) = 3.0
    KRATOS_EXPECT_EQ(calculated_coefficients.size(), 1);
    KRATOS_EXPECT_NEAR(calculated_coefficients[0], 3.0, Defaults::absolute_tolerance);
}

} // namespace Kratos::Testing
