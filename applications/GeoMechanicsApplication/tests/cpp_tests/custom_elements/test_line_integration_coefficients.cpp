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

#include "custom_elements/integration_coefficient_modifier_for_line_element.h"
#include "geometries/line_2d_2.h"
#include "includes/element.h"
#include "includes/properties.h"
#include "structural_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, LineIntegrationCoefficientsCalculator_ReturnsCorrectValue)
{
    // Set
    const auto line_integration_coefficient_calculator =
        IntegrationCoefficientsCalculator{std::make_unique<IntegrationCoefficientModifierForLineElement>()};
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
        line_integration_coefficient_calculator.Run<>(integration_points, detJs, &element);

    // Assert
    // The expected number is calculated as follows:
    // 2.0 (detJ) * 0.5 (weight) * 0.5 (cross area) = 0.5
    EXPECT_EQ(calculated_coefficients.size(), 1);
    EXPECT_NEAR(calculated_coefficients[0], 0.5, Defaults::absolute_tolerance);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, LineIntegrationCoefficientsCalculator_CloneReturnsNotNullptr)
{
    // Set
    const auto line_integration_coefficient_calculator =
        IntegrationCoefficientsCalculator{std::make_unique<IntegrationCoefficientModifierForLineElement>()};

    // Act
    const auto clone_modifier = line_integration_coefficient_calculator.CloneModifier();

    // Assert
    EXPECT_NE(clone_modifier, nullptr);
}
} // namespace Kratos::Testing
