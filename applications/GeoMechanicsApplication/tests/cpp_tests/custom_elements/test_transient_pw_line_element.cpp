// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Gennady Markelov
//

#include "containers/model.h"
#include "custom_elements/transient_Pw_line_element.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(TransientPwLineElement_GetIntegrationMethodForAllRegisteredElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const std::vector<CalculationContribution> contributions;
    PointerVector<Node>                        nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));

    // Act and Assert
    auto p_transient_pw_line_element_2D2N = make_intrusive<TransientPwLineElement<2, 2>>(
        1, std::make_shared<Line2D2<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_pw_line_element_3D2N = make_intrusive<TransientPwLineElement<3, 2>>(
        1, std::make_shared<Line3D2<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_3D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    auto p_transient_pw_line_element_2D3N = make_intrusive<TransientPwLineElement<2, 3>>(
        1, std::make_shared<Line2D3<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    auto p_transient_pw_line_element_3D3N = make_intrusive<TransientPwLineElement<3, 3>>(
        1, std::make_shared<Line3D3<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_3D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    auto p_transient_pw_line_element_2D4N = make_intrusive<TransientPwLineElement<2, 4>>(
        1, std::make_shared<Line2D4<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_3);

    nodes.push_back(make_intrusive<Node>(5, 1.0, 0.5, 0.0));
    auto p_transient_pw_line_element_2D5N = make_intrusive<TransientPwLineElement<2, 5>>(
        1, std::make_shared<Line2D5<Node>>(nodes), contributions,
        std::make_unique<IntegrationCoefficientModifierForPwLineElement>());
    KRATOS_EXPECT_EQ(p_transient_pw_line_element_2D5N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_5);
}

} // namespace Kratos::Testing
