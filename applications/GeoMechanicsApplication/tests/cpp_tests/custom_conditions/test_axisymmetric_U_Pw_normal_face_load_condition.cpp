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
#include "custom_conditions/axisymmetric_U_Pw_normal_face_load_condition.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <geometries/line_2d_4.h>
#include <geometries/line_2d_5.h>

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(AxisymmetricUPwNormalFaceLoadCondition_GetIntegrationMethodForAllRegisteredElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto          p_properties = std::make_shared<Properties>();
    PointerVector<Node> nodes;
    nodes.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    nodes.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));

    // Act and Assert
    auto p_axisymmetric_U_Pw_normal_face_load_condition_2D2N =
        make_intrusive<AxisymmetricUPwNormalFaceLoadCondition<2, 2>>(
            1, std::make_shared<Line2D2<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_axisymmetric_U_Pw_normal_face_load_condition_2D2N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    auto p_axisymmetric_U_Pw_normal_face_load_condition_2D3N =
        make_intrusive<AxisymmetricUPwNormalFaceLoadCondition<2, 3>>(
            1, std::make_shared<Line2D3<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_axisymmetric_U_Pw_normal_face_load_condition_2D3N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_2);

    nodes.push_back(make_intrusive<Node>(4, 0.5, 0.0, 0.0));
    auto p_axisymmetric_U_Pw_normal_face_load_condition_2D4N =
        make_intrusive<AxisymmetricUPwNormalFaceLoadCondition<2, 4>>(
            1, std::make_shared<Line2D4<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_axisymmetric_U_Pw_normal_face_load_condition_2D4N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_3);

    nodes.push_back(make_intrusive<Node>(5, 1.0, 0.5, 0.0));
    auto p_axisymmetric_U_Pw_normal_face_load_condition_2D5N =
        make_intrusive<AxisymmetricUPwNormalFaceLoadCondition<2, 5>>(
            1, std::make_shared<Line2D5<Node>>(nodes), p_properties);
    KRATOS_EXPECT_EQ(p_axisymmetric_U_Pw_normal_face_load_condition_2D5N->GetIntegrationMethod(),
                     GeometryData::IntegrationMethod::GI_GAUSS_5);
}

} // namespace Kratos::Testing
