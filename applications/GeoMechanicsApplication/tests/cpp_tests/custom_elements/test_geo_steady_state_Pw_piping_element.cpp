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

#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoSteadyStatePipingElement_DoesNotLosePipingStateAfterInitializingTwice,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Dummy");
    auto  p_node_1     = r_model_part.CreateNewNode(1, 0.0, -0.1, 0.0);
    auto  p_node_2     = r_model_part.CreateNewNode(2, 1.0, -0.1, 0.0);
    auto  p_element    = make_intrusive<GeoSteadyStatePwPipingElement<2, 2>>(
        1, Kratos::make_shared<Line2D2<Node>>(p_node_1, p_node_2));

    p_element->Initialize(r_model_part.GetProcessInfo());
    KRATOS_EXPECT_FALSE(p_element->GetValue(PIPE_ACTIVE));

    p_element->SetValue(PIPE_ACTIVE, true);
    p_element->Initialize(r_model_part.GetProcessInfo());

    KRATOS_EXPECT_TRUE(p_element->GetValue(PIPE_ACTIVE));
}

} // namespace Kratos::Testing