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

#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

namespace Kratos::Testing
{
KRATOS_TEST_CASE_IN_SUITE(FindNeighbourElementsOfConditions, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(1));
    nodes.push_back(r_model_part.pGetNode(3));
    nodes.push_back(r_model_part.pGetNode(2));
    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);
    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighbourElementsOfConditionsSecondConfiguration, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(2));
    nodes.push_back(r_model_part.pGetNode(1));
    nodes.push_back(r_model_part.pGetNode(3));
    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);
    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighbourElementsOfConditionsThirdConfiguration, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);
    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(3));
    nodes.push_back(r_model_part.pGetNode(2));
    nodes.push_back(r_model_part.pGetNode(1));
    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);
    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}
} // namespace Kratos::Testing
