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

class ParametrizedFindNeighbourElementsOfConditions : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindNeighbourElementsOfConditions, NeighboringTetraElementsAreFoundForDifferentNodeOrderings)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);

    PointerVector<Node> nodes;
    const auto          order = GetParam();
    nodes.push_back(r_model_part.pGetNode(order[0]));
    nodes.push_back(r_model_part.pGetNode(order[1]));
    nodes.push_back(r_model_part.pGetNode(order[2]));

    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindNeighbourElementsOfConditions,
                         ::testing::Values(std::vector<std::size_t>{1, 3, 2},
                                           std::vector<std::size_t>{2, 1, 3},
                                           std::vector<std::size_t>{3, 2, 1}));

class ParametrizedFindHexaNeighbourElementsOfConditions
    : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindHexaNeighbourElementsOfConditions, NeighboringHexaElementsAreFoundForDifferentNodeOrderings)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D8NElement(model);

    PointerVector<Node> nodes;
    const auto          order = GetParam();
    nodes.push_back(r_model_part.pGetNode(order[0]));
    nodes.push_back(r_model_part.pGetNode(order[1]));
    nodes.push_back(r_model_part.pGetNode(order[2]));
    nodes.push_back(r_model_part.pGetNode(order[3]));

    auto p_condition = ElementSetupUtilities::Create3D4NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindHexaNeighbourElementsOfConditions,
                         ::testing::Values(std::vector<std::size_t>{4, 3, 2, 1},
                                           std::vector<std::size_t>{1, 4, 3, 2},
                                           std::vector<std::size_t>{2, 1, 4, 3},
                                           std::vector<std::size_t>{3, 2, 1, 4}));

class ParametrizedFindQuadraticTetraNeighbourElementsOfConditions
    : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindQuadraticTetraNeighbourElementsOfConditions, NeighboringHexaElementsAreFoundForDifferentNodeOrderings)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(model);

    PointerVector<Node> nodes;
    const auto          order = GetParam();
    nodes.push_back(r_model_part.pGetNode(order[0]));
    nodes.push_back(r_model_part.pGetNode(order[1]));
    nodes.push_back(r_model_part.pGetNode(order[2]));
    nodes.push_back(r_model_part.pGetNode(order[3]));
    nodes.push_back(r_model_part.pGetNode(order[4]));
    nodes.push_back(r_model_part.pGetNode(order[5]));

    auto p_condition = ElementSetupUtilities::Create3D6NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindQuadraticTetraNeighbourElementsOfConditions,
                         ::testing::Values(std::vector<std::size_t>{1, 3, 2, 7, 6, 5},
                                           std::vector<std::size_t>{2, 1, 3, 5, 7, 6},
                                           std::vector<std::size_t>{3, 2, 1, 6, 5, 7}));
} // namespace Kratos::Testing
