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

#include "custom_processes/find_neighbour_elements_of_conditions_process.h"
#include "test_setup_utilities/element_setup_utilities.hpp"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

ModelPart& CreateModelPartWithSingleElement(const std::string& rType, Model& rModel)
{
    if (rType == "3D4NElement")
        return ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(rModel);
    if (rType == "3D8NElement")
        return ModelSetupUtilities::CreateModelPartWithASingle3D8NElement(rModel);
    if (rType == "3D10NElement")
        return ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(rModel);
    if (rType == "3D20NElement")
        return ModelSetupUtilities::CreateModelPartWithASingle3D20NElement(rModel);
    if (rType == "3D6NInterfaceElement")
        return ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(rModel);
    if (rType == "2D2NElement")
        return ModelSetupUtilities::CreateModelPartWithASingle2D2NElement(rModel);

    KRATOS_ERROR << "Element type " << rType << " not recognized.";
}

PointerVector<Node> GetNodesFromIds(const std::vector<std::size_t>& rNodeIds, ModelPart& rModelPart)
{
    PointerVector<Node> result(rNodeIds.size());
    std::ranges::transform(rNodeIds, result.ptr_begin(),
                           [&rModelPart](auto Id) { return rModelPart.pGetNode(Id); });

    return result;
}

} // namespace

namespace Kratos::Testing
{

class ParametrizedFindNeighbouringElementsForDifferentTypesAndOrderingsFixture
    : public ::testing::TestWithParam<std::tuple<std::vector<std::size_t>, std::string, std::string>>
{
};

TEST_P(ParametrizedFindNeighbouringElementsForDifferentTypesAndOrderingsFixture, NeighboursAreFoundForDifferentNodeOrderings)
{
    // Arrange
    const auto& [order, condition_type, element_type] = GetParam();

    Model model;
    auto& r_model_part = CreateModelPartWithSingleElement(element_type, model);
    auto  nodes        = GetNodesFromIds(order, r_model_part);
    auto  p_condition  = ElementSetupUtilities::CreateCondition(condition_type, nodes);
    r_model_part.AddCondition(p_condition);
    FindNeighbourElementsOfConditionsProcess process(r_model_part);
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);

    // Act
    process.Execute();

    // Assert
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedFindNeighbouringElementsForDifferentTypesAndOrderingsFixture,
    ::testing::Values(
        std::make_tuple(std::vector<std::size_t>{1, 3, 2}, "3D3NCondition", "3D4NElement"),
        std::make_tuple(std::vector<std::size_t>{2, 1, 3}, "3D3NCondition", "3D4NElement"),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1}, "3D3NCondition", "3D4NElement"),
        std::make_tuple(std::vector<std::size_t>{4, 3, 2, 1}, "3D4NCondition", "3D8NElement"),
        std::make_tuple(std::vector<std::size_t>{1, 4, 3, 2}, "3D4NCondition", "3D8NElement"),
        std::make_tuple(std::vector<std::size_t>{2, 1, 4, 3}, "3D4NCondition", "3D8NElement"),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 4}, "3D4NCondition", "3D8NElement"),
        std::make_tuple(std::vector<std::size_t>{1, 3, 2, 7, 6, 5}, "3D6NCondition", "3D10NElement"),
        std::make_tuple(std::vector<std::size_t>{2, 1, 3, 5, 7, 6}, "3D6NCondition", "3D10NElement"),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 6, 5, 7}, "3D6NCondition", "3D10NElement"),
        std::make_tuple(std::vector<std::size_t>{4, 3, 2, 1, 11, 10, 9, 12}, "3D8NCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{1, 4, 3, 2, 12, 11, 10, 9}, "3D8NCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{2, 1, 4, 3, 9, 12, 11, 10}, "3D8NCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 4, 10, 9, 12, 11}, "3D8NCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{1, 2, 3}, "3D3NCondition", "3D6NInterfaceElement"),
        std::make_tuple(std::vector<std::size_t>{3, 1, 2}, "3D3NCondition", "3D6NInterfaceElement"),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, "3D3NCondition", "3D6NInterfaceElement"),
        std::make_tuple(std::vector<std::size_t>{1}, "3D1NCondition", "3D6NInterfaceElement"),
        std::make_tuple(std::vector<std::size_t>{1, 2, 9}, "3D3NLineCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{7, 8, 19}, "3D3NLineCondition", "3D20NElement"),
        std::make_tuple(std::vector<std::size_t>{1, 2}, "2D2NCondition", "2D2NElement")));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestFindNeighboursForMultipleConditionsOnTheSameNodes)
{
    // Arrange
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(model);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(1));

    auto p_condition1 = ElementSetupUtilities::Create3D1NCondition(nodes);
    auto p_condition2 = ElementSetupUtilities::Create3D1NCondition(nodes);
    r_model_part.AddCondition(p_condition1);
    p_condition2->SetId(2);
    r_model_part.AddCondition(p_condition2);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition1->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    EXPECT_EQ(p_condition2->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);

    // Act
    process.Execute();

    // Assert
    EXPECT_EQ(p_condition1->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_condition1->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);

    EXPECT_EQ(p_condition2->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_condition2->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

class ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture
    : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture, ProcessThrowsWhenNoNeighboringElementsAreFound)
{
    // Arrange
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);

    const auto& order = GetParam();
    auto        nodes = GetNodesFromIds(order, r_model_part);

    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Execute(),
                                      "The condition(s) with the following ID(s) is/are found "
                                      "without any corresponding element: [1]");
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture,
                         ::testing::Values(std::vector<std::size_t>{1, 2, 3},
                                           std::vector<std::size_t>{4, 2, 1},
                                           std::vector<std::size_t>{1, 3, 4}));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, FindNeighboursWithMultipleNeighbouringElements_AddsAllElementsAsNeighbour)
{
    // Arrange
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(1));
    nodes.push_back(r_model_part.pGetNode(2));
    nodes.push_back(r_model_part.pGetNode(3));

    r_model_part.AddElement(r_model_part.GetElement(1).Clone(2, nodes));

    PointerVector<Node> condition_nodes;
    condition_nodes.push_back(r_model_part.pGetNode(1));
    condition_nodes.push_back(r_model_part.pGetNode(2));
    auto p_condition = ElementSetupUtilities::Create2D2NCondition(condition_nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);

    // Act
    process.Execute();

    // Assert
    const auto& r_neighbours = p_condition->GetValue(NEIGHBOUR_ELEMENTS);
    EXPECT_EQ(r_neighbours.size(), 2);

    KRATOS_EXPECT_TRUE(std::any_of(r_neighbours.ptr_begin(), r_neighbours.ptr_end(),
                                   [](const auto& element) { return element->GetId() == 1; }))
    KRATOS_EXPECT_TRUE(std::any_of(r_neighbours.ptr_begin(), r_neighbours.ptr_end(),
                                   [](const auto& element) { return element->GetId() == 2; }))
}

KRATOS_TEST_CASE_IN_SUITE(CheckInfoFindNeighbourElementsOfConditionsProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    Model model;
    auto& r_empty_model_part = model.CreateModelPart("foo");
    const FindNeighbourElementsOfConditionsProcess process(r_empty_model_part);
    // Act & assert
    KRATOS_EXPECT_EQ(process.Info(), "FindNeighbourElementsOfConditionsProcess");
}

} // namespace Kratos::Testing
