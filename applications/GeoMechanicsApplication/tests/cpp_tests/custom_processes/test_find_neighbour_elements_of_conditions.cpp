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
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities/element_setup_utilities.h"
#include "tests/cpp_tests/test_utilities/model_setup_utilities.h"

// namespace
// {

// using namespace Kratos;

// }

namespace Kratos::Testing
{

ModelPart& CreateModelPart3D4N(Model& rModel)
{
    return ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(rModel);
}

ModelPart& CreateModelPart3D8N(Model& rModel)
{
    return ModelSetupUtilities::CreateModelPartWithASingle3D8NElement(rModel);
}

ModelPart& CreateModelPart3D10N(Model& rModel)
{
    return ModelSetupUtilities::CreateModelPartWithASingle3D10NUPwDiffOrderElement(rModel);
}

ModelPart& CreateModelPart3D20N(Model& rModel)
{
    return ModelSetupUtilities::CreateModelPartWithASingle3D20NElement(rModel);
}

ModelPart& CreateModelPart3D6NInterface(Model& rModel)
{
    return ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(rModel);
}

class ParametrizedFindNeighbouring3D4NElementsFixture
    : public ::testing::TestWithParam<
          std::tuple<std::vector<std::size_t>, std::function<Condition::Pointer(const PointerVector<Node>&)>, std::function<ModelPart&(Model&)>>>
{
};

TEST_P(ParametrizedFindNeighbouring3D4NElementsFixture, NeighboursAreFoundForDifferentNodeOrderings)
{
    const auto& [order, condition_creator, model_part_creator] = GetParam();

    Model model;
    auto& r_model_part = model_part_creator(model);

    PointerVector<Node> nodes;
    for (const auto& r_node_id : order) {
        nodes.push_back(r_model_part.pGetNode(r_node_id));
    }

    auto p_condition = condition_creator(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

INSTANTIATE_TEST_SUITE_P(
    KratosGeoMechanicsFastSuiteWithoutKernel,
    ParametrizedFindNeighbouring3D4NElementsFixture,
    ::testing::Values(
        std::make_tuple(std::vector<std::size_t>{1, 3, 2}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D4N),
        std::make_tuple(std::vector<std::size_t>{2, 1, 3}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D4N),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D4N),
        std::make_tuple(std::vector<std::size_t>{4, 3, 2, 1}, ElementSetupUtilities::Create3D4NCondition, CreateModelPart3D8N),
        std::make_tuple(std::vector<std::size_t>{1, 4, 3, 2}, ElementSetupUtilities::Create3D4NCondition, CreateModelPart3D8N),
        std::make_tuple(std::vector<std::size_t>{2, 1, 4, 3}, ElementSetupUtilities::Create3D4NCondition, CreateModelPart3D8N),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 4}, ElementSetupUtilities::Create3D4NCondition, CreateModelPart3D8N),
        std::make_tuple(std::vector<std::size_t>{1, 3, 2, 7, 6, 5}, ElementSetupUtilities::Create3D6NCondition, CreateModelPart3D10N),
        std::make_tuple(std::vector<std::size_t>{2, 1, 3, 5, 7, 6}, ElementSetupUtilities::Create3D6NCondition, CreateModelPart3D10N),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 6, 5, 7}, ElementSetupUtilities::Create3D6NCondition, CreateModelPart3D10N),
        std::make_tuple(std::vector<std::size_t>{4, 3, 2, 1, 11, 10, 9, 12}, ElementSetupUtilities::Create3D8NCondition, CreateModelPart3D20N),
        std::make_tuple(std::vector<std::size_t>{1, 4, 3, 2, 12, 11, 10, 9}, ElementSetupUtilities::Create3D8NCondition, CreateModelPart3D20N),
        std::make_tuple(std::vector<std::size_t>{2, 1, 4, 3, 9, 12, 11, 10}, ElementSetupUtilities::Create3D8NCondition, CreateModelPart3D20N),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1, 4, 10, 9, 12, 11}, ElementSetupUtilities::Create3D8NCondition, CreateModelPart3D20N),
        std::make_tuple(std::vector<std::size_t>{1, 2, 3}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{3, 1, 2}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{2, 3, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface),
        std::make_tuple(std::vector<std::size_t>{3, 2, 1}, ElementSetupUtilities::Create3D3NCondition, CreateModelPart3D6NInterface)));

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, FindNeighbourElementsOfConditionsProcess_TestPointCondition)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(model);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(1));

    auto p_condition = ElementSetupUtilities::Create3D1NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, FindNeighbourElementsOfConditionsProcess_Test1DElement)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);

    PointerVector<Node> nodes;
    nodes.push_back(r_model_part.pGetNode(1));
    nodes.push_back(r_model_part.pGetNode(2));
    auto p_element = ElementSetupUtilities::Create2D2NElement(nodes, std::make_shared<Properties>(0));
    r_model_part.AddElement(p_element);

    auto p_condition = ElementSetupUtilities::Create2D2NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, TestFindNeighboursForMultipleConditionsOnTheSameNodes)
{
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
    process.Execute();
    EXPECT_EQ(p_condition1->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_condition2->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

class ParametrizedFindNeighbouring3D20NElementsOfLineConditionsFixture
    : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindNeighbouring3D20NElementsOfLineConditionsFixture, NeighboursAreFoundForDifferentNodeOrderings)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D20NElement(model);

    PointerVector<Node> nodes;
    const auto&         order = GetParam();
    for (const auto& r_node_id : order) {
        nodes.push_back(r_model_part.pGetNode(r_node_id));
    }

    auto p_condition = ElementSetupUtilities::Create3D3NLineCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
    process.Execute();
    EXPECT_EQ(p_condition->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindNeighbouring3D20NElementsOfLineConditionsFixture,
                         ::testing::Values(std::vector<std::size_t>{1, 2, 9},
                                           std::vector<std::size_t>{7, 8, 19}));

class ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture
    : public ::testing::TestWithParam<std::vector<std::size_t>>
{
};

TEST_P(ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture, ProcessThrowsWhenNoNeighboringElementsAreFound)
{
    Model model;
    auto& r_model_part = ModelSetupUtilities::CreateModelPartWithASingle3D4NElement(model);

    PointerVector<Node> nodes;
    const auto&         order = GetParam();
    for (const auto& r_node_id : order) {
        nodes.push_back(r_model_part.pGetNode(r_node_id));
    }

    auto p_condition = ElementSetupUtilities::Create3D3NCondition(nodes);
    r_model_part.AddCondition(p_condition);

    FindNeighbourElementsOfConditionsProcess process(r_model_part);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.Execute(),
                                      "Some conditions found without any corresponding element");
}

INSTANTIATE_TEST_SUITE_P(KratosGeoMechanicsFastSuiteWithoutKernel,
                         ParametrizedFindNeighbouring3D4NElementsThrowsWhenNotFoundFixture,
                         ::testing::Values(std::vector<std::size_t>{1, 2, 3},
                                           std::vector<std::size_t>{4, 2, 1},
                                           std::vector<std::size_t>{1, 3, 4}));

} // namespace Kratos::Testing
