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

#include "custom_utilities/neighbouring_entity_finder.h"
#include "test_setup_utilities/element_setup_utilities.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_ReturnsEmptyListWhenNoNeighbouringElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("empty");

    NeighbouringEntityFinder finder;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_ReturnsCorrectNeighbouringElementOfElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    NeighbouringEntityFinder finder;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("empty");

    std::vector         neighbour_node_ids = {1, 2};
    PointerVector<Node> neighbour_nodes(neighbour_node_ids.size());
    std::ranges::transform(neighbour_node_ids, neighbour_nodes.ptr_begin(),
                           [&r_model_part_for_entities_for_finding](auto Id) {
        return r_model_part_for_entities_for_finding.pGetNode(Id);
    });
    auto p_element = ElementSetupUtilities::Create2D2NElement(neighbour_nodes, {});
    p_element->SetId(42);
    r_model_part_for_neighbouring_elements.AddElement(p_element);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 42);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_NeverRefersToItselfAsNeighbour, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    NeighbouringEntityFinder finder;
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_entities_for_finding.Elements(), boundary_generators);

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_FindsNeighboursBetweenTwoContinua, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    r_model_part.CreateNewNode(3, 1.0, 1.0, 0.0);
    r_model_part.CreateNewNode(4, 0.0, 1.0, 0.0);

    std::vector<std::size_t> node_ids_element_1 = {1, 2, 3};
    PointerVector<Node>      nodes_element_1(node_ids_element_1.size());
    std::ranges::transform(node_ids_element_1, nodes_element_1.ptr_begin(),
                           [&r_model_part](auto Id) { return r_model_part.pGetNode(Id); });
    auto p_element_1 = ElementSetupUtilities::Create2D3NElement(nodes_element_1, {});
    r_model_part.AddElement(p_element_1);

    std::vector<std::size_t> node_ids_element_2 = {1, 3, 4};
    PointerVector<Node>      nodes_element_2(node_ids_element_2.size());
    std::ranges::transform(node_ids_element_2, nodes_element_2.ptr_begin(),
                           [&r_model_part](auto Id) { return r_model_part.pGetNode(Id); });
    auto p_element_2 = ElementSetupUtilities::Create2D3NElement(nodes_element_2, {});
    p_element_2->SetId(2);
    r_model_part.AddElement(p_element_2);

    constexpr auto           alsoSearchReverse = true;
    NeighbouringEntityFinder finder(alsoSearchReverse);

    auto generate_generic_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateBoundariesEntities();
    };
    finder.FindEntityNeighboursBasedOnBoundaryType(generate_generic_boundaries, r_model_part.Elements());
    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_model_part.Elements(), r_model_part.Elements(), boundary_generators);

    ASSERT_EQ(p_element_1->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_element_1->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 2);

    ASSERT_EQ(p_element_2->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_element_2->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

} // namespace Kratos::Testing