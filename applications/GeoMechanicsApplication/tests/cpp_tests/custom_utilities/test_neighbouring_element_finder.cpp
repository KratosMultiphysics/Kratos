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

#include "custom_utilities/neighbouring_element_finder.hpp"
#include "test_setup_utilities/element_setup_utilities.h"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{
using namespace Kratos;

void CreateNumberOfNewNodes(ModelPart& rModelPart, std::size_t NumberOfNodes)
{
    for (std::size_t i = 0; i < NumberOfNodes; ++i) {
        rModelPart.CreateNewNode(i + 1, 0.0, 0.0, 0.0);
    }
}

PointerVector<Node> GetNodesFromIds(ModelPart& rModelPart, const std::vector<std::size_t>& rNodeIds)
{
    PointerVector<Node> result(rNodeIds.size());
    std::ranges::transform(rNodeIds, result.ptr_begin(),
                           [&rModelPart](auto Id) { return rModelPart.pGetNode(Id); });
    return result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_ReturnsEmptyListWhenNoNeighbouringElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("empty");

    NeighbouringElementFinder finder;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_ReturnsCorrectNeighbouringElementOfElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("Main");

    std::vector<std::size_t> neighbour_node_ids = {1, 2};
    PointerVector<Node>      neighbour_nodes =
        GetNodesFromIds(r_model_part_for_entities_for_finding, neighbour_node_ids);
    auto p_element = ElementSetupUtilities::Create2D2NElement(neighbour_nodes, {});
    p_element->SetId(42);
    r_model_part_for_neighbouring_elements.AddElement(p_element);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    NeighbouringElementFinder finder;
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 42);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_ReturnsCorrectNeighbouringElementOfElement_WhenRunningSearchTwice,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    auto& r_model_part_for_neighbouring_elements = model.CreateModelPart("Main");

    std::vector<std::size_t> neighbour_node_ids = {1, 2};
    PointerVector<Node>      neighbour_nodes =
        GetNodesFromIds(r_model_part_for_entities_for_finding, neighbour_node_ids);
    auto p_element = ElementSetupUtilities::Create2D2NElement(neighbour_nodes, {});
    p_element->SetId(42);
    r_model_part_for_neighbouring_elements.AddElement(p_element);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    NeighbouringElementFinder finder;
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    // The second run should clear the previous results and find the same neighbour again.
    // This would happen when running a multi-stage analysis (where the search is called
    // at the start of each stage).
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_neighbouring_elements.Elements(), boundary_generators);

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 42);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_NeverRefersToItselfAsNeighbour, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    NeighbouringElementFinder finder;
    finder.FindEntityNeighbours(r_model_part_for_entities_for_finding.Elements(),
                                r_model_part_for_entities_for_finding.Elements(), boundary_generators);

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_FindsNeighboursBetweenTwoContinua, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_model_part, 4);

    std::vector<std::size_t> node_ids_element_1 = {1, 2, 3};
    PointerVector<Node>      nodes_element_1    = GetNodesFromIds(r_model_part, node_ids_element_1);
    auto p_element_1 = ElementSetupUtilities::Create2D3NElement(nodes_element_1, {});
    r_model_part.AddElement(p_element_1);

    std::vector<std::size_t> node_ids_element_2 = {1, 3, 4};
    PointerVector<Node>      nodes_element_2    = GetNodesFromIds(r_model_part, node_ids_element_2);
    auto p_element_2 = ElementSetupUtilities::Create2D3NElement(nodes_element_2, {});
    p_element_2->SetId(2);
    r_model_part.AddElement(p_element_2);

    constexpr auto            also_search_reverse = true;
    NeighbouringElementFinder finder(also_search_reverse);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_model_part.Elements(), r_model_part.Elements(), boundary_generators);

    ASSERT_EQ(p_element_1->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_element_1->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 2);

    ASSERT_EQ(p_element_2->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_element_2->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_FindsNeighboursBetweenInterfaceAndContinuum,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_model_part, 5);

    std::vector<std::size_t> node_ids_continuum_element = {3, 4, 5};
    PointerVector<Node> nodes_continuum_element = GetNodesFromIds(r_model_part, node_ids_continuum_element);
    auto p_continuum_element = ElementSetupUtilities::Create2D3NElement(nodes_continuum_element, {});
    r_model_part.AddElement(p_continuum_element);

    std::vector<std::size_t> node_ids_element_2 = {1, 2, 3, 4};
    PointerVector<Node>      nodes_element_2    = GetNodesFromIds(r_model_part, node_ids_element_2);
    auto p_interface_element = ElementSetupUtilities::Create2D4NInterfaceElement(nodes_element_2, {});
    p_interface_element->SetId(2);
    r_model_part.AddElement(p_interface_element);

    auto& r_interface_model_part = model.CreateModelPart("Interfaces");
    r_interface_model_part.AddElement(p_interface_element);

    constexpr auto            also_search_reverse = true;
    NeighbouringElementFinder finder(also_search_reverse);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{1}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_interface_model_part.Elements(), r_model_part.Elements(), boundary_generators);

    ASSERT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_FindsNeighboursBetweenQuadraticInterfaceAndContinuum,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_model_part, 9);

    std::vector<std::size_t> node_ids_continuum_element = {4, 5, 8, 6, 9, 7};
    PointerVector<Node> nodes_continuum_element = GetNodesFromIds(r_model_part, node_ids_continuum_element);
    auto p_continuum_element = ElementSetupUtilities::Create2D6NElement(nodes_continuum_element, {});
    r_model_part.AddElement(p_continuum_element);

    std::vector<std::size_t> node_ids_element_2 = {1, 2, 3, 4, 5, 6};
    PointerVector<Node>      nodes_element_2    = GetNodesFromIds(r_model_part, node_ids_element_2);
    auto p_interface_element = ElementSetupUtilities::Create2D6NInterfaceElement(nodes_element_2, {});
    p_interface_element->SetId(2);
    r_model_part.AddElement(p_interface_element);

    auto& r_interface_model_part = model.CreateModelPart("Interfaces");
    r_interface_model_part.AddElement(p_interface_element);

    constexpr auto            also_search_reverse = true;
    NeighbouringElementFinder finder(also_search_reverse);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{1}] = std::make_unique<EdgesGenerator>();
    finder.FindEntityNeighbours(r_interface_model_part.Elements(), r_model_part.Elements(), boundary_generators);

    ASSERT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringElementFinder_FindsNeighboursBetweenQuadraticSurfaceInterfaceAnd3DContinuum,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_model_part, 16);

    std::vector<std::size_t> node_ids_continuum_element = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    PointerVector<Node> nodes_continuum_element = GetNodesFromIds(r_model_part, node_ids_continuum_element);
    auto p_continuum_element = ElementSetupUtilities::Create3D10NElement(nodes_continuum_element, {});
    r_model_part.AddElement(p_continuum_element);
    p_continuum_element->SetId(1);

    // The node ordering of the interface element is chosen such that it is both reversed and
    // permutated, to test the robustness of the neighbour finding.
    std::vector<std::size_t> node_ids_element_2 = {3, 1, 2, 7, 5, 6, 11, 12, 13, 14, 15, 16};
    PointerVector<Node>      nodes_element_2    = GetNodesFromIds(r_model_part, node_ids_element_2);
    auto p_interface_element = ElementSetupUtilities::Create3D12NInterfaceElement(nodes_element_2, {});

    p_interface_element->SetId(2);
    r_model_part.AddElement(p_interface_element);

    auto& r_interface_model_part = model.CreateModelPart("Interfaces");
    r_interface_model_part.AddElement(p_interface_element);

    constexpr auto            also_search_reverse = true;
    NeighbouringElementFinder finder(also_search_reverse);

    std::map<std::size_t, std::unique_ptr<BoundaryGenerator>> boundary_generators;
    boundary_generators[std::size_t{2}] = std::make_unique<FacesGenerator>();
    finder.FindEntityNeighbours(r_interface_model_part.Elements(), r_model_part.Elements(), boundary_generators);

    ASSERT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

} // namespace Kratos::Testing