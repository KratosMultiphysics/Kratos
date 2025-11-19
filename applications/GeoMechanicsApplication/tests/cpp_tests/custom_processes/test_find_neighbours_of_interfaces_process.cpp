// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "containers/model.h"
#include "custom_processes/find_neighbours_of_interfaces_process.h"
#include "includes/kratos_parameters.h"
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

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_IsAProcess, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto model = Model{};
    model.CreateModelPart("EmptyModelPart");
    auto test_settings = Parameters{};
    test_settings.AddString("model_part_name", "EmptyModelPart");
    test_settings.AddString("model_part_name_for_neighbouring_elements", "EmptyModelPart");
    const auto process = FindNeighboursOfInterfacesProcess{model, test_settings};

    // Act & assert
    EXPECT_NE(dynamic_cast<const Process*>(&process), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_FindsNoNeighboursWhenAdjacentModelPartIsEmpty,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto  model = Model{};
    auto& r_model_part_with_interface_element =
        ModelSetupUtilities::CreateModelPartWithASingle3D6NInterfaceElement(model);
    model.CreateModelPart("EmptyModelPart");
    auto test_settings = Parameters{};
    test_settings.AddString("model_part_name", "Main");
    test_settings.AddString("model_part_name_for_neighbouring_elements", "EmptyModelPart");
    auto process = FindNeighboursOfInterfacesProcess{model, test_settings};

    // Act
    process.ExecuteInitialize();

    // Assert
    EXPECT_EQ(r_model_part_with_interface_element.Elements().front().GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_FindsContinuumNeighbourOfInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_computational_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_computational_model_part, 9);

    const auto node_ids_continuum_element = std::vector<std::size_t>{4, 5, 8, 6, 9, 7};
    const auto nodes_continuum_element = GetNodesFromIds(r_computational_model_part, node_ids_continuum_element);
    auto p_continuum_element = ElementSetupUtilities::Create2D6NElement(nodes_continuum_element, {});
    r_computational_model_part.AddElement(p_continuum_element);

    const auto node_ids_element_2 = std::vector<std::size_t>{1, 2, 3, 4, 5, 6};
    const auto nodes_element_2    = GetNodesFromIds(r_computational_model_part, node_ids_element_2);
    auto p_interface_element = ElementSetupUtilities::Create2D6NInterfaceElement(nodes_element_2, {});
    p_interface_element->SetId(2);
    r_computational_model_part.AddElement(p_interface_element);

    auto& r_interface_model_part = model.CreateModelPart("Interfaces");
    r_interface_model_part.AddElement(p_interface_element);

    FindNeighboursOfInterfacesProcess process(model, Parameters(R"({"model_part_name": "Interfaces",
"model_part_name_for_neighbouring_elements": "Main"})"));

    process.ExecuteInitialize();

    ASSERT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(FindNeighboursOfInterfacesProcess_OnlyFindsNeighbourWhenLocalDimensionOfNeighbourIsHigher,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;
    auto& r_computational_model_part = model.CreateModelPart("Main");
    CreateNumberOfNewNodes(r_computational_model_part, 6);

    const auto node_ids_line_element = std::vector<std::size_t>{1, 2, 3};
    const auto nodes_line_element = GetNodesFromIds(r_computational_model_part, node_ids_line_element);
    auto p_line_element = ElementSetupUtilities::Create2D3NLineElement(nodes_line_element, {});
    r_computational_model_part.AddElement(p_line_element);

    const auto node_ids_element_2 = std::vector<std::size_t>{1, 2, 3, 4, 5, 6};
    const auto nodes_element_2    = GetNodesFromIds(r_computational_model_part, node_ids_element_2);
    auto p_interface_element = ElementSetupUtilities::Create2D6NInterfaceElement(nodes_element_2, {});
    p_interface_element->SetId(2);
    r_computational_model_part.AddElement(p_interface_element);

    auto& r_interface_model_part = model.CreateModelPart("Interfaces");
    r_interface_model_part.AddElement(p_interface_element);

    FindNeighboursOfInterfacesProcess process(model, Parameters(R"({"model_part_name": "Interfaces",
"model_part_name_for_neighbouring_elements": "Main"})"));

    process.ExecuteInitialize();

    EXPECT_EQ(p_interface_element->GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

} // namespace Kratos::Testing
