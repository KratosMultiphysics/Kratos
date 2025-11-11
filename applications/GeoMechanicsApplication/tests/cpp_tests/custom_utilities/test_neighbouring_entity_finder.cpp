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

    NodeIdsToEntitiesHashMap map;
    std::ranges::transform(r_model_part_for_entities_for_finding.Elements(),
                           std::inserter(map, map.end()), [](auto& rCondition) {
        return NodeIdsToEntitiesHashMap::value_type(
            GeometryUtilities::GetNodeIdsFromGeometry(rCondition.GetGeometry()), {&rCondition});
    });
    finder.InitializeBoundaryMaps(map);

    auto generate_generic_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateBoundariesEntities();
    };

    finder.FindConditionNeighboursBasedOnBoundaryType(
        generate_generic_boundaries, r_model_part_for_neighbouring_elements.Elements());

    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(NeighbouringEntityFinder_ReturnsCorrectNeighbouringElementOfElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    Model model;

    NeighbouringEntityFinder finder;

    auto& r_model_part_for_entities_for_finding =
        ModelSetupUtilities::CreateModelPartWithASingle2D3NElement(model);
    NodeIdsToEntitiesHashMap map;
    for (auto& rElement : r_model_part_for_entities_for_finding.Elements()) {
        std::ranges::transform(rElement.GetGeometry().GenerateBoundariesEntities(),
                               std::inserter(map, map.end()), [&rElement](auto& BoundaryGeometry) {
            return NodeIdsToEntitiesHashMap::value_type(
                GeometryUtilities::GetNodeIdsFromGeometry(BoundaryGeometry), {&rElement});
        });
    }

    finder.InitializeBoundaryMaps(map);
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

    auto generate_generic_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateEdges();
    };

    finder.FindConditionNeighboursBasedOnBoundaryType(
        generate_generic_boundaries, r_model_part_for_neighbouring_elements.Elements());

    ASSERT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS).size(), 1);
    EXPECT_EQ(r_model_part_for_entities_for_finding.GetElement(1).GetValue(NEIGHBOUR_ELEMENTS)[0].GetId(), 42);
}

} // namespace Kratos::Testing