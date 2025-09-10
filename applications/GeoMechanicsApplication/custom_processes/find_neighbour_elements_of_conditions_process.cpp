// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#include "custom_processes/find_neighbour_elements_of_conditions_process.hpp"

#include "geometries/geometry.h"
#include "includes/kratos_flags.h"

namespace Kratos
{

void FindNeighbourElementsOfConditionsProcess::Execute()
{
    if (mrModelPart.Conditions().empty()) return;

    NodeIdToConditionsHashMap condition_node_ids_to_condition;
    for (auto& r_condition : mrModelPart.Conditions()) {
        r_condition.Set(VISITED, false);
        auto& r_geometry = r_condition.GetGeometry();

        std::vector<IndexType> Ids(r_geometry.size());
        std::ranges::transform(r_geometry, Ids.begin(), [](const auto& rNode) { return rNode.Id(); });

        condition_node_ids_to_condition.insert(NodeIdToConditionsHashMap::value_type(Ids, {&r_condition}));
    }

    SortedToUnsortedNodeIdsHashMap sorted_to_unsorted_condition_node_ids;
    std::ranges::transform(
        condition_node_ids_to_condition,
        std::inserter(sorted_to_unsorted_condition_node_ids, sorted_to_unsorted_condition_node_ids.end()),
        [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        std::ranges::sort(sorted_ids);
        return std::make_pair(sorted_ids, rPair.first);
    });

    auto generate_boundaries = [](const auto& rGeometry) {
        return rGeometry.GenerateBoundariesEntities();
    };
    auto generate_points   = [](const auto& rGeometry) { return rGeometry.GeneratePoints(); };
    auto generate_edges_3d = [](const auto& rGeometry) {
        return rGeometry.LocalSpaceDimension() == 3 ? rGeometry.GenerateEdges()
                                                    : PointerVector<Geometry<Node>>();
    };
    auto generate_edges_1d = [](const auto& rGeometry) {
        return rGeometry.LocalSpaceDimension() == 1 ? rGeometry.GenerateEdges()
                                                    : PointerVector<Geometry<Node>>();
    };

    // Note the order in the generators: the 1D elements are only added
    // as neighbours when the condition is not neighbouring 2D or 3D elements
    const std::vector<std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)>> boundary_generators = {
        generate_boundaries, generate_points, generate_edges_3d, generate_edges_1d};

    for (const auto& r_generator : boundary_generators) {
        CheckBoundaryTypeForAllElements(r_generator, condition_node_ids_to_condition,
                                        sorted_to_unsorted_condition_node_ids);
        if (AllConditionsAreVisited()) return;
    }

    ReportConditionsWithoutNeighbours();
    KRATOS_ERROR << "Some conditions found without any corresponding element" << std::endl;
}

void FindNeighbourElementsOfConditionsProcess::CheckBoundaryTypeForAllElements(
    auto                            generate_boundaries,
    NodeIdToConditionsHashMap&      condition_node_ids_to_condition,
    SortedToUnsortedNodeIdsHashMap& sorted_condition_node_ids_to_condition)
{
    for (auto& r_element : mrModelPart.Elements()) {
        const auto& rGeometryElement    = r_element.GetGeometry();
        const auto  rBoundaryGeometries = generate_boundaries(rGeometryElement);
        AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
            condition_node_ids_to_condition, sorted_condition_node_ids_to_condition, r_element, rBoundaryGeometries);
    }
}

void FindNeighbourElementsOfConditionsProcess::AddNeighboringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
    NodeIdToConditionsHashMap&                 FacesMap,
    const SortedToUnsortedNodeIdsHashMap&      FacesMapSorted,
    Element&                                   rElement,
    const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries) const
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        std::vector<IndexType> element_boundary_node_ids(r_boundary_geometry.size());
        std::ranges::transform(r_boundary_geometry, element_boundary_node_ids.begin(),
                               [](const Node& rNode) { return rNode.Id(); });
        std::vector<std::size_t> adjacent_condition_node_ids;
        auto                     itFace = FacesMap.find(element_boundary_node_ids);
        if (itFace != FacesMap.end()) {
            adjacent_condition_node_ids = itFace->first;
        }

        if (itFace == FacesMap.end() && r_boundary_geometry.LocalSpaceDimension() == 2) {
            // condition is not found but might be a problem of ordering in 2D boundary geometries!
            std::vector<std::size_t> face_ids_sorted = element_boundary_node_ids;
            std::ranges::sort(face_ids_sorted);

            auto it = FacesMapSorted.find(face_ids_sorted);
            if (it != FacesMapSorted.end()) {
                bool permutations_found = false;
                switch (r_boundary_geometry.GetGeometryOrderType()) {
                    using enum GeometryData::KratosGeometryOrderType;
                case Kratos_Linear_Order:
                    permutations_found = FindPermutations(element_boundary_node_ids, it->second);
                    break;
                case Kratos_Quadratic_Order:
                    permutations_found = FindPermutationsQuadratic(element_boundary_node_ids, it->second);
                    break;
                default:
                    break;
                }
                if (permutations_found) {
                    adjacent_condition_node_ids = it->second;
                }
            }
        }

        if (!adjacent_condition_node_ids.empty()) {
            // condition is found!
            // but check if there are more than one condition on the element
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(FacesMap, adjacent_condition_node_ids, &rElement);
        }
    }
}

bool FindNeighbourElementsOfConditionsProcess::AllConditionsAreVisited() const
{
    return std::ranges::all_of(mrModelPart.Conditions(),
                               [](const auto& rCondition) { return rCondition.Is(VISITED); });
}

void FindNeighbourElementsOfConditionsProcess::SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(
    NodeIdToConditionsHashMap& rFacesMap, const std::vector<std::size_t>& rConditionNodeIds, Element* pElement)
{
    const auto [start, end] = rFacesMap.equal_range(rConditionNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto&                   r_conditions = it->second;
        GlobalPointersVector<Element> vector_of_neighbours;
        vector_of_neighbours.resize(1);
        vector_of_neighbours(0) = Element::WeakPointer(pElement);

        for (auto& p_condition : r_conditions) {
            p_condition->Set(VISITED, true);
            p_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

bool FindNeighbourElementsOfConditionsProcess::FindPermutations(std::vector<std::size_t> elements_boundary_node_ids,
                                                                const std::vector<std::size_t>& condition_node_ids) const
{
    const auto amount_of_needed_rotations =
        std::ranges::find(elements_boundary_node_ids, condition_node_ids[0]) -
        elements_boundary_node_ids.begin();
    std::ranges::rotate(elements_boundary_node_ids, elements_boundary_node_ids.begin() + amount_of_needed_rotations);
    return elements_boundary_node_ids == condition_node_ids;
}

bool FindNeighbourElementsOfConditionsProcess::FindPermutationsQuadratic(
    std::vector<std::size_t> elements_boundary_node_ids, const std::vector<std::size_t>& condition_node_ids) const
{
    const auto position_of_first_condition_node_in_element_boundary =
        std::ranges::find(elements_boundary_node_ids, condition_node_ids[0]);
    const auto amount_of_needed_rotations =
        position_of_first_condition_node_in_element_boundary - elements_boundary_node_ids.begin();

    std::rotate(elements_boundary_node_ids.begin(), elements_boundary_node_ids.begin() + amount_of_needed_rotations,
                elements_boundary_node_ids.begin() + elements_boundary_node_ids.size() / 2);

    std::rotate(elements_boundary_node_ids.begin() + elements_boundary_node_ids.size() / 2,
                elements_boundary_node_ids.begin() + elements_boundary_node_ids.size() / 2 + amount_of_needed_rotations,
                elements_boundary_node_ids.end());
    return elements_boundary_node_ids == condition_node_ids;
}

void FindNeighbourElementsOfConditionsProcess::ReportConditionsWithoutNeighbours() const
{
    std::vector<Condition> unvisited_conditions;
    std::ranges::copy_if(mrModelPart.Conditions(), std::back_inserter(unvisited_conditions),
                         [](const auto& rCondition) { return !rCondition.Is(VISITED); });
    for (const auto& r_condition : unvisited_conditions) {
        KRATOS_INFO("Condition without any corresponding element, ID ") << r_condition.Id() << std::endl;
    }
}

} // namespace Kratos