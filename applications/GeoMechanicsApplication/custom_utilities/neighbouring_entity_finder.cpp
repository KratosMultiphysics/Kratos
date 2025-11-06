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

#include "neighbouring_entity_finder.h"
#include <algorithm>

namespace Kratos
{

void NeighbouringEntityFinder::InitializeConditionMaps(ModelPart::ConditionsContainerType& rConditions)
{
    mConditionNodeIdsToConditions.clear();
    std::ranges::transform(
        rConditions, std::inserter(mConditionNodeIdsToConditions, mConditionNodeIdsToConditions.end()),
        [](auto& rCondition) {
        return NodeIdsToConditionsHashMap::value_type(
            GetNodeIdsFromGeometry(rCondition.GetGeometry()), {&rCondition});
    });

    mSortedToUnsortedConditionNodeIds.clear();
    std::ranges::transform(
        mConditionNodeIdsToConditions,
        std::inserter(mSortedToUnsortedConditionNodeIds, mSortedToUnsortedConditionNodeIds.end()),
        [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        std::ranges::sort(sorted_ids);
        return std::make_pair(sorted_ids, rPair.first);
    });
}

std::vector<std::size_t> NeighbouringEntityFinder::GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry)
{
    std::vector<std::size_t> result;
    result.reserve(rGeometry.size());
    std::ranges::transform(rGeometry, std::back_inserter(result),
                           [](const auto& rNode) { return rNode.Id(); });
    return result;
}

void NeighbouringEntityFinder::FindConditionNeighboursBasedOnBoundaryType(
    std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)> GenerateBoundaries,
    ModelPart::ElementsContainerType&                                   rElements)
{
    for (auto& r_element : rElements) {
        const auto& r_element_geometry  = r_element.GetGeometry();
        const auto  boundary_geometries = GenerateBoundaries(r_element_geometry);
        AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(r_element, boundary_geometries);
    }
}

void NeighbouringEntityFinder::AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
    Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries)
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        const auto element_boundary_node_ids = GetNodeIdsFromGeometry(r_boundary_geometry);

        if (mConditionNodeIdsToConditions.contains(element_boundary_node_ids)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
        } else if (r_boundary_geometry.LocalSpaceDimension() == 2) {
            // No condition is directly found for this boundary, but it might be a rotated equivalent
            SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
                rElement, element_boundary_node_ids, r_boundary_geometry.GetGeometryOrderType());
        }
    }
}

void NeighbouringEntityFinder::SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(
    const std::vector<std::size_t>& rConditionNodeIds, Element* pElement)
{
    const auto [start, end] = mConditionNodeIdsToConditions.equal_range(rConditionNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto& r_conditions  = it->second;
        auto vector_of_neighbours = GlobalPointersVector<Element>{Element::WeakPointer{pElement}};

        for (auto& p_condition : r_conditions) {
            p_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

void NeighbouringEntityFinder::SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
    Element&                                     rElement,
    const std::vector<std::size_t>&              element_boundary_node_ids,
    const GeometryData::KratosGeometryOrderType& r_order_type)
{
    auto sorted_boundary_node_ids = element_boundary_node_ids;
    std::ranges::sort(sorted_boundary_node_ids);

    if (mSortedToUnsortedConditionNodeIds.contains(sorted_boundary_node_ids)) {
        const auto unsorted_condition_node_ids =
            mSortedToUnsortedConditionNodeIds.find(sorted_boundary_node_ids)->second;
        if (AreRotatedEquivalents(element_boundary_node_ids, unsorted_condition_node_ids, r_order_type)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(unsorted_condition_node_ids, &rElement);
        }
    }
}

bool NeighbouringEntityFinder::AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                     const std::vector<std::size_t>& rSecond,
                                                     const GeometryData::KratosGeometryOrderType& rOrderType)
{
    switch (rOrderType) {
        using enum GeometryData::KratosGeometryOrderType;
    case Kratos_Linear_Order:
        return AreLinearRotatedEquivalents(rFirst, rSecond);
    case Kratos_Quadratic_Order:
        return AreQuadraticRotatedEquivalents(rFirst, rSecond);
    default:
        return false;
    }
}

bool NeighbouringEntityFinder::AreLinearRotatedEquivalents(std::vector<std::size_t>        First,
                                                           const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, First.end());
    return First == rSecond;
}

bool NeighbouringEntityFinder::AreQuadraticRotatedEquivalents(std::vector<std::size_t> First,
                                                              const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    auto       first_mid_side_node_id     = First.begin() + First.size() / 2;
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, first_mid_side_node_id);

    std::rotate(first_mid_side_node_id, first_mid_side_node_id + amount_of_needed_rotations, First.end());

    return First == rSecond;
}

} // namespace Kratos