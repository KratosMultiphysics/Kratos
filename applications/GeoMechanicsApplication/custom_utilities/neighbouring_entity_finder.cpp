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
#include "geometry_utilities.h"
#include <algorithm>

namespace Kratos
{

NeighbouringEntityFinder::NeighbouringEntityFinder(bool alsoSearchReverse)
    : mAlsoSearchReverse(alsoSearchReverse)
{
}

void NeighbouringEntityFinder::InitializeBoundaryMaps(NodeIdsToEntitiesHashMap GeometryNodeIdsToEntityMapping)
{
    mGeometryNodeIdsToEntities = std::move(GeometryNodeIdsToEntityMapping);

    mSortedToUnsortedEntityNodeIds.clear();
    std::ranges::transform(
        mGeometryNodeIdsToEntities,
        std::inserter(mSortedToUnsortedEntityNodeIds, mSortedToUnsortedEntityNodeIds.end()),
        [](const auto& rPair) {
        auto sorted_ids = rPair.first;
        std::ranges::sort(sorted_ids);
        return std::make_pair(sorted_ids, rPair.first);
    });
}

void NeighbouringEntityFinder::FindEntityNeighboursBasedOnBoundaryType(const BoundaryGenerator& rBoundaryGenerator,
                                                                       ModelPart::ElementsContainerType& rElements)
{
    for (auto& r_element : rElements) {
        const auto& r_element_geometry  = r_element.GetGeometry();
        const auto  boundary_geometries = rBoundaryGenerator(r_element_geometry);
        AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(r_element, boundary_geometries);
    }
}

void NeighbouringEntityFinder::AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
    Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries)
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        auto element_boundary_node_ids = GeometryUtilities::GetNodeIdsFromGeometry(r_boundary_geometry);
        if (mGeometryNodeIdsToEntities.contains(element_boundary_node_ids)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
        } else if (r_boundary_geometry.LocalSpaceDimension() == 2) {
            // No condition is directly found for this boundary, but it might be a rotated equivalent
            SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
                rElement, element_boundary_node_ids, r_boundary_geometry.GetGeometryOrderType());
        }

        if (mAlsoSearchReverse) {
            GeometryUtilities::ReverseNodes(element_boundary_node_ids,
                                           r_boundary_geometry.GetGeometryFamily(),
                                           r_boundary_geometry.GetGeometryOrderType());
            if (mGeometryNodeIdsToEntities.contains(element_boundary_node_ids)) {
                SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
            }
        }
    }
}

void NeighbouringEntityFinder::SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(
    const std::vector<std::size_t>& rConditionNodeIds, Element* pElement)
{
    const auto [start, end] = mGeometryNodeIdsToEntities.equal_range(rConditionNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto& r_conditions  = it->second;
        auto vector_of_neighbours = GlobalPointersVector<Element>{Element::WeakPointer{pElement}};

        for (auto& rp_condition : r_conditions) {
            if (rp_condition->GetGeometry().Id() == pElement->GetGeometry().Id()) continue;
            rp_condition->SetValue(NEIGHBOUR_ELEMENTS, vector_of_neighbours);
        }
    }
}

void NeighbouringEntityFinder::SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
    Element& rElement, const std::vector<std::size_t>& rElementBoundaryNodeIds, GeometryData::KratosGeometryOrderType OrderType)
{
    auto sorted_boundary_node_ids = rElementBoundaryNodeIds;
    std::ranges::sort(sorted_boundary_node_ids);

    if (mSortedToUnsortedEntityNodeIds.contains(sorted_boundary_node_ids)) {
        const auto unsorted_condition_node_ids =
            mSortedToUnsortedEntityNodeIds.find(sorted_boundary_node_ids)->second;
        if (AreRotatedEquivalents(rElementBoundaryNodeIds, unsorted_condition_node_ids, OrderType)) {
            SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(unsorted_condition_node_ids, &rElement);
        }
    }
}

bool NeighbouringEntityFinder::AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
                                                     const std::vector<std::size_t>& rSecond,
                                                     GeometryData::KratosGeometryOrderType OrderType)
{
    switch (OrderType) {
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
    const auto first_mid_side_node_id     = First.begin() + First.size() / 2;
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, first_mid_side_node_id);

    std::rotate(first_mid_side_node_id, first_mid_side_node_id + amount_of_needed_rotations, First.end());

    return First == rSecond;
}

} // namespace Kratos