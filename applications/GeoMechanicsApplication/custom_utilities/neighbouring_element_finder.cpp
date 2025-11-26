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

#include "neighbouring_element_finder.hpp"
#include "geometry_utilities.h"
#include <algorithm>

namespace Kratos
{

NeighbouringElementFinder::NeighbouringElementFinder(bool EnableReverseSearch)
    : mEnableReverseSearch(EnableReverseSearch)
{
}

void NeighbouringElementFinder::InitializeSortedToUnsortedEntityNodeIds()
{
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

void NeighbouringElementFinder::FindEntityNeighboursBasedOnBoundaryType(const BoundaryGenerator& rBoundaryGenerator,
                                                                        ModelPart::ElementsContainerType& rElements)
{
    for (auto& r_element : rElements) {
        const auto boundary_geometries = rBoundaryGenerator(r_element.GetGeometry());
        AddNeighbouringElementsToEntitiesBasedOnOverlappingBoundaryGeometries(r_element, boundary_geometries);
    }
}

void NeighbouringElementFinder::AddNeighbouringElementsToEntitiesBasedOnOverlappingBoundaryGeometries(
    Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries)
{
    for (const auto& r_boundary_geometry : rBoundaryGeometries) {
        AddNeighbouringElementsBasedOnBoundaryGeometry(rElement, r_boundary_geometry);
        if (mEnableReverseSearch) {
            constexpr auto reverse_search = true;
            AddNeighbouringElementsBasedOnBoundaryGeometry(rElement, r_boundary_geometry, reverse_search);
        }
    }
}

void NeighbouringElementFinder::AddNeighbouringElementsBasedOnBoundaryGeometry(Element& rElement,
                                                                               const Geometry<Node>& rBoundaryGeometry,
                                                                               bool ReverseSearch)
{
    auto element_boundary_node_ids = GenericUtilities::CollectIdsFromEntity(rBoundaryGeometry);
    if (ReverseSearch) {
        GeometryUtilities::ReverseNodes(element_boundary_node_ids, rBoundaryGeometry.GetGeometryFamily(),
                                        rBoundaryGeometry.GetGeometryOrderType());
    }

    if (mGeometryNodeIdsToEntities.contains(element_boundary_node_ids)) {
        SetElementAsNeighbourOfAllEntitiesWithIdenticalNodeIds(element_boundary_node_ids, &rElement);
    } else if (rBoundaryGeometry.LocalSpaceDimension() == 2) {
        // No entity is directly found for this boundary, but it might be a rotated equivalent
        SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
            rElement, element_boundary_node_ids, rBoundaryGeometry.GetGeometryOrderType());
    }
}

void NeighbouringElementFinder::SetElementAsNeighbourOfAllEntitiesWithIdenticalNodeIds(
    const std::vector<std::size_t>& rEntityNodeIds, Element* pElement)
{
    const auto [start, end] = mGeometryNodeIdsToEntities.equal_range(rEntityNodeIds);
    for (auto it = start; it != end; ++it) {
        const auto& r_entities = it->second;
        for (auto& rp_entity : r_entities) {
            if (rp_entity->GetGeometry().Id() == pElement->GetGeometry().Id()) continue;
            rp_entity->GetValue(NEIGHBOUR_ELEMENTS).push_back(Element::WeakPointer{pElement});
        }
    }
}

void NeighbouringElementFinder::SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(
    Element& rElement, const std::vector<std::size_t>& rElementBoundaryNodeIds, GeometryData::KratosGeometryOrderType OrderType)
{
    auto sorted_boundary_node_ids = rElementBoundaryNodeIds;
    std::ranges::sort(sorted_boundary_node_ids);

    if (mSortedToUnsortedEntityNodeIds.contains(sorted_boundary_node_ids)) {
        const auto unsorted_entity_node_ids =
            mSortedToUnsortedEntityNodeIds.find(sorted_boundary_node_ids)->second;
        if (AreRotatedEquivalents(rElementBoundaryNodeIds, unsorted_entity_node_ids, OrderType)) {
            SetElementAsNeighbourOfAllEntitiesWithIdenticalNodeIds(unsorted_entity_node_ids, &rElement);
        }
    }
}

bool NeighbouringElementFinder::AreRotatedEquivalents(const std::vector<std::size_t>& rFirst,
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
        // Rotations are not supported for orders higher than quadratic
        return false;
    }
}

bool NeighbouringElementFinder::AreLinearRotatedEquivalents(std::vector<std::size_t>        First,
                                                            const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, First.end());
    return First == rSecond;
}

bool NeighbouringElementFinder::AreQuadraticRotatedEquivalents(std::vector<std::size_t> First,
                                                               const std::vector<std::size_t>& rSecond)
{
    const auto amount_of_needed_rotations = std::ranges::find(First, rSecond[0]) - First.begin();
    const auto first_mid_side_node_id     = First.begin() + First.size() / 2;
    std::rotate(First.begin(), First.begin() + amount_of_needed_rotations, first_mid_side_node_id);

    std::rotate(first_mid_side_node_id, first_mid_side_node_id + amount_of_needed_rotations, First.end());

    return First == rSecond;
}

} // namespace Kratos