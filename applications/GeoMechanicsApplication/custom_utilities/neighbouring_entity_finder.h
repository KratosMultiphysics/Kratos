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

#pragma once

#include "includes/condition.h"
#include "includes/key_hash.h"
#include "includes/kratos_export_api.h"
#include "includes/model_part.h"

namespace Kratos
{
using NodeIdsToConditionsHashMap     = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<Condition::Pointer>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;
using SortedToUnsortedNodeIdsHashMap = std::unordered_multimap<std::vector<std::size_t>,
                                                               std::vector<std::size_t>,
                                                               KeyHasherRange<std::vector<std::size_t>>,
                                                               KeyComparorRange<std::vector<std::size_t>>>;

class KRATOS_API(GEO_MECHANICS_APPLICATION) NeighbouringEntityFinder
{
public:
    NodeIdsToConditionsHashMap     mConditionNodeIdsToConditions;
    SortedToUnsortedNodeIdsHashMap mSortedToUnsortedConditionNodeIds;
    void InitializeConditionMaps(ModelPart::ConditionsContainerType& rConditions);
    std::vector<std::size_t> GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry);
    void                     FindConditionNeighboursBasedOnBoundaryType(
                            std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)> GenerateBoundaries,
                            ModelPart::ElementsContainerType&                                   rElements);

    void AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries);
    void SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(const std::vector<std::size_t>& rConditionNodeIds,
                                                                  Element* pElement);
    void SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(Element& rElement,
                                                            const std::vector<std::size_t>& element_boundary_node_ids,
                                                            const GeometryData::KratosGeometryOrderType& r_order_type);
    bool AreRotatedEquivalents(const std::vector<std::size_t>&              rFirst,
                               const std::vector<std::size_t>&              rSecond,
                               const GeometryData::KratosGeometryOrderType& rOrderType);
    bool AreLinearRotatedEquivalents(std::vector<std::size_t> First, const std::vector<std::size_t>& rSecond);
    bool AreQuadraticRotatedEquivalents(std::vector<std::size_t> First, const std::vector<std::size_t>& rSecond);
};

} // namespace Kratos