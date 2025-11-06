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

#include <unordered_map>

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
    void InitializeConditionMaps(ModelPart::ConditionsContainerType& rConditions);
    void FindConditionNeighboursBasedOnBoundaryType(
        std::function<PointerVector<Geometry<Node>>(const Geometry<Node>&)> GenerateBoundaries,
        ModelPart::ElementsContainerType&                                   rElements);

private:
    void AddNeighbouringElementsToConditionsBasedOnOverlappingBoundaryGeometries(
        Element& rElement, const Geometry<Node>::GeometriesArrayType& rBoundaryGeometries);
    void SetElementAsNeighbourOfAllConditionsWithIdenticalNodeIds(const std::vector<std::size_t>& rConditionNodeIds,
                                                                  Element* pElement);
    void        SetElementAsNeighbourIfRotatedNodeIdsAreEquivalent(Element& rElement,
                                                                   const std::vector<std::size_t>& element_boundary_node_ids,
                                                                   const GeometryData::KratosGeometryOrderType& r_order_type);
    static bool AreRotatedEquivalents(const std::vector<std::size_t>&              rFirst,
                                      const std::vector<std::size_t>&              rSecond,
                                      const GeometryData::KratosGeometryOrderType& rOrderType);
    static bool AreLinearRotatedEquivalents(std::vector<std::size_t>        First,
                                            const std::vector<std::size_t>& rSecond);
    static bool AreQuadraticRotatedEquivalents(std::vector<std::size_t>        First,
                                               const std::vector<std::size_t>& rSecond);
    static std::vector<std::size_t> GetNodeIdsFromGeometry(const Geometry<Node>& rGeometry);

    NodeIdsToConditionsHashMap     mConditionNodeIdsToConditions;
    SortedToUnsortedNodeIdsHashMap mSortedToUnsortedConditionNodeIds;
};
} // namespace Kratos