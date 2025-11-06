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


void NeighbouringEntityFinder::InitializeConditionMaps(
    ModelPart::ConditionsContainerType& rConditions)
{
    mConditionNodeIdsToConditions.clear();
    std::ranges::transform(
        rConditions,
        std::inserter(mConditionNodeIdsToConditions, mConditionNodeIdsToConditions.end()), [this](auto& rCondition) {
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
}