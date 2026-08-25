// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_utilities/seepage_boundary_utilities.h"

#include <algorithm>

#include "custom_conditions/geo_seepage_condition.h"
#include "includes/variables.h"

namespace Kratos::Geo::SeepageBoundaryUtilities
{

std::vector<Node*> CollectSeepageNodes(ModelPart& rModelPart)
{
    auto result = std::vector<Node*>{};

    for (auto& r_condition : rModelPart.Conditions()) {
        if (dynamic_cast<const GeoSeepageCondition*>(&r_condition) == nullptr) continue;

        for (auto& r_node : r_condition.GetGeometry()) {
            result.push_back(&r_node);
        }
    }

    // Adjacent conditions share end nodes, so the same node can be collected more than once.
    std::sort(result.begin(), result.end(),
              [](const Node* pLeft, const Node* pRight) { return pLeft->Id() < pRight->Id(); });
    result.erase(std::unique(result.begin(), result.end()), result.end());

    return result;
}

} // namespace Kratos::Geo::SeepageBoundaryUtilities

