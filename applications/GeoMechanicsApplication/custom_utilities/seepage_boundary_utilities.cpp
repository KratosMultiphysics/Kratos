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

void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>& rDofs,
                                    const Vector&                    rRightHandSide,
                                    NodalFlowMap&                    rNodalFlows)
{
    KRATOS_ERROR_IF(rDofs.size() != rRightHandSide.size())
        << "Number of degrees of freedom (" << rDofs.size()
        << ") does not match the size of the right hand side (" << rRightHandSide.size() << ")"
        << std::endl;

    for (std::size_t i = 0; i < rDofs.size(); ++i) {
        if (rDofs[i]->GetVariable() != WATER_PRESSURE) continue;

        rNodalFlows[rDofs[i]->Id()] += rRightHandSide[i];
    }
}

NodalFlowMap CalculateNodalWaterFlows(ModelPart& rModelPart, const ProcessInfo& rProcessInfo)
{
    auto result = NodalFlowMap{};

    auto dofs            = std::vector<Dof<double>*>{};
    auto right_hand_side = Vector{};

    for (auto& r_element : rModelPart.Elements()) {
        r_element.GetDofList(dofs, rProcessInfo);
        r_element.CalculateRightHandSide(right_hand_side, rProcessInfo);

        AccumulateWaterPressureEntries(dofs, right_hand_side, result);
    }

    return result;
}

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


