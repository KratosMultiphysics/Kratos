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
#include <limits>

#include "custom_conditions/geo_seepage_condition.h"
#include "geo_mechanics_application_variables.h"
#include "includes/variables.h"

namespace Kratos::Geo::SeepageBoundaryUtilities
{

void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>& rDofs,
                                    const Vector&                    rRightHandSide,
                                    NodalFlowMap&                    rNodalFlows)
{
    KRATOS_ERROR_IF(rDofs.size() != rRightHandSide.size())
        << "Number of degrees of freedom (" << rDofs.size() << ") does not match the size of the right hand side ("
        << rRightHandSide.size() << ")" << std::endl;

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

void AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_WATER_FLOW) = 0.0;
    }

    for (const auto& [node_id, flow] : rNodalFlows) {
        rModelPart.GetNode(node_id).FastGetSolutionStepValue(NODAL_WATER_FLOW) = flow;
    }
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

namespace
{

// Returns the node maximising the given score among the candidates, or nullptr when there are none.
// Candidates are visited in ascending node id order, and a strict comparison keeps the first of any
// tie, which makes the choice reproducible.
template <typename PredicateType, typename ScoreType>
Node* SelectBestCandidate(const std::vector<Node*>& rNodes, PredicateType IsCandidate, ScoreType Score)
{
    Node* p_best     = nullptr;
    auto  best_score = 0.0;

    for (auto* p_node : rNodes) {
        if (!IsCandidate(*p_node)) continue;

        const auto score = Score(*p_node);
        if (p_best == nullptr || score > best_score) {
            p_best     = p_node;
            best_score = score;
        }
    }

    return p_best;
}

} // namespace

bool SwitchOneSeepageNode(const std::vector<Node*>& rSeepageNodes, const NodalFlowMap& rNodalFlows)
{
    KRATOS_TRY

    for (auto* p_node : rSeepageNodes) {
        KRATOS_INFO("Node") << p_node->Id()
                            << " pressure = " << p_node->FastGetSolutionStepValue(WATER_PRESSURE) << "\n";
        KRATOS_INFO("Node") << p_node->Id() << " fixed = " << p_node->IsFixed(WATER_PRESSURE) << "\n";
    }
    // A free node under positive pressure is unsaturated, so it cannot be a draining face. Fixing
    // takes precedence over releasing.
    if (auto* p_node = SelectBestCandidate(rSeepageNodes, [](const Node& rNode) {
        return !rNode.IsFixed(WATER_PRESSURE) && rNode.FastGetSolutionStepValue(WATER_PRESSURE) < 0.0;
    }, [](const Node& rNode) { return -1.0 * rNode.FastGetSolutionStepValue(WATER_PRESSURE); })) {
        KRATOS_INFO("Switch") << "Node " << p_node->Id() << " switched to Dirichlet, because pressure was "
                              << p_node->FastGetSolutionStepValue(WATER_PRESSURE) << "\n";
        p_node->FastGetSolutionStepValue(WATER_PRESSURE) = 0.0;
        p_node->Fix(WATER_PRESSURE);
        return true;
    }

    // Otherwise release the prescribed node carrying the largest inflow.
    const auto flow_of = [&rNodalFlows](const Node& rNode) {
        const auto it = rNodalFlows.find(rNode.Id());
        return it == rNodalFlows.end() ? 0.0 : it->second;
    };

    if (auto* p_node = SelectBestCandidate(rSeepageNodes, [&flow_of](const Node& rNode) {
        return rNode.IsFixed(WATER_PRESSURE) && flow_of(rNode) < -1e-11;
    }, [&flow_of](const Node& rNode) { return -1.0 * flow_of(rNode); })) {
        KRATOS_INFO("Switch") << "Node " << p_node->Id() << " switched to Neumann, because flow was "
                              << flow_of(*p_node) << "\n";
        p_node->Free(WATER_PRESSURE);
        return true;
    }

    return false;

    KRATOS_CATCH("")
}

} // namespace Kratos::Geo::SeepageBoundaryUtilities
