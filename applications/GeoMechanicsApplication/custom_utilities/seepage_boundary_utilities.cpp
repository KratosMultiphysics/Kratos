// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Wijtze Pieter Kikstra

#include <algorithm>
#include <functional>
#include <set>

#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/seepage_boundary_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/variables.h"

namespace Kratos::Geo
{

namespace
{

void AccumulateWaterPressureEntries(const std::vector<Dof<double>*>&        rElementDofs,
                                    const Vector&                           rElementRightHandSide,
                                    SeepageBoundaryUtilities::NodalFlowMap& rNodalFlows)
{
    for (auto i = std::size_t{0}; i < rElementDofs.size(); ++i) {
        if (rElementDofs[i]->GetVariable() != WATER_PRESSURE) continue;

        rNodalFlows[rElementDofs[i]->Id()] += rElementRightHandSide[i];
    }
}

} // namespace

SeepageBoundaryUtilities::NodalFlowMap SeepageBoundaryUtilities::CalculateNodalWaterFlows(
    ModelPart::ElementsContainerType& rElements, const ProcessInfo& rProcessInfo)
{
    auto result = NodalFlowMap{};

    for (auto& r_element : rElements) {
        if (!r_element.IsActive()) continue;

        auto dofs = std::vector<Dof<double>*>{};
        r_element.GetDofList(dofs, rProcessInfo);
        auto right_hand_side = Vector{};
        r_element.CalculateRightHandSide(right_hand_side, rProcessInfo);

        AccumulateWaterPressureEntries(dofs, right_hand_side, result);
    }

    return result;
}

void SeepageBoundaryUtilities::AssignNodalWaterFlows(ModelPart& rModelPart, const NodalFlowMap& rNodalFlows)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(NODAL_WATER_FLOW) = 0.0;
    }

    for (const auto& [node_id, flow] : rNodalFlows) {
        rModelPart.GetNode(node_id).FastGetSolutionStepValue(NODAL_WATER_FLOW) = flow;
    }
}

namespace
{

struct NodeComparator {
    bool operator()(const Node* pLeft, const Node* pRight) const
    {
        return pLeft->Id() < pRight->Id();
    }
};

} // namespace

std::vector<Node*> SeepageBoundaryUtilities::CollectSeepageNodes(ModelPart& rModelPart)
{
    auto result = std::set<Node*, NodeComparator>{};

    for (auto& r_condition : rModelPart.Conditions()) {
        if (!dynamic_cast<const GeoSeepageCondition*>(&r_condition)) continue;

        std::ranges::transform(r_condition.GetGeometry(), std::inserter(result, result.end()),
                               [](Node& r_node) { return &r_node; });
    }

    return {result.begin(), result.end()};
}

namespace
{

using CandidatePredicateType = std::function<bool(const Node*)>;
using ScoreCalculatorType    = std::function<double(const Node*)>;

// Returns the node maximising the given score among the candidates, or nullptr when there are none.
// Candidates are visited in ascending node id order, and a strict comparison keeps the first of any
// tie, which makes the choice reproducible.
Node* SelectBestCandidate(const std::vector<Node*>&     rNodes,
                          const CandidatePredicateType& rIsCandidate,
                          const ScoreCalculatorType&    rScoreCalculator)
{
    // For filtering the candidates, we'd prefer to use std::views::filter, but unfortunately not
    // all compilers on GitHub support it yet. Therefore, we copy the candidates into a separate
    // vector.
    auto candidates = std::vector<Node*>{};
    std::ranges::copy_if(rNodes, std::back_inserter(candidates), rIsCandidate);
    auto first_score_less_than_second_score = [&rScoreCalculator](auto* pFirstNode, auto* pSecondNode) {
        return rScoreCalculator(pFirstNode) < rScoreCalculator(pSecondNode);
    };
    auto iter = std::ranges::max_element(candidates, first_score_less_than_second_score);
    return iter != candidates.end() ? *iter : nullptr;
}

} // namespace

bool SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(const std::vector<Node*>& rSeepageNodes,
                                                            const NodalFlowMap&       rNodalFlows,
                                                            int                       EchoLevel)
{
    if (EchoLevel > 1) {
        for (auto* p_node : rSeepageNodes) {
            KRATOS_INFO("Node") << p_node->Id()
                                << " pressure = " << p_node->FastGetSolutionStepValue(WATER_PRESSURE)
                                << "\n";
            KRATOS_INFO("Node") << p_node->Id() << " fixed = " << p_node->IsFixed(WATER_PRESSURE) << "\n";
        }
    }

    // A free node with negative water pressure is a candidate for fixing.
    auto is_candidate     = CandidatePredicateType{[](const auto* pNode) {
        return pNode && !pNode->IsFixed(WATER_PRESSURE) &&
               pNode->FastGetSolutionStepValue(WATER_PRESSURE) < 0.0;
    }};
    auto score_calculator = ScoreCalculatorType{[](const auto* pNode) {
        return pNode ? -1.0 * pNode->FastGetSolutionStepValue(WATER_PRESSURE) : 0.0;
    }};
    if (auto* p_node = SelectBestCandidate(rSeepageNodes, is_candidate, score_calculator)) {
        KRATOS_INFO_IF("Switch", EchoLevel > 1)
            << "Node " << p_node->Id() << " switched to Dirichlet, because pressure was "
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
    is_candidate     = CandidatePredicateType{[&flow_of](const auto* pNode) {
        return pNode && pNode->IsFixed(WATER_PRESSURE) && flow_of(*pNode) < -1e-11;
    }};
    score_calculator = ScoreCalculatorType{
        [&flow_of](const auto* pNode) { return pNode ? -1.0 * flow_of(*pNode) : 0.0; }};
    if (auto* p_node = SelectBestCandidate(rSeepageNodes, is_candidate, score_calculator)) {
        KRATOS_INFO_IF("Switch", EchoLevel > 1) << "Node " << p_node->Id() << " switched to Neumann, because flow was "
                                                << flow_of(*p_node) << "\n";
        p_node->Free(WATER_PRESSURE);
        return true;
    }

    return false;
}

} // namespace Kratos::Geo
