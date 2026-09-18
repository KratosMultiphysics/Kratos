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
//                   Wijtze Pieter Kikstra

#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/seepage_boundary_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "includes/variables.h"
#include "test_setup_utilities/element_setup_utilities.hpp"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;
using namespace std::string_literals;

namespace
{

// Builds a model part with the given number of nodes, each owning a WATER_PRESSURE DoF.
ModelPart& CreateModelPartWithNodes(Model& rModel, std::size_t NumberOfNodes)
{
    auto& r_model_part = rModel.CreateModelPart("Main"s);
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    Testing::ModelSetupUtilities::CreateNumberOfNewNodes(r_model_part, NumberOfNodes);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return r_model_part;
}

// Adds a two-noded seepage condition spanning the two given node ids. All conditions share one
// Properties object, since creating the same property id twice is an error.
void AddSeepageCondition(ModelPart& rModelPart, std::size_t Id, std::size_t FirstNodeId, std::size_t SecondNodeId)
{
    auto p_geometry =
        std::make_shared<Line2D2<Node>>(rModelPart.pGetNode(static_cast<int>(FirstNodeId)),
                                        rModelPart.pGetNode(static_cast<int>(SecondNodeId)));
    auto p_properties = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0)
                                                    : rModelPart.CreateNewProperties(0);
    rModelPart.AddCondition(make_intrusive<GeoSeepageCondition>(static_cast<int>(Id), p_geometry, p_properties));
}

// Test-only shortcut: treats every node of the model part as a seepage node, so the decision logic
// can be exercised without constructing conditions. Production code uses CollectSeepageNodes.
std::vector<Node*> AllNodesOf(ModelPart& rModelPart)
{
    auto result = std::vector<Node*>{};
    std::ranges::transform(rModelPart.Nodes(), std::back_inserter(result),
                           [](auto& rNode) { return &rNode; });
    return result;
}

class MockPwElementForSeepageTests : public Element
{
public:
    explicit MockPwElementForSeepageTests(const std::vector<intrusive_ptr<Node>>& rNodes);

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo&) override;

private:
    std::vector<intrusive_ptr<Node>> mNodes;
};

MockPwElementForSeepageTests::MockPwElementForSeepageTests(const std::vector<intrusive_ptr<Node>>& rNodes)
    : mNodes{rNodes}
{
}

void MockPwElementForSeepageTests::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList.clear();
    for (const auto& rp_node : mNodes) {
        std::ranges::transform(rp_node->GetDofs(), std::back_inserter(rElementalDofList),
                               [](const auto& rpDof) { return rpDof.get(); });
    }
}

void MockPwElementForSeepageTests::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo&)
{
    auto dofs = DofsVectorType{};
    this->GetDofList(dofs, ProcessInfo{});

    rRightHandSideVector.resize(dofs.size());
    // For each water pressure DoF, use the node ID as the right hand side value. For every other
    // DoF, use twice the node ID.
    auto calculate_dof_value = [](const auto* pDof) {
        return pDof->GetVariable() == WATER_PRESSURE ? static_cast<double>(pDof->Id())
                                                     : 2.0 * static_cast<double>(pDof->Id());
    };
    std::ranges::transform(dofs, rRightHandSideVector.begin(), calculate_dof_value);
}

// Creates two triangular mock elements sharing an edge, with WATER_PRESSURE DoFs on each node.

//   3      4
//   *------*
//   | \    |
//   |   \  |
//   |     \|
//   *------*
//   1      2
auto CreateTwoConnectedMockElements(const Geo::ConstVariableDataRefs& rSolutionStepVariables,
                                    const Geo::ConstVariableRefs&     rDegreesOfFreedom)
{
    const auto nodal_positions = std::vector{Point{0.0, 0.0, 0.0}, Point{1.0, 0.0, 0.0},
                                             Point{0.0, 1.0, 0.0}, Point{1.0, 1.0, 0.0}};
    auto       nodes           = Testing::ElementSetupUtilities::GenerateNodes(nodal_positions);

    Testing::ElementSetupUtilities::AddVariablesToNodes(nodes, rSolutionStepVariables, rDegreesOfFreedom);

    const auto nodes_of_element_1 =
        std::vector{nodes.GetContainer()[0], nodes.GetContainer()[1], nodes.GetContainer()[2]};
    const auto nodes_of_element_2 =
        std::vector{nodes.GetContainer()[1], nodes.GetContainer()[3], nodes.GetContainer()[2]};

    auto result = ModelPart::ElementsContainerType{};
    result.push_back(make_intrusive<MockPwElementForSeepageTests>(nodes_of_element_1));
    result.push_back(make_intrusive<MockPwElementForSeepageTests>(nodes_of_element_2));
    return result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CollectSeepageNodesReturnsNothingWhenThereAreNoSeepageConditions,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::CollectSeepageNodes(r_model_part).empty())
}

KRATOS_TEST_CASE_IN_SUITE(CollectSeepageNodesReturnsSharedNodesOnlyOnce, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    // Two adjacent conditions sharing node 2.
    AddSeepageCondition(r_model_part, 1, 1, 2);
    AddSeepageCondition(r_model_part, 2, 2, 3);

    const auto nodes = Geo::SeepageBoundaryUtilities::CollectSeepageNodes(r_model_part);

    ASSERT_EQ(nodes.size(), 3);
    KRATOS_EXPECT_EQ(nodes[0]->Id(), 1);
    KRATOS_EXPECT_EQ(nodes[1]->Id(), 2);
    KRATOS_EXPECT_EQ(nodes[2]->Id(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(CalculateNodalWaterFlowsSumsRHSContributionsFromSeveralElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto elements = CreateTwoConnectedMockElements({std::cref(WATER_PRESSURE)}, {std::cref(WATER_PRESSURE)});

    const auto nodal_flow_map =
        Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(elements, ProcessInfo{});

    ASSERT_EQ(nodal_flow_map.size(), 4);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(1), 1 * 1.0); // only element 1 contributes
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(2), 2 * 2.0); // both elements contribute
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(3), 2 * 3.0); // both elements contribute
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(4), 1 * 4.0); // only element 2 contributes
}

KRATOS_TEST_CASE_IN_SUITE(CalculateNodalWaterFlowsSkipsInactiveElements, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto elements = CreateTwoConnectedMockElements({std::cref(WATER_PRESSURE)}, {std::cref(WATER_PRESSURE)});
    elements.back().Set(ACTIVE, false); // deactivate element 2

    const auto nodal_flow_map =
        Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(elements, ProcessInfo{});

    ASSERT_EQ(nodal_flow_map.size(), 3);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(1), 1 * 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(2), 1 * 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(3), 1 * 3.0);
    // Node 4 is not part of any active element, so it should not appear in the map
}

KRATOS_TEST_CASE_IN_SUITE(CalculateNodalWaterFlowsConsidersPwDoFsOnly, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto elements = CreateTwoConnectedMockElements(
        {std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT)},
        {std::cref(WATER_PRESSURE), std::cref(DISPLACEMENT_X), std::cref(DISPLACEMENT_Y)});

    const auto nodal_flow_map =
        Geo::SeepageBoundaryUtilities::CalculateNodalWaterFlows(elements, ProcessInfo{});

    ASSERT_EQ(nodal_flow_map.size(), 4);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(1), 1 * 1.0); // only element 1 contributes
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(2), 2 * 2.0); // both elements contribute
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(3), 2 * 3.0); // both elements contribute
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flow_map.at(4), 1 * 4.0); // only element 2 contributes
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeDoesNothingWhenNoNodeViolatesItsCondition,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 fixed with no inflow, node 2 free and under suction: both are consistent.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = 5.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 1.0}, {2, 0.0}};

    KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(nodes, nodal_flows))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeFixesTheHighestPressureFreeNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Free(WATER_PRESSURE);
    }
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(WATER_PRESSURE) = -2.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = -9.0; // highest
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;

    const auto nodes = AllNodesOf(r_model_part);

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(
        nodes, Geo::SeepageBoundaryUtilities::NodalFlowMap{}))

    // Only node 2 switches, and it is prescribed at zero pressure.
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeReleasesTheLargestInflowFixedNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    const auto nodes = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, -4.0}, {2, -11.0}, {3, 2.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(nodes, nodal_flows))

    // Only node 2, which has the largest outflow, is released.
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodePrefersFixingOverReleasing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 is a fixed node with a large inflow, node 2 is a free node under negative pressure.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = -1.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, -100.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(nodes, nodal_flows))

    // The Neumann to Dirichlet switch wins, and node 1 is left alone this iteration.
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeBreaksTiesByLowestNodeId, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, -5.0}, {2, -5.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNodeIfNeeded(nodes, nodal_flows))

    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main"s);
    r_model_part.AddNodalSolutionStepVariable(NODAL_WATER_FLOW);
    for (auto i = std::size_t{1}; i <= 3; ++i) {
        r_model_part.CreateNewNode(static_cast<int>(i), static_cast<double>(i), 0.0, 0.0);
    }
    // Pre-seed a stale value on node 3 to prove it gets overwritten.
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW) = 99.0;

    // Node 2 is deliberately absent from the map and must end up at 0.0.
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 4.0}, {3, -2.0}};

    Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(r_model_part, nodal_flows);

    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW), -2.0);
}

} // namespace Kratos::Testing
