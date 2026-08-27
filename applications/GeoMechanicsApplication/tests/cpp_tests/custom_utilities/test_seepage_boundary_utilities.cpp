// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "containers/model.h"
#include "custom_conditions/geo_seepage_condition.h"
#include "custom_utilities/seepage_boundary_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/line_2d_2.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
using namespace Kratos;

namespace
{

// Builds a model part with a chain of nodes along the x axis, each owning a WATER_PRESSURE dof.
Kratos::ModelPart& CreateModelPartWithNodes(Kratos::Model& rModel, std::size_t NumberOfNodes)
{
    auto& r_model_part = rModel.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    for (std::size_t i = 1; i <= NumberOfNodes; ++i) {
        r_model_part.CreateNewNode(static_cast<int>(i), static_cast<double>(i), 0.0, 0.0);
    }
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    return r_model_part;
}

// Adds a two-noded seepage condition spanning the two given node ids. All conditions share one
// Properties object, since creating the same property id twice is an error.
void AddSeepageCondition(Kratos::ModelPart& rModelPart, std::size_t Id, std::size_t FirstNodeId, std::size_t SecondNodeId)
{
    auto p_geometry =
        Kratos::make_shared<Line2D2<Node>>(rModelPart.pGetNode(static_cast<int>(FirstNodeId)),
                                           rModelPart.pGetNode(static_cast<int>(SecondNodeId)));
    auto p_properties = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0)
                                                    : rModelPart.CreateNewProperties(0);
    rModelPart.AddCondition(
        Kratos::make_intrusive<GeoSeepageCondition>(static_cast<int>(Id), p_geometry, p_properties));
}

// Test-only shortcut: treats every node of the model part as a seepage node, so the decision logic
// can be exercised without constructing conditions. Production code uses CollectSeepageNodes.
std::vector<Kratos::Node*> AllNodesOf(Kratos::ModelPart& rModelPart)
{
    auto result = std::vector<Kratos::Node*>{};
    for (auto& r_node : rModelPart.Nodes()) {
        result.push_back(&r_node);
    }
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

    KRATOS_EXPECT_EQ(nodes.size(), 3);
    KRATOS_EXPECT_EQ(nodes[0]->Id(), 1);
    KRATOS_EXPECT_EQ(nodes[1]->Id(), 2);
    KRATOS_EXPECT_EQ(nodes[2]->Id(), 3);
}

KRATOS_TEST_CASE_IN_SUITE(AccumulateWaterPressureEntriesPicksOnlyWaterPressureDofs, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X);
        r_node.AddDof(DISPLACEMENT_Y);
        r_node.AddDof(WATER_PRESSURE);
    }

    // Mimic a U-Pw ordering: ux, uy, p for node 1, then ux, uy, p for node 2.
    auto dofs = std::vector<Dof<double>*>{r_model_part.pGetNode(1)->pGetDof(DISPLACEMENT_X),
                                          r_model_part.pGetNode(1)->pGetDof(DISPLACEMENT_Y),
                                          r_model_part.pGetNode(1)->pGetDof(WATER_PRESSURE),
                                          r_model_part.pGetNode(2)->pGetDof(DISPLACEMENT_X),
                                          r_model_part.pGetNode(2)->pGetDof(DISPLACEMENT_Y),
                                          r_model_part.pGetNode(2)->pGetDof(WATER_PRESSURE)};

    auto right_hand_side = Vector{6};
    right_hand_side[0]   = 10.0; // ux node 1, must be ignored
    right_hand_side[1]   = 20.0; // uy node 1, must be ignored
    right_hand_side[2]   = 3.0;  // p  node 1
    right_hand_side[3]   = 30.0; // ux node 2, must be ignored
    right_hand_side[4]   = 40.0; // uy node 2, must be ignored
    right_hand_side[5]   = 7.0;  // p  node 2

    auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{};
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, right_hand_side, nodal_flows);

    KRATOS_EXPECT_EQ(nodal_flows.size(), 2);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(1), 3.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(2), 7.0);
}

KRATOS_TEST_CASE_IN_SUITE(AccumulateWaterPressureEntriesSumsContributionsFromSeveralElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);

    auto dofs = std::vector<Dof<double>*>{r_model_part.pGetNode(1)->pGetDof(WATER_PRESSURE),
                                          r_model_part.pGetNode(2)->pGetDof(WATER_PRESSURE)};

    auto first_contribution = Vector{2};
    first_contribution[0]   = 1.5;
    first_contribution[1]   = 2.5;

    auto second_contribution = Vector{2};
    second_contribution[0]   = 0.5;
    second_contribution[1]   = -2.5;

    auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{};
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, first_contribution, nodal_flows);
    Geo::SeepageBoundaryUtilities::AccumulateWaterPressureEntries(dofs, second_contribution, nodal_flows);

    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(nodal_flows.at(2), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeDoesNothingWhenNoNodeViolatesItsCondition,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 fixed with no outflow, node 2 free and under suction: both are consistent.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = -5.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, -1.0}, {2, 0.0}};

    KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))
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
    r_model_part.pGetNode(1)->FastGetSolutionStepValue(WATER_PRESSURE) = 2.0;
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = 9.0; // highest
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(WATER_PRESSURE) = -1.0;

    const auto nodes = AllNodesOf(r_model_part);

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(
        nodes, Geo::SeepageBoundaryUtilities::NodalFlowMap{}))

    // Only node 2 switches, and it is prescribed at zero pressure.
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE), 0.0);
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodeReleasesTheLargestOutflowFixedNode, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 3);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.Fix(WATER_PRESSURE);
    }

    const auto nodes = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 4.0}, {2, 11.0}, {3, -2.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

    // Only node 2, which has the largest outflow, is released.
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(3)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(SwitchOneSeepageNodePrefersFixingOverReleasing, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = CreateModelPartWithNodes(model, 2);
    // Node 1 is a fixed node with a large outflow, node 2 is a free node under positive pressure.
    r_model_part.pGetNode(1)->Fix(WATER_PRESSURE);
    r_model_part.pGetNode(2)->Free(WATER_PRESSURE);
    r_model_part.pGetNode(2)->FastGetSolutionStepValue(WATER_PRESSURE) = 1.0;

    const auto nodes       = AllNodesOf(r_model_part);
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 100.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

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
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 5.0}, {2, 5.0}};

    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::SwitchOneSeepageNode(nodes, nodal_flows))

    KRATOS_EXPECT_FALSE(r_model_part.pGetNode(1)->IsFixed(WATER_PRESSURE))
    KRATOS_EXPECT_TRUE(r_model_part.pGetNode(2)->IsFixed(WATER_PRESSURE))
}

KRATOS_TEST_CASE_IN_SUITE(ShouldReleaseToNeumannIsTrueOnlyForPositiveFlow, KratosGeoMechanicsFastSuiteWithoutKernel){
    KRATOS_EXPECT_TRUE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(1.0e-9))
        KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(0.0))
            KRATOS_EXPECT_FALSE(Geo::SeepageBoundaryUtilities::ShouldReleaseToNeumann(-1.0))}

KRATOS_TEST_CASE_IN_SUITE(AssignNodalWaterFlowsWritesMappedValuesAndZeroesTheRest, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(NODAL_WATER_FLOW);
    for (std::size_t i = 1; i <= 3; ++i) {
        r_model_part.CreateNewNode(static_cast<int>(i), static_cast<double>(i), 0.0, 0.0);
    }
    // Pre-seed a stale value on node 3 to prove it gets zeroed.
    r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW) = 99.0;

    // Node 2 is deliberately absent from the map and must end up at 0.0.
    const auto nodal_flows = Geo::SeepageBoundaryUtilities::NodalFlowMap{{1, 4.0}, {3, -2.0}};

    Geo::SeepageBoundaryUtilities::AssignNodalWaterFlows(r_model_part, nodal_flows);

    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(1)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 4.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(2)->FastGetSolutionStepValue(NODAL_WATER_FLOW), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_model_part.pGetNode(3)->FastGetSolutionStepValue(NODAL_WATER_FLOW), -2.0);
}

} // namespace Kratos::Testing
