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
    auto p_geometry = Kratos::make_shared<Line2D2<Node>>(
        rModelPart.pGetNode(static_cast<int>(FirstNodeId)), rModelPart.pGetNode(static_cast<int>(SecondNodeId)));
    auto p_properties = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0)
                                                    : rModelPart.CreateNewProperties(0);
    rModelPart.AddCondition(
        Kratos::make_intrusive<GeoSeepageCondition>(static_cast<int>(Id), p_geometry, p_properties));
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CollectSeepageNodesReturnsNothingWhenThereAreNoSeepageConditions, KratosGeoMechanicsFastSuiteWithoutKernel)
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

} // namespace Kratos::Testing


