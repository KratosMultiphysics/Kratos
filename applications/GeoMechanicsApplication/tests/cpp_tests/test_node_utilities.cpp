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
#include "custom_utilities/node_utilities.h"
#include "geo_mechanics_fast_suite.h"
#include "includes/node.h"
#include "includes/variables.h"

using namespace Kratos;

namespace
{

void AddAcceleration(Node& rNode)
{
    auto variables_list = make_intrusive<VariablesList>();
    variables_list->Add(ACCELERATION);
    rNode.SetSolutionStepVariablesList(variables_list);
    rNode.AddDof(ACCELERATION_X);
    rNode.AddDof(ACCELERATION_Y);
    rNode.AddDof(ACCELERATION_Z);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(NodeUtilities_DoesNotUpdateFixedComponent, KratosGeoMechanicsFastSuite)
{
    Node node(1, 0.0, 0.0, 0.0);
    AddAcceleration(node);

    node.Fix(ACCELERATION_X);
    node.Fix(ACCELERATION_Y);

    const array_1d<double, 3> new_values{1.0, 2.0, 3.0};
    NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(node, ACCELERATION, new_values);

    const array_1d<double, 3> expected_vector{0.0, 0.0, 3.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_vector, node.FastGetSolutionStepValue(ACCELERATION, 0), 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(NodeUtilities_DoesNotUpdateAnythingWhenAllComponentsAreFixed, KratosGeoMechanicsFastSuite)
{
    Node node(1, 0.0, 0.0, 0.0);
    AddAcceleration(node);

    node.Fix(ACCELERATION_X);
    node.Fix(ACCELERATION_Y);
    node.Fix(ACCELERATION_Z);

    const array_1d<double, 3> new_values{1.0, 2.0, 3.0};
    NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(node, ACCELERATION, new_values);

    const array_1d<double, 3> expected_vector{0.0, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_vector, node.FastGetSolutionStepValue(ACCELERATION, 0), 1e-12)
}

KRATOS_TEST_CASE_IN_SUITE(NodeUtilities_UpdatesEverythingWhenNoComponentIsFixed, KratosGeoMechanicsFastSuite)
{
    Node node(1, 0.0, 0.0, 0.0);
    AddAcceleration(node);

    const array_1d<double, 3> new_values{1.0, 2.0, 3.0};
    NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(node, ACCELERATION, new_values);

    const array_1d<double, 3> expected_vector{1.0, 2.0, 3.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_vector, node.FastGetSolutionStepValue(ACCELERATION, 0), 1e-12)
}

} // namespace Kratos::Testing
