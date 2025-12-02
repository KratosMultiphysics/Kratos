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
#include "includes/expect.h"
#include "includes/node.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"
#include "tests/cpp_tests/test_utilities.h"

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

PointerVectorSet<Node, IndexedObject> NodeContainerWithFixedAccelerationY(const intrusive_ptr<Node>& pNode)
{
    AddAcceleration(*pNode);
    pNode->Fix(ACCELERATION_Y);
    // Setting the buffer size only has an effect when the variables list has been set
    pNode->SetBufferSize(2);

    auto result = PointerVectorSet<Node, IndexedObject>{};
    result.SetMaxBufferSize(2);
    result.insert(pNode);
    return result;
}

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AssignUpdatedVectorVariableToNonFixedComponents_DoesNotUpdateFixedComponent)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       AssignUpdatedVectorVariableToNonFixedComponents_DoesNotUpdateAnythingWhenAllComponentsAreFixed)
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

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel,
       AssignUpdatedVectorVariableToNonFixedComponents_UpdatesEverythingWhenNoComponentIsFixed)
{
    Node node(1, 0.0, 0.0, 0.0);
    AddAcceleration(node);

    const array_1d<double, 3> new_values{1.0, 2.0, 3.0};
    NodeUtilities::AssignUpdatedVectorVariableToNonFixedComponents(node, ACCELERATION, new_values);

    const array_1d<double, 3> expected_vector{1.0, 2.0, 3.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_vector, node.FastGetSolutionStepValue(ACCELERATION, 0), 1e-12)
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AssignUpdatedVectorVariableToNodes_UpdatesEverythingRegardlessOfFixities)
{
    // Arrange
    auto p_node          = make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto nodes_container = NodeContainerWithFixedAccelerationY(p_node);

    // Act
    const array_1d<double, 3> new_acceleration_vector{1.0, 2.0, 4.0};
    NodeUtilities::AssignUpdatedVectorVariableToNodes(nodes_container, ACCELERATION, new_acceleration_vector, 0);
    NodeUtilities::AssignUpdatedVectorVariableToNodes(nodes_container, ACCELERATION, new_acceleration_vector, 1);

    // Assert
    const array_1d<double, 3> expected_acceleration_vector{1.0, 2.0, 4.0};
    KRATOS_EXPECT_VECTOR_NEAR(expected_acceleration_vector,
                              p_node->FastGetSolutionStepValue(ACCELERATION, 0), Defaults::absolute_tolerance)
    KRATOS_EXPECT_VECTOR_NEAR(expected_acceleration_vector,
                              p_node->FastGetSolutionStepValue(ACCELERATION, 1), Defaults::absolute_tolerance)
}

} // namespace Kratos::Testing
