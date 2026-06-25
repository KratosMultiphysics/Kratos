//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "processes/structured_mesh_generator_process.h"
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "testing/testing.h"

namespace Kratos::Testing
{

namespace
{

void SetTestDofArrayModelPart(ModelPart& rModelPart)
{
    // Set the background mesh model part
    rModelPart.GetProcessInfo().SetValue(DOMAIN_SIZE, 2);
    rModelPart.AddNodalSolutionStepVariable(DISTANCE);
    const Variable<double>* p_distance_var = &DISTANCE;
    rModelPart.GetNodalSolutionStepVariablesList().AddDof(p_distance_var);

    // Generate the background mesh (done with the StructuredMeshGeneratorProcess)
    auto p_point_1 = Kratos::make_intrusive<Node>(1, -0.5, -0.5,  0.0);
    auto p_point_2 = Kratos::make_intrusive<Node>(2, -0.5,  0.5,  0.0);
    auto p_point_3 = Kratos::make_intrusive<Node>(3,  0.5,  0.5,  0.0);
    auto p_point_4 = Kratos::make_intrusive<Node>(4,  0.5, -0.5,  0.0);
    Quadrilateral2D4<Node> geometry(p_point_1, p_point_2, p_point_3, p_point_4);
    Parameters mesher_parameters(R"({
        "number_of_divisions": 3,
        "element_name": "DistanceCalculationElementSimplex2D3N",
        "condition_name": "LineCondition"
    })");
    StructuredMeshGeneratorProcess(geometry, rModelPart, mesher_parameters).Execute();

    // Add DISTANCE DOFs to mesh nodes
    for (auto& rNode : rModelPart.Nodes()) {
        rNode.AddDof(*p_distance_var);
    }
};

}

KRATOS_TEST_CASE_IN_SUITE(SetUpDofArray, KratosCoreFastSuite)
{
    // Set a mesh with DISTANCE DOFs elements
    Model model;
    auto& r_model_part = model.CreateModelPart("TestModelPart");
    SetTestDofArrayModelPart(r_model_part);

    // Set up the DOFs array
    DofArrayUtilities::DofsArrayType dofs_array;
    DofArrayUtilities::SetUpDofArray(r_model_part, dofs_array);

    // Check results
    KRATOS_CHECK_EQUAL(dofs_array.size(), 16);
    KRATOS_CHECK_EQUAL(dofs_array.begin()->GetVariable(), DISTANCE);
}

KRATOS_TEST_CASE_IN_SUITE(SetUpEffectiveDofArray, KratosCoreFastSuite)
{
    // Set a mesh with DISTANCE DOFs elements
    Model model;
    auto& r_model_part = model.CreateModelPart("TestModelPart");
    SetTestDofArrayModelPart(r_model_part);

    // Create some constraints
    r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_model_part.GetNode(1), DISTANCE, r_model_part.GetNode(3), DISTANCE, 1.0, 1.0);
    r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_model_part.GetNode(1), DISTANCE, r_model_part.GetNode(4), DISTANCE, 1.0, 1.0);

    // Set up the DOFs array
    DofArrayUtilities::DofsArrayType dofs_array;
    DofArrayUtilities::SetUpDofArray(r_model_part, dofs_array);

    // Set up the effective DOFs array
    DofArrayUtilities::DofsArrayType eff_dofs_array;
    DofArrayUtilities::SetUpEffectiveDofArray(r_model_part, dofs_array, eff_dofs_array);

    // Check results
    KRATOS_CHECK_EQUAL(eff_dofs_array.size(), 14);
}

KRATOS_TEST_CASE_IN_SUITE(SetDofEquationIds, KratosCoreFastSuite)
{
    // Set a mesh with DISTANCE DOFs elements
    Model model;
    auto& r_model_part = model.CreateModelPart("TestModelPart");
    SetTestDofArrayModelPart(r_model_part);

    // Set up the DOFs array
    DofArrayUtilities::DofsArrayType dofs_array;
    DofArrayUtilities::SetUpDofArray(r_model_part, dofs_array);

    // Set up the DOF equation ids
    DofArrayUtilities::SetDofEquationIds(dofs_array);

    // Set up the effective DOF equation ids
    // As the effective DOF array is the same of the "standard" the effective DOFs will be the same
    DofArrayUtilities::SetEffectiveDofEquationIds(dofs_array, dofs_array);

    // Check results
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(1).GetDof(DISTANCE).EquationId(), 0);
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(1).GetDof(DISTANCE).EffectiveEquationId(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(SetDofEquationIdsWithConstraints, KratosCoreFastSuite)
{
    // Set a mesh with DISTANCE DOFs elements
    Model model;
    auto& r_model_part = model.CreateModelPart("TestModelPart");
    SetTestDofArrayModelPart(r_model_part);

    // Create some constraints
    r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_model_part.GetNode(1), DISTANCE, r_model_part.GetNode(3), DISTANCE, 1.0, 1.0);
    r_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_model_part.GetNode(1), DISTANCE, r_model_part.GetNode(4), DISTANCE, 1.0, 1.0);

    // Set up the DOFs array
    DofArrayUtilities::DofsArrayType dofs_array;
    DofArrayUtilities::SetUpDofArray(r_model_part, dofs_array);

    // Set up the effective DOFs array
    DofArrayUtilities::DofsArrayType eff_dofs_array;
    DofArrayUtilities::SetUpEffectiveDofArray(r_model_part, dofs_array, eff_dofs_array);

    // Set up the DOF equation ids
    DofArrayUtilities::SetDofEquationIds(dofs_array);

    // Set up the effective DOF equation ids
    // As the effective DOF array is the same of the "standard" the effective DOFs will be the same
    DofArrayUtilities::SetEffectiveDofEquationIds(dofs_array, eff_dofs_array);

    // Check results
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(1).GetDof(DISTANCE).EquationId(), 0);
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(5).GetDof(DISTANCE).EquationId(), 4);
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(1).GetDof(DISTANCE).EffectiveEquationId(), 0);
    KRATOS_CHECK_EQUAL(r_model_part.GetNode(5).GetDof(DISTANCE).EffectiveEquationId(), 2);
}

} // namespace Kratos::Testing
