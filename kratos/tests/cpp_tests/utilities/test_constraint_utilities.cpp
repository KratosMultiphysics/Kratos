//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "testing/testing.h"
#include "includes/model_part.h"

// Utilities
#include "utilities/constraint_utilities.h"
#include "utilities/cpp_tests_utilities.h"

namespace Kratos::Testing
{

using DofsArrayType = ModelPart::DofsArrayType;

KRATOS_TEST_CASE_IN_SUITE(ConstraintUtilitiesComputeActiveDofs, KratosCoreFastSuite)
{
    // Create the test model part
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("Main");
    r_model_part.AddNodalSolutionStepVariable(TEMPERATURE);
    CppTestsUtilities::Create2DGeometry(r_model_part, "Element2D3N", true, true);

    // Define the active dofs vector
    std::vector<int> active_dofs;

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&TEMPERATURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(TEMPERATURE);
    }

    // Set the auxiliary arrays
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.pGetDof(TEMPERATURE)->SetEquationId(r_node.Id() - 1);
    }
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(r_model_part.NumberOfNodes());
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(TEMPERATURE));
    }
    aux_dof_set.Sort();

    // Adding MPC
    auto p_mpc_1 = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, r_model_part.GetNode(1), TEMPERATURE, r_model_part.GetNode(2), TEMPERATURE, 1.0, 0.0);
    auto p_mpc_2 = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, r_model_part.GetNode(3), TEMPERATURE, r_model_part.GetNode(4), TEMPERATURE, 1.0, 0.0);

    // Fix nodes
    r_model_part.GetNode(1).Fix(TEMPERATURE);
    r_model_part.GetNode(3).Fix(TEMPERATURE);

    // Compute the active dofs
    ConstraintUtilities::ComputeActiveDofs(r_model_part, active_dofs, aux_dof_set);

    // Check the results
    KRATOS_EXPECT_EQ(active_dofs.size(), 6);
    unsigned int counter = 0;
    for (int i : active_dofs) {
        if (counter == 0 || counter == 1 || counter == 2 || counter == 3) {
            KRATOS_EXPECT_EQ(i, 0);
        } else {
            KRATOS_EXPECT_EQ(i, 1);
        }
        ++counter;
    }
}

} // namespace Kratos::Testing.