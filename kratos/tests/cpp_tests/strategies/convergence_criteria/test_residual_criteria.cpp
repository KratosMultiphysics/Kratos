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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "spaces/ublas_space.h"

namespace Kratos::Testing
{
using GeometryType = Geometry<Node>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

using DofsArrayType = ModelPart::DofsArrayType;

using ResidualCriteriaType = ResidualCriteria<SparseSpaceType, LocalSpaceType>;

void GenerateTestResidualCriteriaModelPart(ModelPart& rModelPart)
{
    // Model part settings
    rModelPart.SetBufferSize(1);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);

    // Create the auxiliary set of nodes
    const unsigned int n_nodes = 10;
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto p_node = rModelPart.CreateNewNode(i_node + 1, 0.0, 0.0, 0.0); // Coordinates do not matter in this test
        rModelPart.AddNode(p_node);
    }
}

/**
 * Checks the residual criteria
 */
KRATOS_TEST_CASE_IN_SUITE(ResidualCriteria, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestResidualCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
    }

    // Create the residual criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    auto residual_criteria = ResidualCriteriaType(rel_tol, abs_tol);

    // Set the auxiliary arrays
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.pGetDof(PRESSURE)->SetEquationId(r_node.Id() - 1);
    }
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(10);
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
    }
    aux_dof_set.Sort();
    typename ResidualCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename ResidualCriteriaType::TSystemVectorType b(10);
    typename ResidualCriteriaType::TSystemVectorType Dx; // Only required to match the API

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;

    // Set the auxiliary fake data to check the convergence for (initialization)
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        b[i] = aux_val;
        i++;
    }

    // Initialize the solution step
    residual_criteria.InitializeSolutionStep(r_model_part, aux_dof_set, A, Dx, b);

    // Check convergence (failing)
    bool convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_FALSE(convergence)

    // Set the auxiliary fake data to check the convergence for (passing)
    i = 0;
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        b[i] = aux_val / 200.0;
        i++;
    }

    // Check convergence (passing)
    convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

/**
 * Checks the residual criteria with MPC
 */
KRATOS_TEST_CASE_IN_SUITE(ResidualCriteriaWithMPC, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestResidualCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
    }

    // Create the residual criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    auto residual_criteria = ResidualCriteriaType(rel_tol, abs_tol);

    // Set the auxiliary arrays
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.pGetDof(PRESSURE)->SetEquationId(r_node.Id() - 1);
    }
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(10);
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
    }
    aux_dof_set.Sort();

    // Adding MPC
    auto p_mpc_1 = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 1, r_model_part.GetNode(1), PRESSURE, r_model_part.GetNode(2), PRESSURE, 1.0, 0.0);
    auto p_mpc_2 = r_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", 2, r_model_part.GetNode(3), PRESSURE, r_model_part.GetNode(4), PRESSURE, 1.0, 0.0);

    // Fix nodes
    r_model_part.GetNode(1).Fix(PRESSURE);
    r_model_part.GetNode(3).Fix(PRESSURE);

    typename ResidualCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename ResidualCriteriaType::TSystemVectorType b(10);
    typename ResidualCriteriaType::TSystemVectorType Dx; // Only required to match the API

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;

    // Set the auxiliary fake data to check the convergence for (initialization)
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        b[i] = aux_val;
        i++;
    }

    // Initialize the solution step
    residual_criteria.InitializeSolutionStep(r_model_part, aux_dof_set, A, Dx, b);

    // Check convergence (failing)
    bool convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_FALSE(convergence)

    // Set the auxiliary fake data to check the convergence for (passing)
    i = 0;
    for (auto& r_node : r_model_part.Nodes()) {
        // Non fixed or MPC nodes
        if (i > 3) {
            const double aux_val = r_node.Id() * aux_constant;
            b[i] = aux_val / 1000.0;
        }
        i++;
    }

    // Check convergence (passing)
    convergence = residual_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

} // namespace Kratos::Testing
