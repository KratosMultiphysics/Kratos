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
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "spaces/ublas_space.h"

namespace Kratos::Testing
{
using GeometryType = Geometry<Node>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

using DofsArrayType = ModelPart::DofsArrayType;

using DisplacementCriteriaType = DisplacementCriteria<SparseSpaceType, LocalSpaceType>;

void GenerateTestDisplacementCriteriaModelPart(ModelPart& rModelPart)
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
 * Checks the displacement criteria
 */
KRATOS_TEST_CASE_IN_SUITE(DisplacementCriteria, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestDisplacementCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
    }

    // Create the displacement criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    auto displacement_criteria = DisplacementCriteriaType(rel_tol, abs_tol);

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
    typename DisplacementCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename DisplacementCriteriaType::TSystemVectorType b; // Only required to match the API
    typename DisplacementCriteriaType::TSystemVectorType Dx(10);

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;

    // Set the auxiliary fake data to check the convergence for (failing)
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_val;
        Dx[i] = aux_val;
        i++;
    }

    // Check convergence (failing)
    bool convergence = displacement_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_FALSE(convergence)

    // Set the auxiliary fake data to check the convergence for (passing)
    i = 0;
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_val;
        Dx[i] = aux_val / 1000.0;
        i++;
    }

    // Check convergence (passing)
    convergence = displacement_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

} // namespace Kratos::Testing
