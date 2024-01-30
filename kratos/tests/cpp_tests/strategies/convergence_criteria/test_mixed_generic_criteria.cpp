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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "geometries/point_3d.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "spaces/ublas_space.h"

namespace Kratos::Testing
{
using NodeType = Node;
using GeometryType = Geometry<NodeType>;
using LocalSpaceType = UblasSpace<double, Matrix, Vector>;
using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

using DofsArrayType = ModelPart::DofsArrayType;

using MixedGenericCriteriaType = MixedGenericCriteria<SparseSpaceType, LocalSpaceType>;
using ConvergenceVariableListType = typename MixedGenericCriteriaType::ConvergenceVariableListType;

void GenerateTestConvergenceCriteriaModelPart(ModelPart& rModelPart)
{
    // Model part settings
    rModelPart.SetBufferSize(1);
    rModelPart.AddNodalSolutionStepVariable(PRESSURE);
    rModelPart.AddNodalSolutionStepVariable(VELOCITY);
    rModelPart.AddNodalSolutionStepVariable(TEMPERATURE);

    // Create the auxiliary set of nodes
    const unsigned int n_nodes = 10;
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        auto p_node = rModelPart.CreateNewNode(i_node + 1, 0.0, 0.0, 0.0); // Coordinates do not matter in this test
        rModelPart.AddNode(p_node);
    }
}

/**
 * Checks the mixed generic criteria with two double variables
 */
KRATOS_TEST_CASE_IN_SUITE(MixedGenericCriteriaDoubleDouble, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestConvergenceCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&TEMPERATURE);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
        r_node.AddDof(TEMPERATURE);
    }

    // Create the mixed generic criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    VariableData* p_pres = &PRESSURE;
    VariableData* p_temp = &TEMPERATURE;
    ConvergenceVariableListType convergence_settings;
    convergence_settings.push_back(std::make_tuple(p_pres, rel_tol, abs_tol));
    convergence_settings.push_back(std::make_tuple(p_temp, rel_tol, abs_tol));
    auto mixed_generic_criteria = MixedGenericCriteriaType(convergence_settings);

    // Set the auxiliary arrays
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(20);
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
        aux_dof_set.push_back(r_node.pGetDof(TEMPERATURE));
    }
    aux_dof_set.Sort();
    typename MixedGenericCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType b; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType Dx(20);

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_val;
        r_node.FastGetSolutionStepValue(TEMPERATURE) = aux_val;
        Dx[2 * i] = aux_val / 100.0;
        Dx[2 * i + 1] = 2.0 * aux_val / 100.0;
        i++;
    }

    // Check convergence
    const bool convergence = mixed_generic_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

/**
 * Checks the mixed generic criteria with a double and a vector variable
 */
KRATOS_TEST_CASE_IN_SUITE(MixedGenericCriteriaDoubleVector, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestConvergenceCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&VELOCITY);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
    }

    // Create the mixed generic criteria
    const double rel_tol = 1.0e-3;
    const double abs_tol = 1.0e-5;
    VariableData* p_pres = &PRESSURE;
    VariableData* p_temp = &VELOCITY;
    ConvergenceVariableListType convergence_settings;
    convergence_settings.push_back(std::make_tuple(p_pres, rel_tol, abs_tol));
    convergence_settings.push_back(std::make_tuple(p_temp, rel_tol, abs_tol));
    auto mixed_generic_criteria = MixedGenericCriteriaType(convergence_settings);

    // Set the auxiliary arrays
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(30);
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
        aux_dof_set.push_back(r_node.pGetDof(VELOCITY_X));
        aux_dof_set.push_back(r_node.pGetDof(VELOCITY_Y));
    }
    aux_dof_set.Sort();
    typename MixedGenericCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType b; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType Dx(30);

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;
    array_1d<double, 3> aux_vel = ZeroVector(3);
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        aux_vel[0] = 2.0 * aux_val;
        aux_vel[1] = 3.0 * aux_val;
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_val;
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        Dx[3 * i] = aux_val / 100.0;
        Dx[3 * i + 1] = 2.0 * aux_val / 100.0;
        Dx[3 * i + 2] = 3.0 * aux_val / 100.0;
        i++;
    }

    // Check convergence
    const bool convergence = mixed_generic_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

/**
 * Checks the mixed generic criteria with a double and a vector variable
 */
KRATOS_TEST_CASE_IN_SUITE(MixedGenericCriteriaDoubleVectorWithParameters, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

    // Set the test model part
    GenerateTestConvergenceCriteriaModelPart(r_model_part); // Create the geometry

    // Add the corresponding test DOFs
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
    r_model_part.GetNodalSolutionStepVariablesList().AddDof(&VELOCITY);
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.AddDof(PRESSURE);
        r_node.AddDof(VELOCITY_X);
        r_node.AddDof(VELOCITY_Y);
    }

    // Create the mixed generic criteria
    Parameters parameters = Parameters(R"(
    {
        "convergence_variables_list" : {
            "pressure" :
            {
                "variable"           : "PRESSURE",
                "relative_tolerance" : 1.0e-3,
                "absolute_tolerance" : 1.0e-5
            },
            "velocity" :
            {
                "variable"           : "VELOCITY",
                "relative_tolerance" : 1.0e-3,
                "absolute_tolerance" : 1.0e-5
            }
        }
    })" );
    auto mixed_generic_criteria = MixedGenericCriteriaType(parameters);

    // Set the auxiliary arrays
    DofsArrayType aux_dof_set;
    aux_dof_set.reserve(30);
    for (auto& r_node : r_model_part.Nodes()) {
        aux_dof_set.push_back(r_node.pGetDof(PRESSURE));
        aux_dof_set.push_back(r_node.pGetDof(VELOCITY_X));
        aux_dof_set.push_back(r_node.pGetDof(VELOCITY_Y));
    }
    aux_dof_set.Sort();
    typename MixedGenericCriteriaType::TSystemMatrixType A; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType b; // Only required to match the API
    typename MixedGenericCriteriaType::TSystemVectorType Dx(30);

    // Set the auxiliary fake data to check the convergence for
    unsigned int i = 0;
    const double aux_constant = 1.0e-3;
    array_1d<double, 3> aux_vel = ZeroVector(3);
    for (auto& r_node : r_model_part.Nodes()) {
        const double aux_val = r_node.Id() * aux_constant;
        aux_vel[0] = 2.0 * aux_val;
        aux_vel[1] = 3.0 * aux_val;
        r_node.FastGetSolutionStepValue(PRESSURE) = aux_val;
        r_node.FastGetSolutionStepValue(VELOCITY) = aux_vel;
        Dx[3 * i] = aux_val / 100.0;
        Dx[3 * i + 1] = 2.0 * aux_val / 100.0;
        Dx[3 * i + 2] = 3.0 * aux_val / 100.0;
        i++;
    }

    // Check convergence
    const bool convergence = mixed_generic_criteria.PostCriteria(r_model_part, aux_dof_set, A, Dx, b);
    KRATOS_EXPECT_TRUE(convergence)
}

} // namespace Kratos::Testing
