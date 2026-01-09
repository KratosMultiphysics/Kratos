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
#include "constraints/linear_master_slave_constraint.h"
#include "containers/model.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
#include "future/containers/define_linear_algebra_serial.h"
#include "future/linear_solvers/amgcl_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "future/solving_strategies/strategies/linear_strategy.h"
#include "test_utilities/solving_strategies_test_utilities.h"
#endif

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto& r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 2);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 3);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithJumpConstraintEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3
    // Create a constraint with jump to impose u_3 = u_1 + 2
    // Note that as u_1 is fixed to 1.0 this actually imposes the Laplacian solution
    // This checks the MPC with fixed DOF but ensuring a non-overconstrained and compatible system of equations
    const double jump = 2.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(3);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 1);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 1);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 1);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithJumpConstraintBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3
    // Create a constraint with jump to impose u_3 = u_1 + 2
    // Note that as u_1 is fixed to 1.0 this actually imposes the Laplacian solution
    // This checks the MPC with fixed DOF but ensuring a non-overconstrained and compatible system of equations
    const double jump = 2.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(3);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 2);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithPeriodicityConstraintEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a periodicity constraint with jump to impose u_4 = u_1 + 1
    const double jump = 1.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(4);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 0.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 2);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 4.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 5.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithPeriodicityConstraintBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a periodicity constraint with jump to impose u_4 = u_1 + 1
    const double jump = 1.0;
    const double weight = 1.0;
    auto& r_master_node = r_test_model_part.GetNode(1);
    auto& r_slave_node = r_test_model_part.GetNode(4);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node, DISTANCE, r_slave_node, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 0.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 3);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 0.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 4.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 5.0/3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithMultipleDofsConstraintsEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a constraint involving multiple master DOFs such that u_4 = 0.5*u_1 + 0.5*u_2 + 1.0
    std::vector<typename Dof<double>::Pointer> slave_dofs(1);
    std::vector<typename Dof<double>::Pointer> master_dofs(2);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);

    Vector jump_vector(1);
    jump_vector[0] = 1.0;

    Matrix relation_matrix(1, 2);
    relation_matrix(0,0) = 0.5;
    relation_matrix(0,1) = 0.5;

    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, jump_vector);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 2);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithMultipleDofsConstraintsBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create a constraint involving multiple master DOFs such that u_4 = 0.5*u_1 + 0.5*u_2 + 1.0
    std::vector<typename Dof<double>::Pointer> slave_dofs(1);
    std::vector<typename Dof<double>::Pointer> master_dofs(2);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);

    Vector jump_vector(1);
    jump_vector[0] = 1.0;

    Matrix relation_matrix(1, 2);
    relation_matrix(0,0) = 0.5;
    relation_matrix(0,1) = 0.5;

    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, jump_vector);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs
    auto p_node_1 = r_test_model_part.pGetNode(1);
    p_node_1->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_1->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 3);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 3);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 3.5, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 3.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithTieConstraintsEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the nodes to apply the tie constraints
    auto p_node_5 = r_test_model_part.CreateNewNode(5, 0.0, 0.0, 0.0);
    auto p_node_6 = r_test_model_part.CreateNewNode(6, 0.0, 0.0, 0.0);
    p_node_5->AddDof(DISTANCE);
    p_node_6->AddDof(DISTANCE);

    // Mesh:  x    x ---- x ---- x ---- x    x
    // DOFs: u_5  u_1 -- u_2 -- u_3 -- u_4  u_6
    // Create a tie constraint at each side of the domain s.t. u_1 = u_5 and u_4 = u_6
    // Note that these nodes are not connected to the domain (i.e., are somehow "flying")
    const double jump = 0.0;
    const double weight = 1.0;

    auto& r_slave_node_a = r_test_model_part.GetNode(1);
    auto& r_master_node_a = r_test_model_part.GetNode(5);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node_a, DISTANCE, r_slave_node_a, DISTANCE, weight, jump);

    auto& r_slave_node_b = r_test_model_part.GetNode(4);
    auto& r_master_node_b = r_test_model_part.GetNode(6);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_master_node_b, DISTANCE, r_slave_node_b, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs to the tying ("flying") nodes
    p_node_5->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_5->Fix(DISTANCE);

    p_node_6->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_6->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 2);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 2);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithTieConstraintsBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create the nodes to apply the tie constraints
    auto p_node_5 = r_test_model_part.CreateNewNode(5, 0.0, 0.0, 0.0);
    auto p_node_6 = r_test_model_part.CreateNewNode(6, 0.0, 0.0, 0.0);
    p_node_5->AddDof(DISTANCE);
    p_node_6->AddDof(DISTANCE);

    // Mesh:  x    x ---- x ---- x ---- x    x
    // DOFs: u_5  u_1 -- u_2 -- u_3 -- u_4  u_6
    // Create a tie constraint at each side of the domain s.t. u_1 = u_5 and u_4 = u_6
    // Note that these nodes are not connected to the domain (i.e., are somehow "flying")
    const double jump = 0.0;
    const double weight = 1.0;

    auto& r_slave_node_a = r_test_model_part.GetNode(1);
    auto& r_master_node_a = r_test_model_part.GetNode(5);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, r_master_node_a, DISTANCE, r_slave_node_a, DISTANCE, weight, jump);

    auto& r_slave_node_b = r_test_model_part.GetNode(4);
    auto& r_master_node_b = r_test_model_part.GetNode(6);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 2, r_master_node_b, DISTANCE, r_slave_node_b, DISTANCE, weight, jump);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Apply Dirichlet BCs to the tying ("flying") nodes
    p_node_5->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_5->Fix(DISTANCE);

    p_node_6->FastGetSolutionStepValue(DISTANCE, 0) = 1.0;
    p_node_6->Fix(DISTANCE);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();

    // Check array sizes
    const auto &r_linear_system_container = p_strategy->GetImplicitStrategyDataContainer();
    KRATOS_CHECK_EQUAL(r_linear_system_container.pRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pLhs->size2(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveRhs->size(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size1(), 4);
    KRATOS_CHECK_EQUAL(r_linear_system_container.pEffectiveLhs->size2(), 4);

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 2.0, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), 1.0, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithRigidBodyMotionConstraintEliminationBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create an average constraint to enforce the average solution to be zero (i.e., no rigid body motion)
    // We do this as u_4 = - u_1 - u_2 - u_3 as u_1 + u_2 + u_3 + u_4 = u_1 + u_2 + u_3 + (- u_1 - u_2 - u_3) = 0
    typename LinearMasterSlaveConstraint::DofPointerVectorType slave_dofs(1);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    typename LinearMasterSlaveConstraint::DofPointerVectorType master_dofs(3);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);
    master_dofs[2] = r_test_model_part.GetNode(3).pGetDof(DISTANCE);
    Matrix relation_matrix(1,3);
    relation_matrix(0,0) = -1.0;
    relation_matrix(0,1) = -1.0;
    relation_matrix(0,2) = -1.0;
    Vector constant_vector = ZeroVector(1);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, constant_vector);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
#endif
}

KRATOS_TEST_CASE_IN_SUITE(LinearStrategyWithRigidBodyMotionConstraintBlockBuild, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 3;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Mesh:  x ---- x ---- x ---- x
    // DOFs: u_1 -- u_2 -- u_3 -- u_4
    // Create an average constraint to enforce the average solution to be zero (i.e., no rigid body motion)
    // We do this as u_4 = - u_1 - u_2 - u_3 as u_1 + u_2 + u_3 + u_4 = u_1 + u_2 + u_3 + (- u_1 - u_2 - u_3) = 0
    typename LinearMasterSlaveConstraint::DofPointerVectorType slave_dofs(1);
    slave_dofs[0] = r_test_model_part.GetNode(4).pGetDof(DISTANCE);
    typename LinearMasterSlaveConstraint::DofPointerVectorType master_dofs(3);
    master_dofs[0] = r_test_model_part.GetNode(1).pGetDof(DISTANCE);
    master_dofs[1] = r_test_model_part.GetNode(2).pGetDof(DISTANCE);
    master_dofs[2] = r_test_model_part.GetNode(3).pGetDof(DISTANCE);
    Matrix relation_matrix(1,3);
    relation_matrix(0,0) = -1.0;
    relation_matrix(0,1) = -1.0;
    relation_matrix(0,2) = -1.0;
    Vector constant_vector = ZeroVector(1);
    r_test_model_part.CreateNewMasterSlaveConstraint(
        "LinearMasterSlaveConstraint", 1, master_dofs, slave_dofs, relation_matrix, constant_vector);

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);

    // Create the linear solver
    Parameters amgcl_settings = Parameters(R"({
    })");
    using AMGCLSolverType = Future::AMGCLSolver<Future::SerialLinearAlgebraTraits>;
    using LinearSolverType = Future::LinearSolver<Future::SerialLinearAlgebraTraits>;
    typename LinearSolverType::Pointer p_amgcl_solver = Kratos::make_shared<AMGCLSolverType>(amgcl_settings);

    // Create the strategy
    Parameters strategy_settings = Parameters(R"({
    })");
    auto p_strategy = Kratos::make_unique<Future::LinearStrategy<Future::SerialLinearAlgebraTraits>>(r_test_model_part, p_scheme, p_amgcl_solver);

    // Solve the problem
    p_strategy->Initialize();
    p_strategy->Check();
    p_strategy->Predict();
    p_strategy->InitializeSolutionStep();
    p_strategy->SolveSolutionStep();
    p_strategy->FinalizeSolutionStep();
    p_strategy->Clear();

    // Check results
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(1).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(2).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(3).FastGetSolutionStepValue(DISTANCE), 0.125, 1.0e-12);
    KRATOS_CHECK_NEAR(r_test_model_part.GetNode(4).FastGetSolutionStepValue(DISTANCE), -0.125, 1.0e-12);
#endif
}

}  // namespace Kratos::Testing.

