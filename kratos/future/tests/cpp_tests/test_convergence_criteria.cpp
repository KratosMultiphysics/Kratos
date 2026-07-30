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
#include "includes/element.h"
#include "includes/model_part.h"
#include "testing/testing.h"

#ifdef KRATOS_USE_FUTURE
#include "future/containers/define_linear_algebra_serial.h"
#include "future/solving_strategies/convergence_criteria/solution_criteria.h"
#include "future/solving_strategies/convergence_criteria/residual_criteria.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "test_utilities/solving_strategies_test_utilities.h"
#endif

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FutureSolutionCriteria, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create an implicit strategy data container to hold DOFs sets and system vectors
    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;

    // Create an auxiliary scheme to initialize the DOFs sets
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);
    p_scheme->Initialize(strategy_data_container);

    // Create an effective solution increment vector
    DenseVector<double> aux_data(3);
    aux_data[0] = 3.0;
    aux_data[1] = 7.0;
    aux_data[2] = 2.0;
    Future::SerialLinearAlgebraTraits::VectorType eff_dx(aux_data);
    auto p_eff_dx = Kratos::make_shared<Future::SerialLinearAlgebraTraits::VectorType>(eff_dx);
    auto p_eff_lin_sys = strategy_data_container.pGetEffectiveLinearSystem();
    p_eff_lin_sys->pSetVector(p_eff_dx, Future::LinearSystemTags::DenseVectorTag::Dx);

    // Set solution values and fix one node to check that the convergence criteria works with fixed DOFs
    for (auto& r_node : r_test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISTANCE) = 2.0 * static_cast<double>(r_node.Id());
    }
    r_test_model_part.pGetNode(1)->Fix(DISTANCE);

    // Create solution criteria
    Parameters solution_criteria_settings = Parameters(R"({
        "echo_level" : 0,
        "variable_name" : "DISTANCE",
        "relative_tolerance" : 1.0e-4,
        "absolute_tolerance" : 1.0e-6
    })");
    auto p_convergence_criteria = Kratos::make_unique<Future::SolutionCriteria<Future::SerialLinearAlgebraTraits>>(r_test_model_part, solution_criteria_settings);

    // Call convergence criteria check
    const bool is_converged = p_convergence_criteria->IsConverged(strategy_data_container);
    const double res_norm = r_test_model_part.GetProcessInfo()[RESIDUAL_NORM];
    const double conv_ratio = r_test_model_part.GetProcessInfo()[CONVERGENCE_RATIO];
    KRATOS_EXPECT_FALSE(is_converged);
    KRATOS_EXPECT_NEAR(res_norm, 5.14781507049, 1e-10);
    KRATOS_EXPECT_NEAR(conv_ratio, 1.00956959603, 1e-10);
    KRATOS_EXPECT_FALSE(p_convergence_criteria->RequiresResidual());
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(FutureResidualCriteria, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create an implicit strategy data container to hold DOFs sets and system vectors
    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;

    // Create an auxiliary scheme to initialize the DOFs sets
    Parameters scheme_settings = Parameters(R"({
        "build_settings" : {
            "name" : "elimination_builder"
        }
    })");
    auto p_scheme = Kratos::make_shared<Future::StaticScheme<Future::SerialLinearAlgebraTraits>>(r_test_model_part, scheme_settings);
    p_scheme->Initialize(strategy_data_container);

    // Create an implicit strategy data container with a fake effective residual vector
    DenseVector<double> aux_data(3);
    aux_data[0] = 3.0;
    aux_data[1] = 7.0;
    aux_data[2] = 2.0;
    Future::SerialLinearAlgebraTraits::VectorType eff_rhs(aux_data);
    auto p_eff_rhs = Kratos::make_shared<Future::SerialLinearAlgebraTraits::VectorType>(eff_rhs);
    auto p_eff_lin_sys = strategy_data_container.pGetEffectiveLinearSystem();
    p_eff_lin_sys->pSetVector(p_eff_rhs, Future::LinearSystemTags::DenseVectorTag::RHS);

    // Create solution criteria
    Parameters residual_criteria_settings = Parameters(R"({
        "echo_level" : 0,
        "relative_tolerance" : 1.0e-4,
        "absolute_tolerance" : 1.0e-6
    })");
    auto p_convergence_criteria = Kratos::make_unique<Future::ResidualCriteria<Future::SerialLinearAlgebraTraits>>(r_test_model_part, residual_criteria_settings);

    // Call InitializeSolutionStep (sets the initial residual norm)
    p_convergence_criteria->InitializeSolutionStep(strategy_data_container);

    // Change the effective residual vector
    eff_rhs[0] = 10.0;
    eff_rhs[1] = 11.0;
    eff_rhs[2] = 12.0;
    p_eff_rhs = Kratos::make_shared<Future::SerialLinearAlgebraTraits::VectorType>(eff_rhs);
    p_eff_lin_sys->pSetVector(p_eff_rhs, Future::LinearSystemTags::DenseVectorTag::RHS);

    // Call convergence criteria check
    const bool is_converged = p_convergence_criteria->IsConverged(strategy_data_container);
    const double res_norm = r_test_model_part.GetProcessInfo()[RESIDUAL_NORM];
    const double conv_ratio = r_test_model_part.GetProcessInfo()[CONVERGENCE_RATIO];
    KRATOS_EXPECT_FALSE(is_converged);
    KRATOS_EXPECT_NEAR(res_norm, 11.0302614052, 1e-10);
    KRATOS_EXPECT_NEAR(conv_ratio, 2.4263340195, 1e-10);
    KRATOS_EXPECT_TRUE(p_convergence_criteria->RequiresResidual());
#else
    true;
#endif
}

}  // namespace Kratos::Testing.

