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
#include "future/solving_strategies/convergence_criteria/convergence_criteria.h"
#include "future/solving_strategies/convergence_criteria/solution_criteria.h"
#include "test_utilities/solving_strategies_test_utilities.h"
#endif

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(ConvergenceCriteria, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems = 2;
    const double elem_size = 1.0;
    SolvingStrategiesTestUtilities::SetUpTestModelPart1D(num_elems, elem_size, r_test_model_part);

    // Create an implicit strategy data container with a fake effective solution increment vector
    DenseVector<double> aux_data(3);
    aux_data[0] = 3.0;
    aux_data[1] = 7.0;
    aux_data[2] = 2.0;
    Future::SerialLinearAlgebraTraits::VectorType eff_dx(aux_data);
    auto p_eff_dx = Kratos::make_shared<Future::SerialLinearAlgebraTraits::VectorType>(eff_dx);

    auto p_eff_lin_sys = Kratos::make_shared<Future::LinearSystem<Future::SerialLinearAlgebraTraits>>();
    p_eff_lin_sys->pSetVector(p_eff_dx, Future::LinearSystemTags::DenseVectorTag::Dx);

    Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;
    strategy_data_container.pSetEffectiveLinearSystem(p_eff_lin_sys);

    // Set solution values and fix one node to check that the convergence criteria works with fixed DOFs
    for (auto& r_node : r_test_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISTANCE) = 2.0 * static_cast<double>(r_node.Id());
    }
    r_test_model_part.pGetNode(1)->Fix(DISTANCE);

    // Create solution criteria
    Parameters solution_criteria_parameters = Parameters(R"({
        "variable_name" : "DISTANCE",
        "relative_tolerance" : 1.0e-4,
        "absolute_tolerance" : 1.0e-6
    })");
    auto p_convergence_criteria = Kratos::make_unique<Future::SolutionCriteria<Future::SerialLinearAlgebraTraits>>(solution_criteria_parameters);

    // Call convergence criteria check
    const bool is_converged = p_convergence_criteria->IsConverged(r_test_model_part, strategy_data_container);
    const double res_norm = r_test_model_part.GetProcessInfo()[RESIDUAL_NORM];
    const double conv_ratio = r_test_model_part.GetProcessInfo()[CONVERGENCE_RATIO];

    std::cout << std::setprecision(12) << res_norm << std::endl;
    std::cout << std::setprecision(12) << conv_ratio << std::endl;

    // Parameters scheme_settings = Parameters(R"({
    //     "build_settings" : {
    //         "name" : "block_builder"
    //     }
    // })");
    // using SchemeType = Future::StaticScheme<Future::SerialLinearAlgebraTraits>;
    // auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // // Set up the matrix graph and arrays
    // // Note that in a standard case this happens at the strategy level
    // Future::ImplicitStrategyData<Future::SerialLinearAlgebraTraits> strategy_data_container;

    // // Call the initialize solution step (note that this sets all the arrays above)
    // p_scheme->Initialize(strategy_data_container);
    // p_scheme->InitializeSolutionStep(strategy_data_container);

    // // Call the build
    // const auto p_linear_system = strategy_data_container.pGetLinearSystem();
    // auto& r_lhs = *(p_linear_system->pGetMatrix(Future::LinearSystemTags::SparseMatrixTag::LHS));
    // auto& r_rhs = *(p_linear_system->pGetVector(Future::LinearSystemTags::DenseVectorTag::RHS));
    // p_scheme->Build(r_lhs, r_rhs);

    // // Check resultant matrices
    // const double tol = 1.0e-12;
    // std::vector<double> expected_rhs = {0.5,1.0,0.5};
    // BoundedMatrix<double,3,3> expected_lhs;
    // expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    // expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    // expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    // KRATOS_CHECK_VECTOR_NEAR(r_rhs, expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    // for (unsigned int i = 0; i < r_lhs.size1(); ++i) {
    //     for (unsigned int j = 0; j < r_lhs.size2(); ++j) {
    //         const double expected_val = expected_lhs(i,j);
    //         if (std::abs(expected_val) > tol) {
    //             KRATOS_CHECK_NEAR(r_lhs(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
    //         }
    //     }
    // }
#else
    true;
#endif
}

}  // namespace Kratos::Testing.

