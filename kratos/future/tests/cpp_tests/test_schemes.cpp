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
#include "future/linear_solvers/amgcl_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"
#include "future/solving_strategies/schemes/static_scheme.h"
#include "test_utilities/solving_strategies_test_utilities.h"
#endif

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(StaticSchemeBuild1D, KratosCoreFastSuite)
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
    using SchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    Future::LinearSystemContainer<CsrMatrix<>, SystemVector<>> linear_system_container;

    // Call the initialize solution step (note that this sets all the arrays above)
    p_scheme->Initialize(linear_system_container);
    p_scheme->InitializeSolutionStep(linear_system_container);

    // Call the build
    auto p_lhs = linear_system_container.pLhs;
    auto p_rhs = linear_system_container.pRhs;
    p_scheme->Build(*p_lhs, *p_rhs);

    // Check resultant matrices
    const double tol = 1.0e-12;
    std::vector<double> expected_rhs = {0.5,1.0,0.5};
    BoundedMatrix<double,3,3> expected_lhs;
    expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    for (unsigned int i = 0; i < p_lhs->size1(); ++i) {
        for (unsigned int j = 0; j < p_lhs->size2(); ++j) {
            const double expected_val = expected_lhs(i,j);
            if (std::abs(expected_val) > tol) {
                KRATOS_CHECK_NEAR(p_lhs->operator()(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
            }
        }
    }
#else
    true;
#endif
}

KRATOS_TEST_CASE_IN_SUITE(StaticSchemeBuild2D, KratosCoreFastSuite)
{
#ifdef KRATOS_USE_FUTURE
    // Set up the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");
    const std::size_t num_elems_x = 25;
    const std::size_t num_elems_y = 25;
    const double elem_size_x = 1.0;
    const double elem_size_y = 1.0;

    BuiltinTimer timer_mesh;
    SolvingStrategiesTestUtilities::SetUpTestModelPart2D(num_elems_x, num_elems_y, elem_size_x, elem_size_y, r_test_model_part);
    std::cout << "Mesh generation time: " << timer_mesh << std::endl;
    std::cout << "Number of elements: " << num_elems_x * num_elems_y << std::endl;

    // Create the scheme
    Parameters scheme_settings = Parameters(R"({
        "echo_level" : 1,
        "build_settings" : {
            "name" : "block_builder"
        }
    })");
    using SchemeType = Future::StaticScheme<CsrMatrix<>, SystemVector<>, SparseContiguousRowGraph<>>;
    auto p_scheme = Kratos::make_unique<SchemeType>(r_test_model_part, scheme_settings);

    // Set up the matrix graph and arrays
    // Note that in a standard case this happens at the strategy level
    Future::LinearSystemContainer<CsrMatrix<>, SystemVector<>> linear_system_container;

    // Call the initialize solution step (note that this sets all the arrays above)
    p_scheme->Initialize(linear_system_container);
    p_scheme->InitializeSolutionStep(linear_system_container);

    // Call the build
    auto p_lhs = linear_system_container.pLhs;
    auto p_rhs = linear_system_container.pRhs;

    sleep(30);

    BuiltinTimer timer_build;
    for(unsigned int i=0; i<20; ++i)
        p_scheme->Build(*p_lhs, *p_rhs);
    std::cout << "Build time: " << timer_build << std::endl;

    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_safe;
    for(unsigned int i=0; i<20; ++i)
        p_scheme->BuildWithSafeAssemble(*p_lhs, *p_rhs);
    std::cout << "Build w/ safe assemble time: " << timer_build_safe << std::endl;

#ifdef KRATOS_USE_TBB
    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_thread_local;
    for(unsigned int i=0; i<20; ++i)
        p_scheme->BuildWithThreadLocal(*p_lhs, *p_rhs);
    std::cout << "Build w/ thread local: " << timer_build_thread_local << std::endl;

    p_lhs->SetValue(0.0);
    p_rhs->SetValue(0.0);

    sleep(30);

    BuiltinTimer timer_build_local_allocation;
    for(unsigned int i=0; i<20; ++i)
    p_scheme->BuildWithLocalAllocation(*p_lhs, *p_rhs);
    std::cout << "Build w/ local allocation: " << timer_build_local_allocation << std::endl;
#endif

    // // Check resultant matrices
    // const double tol = 1.0e-12;
    // std::vector<double> expected_rhs = {0.5,1.0,0.5};
    // BoundedMatrix<double,3,3> expected_lhs;
    // expected_lhs(0,0) = 1.0; expected_lhs(0,1) = -1.0; expected_lhs(0,2) = 0.0;
    // expected_lhs(1,0) = -1.0; expected_lhs(1,1) = 2.0; expected_lhs(1,2) = -1.0;
    // expected_lhs(2,0) = 0.0; expected_lhs(2,1) = -1.0; expected_lhs(2,2) = 1.0;
    // KRATOS_CHECK_VECTOR_NEAR((*p_rhs), expected_rhs, tol); // Note that as there are not non-zero entries in the sparse vector we can use the standard macro
    // for (unsigned int i = 0; i < p_lhs->size1(); ++i) {
    //     for (unsigned int j = 0; j < p_lhs->size2(); ++j) {
    //         const double expected_val = expected_lhs(i,j);
    //         if (std::abs(expected_val) > tol) {
    //             KRATOS_CHECK_NEAR(p_lhs->operator()(i,j), expected_val, tol); // Note that we check if the expected value is non-zero as this is a CSR matrix
    //         }
    //     }
    // }
#else
    true;
#endif
}

}  // namespace Kratos::Testing.

