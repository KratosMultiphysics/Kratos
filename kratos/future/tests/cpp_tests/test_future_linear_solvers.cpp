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

// Project includes
#include "testing/testing.h"
#include "future/linear_solvers/amgcl_solver.h"
#include "future/linear_solvers/skyline_lu_factorization_solver.h"

namespace Kratos::Testing
{

namespace
{
    std::size_t SetLinearSystem(
        CsrMatrix<>& rLHS,
        SystemVector<>& rRHS)
    {
        // Set the matrix map and assign it to the provided LHS matrix
        typename CsrMatrix<double>::MatrixMapType matrix_map;
        matrix_map[{0, 0}] = 10; matrix_map[{0, 3}] = -2;
        matrix_map[{1, 1}] = 8; matrix_map[{1, 2}] = -1;
        matrix_map[{2, 2}] = 5; matrix_map[{2, 1}] = -1; matrix_map[{2, 4}] = -2;
        matrix_map[{3, 3}] = 7; matrix_map[{3, 0}] = -2; matrix_map[{3, 4}] = -1;
        matrix_map[{4, 4}] = 6; matrix_map[{4, 2}] = -2; matrix_map[{4, 3}] = -1;
        rLHS = CsrMatrix<double>(matrix_map);

        // Set the RHS vector
        DenseVector<double> aux_data(5);
        aux_data[0] = 3.0;
        aux_data[1] = 7.0;
        aux_data[2] = 2.0;
        aux_data[3] = 5.0;
        aux_data[4] = 1.0;
        rRHS = SystemVector<>(aux_data);

        // Return the system size
        return rRHS.size();
    }
}

KRATOS_TEST_CASE_IN_SUITE(FutureLinearSolversSkylineLUFactorizationSolver, KratosCoreFutureSuite)
{
    // Set up the system to be solved
    CsrMatrix<> LHS;
    SystemVector<> RHS;
    const std::size_t system_size = SetLinearSystem(LHS, RHS);

    // Set the linear solver to be tested
    Parameters skyline_lu_settings(R"({
    })");
    auto p_linear_solver = Kratos::make_unique<Future::SkylineLUFactorizationSolver<CsrMatrix<>, SystemVector<>>>(skyline_lu_settings);

    // Solve the problem
    SystemVector<> sol(system_size);
    p_linear_solver->Solve(LHS, sol, RHS);

    // Check the obtained results
    std::vector<double> ref_sol = {0.487946221604, 0.979601298099, 0.836810384794, 0.93973110802, 0.602225312935};
    KRATOS_EXPECT_VECTOR_NEAR(sol.data(), ref_sol, 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(FutureLinearSolversAmgcl, KratosCoreFutureSuite)
{
    // Set up the system to be solved
    CsrMatrix<> LHS;
    SystemVector<> RHS;
    const std::size_t system_size = SetLinearSystem(LHS, RHS);

    // Set the linear solver to be tested
    Parameters amgcl_settings(R"({
        "solver_type"                    : "AMGCL",
        "max_iteration"                  : 100,
        "verbosity"                      : 1,
        "tolerance"                      : 1e-6
    })");
    auto p_linear_solver = Kratos::make_unique<Future::AMGCLSolver<CsrMatrix<>, SystemVector<>>>(amgcl_settings);

    // Solve the problem
    SystemVector<> sol(system_size);
    p_linear_solver->Solve(LHS, sol, RHS);

    // Check the obtained results
    std::vector<double> ref_sol = {0.487946221604, 0.979601298099, 0.836810384794, 0.93973110802, 0.602225312935};
    KRATOS_EXPECT_VECTOR_NEAR(sol.data(), ref_sol, 1e-12);
}

}