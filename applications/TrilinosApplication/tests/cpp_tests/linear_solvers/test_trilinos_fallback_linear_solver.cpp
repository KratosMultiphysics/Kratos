//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "trilinos_space.h"
#include "tests/cpp_tests/trilinos_cpp_test_utilities.h"
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "linear_solvers/fallback_linear_solver.h"

namespace Kratos::Testing
{

template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DummyLinearSolver
    : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>;
    using SparseMatrixType = typename BaseType::SparseMatrixType;
    using VectorType = typename BaseType::VectorType;

    /// Pointer definition of DummyLinearSolver
    KRATOS_CLASS_POINTER_DEFINITION(DummyLinearSolver);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    DummyLinearSolver() : BaseType() {}

    /// Destructor.
    ~DummyLinearSolver() override = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Override Solve method to always return false
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        return false; // Indicate that the solver didn't actually solve anything
    }

    ///@}

}; // Class DummyLinearSolver

// Define the types of spaces
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;

// Define the vector and matrix types
using SparseMatrixType = typename TrilinosSparseSpaceType::MatrixType;
using VectorType = typename TrilinosSparseSpaceType::VectorType;
using DenseMatrixType = typename TrilinosLocalSpaceType::MatrixType;
using DenseVectorType = typename TrilinosLocalSpaceType::VectorType;

// Define the types of solvers
using TrilinosLinearSolverType = LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
using TrilinosDummyLinearSolverType = DummyLinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
using TrilinosFallbackLinearSolverType = FallbackLinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
using TrilinosLinearSolverFactoryType = LinearSolverFactory<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosFallbackLinearSolverConstructorSolvers, KratosTrilinosApplicationMPITestSuite)
{
    // Create the solvers
    auto p_solver1 = Kratos::make_shared<TrilinosDummyLinearSolverType>();
    Parameters amgcl_parameters = Parameters(R"({
        "solver_type": "amgcl"
    })");
    auto p_solver2 = TrilinosLinearSolverFactoryType().Create(amgcl_parameters);

    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // Create the matrix and vectors
    const int size = 2 * r_comm.Size();
    auto A = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    auto b = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto x = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Create a simple fallback solver
    std::vector<TrilinosLinearSolverType::Pointer> solvers = {p_solver1, p_solver2};
    TrilinosFallbackLinearSolverType simple_fallback_solver(solvers);

    // Solve the system
    const auto solved = simple_fallback_solver.Solve(A, x, b);

    // Check that the system was solved
    KRATOS_EXPECT_TRUE(solved);

    // Check index is 1 (solve 0 failed)
    KRATOS_EXPECT_EQ(simple_fallback_solver.GetCurrentSolverIndex(), 1);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(TrilinosFallbackLinearSolverConstructorParameters, KratosTrilinosApplicationMPITestSuite)
{
    // The data communicator
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    // Create the matrix and vectors
    const int size = 2 * r_comm.Size();
    auto A = TrilinosCPPTestUtilities::GenerateDummySparseMatrix(r_comm, size, 1.0);
    auto b = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);
    auto x = TrilinosCPPTestUtilities::GenerateDummySparseVector(r_comm, size);

    // Create a simple fallback solver
    Parameters parameters = Parameters(R"({
        "solver_type": "fallback_linear_solver",
        "solvers"    : [
            {
                "solver_type": "amgcl"
            },
            {
                "solver_type": "amgcl"
            }
        ],
        "reset_solver_each_try": false
    })");
    TrilinosFallbackLinearSolverType simple_fallback_solver(parameters);

    // Solve the system
    bool solved = simple_fallback_solver.Solve(A, x, b);

    // Check that the system was solved
    KRATOS_EXPECT_TRUE(solved);

    // Check index is 0 (solve 0 succeeded)
    KRATOS_EXPECT_EQ(simple_fallback_solver.GetCurrentSolverIndex(), 0);
}

} // namespace Kratos::Testing