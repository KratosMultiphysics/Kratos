//
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
#include "testing/testing.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
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
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "Dummy linear solver";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Dummy linear solver";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

}; // Class DummyLinearSolver

// Define the types of spaces
using SpaceType = UblasSpace<double, boost::numeric::ublas::compressed_matrix<double>, boost::numeric::ublas::vector<double>>;
using LocalSpaceType = UblasSpace<double, DenseMatrix<double>, DenseVector<double>>;

// Define the vector and matrix types
using SparseMatrixType = typename SpaceType::MatrixType;
using VectorType = typename SpaceType::VectorType;
using DenseMatrixType = typename LocalSpaceType::MatrixType;
using DenseVectorType = typename LocalSpaceType::VectorType;

// Define the types of solvers
using LinearSolverType = LinearSolver<SpaceType, LocalSpaceType>;
using DummyLinearSolverType = DummyLinearSolver<SpaceType, LocalSpaceType>;
using SkylineLUFactorizationSolverType = SkylineLUFactorizationSolver<SpaceType, LocalSpaceType>;
using FallbackLinearSolverType = FallbackLinearSolver<SpaceType, LocalSpaceType>;
using LinearSolverFactoryType = LinearSolverFactory<SpaceType, LocalSpaceType>;

KRATOS_TEST_CASE_IN_SUITE(FallbackLinearSolverConstructorSolvers, KratosCoreFastSuite)
{
    // Create the solvers
    auto p_solver1 = Kratos::make_shared<DummyLinearSolverType>();
    auto p_solver2 = Kratos::make_shared<SkylineLUFactorizationSolverType>();

    // Create the matrix and vectors
    const std::size_t size = 3;
    SparseMatrixType A(size, size);
    // Push back 1.0 in the diagonal
    for (std::size_t i = 0; i < size; ++i) {
        A.push_back(i, i, 1.0);
    }
    A.set_filled(size + 1, size);
    VectorType b(size);
    VectorType x(size);

    // Create a simple fallback solver
    std::vector<LinearSolverType::Pointer> solvers = {p_solver1, p_solver2};
    FallbackLinearSolverType simple_fallback_solver(solvers);

    // Solve the system
    const auto solved = simple_fallback_solver.Solve(A, x, b);

    // Check that the system was solved
    KRATOS_EXPECT_TRUE(solved);

    // Check index is 1 (solve 0 failed)
    KRATOS_EXPECT_EQ(simple_fallback_solver.GetCurrentSolverIndex(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(FallbackLinearSolverConstructorParameters, KratosCoreFastSuite)
{
    // Create the matrix and vectors
    const std::size_t size = 3;
    SparseMatrixType A(size, size);
    // Push back 1.0 in the diagonal
    for (std::size_t i = 0; i < size; ++i) {
        A.push_back(i, i, 1.0);
    }
    A.set_filled(size + 1, size);
    VectorType b(size);
    VectorType x(size);

    // Create a simple fallback solver
    Parameters parameters = Parameters(R"({
        "solver_type": "fallback_linear_solver",
        "solvers"    : [
            {
                "solver_type": "skyline_lu_factorization"
            },
            {
                "solver_type": "cg"
            }
        ],
        "reset_solver_each_try": false
    })");
    FallbackLinearSolverType simple_fallback_solver(parameters);

    // Solve the system
    bool solved = simple_fallback_solver.Solve(A, x, b);

    // Check that the system was solved
    KRATOS_EXPECT_TRUE(solved);

    // Check index is 0 (solve 0 succeeded)
    KRATOS_EXPECT_EQ(simple_fallback_solver.GetCurrentSolverIndex(), 0);
}

} // namespace Kratos::Testing