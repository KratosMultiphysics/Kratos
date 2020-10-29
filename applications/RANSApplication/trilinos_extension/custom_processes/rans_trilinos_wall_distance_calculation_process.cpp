//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes
#include "Epetra_FEVector.h"
#include "Epetra_MpiComm.h"
#include "custom_factories/trilinos_linear_solver_factory.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"
#include "processes/variational_distance_calculation_process.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

// Application includes

// Include base h
#include "rans_trilinos_wall_distance_calculation_process.h"

namespace Kratos
{
std::string TrilinosRansWallDistanceCalculationProcess::Info() const
{
    return std::string("TrilinosRansWallDistanceCalculationProcess");
}

Process::Pointer TrilinosRansWallDistanceCalculationProcess::GetWallDistanceCalculationProcess(
    ModelPart& rModelPart,
    Parameters LinearSolverParameters,
    const int MaxIterations)
{
    KRATOS_TRY

    using SparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;

    auto p_linear_solver =
        LinearSolverFactory<SparseSpaceType, DenseSpaceType>().Create(LinearSolverParameters);

    const int domain_size = rModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const int row_size_guess = (domain_size == 2 ? 15 : 40);

    auto p_builder_and_solver =
        Kratos::make_shared<TrilinosBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>>(
            mMPIComm, row_size_guess, p_linear_solver);

    Process::Pointer p_process = nullptr;

    if (domain_size == 2) {
        p_process =
            Kratos::make_shared<VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType>>(
                rModelPart, p_linear_solver, p_builder_and_solver, MaxIterations);
    } else if (domain_size == 3) {
        p_process =
            Kratos::make_shared<VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType>>(
                rModelPart, p_linear_solver, p_builder_and_solver, MaxIterations);
    } else {
        KRATOS_ERROR << "Unknown domain size = " << domain_size;
    }

    return p_process;

    KRATOS_CATCH("");
}

} // namespace Kratos.
