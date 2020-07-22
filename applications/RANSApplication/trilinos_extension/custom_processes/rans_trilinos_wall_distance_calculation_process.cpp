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
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

// Application includes

// Include base h
#include "rans_trilinos_wall_distance_calculation_process.h"

namespace Kratos
{
template <>
void TrilinosRansWallDistanceCalculationProcess<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>>::CreateLinearSolver()
{
    mpLinearSolver = LinearSolverFactory<SparseSpaceType, DenseSpaceType>().Create(
        mrParameters["linear_solver_settings"]);
}

template <>
void TrilinosRansWallDistanceCalculationProcess<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>>::CreateBuilderAndSolver()
{
    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];
    const int row_size_guess = (domain_size == 2 ? 15 : 40);

    mpBuilderAndSolver =
        Kratos::make_shared<TrilinosBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>>(
            mMPIComm, row_size_guess, mpLinearSolver);
}

template <>
std::string TrilinosRansWallDistanceCalculationProcess<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>>::Info() const
{
    return std::string("TrilinosRansWallDistanceCalculationProcess");
}

} // namespace Kratos.
