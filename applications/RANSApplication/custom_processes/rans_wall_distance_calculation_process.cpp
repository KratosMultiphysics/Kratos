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
#include "factories/linear_solver_factory.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

// Application includes

// Include base h
#include "rans_wall_distance_calculation_process.h"

namespace Kratos
{
template <>
void RansWallDistanceCalculationProcess<
    UblasSpace<double, CompressedMatrix, Vector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>>::CreateLinearSolver()
{
    mpLinearSolver = LinearSolverFactory<SparseSpaceType, DenseSpaceType>().Create(
        mrParameters["linear_solver_settings"]);
}

template <>
void RansWallDistanceCalculationProcess<
    UblasSpace<double, CompressedMatrix, Vector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>>::CreateBuilderAndSolver()
{
    mpBuilderAndSolver =
        Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>>(
            mpLinearSolver);
}

template <>
std::string RansWallDistanceCalculationProcess<
    UblasSpace<double, CompressedMatrix, Vector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>>::Info() const
{
    return std::string("RansWallDistanceCalculationProcess");
}

} // namespace Kratos.
