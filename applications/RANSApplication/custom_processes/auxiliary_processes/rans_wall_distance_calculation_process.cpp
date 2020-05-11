//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// External includes

// Project includes
#include "processes/variational_distance_calculation_process.h"

// Application includes

// Include base h
#include "rans_wall_distance_calculation_process.h"

namespace Kratos
{
template <>
void RansWallDistanceCalculationProcess<
    UblasSpace<double, CompressedMatrix, Vector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>>::ExecuteVariationalDistanceCalculationProcess()
{
    KRATOS_TRY

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType =
        LinearSolver<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>;

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    if (domain_size == 2)
    {
        if (mpBuilderAndSolver)
        {
            VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
                r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else if (!r_model_part.IsDistributed())
        {
            VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
                r_model_part, mpLinearSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else
        {
            KRATOS_ERROR << "Distributed run requires to use "
                            "SetBuilderAndSolver method to provide "
                            "appropriate distributed builder and solver.\n";
        }
    }
    else if (domain_size == 3)
    {
        if (mpBuilderAndSolver)
        {
            VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
                r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else if (!r_model_part.IsDistributed())
        {
            VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
                r_model_part, mpLinearSolver, mMaxIterations);
            distance_calculation_process.Execute();
        }
        else
        {
            KRATOS_ERROR << "Distributed run requires to use "
                            "SetBuilderAndSolver method to provide "
                            "appropriate distributed builder and solver.\n";
        }
    }
    else
    {
        KRATOS_ERROR << "Unknown domain size = " << domain_size;
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
