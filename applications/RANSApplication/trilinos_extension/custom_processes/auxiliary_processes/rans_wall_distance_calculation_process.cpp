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
#include "Epetra_FEVector.h"
#include "processes/variational_distance_calculation_process.h"
#include "trilinos_space.h"

// Application includes

// Include base h
#include "custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.h"

namespace Kratos
{
template <>
void RansWallDistanceCalculationProcess<
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
    UblasSpace<double, Matrix, Vector>,
    LinearSolver<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>>::ExecuteVariationalDistanceCalculationProcess()
{
    KRATOS_TRY

    using SparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
    using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;

    ModelPart& r_model_part = mrModel.GetModelPart(mModelPartName);

    const int domain_size = r_model_part.GetProcessInfo()[DOMAIN_SIZE];

    KRATOS_ERROR_IF(!mpBuilderAndSolver)
        << "Distributed run requires to use "
           "SetBuilderAndSolver method to provide "
           "appropriate distributed builder and solver.\n";

    if (domain_size == 2)
    {
        VariationalDistanceCalculationProcess<2, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
            r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
        distance_calculation_process.Execute();
    }
    else if (domain_size == 3)
    {
        VariationalDistanceCalculationProcess<3, SparseSpaceType, DenseSpaceType, LinearSolverType> distance_calculation_process(
            r_model_part, mpLinearSolver, mpBuilderAndSolver, mMaxIterations);
        distance_calculation_process.Execute();
    }
    else
    {
        KRATOS_ERROR << "Unknown domain size = " << domain_size;
    }

    KRATOS_CATCH("");
}

} // namespace Kratos.
