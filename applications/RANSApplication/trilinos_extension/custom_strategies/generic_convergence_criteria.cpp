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
#include <unordered_map>
#include <cmath>

/* Project includes */
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/openmp_utils.h"

#include "Epetra_FEVector.h"
#include "trilinos_space.h"

// Include base h
#include "custom_strategies/generic_convergence_criteria.h"

namespace Kratos
{
template <>
void GenericConvergenceCriteria<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>, UblasSpace<double, Matrix, Vector>>::CalculateConvergenceCheckNorms(
    double& rSolutionNorm,
    double& rIncreaseNorm,
    double& rDofSize,
    ModelPart& rModelPart,
    DofsArrayType& rDofSet,
    const TSystemMatrixType& A,
    const TSystemVectorType& Dx,
    const TSystemVectorType& b)
{
    KRATOS_TRY

    int NumDofs = rDofSet.size();

    double solution_norm{0.0}, increase_norm{0.0};
    int dof_num{0};

    // Set a partition for OpenMP
    PartitionVector DofPartition;
    const int NumThreads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::DivideInPartitions(NumDofs, NumThreads, DofPartition);

    const Communicator& r_communicator = rModelPart.GetCommunicator();
    const int my_pid = r_communicator.MyPID();

    // Loop over Dofs
#pragma omp parallel reduction(+ : solution_norm, increase_norm, dof_num)
    {
        const int k = OpenMPUtils::ThisThread();
        typename DofsArrayType::iterator DofBegin = rDofSet.begin() + DofPartition[k];
        typename DofsArrayType::iterator DofEnd = rDofSet.begin() + DofPartition[k + 1];

        std::size_t DofId;
        TDataType DofValue;
        TDataType DofIncr;

        for (typename DofsArrayType::iterator itDof = DofBegin; itDof != DofEnd; ++itDof)
        {
            if (itDof->IsFree() && itDof->GetSolutionStepValue(PARTITION_INDEX) == my_pid)
            {
                DofId = itDof->EquationId();
                DofValue = itDof->GetSolutionStepValue(0);
                DofIncr = SparseSpaceType::GetValue(Dx, DofId);

                solution_norm += DofValue * DofValue;
                increase_norm += DofIncr * DofIncr;
                dof_num += 1;
            }
        }
    }

    std::vector<double> residual_norms = {increase_norm, solution_norm,
                                          static_cast<double>(dof_num)};
    const std::vector<double>& total_residual_norms =
        r_communicator.GetDataCommunicator().SumAll(residual_norms);

    rSolutionNorm = std::sqrt(total_residual_norms[1]);
    rIncreaseNorm = std::sqrt(total_residual_norms[0]);
    rDofSize = total_residual_norms[2];

    KRATOS_CATCH("");
}

// template instantiations

template class GenericConvergenceCriteria<TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>,
                                          UblasSpace<double, Matrix, Vector>>;

} // namespace Kratos
