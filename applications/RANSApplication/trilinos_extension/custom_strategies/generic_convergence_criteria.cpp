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
#include <cmath>
#include <unordered_map>

/* Project includes */
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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
    const TSystemMatrixType& rA,
    const TSystemVectorType& rDx,
    const TSystemVectorType& rb)
{
    KRATOS_TRY

    const Communicator& r_communicator = rModelPart.GetCommunicator();
    const int my_pid = r_communicator.MyPID();

    double solution_norm, increase_norm;
    int dof_num;

    std::tie(solution_norm, increase_norm, dof_num) =
        BlockPartition<DofsArrayType>(rDofSet)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<int>>>(
                [&](const DofType& rDof) -> std::tuple<double, double, int> {
                    if (rDof.IsFree() && rDof.GetSolutionStepValue(PARTITION_INDEX) == my_pid) {
                        const auto dof_value = rDof.GetSolutionStepValue(0);
                        const auto dof_increment =
                            SparseSpaceType::GetValue(rDx, rDof.EquationId());

                        return std::make_tuple<double, double, int>(
                            dof_value * dof_value, dof_increment * dof_increment, 1);
                    } else {
                        return std::make_tuple<double, double, int>(0.0, 0.0, 0);
                    }
                });

    auto residual_norms = {increase_norm, solution_norm, static_cast<double>(dof_num)};
    const auto& total_residual_norms =
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
