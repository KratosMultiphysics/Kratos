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

/* Project includes */
#include "includes/define.h"
#include "spaces/ublas_space.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/openmp_utils.h"

// Include base h
#include "generic_convergence_criteria.h"

namespace Kratos
{
template <>
void GenericConvergenceCriteria<UblasSpace<double, CompressedMatrix, Vector>, UblasSpace<double, Matrix, Vector>>::CalculateConvergenceCheckNorms(
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

    double solution_norm, increase_norm;
    int dof_num;

    std::tie(solution_norm, increase_norm, dof_num) =
        BlockPartition<DofsArrayType>(rDofSet)
            .for_each<CombinedReduction<SumReduction<double>, SumReduction<double>, SumReduction<int>>>(
                [&](const DofType& rDof) -> std::tuple<double, double, int> {
                    if (rDof.IsFree()) {
                        const auto dof_value = rDof.GetSolutionStepValue(0);
                        const auto dof_increment = rDx[rDof.EquationId()];

                        return std::make_tuple<double, double, int>(
                            dof_value * dof_value, dof_increment * dof_increment, 1);
                    } else {
                        return std::make_tuple<double, double, int>(0.0, 0.0, 0);
                    }
                });

    rSolutionNorm = std::sqrt(solution_norm);
    rIncreaseNorm = std::sqrt(increase_norm);
    rDofSize = static_cast<double>(dof_num);

    KRATOS_CATCH("");
}

// template instantiations

template class GenericConvergenceCriteria<UblasSpace<double, CompressedMatrix, Vector>,
                                          UblasSpace<double, Matrix, Vector>>;

} // namespace Kratos
