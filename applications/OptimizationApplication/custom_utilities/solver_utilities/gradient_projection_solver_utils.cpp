//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <limits>

// Project includes
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "gradient_projection_solver_utils.h"

namespace Kratos {

void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(
    Vector& rSearchDirection,
    Vector& rSearchCorrection,
    LinearSolver<DenseSpace, DenseSpace>& rSolver,
    const Vector& rConstraintValues,
    const Vector& rObjectiveGradient,
    const std::vector<Vector>& rConstraintsGradients)
{
    KRATOS_TRY

    const IndexType dofs_size = rObjectiveGradient.size();
    const IndexType constraints_size = rConstraintsGradients.size();

    if (rSearchDirection.size() != dofs_size) {
        rSearchDirection.resize(dofs_size, false);
        rSearchCorrection.resize(dofs_size, false);
    }

    KRATOS_ERROR_IF_NOT(constraints_size == rConstraintValues.size())
        << "Number of constraints values provided does not match with number "
           "of constraint gradient vectors provided. [ number of constraint "
           "values = "
        << rConstraintValues.size()
        << ", number of constraints gradient vectors = " << constraints_size << " ].\n";

    for (const auto& r_vector : rConstraintsGradients) {
        KRATOS_ERROR_IF_NOT(r_vector.size() == dofs_size)
            << "Objective gradient and constraint gradients degree of freedom "
               "does not match. [ "
            << "rObjectiveGradient.size() = " << dofs_size
            << ", rConstraintGradient.size() = " << r_vector.size() << " ].\n";
    }

    // build the NTN matrix where N matrix has rows with the derivatives,
    // columns with the constraints
    Matrix NTN(constraints_size, constraints_size);
    // since it is symmetric
    for (IndexType i = 0; i < constraints_size; ++i) {
        for (IndexType j = i; j < constraints_size; ++j) {
            NTN(i, j) = IndexPartition<IndexType>(dofs_size).for_each<SumReduction<double>>(
                [&](const IndexType DofIndex) {
                    return rConstraintsGradients[i][DofIndex] *
                           rConstraintsGradients[j][DofIndex];
                });
            NTN(j, i) = NTN(i, j);
        }
    }

    // Now computing the inverse of NTN
    Matrix NTN_inv(NTN.size1(), NTN.size2());
    Matrix I = IdentityMatrix(NTN.size2());
    rSolver.Solve(NTN, NTN_inv, I); // solve with identity to get the inverse

    // calculates (N^T).nsbla_f
    Vector NT_nabla_f(constraints_size);
    for (IndexType i = 0; i < constraints_size; ++i) {
        NT_nabla_f[i] = IndexPartition<IndexType>(dofs_size).for_each<SumReduction<double>>(
            [&](const IndexType DofIndex) {
                return rConstraintsGradients[i][DofIndex] * rObjectiveGradient[DofIndex];
            });
    }

    // NTN_inv x ((N^T).nsbla_f)
    const Vector& NTN_inv_NT_nabla_f = prod(NTN_inv, NT_nabla_f);
    const Vector& NTN_inv_constraint_values = prod(NTN_inv, rConstraintValues);

    IndexPartition<IndexType>(dofs_size).for_each([&](const IndexType Index) {
        double search_value = -rObjectiveGradient[Index];
        double correction_value = 0.0;

        for (IndexType i = 0; i < constraints_size; ++ i) {
            search_value += rConstraintsGradients[i][Index] * NTN_inv_NT_nabla_f[i];
            correction_value -= rConstraintsGradients[i][Index] * NTN_inv_constraint_values[i];
        }

        rSearchDirection[Index] = search_value;
        rSearchCorrection[Index] = correction_value;
    });

    KRATOS_CATCH("");
}

void GradientProjectionSolverUtils::CalculateControlUpdate(
    Vector& rControlUpdate,
    Vector& rSearchDirection,
    Vector& rSearchCorrection,
    const double StepSize,
    const double MaxCorrectionShare)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rSearchDirection.size() == rSearchCorrection.size())
        << "Search direction and search correction vector size mismatch. [ Search direction vector size = "
        << rSearchDirection.size() << ", search correction vector size = " << rSearchCorrection.size()
        << " ].\n";

    if (rControlUpdate.size() != rSearchCorrection.size()) {
        rControlUpdate.resize(rSearchCorrection.size(), false);
    }

    const double search_direction_norm_inf = OptimizationUtils::NormInf(rSearchDirection);
    const double search_correction_norm_inf = OptimizationUtils::NormInf(rSearchCorrection);

    if (std::abs(search_correction_norm_inf) > std::numeric_limits<double>::epsilon()) {
        if (search_correction_norm_inf <= MaxCorrectionShare * StepSize) {
            OptimizationUtils::MultiplyVector(rSearchDirection, rSearchDirection, (StepSize - search_correction_norm_inf) / search_direction_norm_inf);
        } else {
            KRATOS_WARNING("GradientProjectionSolverUtils")
                << "Correction is scaled down from " << search_correction_norm_inf
                << " to " << MaxCorrectionShare * StepSize << ".\n";

            OptimizationUtils::MultiplyVector(rSearchCorrection, rSearchCorrection, MaxCorrectionShare * StepSize / search_correction_norm_inf);
            OptimizationUtils::MultiplyVector(rSearchDirection, rSearchDirection, (1.0 - MaxCorrectionShare) * StepSize / search_direction_norm_inf);
        }
    } else {
        OptimizationUtils::MultiplyVector(rSearchDirection, rSearchDirection, StepSize / search_direction_norm_inf);
    }

    OptimizationUtils::AddVectors(rControlUpdate, rSearchDirection, rSearchCorrection);

    KRATOS_CATCH("");
}

} // namespace Kratos