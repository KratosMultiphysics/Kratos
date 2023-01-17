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

template<>
void GradientProjectionSolverUtils::AddValue(
    double& rOutput,
    const IndexType VariableComponentIndex,
    const IndexType DofStartingIndex,
    const Vector& rValues)
{
    rOutput += rValues[DofStartingIndex];
}


template<>
void GradientProjectionSolverUtils::AddValue(
    array_1d<double, 3>& rOutput,
    const IndexType VariableComponentIndex,
    const IndexType DofStartingIndex,
    const Vector& rValues)
{
    rOutput[VariableComponentIndex] += rValues[DofStartingIndex + VariableComponentIndex];
}

template<>
double GradientProjectionSolverUtils::CalculateValueNormSquare(
    const double& rValue)
{
    return std::pow(rValue, 2);
}

template<>
double GradientProjectionSolverUtils::CalculateValueNormSquare(
    const array_1d<double, 3>& rValue)
{
    return std::pow(rValue[0], 2) + std::pow(rValue[1], 2) + std::pow(rValue[2], 2);
}

template<class TContainerType, class TDataType>
void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(
    TContainerType& rContainer,
    const IndexType DomainSize,
    LinearSolver<DenseSpace, DenseSpace>& rSolver,
    const Variable<TDataType>& rSearchDirectionVariable,
    const Variable<TDataType>& rSearchDirectionCorrectionVariable,
    const Vector& rConstraintValues,
    const Vector& rObjectiveGradient,
    const std::vector<Vector>& rConstraintsGradients)
{
    KRATOS_TRY

    const IndexType dofs_size = rObjectiveGradient.size();
    const IndexType constraints_size = rConstraintsGradients.size();
    const IndexType local_size = OptimizationUtils::GetLocalSize<TDataType>(DomainSize);

    KRATOS_ERROR_IF_NOT(rContainer.size() * local_size == dofs_size)
        << "Number of dofs in the contaienr does not match with number of dofs "
           "given in the objective gradient vector. [ Dofs in contaienr = "
        << rContainer.size() * local_size
        << ", number dofs in the objective gradient vector = " << dofs_size << " ].\n";

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
    const Vector& negative_NTN_inv_constraint_values = -prod(NTN_inv, rConstraintValues);
    const Vector& negative_objective_gradients = -rObjectiveGradient;

    IndexPartition<IndexType>(rContainer.size()).for_each([&](const IndexType rEntityIndex) {
        TDataType search_direction = rSearchDirectionVariable.Zero();
        TDataType correction = rSearchDirectionVariable.Zero();

        for (IndexType i = 0; i < local_size; ++i) {
            const IndexType dof_starting_index = rEntityIndex * local_size;

            AddValue(search_direction, i, dof_starting_index, negative_objective_gradients);
            for (IndexType j = 0; j < constraints_size; ++j) {
                AddValue(search_direction, i, dof_starting_index, rConstraintsGradients[j] * NTN_inv_NT_nabla_f[j]);
                AddValue(correction, i, dof_starting_index, rConstraintsGradients[j] * negative_NTN_inv_constraint_values[j]);
            }
        }

        auto& r_entity = *(rContainer.begin() + rEntityIndex);
        r_entity.SetValue(rSearchDirectionVariable, search_direction);
        r_entity.SetValue(rSearchDirectionCorrectionVariable, correction);
    });



    KRATOS_CATCH("");
}

template<class TContainerType, class TDataType>
void GradientProjectionSolverUtils::CalculateControlChange(
    TContainerType& rContainer,
    const DataCommunicator& rDataCommunicator,
    const Variable<TDataType>& rSearchDirectionVariable,
    const Variable<TDataType>& rSearchDirectionCorrectionVariable,
    const Variable<TDataType>& rControlChangeVariable,
    const double StepSize,
    const double MaxCorrectionShare)
{
    KRATOS_TRY

    double local_search_direction_sum_square, local_search_correction_sum_square;
    std::tie(local_search_direction_sum_square, local_search_correction_sum_square) = block_for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(rContainer, [&](const auto& rEntity) {
        return std::make_tuple(
            CalculateValueNormSquare(rEntity.GetValue(rSearchDirectionVariable)),
            CalculateValueNormSquare(rEntity.GetValue(rSearchDirectionCorrectionVariable)));
    });
    const auto& norms = rDataCommunicator.SumAll(array_1d<double, 3>({local_search_direction_sum_square, local_search_correction_sum_square, 0.0}));
    const double search_direction_norm = std::sqrt(norms[0]);
    const double correction_norm = std::sqrt(norms[1]);

    if (search_direction_norm > 0.0) {
        double search_direction_multiplyer;
        double correction_multiplyer;
        if (correction_norm > 0.0) {
            if (correction_norm <= MaxCorrectionShare * StepSize) {
                const double delta = StepSize - correction_norm;
                search_direction_multiplyer = delta / search_direction_norm;
            } else {
                KRATOS_WARNING("GradientProjectionSolverUtils")
                    << "Correction is scaled down from " << correction_norm
                    << " to " << MaxCorrectionShare * StepSize << ".\n";
                correction_multiplyer = MaxCorrectionShare * StepSize / correction_norm;
                search_direction_multiplyer = (1.0 - MaxCorrectionShare) * StepSize / search_direction_norm;
            }
        } else {
            search_direction_multiplyer = StepSize / search_direction_norm;
        }

        block_for_each(rContainer, [&](auto& rEntity) {
            auto& r_search_direction = rEntity.GetValue(rSearchDirectionVariable);
            r_search_direction *= search_direction_multiplyer;
            auto& r_search_direction_correction = rEntity.GetValue(rSearchDirectionCorrectionVariable);
            r_search_direction_correction *= correction_multiplyer;
            rEntity.SetValue(rControlChangeVariable, r_search_direction + r_search_direction_correction);
        });
    } else {
        KRATOS_WARNING("GradientProjectionSolverUtils") << "Norm of " << rSearchDirectionVariable.Name() << " is zero. Hence skipping control change calculation.\n";
        VariableUtils().SetNonHistoricalVariableToZero(rControlChangeVariable, rContainer);
    }

    KRATOS_CATCH("");
}

// template instantiations
template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::NodesContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<double>&, const Variable<double>&, const Vector&, const Vector&, const std::vector<Vector>&);
template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::ConditionsContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<double>&, const Variable<double>&, const Vector&, const Vector&, const std::vector<Vector>&);
template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::ElementsContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<double>&, const Variable<double>&, const Vector&, const Vector&, const std::vector<Vector>&);

template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::NodesContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Vector&, const Vector&, const std::vector<Vector>&);
template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::ConditionsContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Vector&, const Vector&, const std::vector<Vector>&);
template void GradientProjectionSolverUtils::CalculateProjectedSearchDirectionAndCorrection(ModelPart::ElementsContainerType&, const IndexType, LinearSolver<DenseSpace, DenseSpace>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Vector&, const Vector&, const std::vector<Vector>&);

template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::NodesContainerType&, const DataCommunicator&, const Variable<double>&, const Variable<double>&, const Variable<double>&, const double, const double);
template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::ConditionsContainerType&, const DataCommunicator&, const Variable<double>&, const Variable<double>&, const Variable<double>&,  const double, const double);
template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::ElementsContainerType&, const DataCommunicator&, const Variable<double>&, const Variable<double>&, const Variable<double>&, const double, const double);

template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::NodesContainerType&, const DataCommunicator&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const double, const double);
template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::ConditionsContainerType&, const DataCommunicator&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const double, const double);
template void GradientProjectionSolverUtils::CalculateControlChange(ModelPart::ElementsContainerType&, const DataCommunicator&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const Variable<array_1d<double, 3>>&, const double, const double);

} // namespace Kratos