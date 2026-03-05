// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include "includes/variables.h"
#include "sparse_system_utilities.h"

using namespace Kratos;

namespace Kratos::Geo
{

void SparseSystemUtilities::GetTotalSolutionStepValueVector(SystemVectorType& rTotalSolutionStepValues,
                                                            const DofsArrayType& rDofSet,
                                                            const ModelPart&     rModelPart,
                                                            const IndexType      BufferIndex)
{
    rTotalSolutionStepValues = SystemVectorType(rDofSet.size(), 0.0);

    // Firstly initialize the rTotalSolutionStepValues with all DOFs in the DofSet, also the Dofs which do not have a "TOTAL" counterpart
    block_for_each(rDofSet, [&rTotalSolutionStepValues](Dof<double>& r_dof) {
        rTotalSolutionStepValues[r_dof.EquationId()] = r_dof.GetSolutionStepValue(0);
    });

    // Mandatory variables (always present)
    const std::vector<std::pair<const Variable<double>*, const Variable<double>*>> mandatory_variables = {
        {&DISPLACEMENT_X, &TOTAL_DISPLACEMENT_X}, {&DISPLACEMENT_Y, &TOTAL_DISPLACEMENT_Y}};

    // Optional variables (only if DOF exists)
    const std::vector<std::pair<const Variable<double>*, const Variable<double>*>> optional_variables = {
        {&DISPLACEMENT_Z, &TOTAL_DISPLACEMENT_Z},
        {&ROTATION_X, &TOTAL_ROTATION_X},
        {&ROTATION_Y, &TOTAL_ROTATION_Y},
        {&ROTATION_Z, &TOTAL_ROTATION_Z}};

    // Overwrite the values in rTotalSolutionStepValues with the TOTAL values from the nodes
    block_for_each(rModelPart.Nodes(), [&rTotalSolutionStepValues, BufferIndex,
                                        &mandatory_variables, &optional_variables](Node& rNode) {
        if (rNode.IsActive()) {
            for (const auto& var_pair : mandatory_variables) {
                const auto equation_id = rNode.GetDof(*var_pair.first).EquationId();
                rTotalSolutionStepValues[equation_id] =
                    rNode.FastGetSolutionStepValue(*var_pair.second, BufferIndex);
            }

            for (const auto& var_pair : optional_variables) {
                if (rNode.HasDofFor(*var_pair.first)) {
                    const auto equation_id = rNode.GetDof(*var_pair.first).EquationId();
                    rTotalSolutionStepValues[equation_id] =
                        rNode.FastGetSolutionStepValue(*var_pair.second, BufferIndex);
                }
            }
        }
    });
}

void SparseSystemUtilities::GetUFirstAndSecondDerivativeVector(SystemVectorType& rFirstDerivativeVector,
                                                               SystemVectorType& rSecondDerivativeVector,
                                                               const DofsArrayType& rDofSet,
                                                               const ModelPart&     rModelPart,
                                                               const IndexType      BufferIndex)
{
    KRATOS_TRY
    rFirstDerivativeVector  = SystemVectorType(rDofSet.size(), 0.0);
    rSecondDerivativeVector = SystemVectorType(rDofSet.size(), 0.0);

    block_for_each(rModelPart.Nodes(),
                   [&rFirstDerivativeVector, &rSecondDerivativeVector, BufferIndex](const Node& rNode) {
        if (rNode.IsActive()) {
            GetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector,
                                      rSecondDerivativeVector, BufferIndex);
            GetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector,
                                      rSecondDerivativeVector, BufferIndex);

            const std::vector<const Variable<double>*> optional_variables = {
                &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

            for (const auto p_variable : optional_variables) {
                GetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector,
                                                  rSecondDerivativeVector, BufferIndex);
            }
        }
    });
    KRATOS_CATCH("")
}

void SparseSystemUtilities::SetUFirstAndSecondDerivativeVector(const SystemVectorType& rFirstDerivativeVector,
                                                               const SystemVectorType& rSecondDerivativeVector,
                                                               ModelPart& rModelPart)
{
    KRATOS_TRY
    block_for_each(rModelPart.Nodes(), [&rFirstDerivativeVector, &rSecondDerivativeVector](Node& rNode) {
        if (rNode.IsActive()) {
            SetDerivativesForVariable(DISPLACEMENT_X, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
            SetDerivativesForVariable(DISPLACEMENT_Y, rNode, rFirstDerivativeVector, rSecondDerivativeVector);

            const std::vector<const Variable<double>*> optional_variables = {
                &ROTATION_X, &ROTATION_Y, &ROTATION_Z, &DISPLACEMENT_Z};

            for (const auto p_variable : optional_variables) {
                SetDerivativesForOptionalVariable(*p_variable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
            }
        }
    });
    KRATOS_CATCH("")
}

void SparseSystemUtilities::ApplyDirichletConditionsSecondaryMatrix(const DofsArrayType& rDofSet,
                                                                    SystemMatrixType& rSecondaryMatrix)
{
    const std::size_t system_size = rSecondaryMatrix.size1();
    Vector            scaling_factors(system_size, 1.0);

    const auto it_dof_iterator_begin = rDofSet.begin();

    // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
    IndexPartition<std::size_t>(rDofSet.size()).for_each([&](std::size_t Index) {
        auto it_dof_iterator = it_dof_iterator_begin + Index;
        if (it_dof_iterator->IsFixed()) {
            scaling_factors[Index] = 0.0;
        }
    });

    double*      Avalues      = rSecondaryMatrix.value_data().begin();
    std::size_t* Arow_indices = rSecondaryMatrix.index1_data().begin();
    std::size_t* Acol_indices = rSecondaryMatrix.index2_data().begin();

    IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index) {
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end   = Arow_indices[Index + 1];
        const double      k_factor  = scaling_factors[Index];
        if (k_factor == 0.0) {
            // Zero out the whole row, except the diagonal
            for (std::size_t j = col_begin; j < col_end; ++j)
                if (Acol_indices[j] != Index) Avalues[j] = 0.0;
        } else {
            // Zero out the column which is associated with the zero'ed row
            for (std::size_t j = col_begin; j < col_end; ++j)
                if (scaling_factors[Acol_indices[j]] == 0) Avalues[j] = 0.0;
        }
    });
}

bool SparseSystemUtilities::HasNonZeroDiagonalEntryOnCurrentRow(const std::size_t RowIndex,
                                                                const unbounded_array<std::size_t>& rCsrIndices1,
                                                                const unbounded_array<std::size_t>& rCsrIndices2,
                                                                const unbounded_array<double>& rCsrValues)
{
    // Check if there is a non-zero entry on the diagonal for the given row index
    for (std::size_t k = rCsrIndices1[RowIndex]; k < rCsrIndices1[RowIndex + 1]; ++k) {
        const std::size_t col = rCsrIndices2[k];
        if (col == RowIndex) {
            return (rCsrValues[k] != 0.0);
        }
    }
    return false;
}

bool SparseSystemUtilities::MatricesHaveSameDiagonalSignature(const SystemMatrixType& rPrimaryMatrix,
                                                              const SystemMatrixType& rSecondaryMatrix,
                                                              const DofsArrayType& rDofSet)
{
    const std::size_t size1 = rPrimaryMatrix.size1();
    if (size1 != rSecondaryMatrix.size1() || rPrimaryMatrix.size2() != rSecondaryMatrix.size2()) {
        return false;
    }

    const auto& r_primary_index1 = rPrimaryMatrix.index1_data();
    const auto& r_primary_index2 = rPrimaryMatrix.index2_data();
    const auto& r_primary_values = rPrimaryMatrix.value_data();

    const auto& r_secondary_index1 = rSecondaryMatrix.index1_data();
    const auto& r_secondary_index2 = rSecondaryMatrix.index2_data();
    const auto& r_secondary_values = rSecondaryMatrix.value_data();

    auto dof_iterator_begin = rDofSet.begin();

    for (std::size_t row_index = 0; row_index < size1; ++row_index) {
        // Find diagonal for Matrix rPrimaryMatrix
        bool has_primary_diag = SparseSystemUtilities::HasNonZeroDiagonalEntryOnCurrentRow(
            row_index, r_primary_index1, r_primary_index2, r_primary_values);

        // Find diagonal for Matrix rSecondaryMatrix
        bool has_secondary_diag = SparseSystemUtilities::HasNonZeroDiagonalEntryOnCurrentRow(
            row_index, r_secondary_index1, r_secondary_index2, r_secondary_values);

        if (has_primary_diag != has_secondary_diag && !(dof_iterator_begin + row_index)->IsFixed())
            return false;
    }

    return true;
}

void SparseSystemUtilities::GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                              const Node&             rNode,
                                                              SystemVectorType& rFirstDerivativeVector,
                                                              SystemVectorType& rSecondDerivativeVector,
                                                              const IndexType BufferIndex)
{
    KRATOS_TRY
    if (rNode.HasDofFor(rVariable)) {
        GetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector, BufferIndex);
    }
    KRATOS_CATCH("")
}

void SparseSystemUtilities::SetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                              Node&                   rNode,
                                                              const SystemVectorType& rFirstDerivativeVector,
                                                              const SystemVectorType& rSecondDerivativeVector)
{
    KRATOS_TRY
    if (rNode.HasDofFor(rVariable)) {
        SetDerivativesForVariable(rVariable, rNode, rFirstDerivativeVector, rSecondDerivativeVector);
    }
    KRATOS_CATCH("")
}

void SparseSystemUtilities::GetDerivativesForVariable(const Variable<double>& rVariable,
                                                      const Node&             rNode,
                                                      SystemVectorType& rFirstDerivativeVector,
                                                      SystemVectorType& rSecondDerivativeVector,
                                                      const IndexType   BufferIndex)
{
    KRATOS_TRY
    const auto& r_first_derivative  = rVariable.GetTimeDerivative();
    const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

    const auto equation_id = rNode.GetDof(rVariable).EquationId();
    rFirstDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_first_derivative, BufferIndex);
    rSecondDerivativeVector[equation_id] = rNode.FastGetSolutionStepValue(r_second_derivative, BufferIndex);

    KRATOS_CATCH("")
}

void SparseSystemUtilities::SetDerivativesForVariable(const Variable<double>& rVariable,
                                                      Node&                   rNode,
                                                      const SystemVectorType& rFirstDerivativeVector,
                                                      const SystemVectorType& rSecondDerivativeVector)
{
    KRATOS_TRY
    const auto& r_first_derivative  = rVariable.GetTimeDerivative();
    const auto& r_second_derivative = r_first_derivative.GetTimeDerivative();

    const auto equation_id = rNode.GetDof(rVariable).EquationId();

    if (!rNode.IsFixed(r_first_derivative)) {
        rNode.FastGetSolutionStepValue(r_first_derivative) = rFirstDerivativeVector[equation_id];
    }

    if (!rNode.IsFixed(r_second_derivative)) {
        rNode.FastGetSolutionStepValue(r_second_derivative) = rSecondDerivativeVector[equation_id];
    }

    KRATOS_CATCH("")
}

} // namespace Kratos::Geo
