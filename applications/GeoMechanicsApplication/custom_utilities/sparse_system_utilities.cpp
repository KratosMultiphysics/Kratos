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
