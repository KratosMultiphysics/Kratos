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

#if !defined(KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED)
#define KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <functional>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/kratos_parameters.h"
#include "solving_strategies/schemes/scheme.h"
#include "response_functions/adjoint_response_function.h"
#include "utilities/variable_utils.h"
#include "utilities/indirect_scalar.h"
#include "utilities/adjoint_extensions.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/schemes/residual_based_adjoint_bossak_scheme.h"

// Application includes
#include "custom_processes/element_refinement_process.h"
#include "custom_utilities/fluid_adjoint_slip_utilities.h"

namespace Kratos
{
///@name Kratos Classes
///@{

template <class TSparseSpace, class TDenseSpace>
class VelocityBossakAdjointScheme : public ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VelocityBossakAdjointScheme);

    using BaseType = ResidualBasedAdjointBossakScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using DofsArrayType = typename BaseType::DofsArrayType;

    using BossakConstants = typename BaseType::BossakConstants;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using EquationIdVectorType = std::vector<IndexType>;

    using AssembledVectorType = std::unordered_map<IndexType, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    VelocityBossakAdjointScheme(
        Parameters Settings,
        AdjointResponseFunction::Pointer pResponseFunction,
        const IndexType Dimension,
        const IndexType BlockSize,
        ElementRefinementProcess::Pointer pElementRefinementProcess = nullptr)
        : BaseType(Settings, pResponseFunction),
          mAdjointSlipUtilities(Dimension, BlockSize),
          mpElementRefinementProcess(pElementRefinementProcess)
    {
        KRATOS_INFO(this->Info()) << this->Info() << " created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";
    }

    /// Destructor.
    ~VelocityBossakAdjointScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        // Allocate auxiliary memory.
        int num_threads = ParallelUtilities::GetNumThreads();
        mAuxiliaryMatrix.resize(num_threads);
        mRotatedMatrix.resize(num_threads);

        if (mpElementRefinementProcess) {
            mpElementRefinementProcess->ExecuteInitialize();

            const auto& dummy_element = *(rModelPart.ElementsBegin());
            const auto& dummy_geometry = dummy_element.GetGeometry();

            EquationIdVectorType equation_ids;
            dummy_element.EquationIdVector(equation_ids, rModelPart.GetProcessInfo());
            const IndexType number_of_dofs_per_node = equation_ids.size() / dummy_geometry.size();

            // reset the nodal interpolation error value holder
            const Vector zero_values(number_of_dofs_per_node, 0.0);
            const Vector max_values(number_of_dofs_per_node, std::numeric_limits<double>::max());
            block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
                rNode.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR, max_values);
                rNode.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY, zero_values);
            });
        }

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, A, Dx, b);

        if (mpElementRefinementProcess) {
            CalculateResponseFunctionInterpolationError(rModelPart);
        }

        KRATOS_CATCH("")
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "VelocityBossakAdjointScheme";
    }

    ///@}

protected:
    ///@name Protected Operations
    ///@{

    void CalculateGradientContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityGradientContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateGradientContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityGradientContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateFirstDerivativeContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityFirstDerivativeContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateFirstDerivativeContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityFirstDerivativeContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateSecondDerivativeContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySecondDerivativeContributions(
            rElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateSecondDerivativeContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySecondDerivativeContributions(
            rCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateTimeSchemeContributions(
        Element& rElement,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityTimeSchemeContributions(rElement, rAdjointTimeSchemeValues2,
                                               rAdjointTimeSchemeValues3,
                                               rCurrentProcessInfo);
    }

    void CalculateTimeSchemeContributions(
        Condition& rCondition,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityTimeSchemeContributions(rCondition, rAdjointTimeSchemeValues2,
                                               rAdjointTimeSchemeValues3,
                                               rCurrentProcessInfo);
    }

    void CalculateAuxiliaryVariableContributions(
        Element& rElement,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityAuxiliaryVariableContributions(
            rElement, rAdjointAuxiliaryValues, rCurrentProcessInfo);
    }

    void CalculateAuxiliaryVariableContributions(
        Condition& rCondition,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        AdjointResponseFunction& rAdjointResponseFunction,
        const BossakConstants& rBossakConstants,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntityAuxiliaryVariableContributions(
            rCondition, rAdjointAuxiliaryValues, rCurrentProcessInfo);
    }

    void CheckAndResizeThreadStorage(unsigned SystemSize) override
    {
        const int k = OpenMPUtils::ThisThread();

        BaseType::CheckAndResizeThreadStorage(SystemSize);

        if (mAuxiliaryMatrix[k].size1() != SystemSize ||
            mAuxiliaryMatrix[k].size2() != SystemSize) {
            mAuxiliaryMatrix[k].resize(SystemSize, SystemSize, false);
        }

        if (mRotatedMatrix[k].size1() != SystemSize || mRotatedMatrix[k].size2() != SystemSize) {
            mRotatedMatrix[k].resize(SystemSize, SystemSize, false);
        }
    }

    ///@}

private:
    ///@name Private Members
    ///@{

    std::vector<Matrix> mAuxiliaryMatrix;
    std::vector<Matrix> mRotatedMatrix;

    const FluidAdjointSlipUtilities mAdjointSlipUtilities;

    ElementRefinementProcess::Pointer mpElementRefinementProcess;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntityGradientContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();

        auto& aux_matrix = this->mAuxiliaryMatrix[k];

        rCurrentEntity.CalculateLeftHandSide(aux_matrix, rCurrentProcessInfo);
        this->mpResponseFunction->CalculateGradient(
            rCurrentEntity, aux_matrix, rRHS_Contribution, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            rLHS_Contribution, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rRHS_Contribution) = -rRHS_Contribution;
    }

    template<class TEntityType>
    void CalculateEntityFirstDerivativeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& rotated_matrix = this->mRotatedMatrix[k];
        auto& response_first_derivatives = this->mFirstDerivsResponseGradient[k];

        rCurrentEntity.CalculateFirstDerivativesLHS(aux_matrix, rCurrentProcessInfo);
        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentEntity, aux_matrix, response_first_derivatives, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            rotated_matrix, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rLHS_Contribution) += this->mBossak.C6 * rotated_matrix;
        noalias(rRHS_Contribution) -= this->mBossak.C6 * response_first_derivatives;
    }

    template<class TEntityType>
    void CalculateEntitySecondDerivativeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo)
    {
        const int k = OpenMPUtils::ThisThread();
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& rotated_matrix = this->mRotatedMatrix[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rCurrentProcessInfo);
        aux_matrix *= (1.0 - this->mBossak.Alpha);
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            rotated_matrix, aux_matrix, rCurrentEntity.GetGeometry());

        noalias(rLHS_Contribution) += this->mBossak.C7 * rotated_matrix;
        noalias(rRHS_Contribution) -= this->mBossak.C7 * response_second_derivatives;
    }

    template<class TEntityType>
    void CalculateEntityTimeSchemeContributions(
        TEntityType& rCurrentEntity,
        LocalSystemVectorType& rAdjointTimeSchemeValues2,
        LocalSystemVectorType& rAdjointTimeSchemeValues3,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        const int k = OpenMPUtils::ThisThread();
        auto& adjoint_values = this->mAdjointValuesVector[k];
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& entity_first_derivatives = this->mFirstDerivsLHS[k];
        auto& entity_second_derivatives = this->mSecondDerivsLHS[k];
        auto& response_first_derivatives = this->mFirstDerivsResponseGradient[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        const auto& r_const_entity_ref = rCurrentEntity;
        r_const_entity_ref.GetValuesVector(adjoint_values);
        this->CheckAndResizeThreadStorage(adjoint_values.size());

        /// starting to build residual for next time step calculations
        rCurrentEntity.CalculateFirstDerivativesLHS(aux_matrix, rProcessInfo);
        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentEntity, aux_matrix, response_first_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            entity_first_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rProcessInfo);
        aux_matrix *= (1.0 - this->mBossak.Alpha);
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            entity_second_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        if (rAdjointTimeSchemeValues2.size() != response_first_derivatives.size()) {
            rAdjointTimeSchemeValues2.resize(response_first_derivatives.size(), false);
        }

        noalias(rAdjointTimeSchemeValues2) =
            -response_first_derivatives - prod(entity_first_derivatives, adjoint_values);

        if (rAdjointTimeSchemeValues3.size() != response_second_derivatives.size()) {
            rAdjointTimeSchemeValues3.resize(response_second_derivatives.size(), false);
        }

        noalias(rAdjointTimeSchemeValues3) =
            -response_second_derivatives - prod(entity_second_derivatives, adjoint_values);

        KRATOS_CATCH("");
    }

    template <class TEntityType>
    void CalculateEntityAuxiliaryVariableContributions(
        TEntityType& rCurrentEntity,
        LocalSystemVectorType& rAdjointAuxiliaryValues,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        const int k = OpenMPUtils::ThisThread();
        auto& adjoint_values = this->mAdjointValuesVector[k];
        auto& aux_matrix = this->mAuxiliaryMatrix[k];
        auto& entity_second_derivatives = this->mSecondDerivsLHS[k];
        auto& response_second_derivatives = this->mSecondDerivsResponseGradient[k];

        const auto& r_const_entity_ref = rCurrentEntity;
        r_const_entity_ref.GetValuesVector(adjoint_values);
        this->CheckAndResizeThreadStorage(adjoint_values.size());

        rCurrentEntity.CalculateSecondDerivativesLHS(aux_matrix, rProcessInfo);
        aux_matrix *= this->mBossak.Alpha;
        this->mpResponseFunction->CalculateSecondDerivativesGradient(
            rCurrentEntity, aux_matrix, response_second_derivatives, rProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedNonSlipVariableDerivatives(
            entity_second_derivatives, aux_matrix, rCurrentEntity.GetGeometry());

        if (rAdjointAuxiliaryValues.size() != entity_second_derivatives.size1()) {
            rAdjointAuxiliaryValues.resize(entity_second_derivatives.size1(), false);
        }
        noalias(rAdjointAuxiliaryValues) =
            prod(entity_second_derivatives, adjoint_values) + response_second_derivatives;

        KRATOS_CATCH("");
    }

    void CalculateResponseFunctionInterpolationError(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const auto& dummy_element = *(rModelPart.ElementsBegin());
        const auto& dummy_geometry = dummy_element.GetGeometry();

        EquationIdVectorType equation_ids;
        dummy_element.EquationIdVector(equation_ids, rModelPart.GetProcessInfo());
        const IndexType number_of_dofs_per_node = equation_ids.size() / dummy_geometry.size();

        // acceptable error for an element
        const double acceptable_element_error = 1.0 / (8 * rModelPart.GetCommunicator().GlobalNumberOfElements());

        // thread local storage
        struct TLS
        {
            Matrix LHS;
            Vector RHS;
            Vector Values;
            Vector Residuals;
            Vector ErrorValuesList;
            EquationIdVectorType EquationIds;
            bool IsInitialized = false;
            std::unordered_map<IndexType, double> AssembledValues;
            std::unordered_map<IndexType, IndexType> ErrorIds;
            ModelPart* pRefinedModelPart;
        };

        const Vector zero_values(number_of_dofs_per_node, 0.0);
        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            noalias(rNode.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY)) = zero_values;
        });

        block_for_each(rModelPart.Elements(), TLS(), [&](Element& rElement, TLS& rTLS) {
            if (!rTLS.IsInitialized) {
                rTLS.pRefinedModelPart = &mpElementRefinementProcess->GetThreadLocalModelPart();

                for (auto& r_element : rTLS.pRefinedModelPart->Elements()) {
                    r_element.EquationIdVector(rTLS.EquationIds, rTLS.pRefinedModelPart->GetProcessInfo());
                    for (IndexType i = 0; i < rTLS.EquationIds.size(); ++i) {
                        const IndexType equation_id = rTLS.EquationIds[i];
                        rTLS.AssembledValues[equation_id];
                        rTLS.ErrorIds[equation_id] = (i % number_of_dofs_per_node);
                    }
                }

                rTLS.IsInitialized = true;
            }

            auto& r_process_info = rTLS.pRefinedModelPart->GetProcessInfo();

            for (auto& p : rTLS.AssembledValues) {
                p.second = 0.0;
            }

            mpElementRefinementProcess->InterpolateThreadLocalRefinedMeshFromCoarseElement(rElement);

            for (auto& r_element : rTLS.pRefinedModelPart->Elements()) {
                r_element.Initialize(r_process_info);
                r_element.InitializeSolutionStep(r_process_info);
            }

            for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
                r_condition.Initialize(r_process_info);
                r_condition.InitializeSolutionStep(r_process_info);
            }

            for (auto& r_element : rTLS.pRefinedModelPart->Elements()) {
                CalculateAndAssembleLocalEntityContributions(
                    r_element, rTLS.LHS, rTLS.Values, rTLS.Residuals,
                    rTLS.EquationIds, rTLS.AssembledValues, r_process_info);
            }

            for (auto& r_condition : rTLS.pRefinedModelPart->Conditions()) {
                if (r_condition.Is(ACTIVE)) {
                    CalculateAndAssembleLocalEntityContributions(
                        r_condition, rTLS.LHS, rTLS.Values, rTLS.Residuals,
                        rTLS.EquationIds, rTLS.AssembledValues, r_process_info);
                }
            }

            if (rTLS.ErrorValuesList.size() != number_of_dofs_per_node) {
                rTLS.ErrorValuesList.resize(number_of_dofs_per_node, false);
            }

            noalias(rTLS.ErrorValuesList) = ZeroVector(number_of_dofs_per_node);

            for (auto& p : rTLS.AssembledValues) {
                rTLS.ErrorValuesList[rTLS.ErrorIds.find(p.first)->second] += 1.0 / std::abs(p.second);
            }

            for (auto& r_node : rElement.GetGeometry()) {
                r_node.SetLock();
                r_node.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY) += rTLS.ErrorValuesList * (acceptable_element_error / rElement.GetGeometry().size());
                r_node.UnSetLock();
            }

        });

        rModelPart.GetCommunicator().AssembleNonHistoricalData(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY);

        // now take the minimum acceptable error out of all the time steps
        block_for_each(rModelPart.Nodes(), [](ModelPart::NodeType& rNode){
            Vector& r_values = rNode.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR);
            const Vector& r_current_time_step_values = rNode.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR_AUXILIARY);

            for (IndexType i = 0; i < r_values.size(); ++i) {
                double& r_value = r_values[i];
                r_value = std::min(r_value, r_current_time_step_values[i]);
            }
        });

        // no need to synchronize here for variable RESPONSE_FUNCTION_INTERPOLATION_ERROR

        KRATOS_CATCH("");
    }

    // This method can be called in parallel regions since this is used inside open mp loops of original model part
    // and iterated over refined elements of the original model part elements
    template<class TEntityType>
    void CalculateAndAssembleLocalEntityContributions(
        TEntityType& rEntity,
        Matrix& rLHS,
        Vector& rValues,
        Vector& rResiduals,
        EquationIdVectorType& rEquationIdsVector,
        AssembledVectorType& rAssembledValues,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        this->CalculateLHSContribution(rEntity, rLHS, rEquationIdsVector, rProcessInfo);
        rEntity.GetValuesVector(rValues);

        if (rResiduals.size() != rValues.size()) {
            rResiduals.resize(rValues.size());
        }

        noalias(rResiduals) = prod(rLHS, rValues);
        for (IndexType i = 0; i < rResiduals.size(); ++i) {
            rAssembledValues.find(rEquationIdsVector[i])->second += rResiduals[i];
        }

        KRATOS_CATCH("");
    }

    ///@}

}; /* Class VelocityBossakAdjointScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_VELOCITY_BOSSAK_ADJOINT_SCHEME_H_INCLUDED defined */
