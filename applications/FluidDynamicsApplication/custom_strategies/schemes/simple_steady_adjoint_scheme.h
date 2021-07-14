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

#if !defined(KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED)
#define KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED

// System includes
#include <string>
#include <unordered_map>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_processes/element_refinement_process.h"
#include "custom_utilities/fluid_adjoint_slip_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for steady adjoint equations.
/**
 * Solves the adjoint equations:
 * \f[
 *  \partial_{\mathbf{w}}\mathbf{r}^T \mathbf{c}\lambda = -\partial_{\mathbf{w}}J^{T}
 * \f]
 *
 * \f$\lambda\f$ is returned by Element::GetFirstDerivativesVector.
 * \f$\partial_{\mathbf{w}}\mathbf{r}^T\f$ is returned by Element::CalculateFirstDerivativesLHS.
 * \f$\partial_{\mathbf{w}}J^{T}\f$ is returned by ResponseFunction::CalculateFirstDerivativesGradient.
 * \f$\mathbf{c}$ is the rotation and slip condition matrix
 */
template <class TSparseSpace, class TDenseSpace>
class SimpleSteadyAdjointScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SimpleSteadyAdjointScheme);

    using IndexType = std::size_t;

    using BaseType = ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    using TSystemVectorType = typename BaseType::TSystemVectorType;

    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    using EquationIdVectorType = std::vector<IndexType>;

    using AssembledVectorType = std::unordered_map<IndexType, double>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit SimpleSteadyAdjointScheme(
        AdjointResponseFunction::Pointer pResponseFunction,
        const IndexType Dimension,
        const IndexType BlockSize,
        ElementRefinementProcess::Pointer pElementRefinementProcess = nullptr)
        : BaseType(pResponseFunction),
          mAdjointSlipUtilities(Dimension, BlockSize),
          mpElementRefinementProcess(pElementRefinementProcess)
    {
        // Allocate auxiliary memory.
        const int number_of_threads = ParallelUtilities::GetNumThreads();
        mAuxMatrices.resize(number_of_threads);

        KRATOS_INFO(this->Info()) << this->Info() << " created [ Dimensionality = " << Dimension << ", BlockSize = " << BlockSize << " ].\n";
    }

    /// Destructor.
    ~SimpleSteadyAdjointScheme() override = default;

    ///@}
    ///@name Operations
    ///@{

    void Initialize(ModelPart& rModelPart) override
    {
        BaseType::Initialize(rModelPart);

        if (mpElementRefinementProcess) {
            mpElementRefinementProcess->ExecuteInitialize();
        }
    }

    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentElement, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        Matrix& aux_matrix = mAuxMatrices[thread_id];
        CalculateEntityLHSContribution(rCurrentElement, aux_matrix, rLHS_Contribution,
                                       rEquationId, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        CalculateEntitySystemContributions(rCurrentCondition, rLHS_Contribution, rRHS_Contribution,
                                           rEquationId, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        const int thread_id = OpenMPUtils::ThisThread();
        Matrix& aux_matrix = mAuxMatrices[thread_id];
        CalculateEntityLHSContribution(rCurrentCondition, aux_matrix, rLHS_Contribution,
                                       rEquationId, rCurrentProcessInfo);

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
        return "SimpleSteadyAdjointScheme";
    }

    ///@}

protected:
    ///@name Protected Member Variables
    ///@{

    std::vector<Matrix> mAuxMatrices;

    const FluidAdjointSlipUtilities mAdjointSlipUtilities;

    ElementRefinementProcess::Pointer mpElementRefinementProcess;

    ///@}
    ///@name Private Operations
    ///@{

    template<class TEntityType>
    void CalculateEntitySystemContributions(
        TEntityType& rEntity,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        typename TEntityType::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int thread_id = OpenMPUtils::ThisThread();
        Matrix& residual_derivatives = mAuxMatrices[thread_id];

        CalculateEntityLHSContribution<TEntityType>(
            rEntity, residual_derivatives, rLHS_Contribution, rEquationId, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rEntity, residual_derivatives, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        if (rLHS_Contribution.size1() != 0) {
            auto& adjoint_values = this->mAdjointValues[thread_id];
            rEntity.GetValuesVector(adjoint_values);
            noalias(rRHS_Contribution) -= prod(rLHS_Contribution, adjoint_values);
        }

        KRATOS_CATCH("");
    }

    template<class TEntityType>
    void CalculateEntityLHSContribution(
        TEntityType& rEntity,
        Matrix& rEntityResidualFirstDerivatives,
        Matrix& rEntityRotatedResidualFirstDerivatives,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        rEntity.CalculateFirstDerivativesLHS(rEntityResidualFirstDerivatives, rCurrentProcessInfo);
        rEntity.EquationIdVector(rEquationId, rCurrentProcessInfo);

        mAdjointSlipUtilities.CalculateRotatedSlipConditionAppliedSlipVariableDerivatives(
            rEntityRotatedResidualFirstDerivatives, rEntityResidualFirstDerivatives, rEntity.GetGeometry());

        KRATOS_CATCH("");
    }

    ///@}
private:
    ///@name Private operations
    ///@{

    void CalculateResponseFunctionInterpolationError(ModelPart& rModelPart)
    {
        KRATOS_TRY

        const auto& dummy_element = *(rModelPart.ElementsBegin());
        const auto& dummy_geometry = dummy_element.GetGeometry();

        EquationIdVectorType equation_ids;
        dummy_element.EquationIdVector(equation_ids, rModelPart.GetProcessInfo());
        const IndexType number_of_dofs_per_node = equation_ids.size() / dummy_geometry.size();

        // reset the nodal interpolation error value holder
        block_for_each(rModelPart.Nodes(), [&](ModelPart::NodeType& rNode) {
            rNode.SetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR, Vector(number_of_dofs_per_node, 0.0));
        });

        // acceptable error for an element
        const double acceptable_element_error = 1.0 / rModelPart.GetCommunicator().GlobalNumberOfElements();

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

            Vector LocalErrors;
        };

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

                rTLS.LocalErrors.resize(number_of_dofs_per_node);

                rTLS.IsInitialized = true;
            }

            // using the process info of the original model part because
            // the thread local model part will not have up-to-date data if CloneTimeStep
            // is used in the  original model part. CloneTimeStep clones existing process info data to new process info data.
            const auto& r_process_info = rModelPart.GetProcessInfo();

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

            for (const auto p : rTLS.AssembledValues) {
                rTLS.ErrorValuesList[rTLS.ErrorIds.find(p.first)->second] += std::abs(p.second);
            }

            const double coeff = acceptable_element_error * rElement.GetGeometry().DomainSize() / (rElement.GetGeometry().size());

            for (auto& r_node : rElement.GetGeometry()) {
                r_node.SetLock();
                auto& r_values = r_node.GetValue(RESPONSE_FUNCTION_INTERPOLATION_ERROR);
                for (IndexType i = 0; i < number_of_dofs_per_node; ++i) {
                    r_values[i] += coeff / rTLS.ErrorValuesList[i];
                }
                r_node.UnSetLock();
            }

        });

        rModelPart.GetCommunicator().AssembleNonHistoricalData(RESPONSE_FUNCTION_INTERPOLATION_ERROR);

        KRATOS_CATCH("");
    }

    ///@}


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

}; /* Class SimpleSteadyAdjointScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED defined */
