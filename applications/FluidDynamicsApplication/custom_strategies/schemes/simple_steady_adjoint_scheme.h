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

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"
#include "utilities/parallel_utilities.h"

// Application includes
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

    using BaseType = ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>;

    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;

    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit SimpleSteadyAdjointScheme(
        AdjointResponseFunction::Pointer pResponseFunction,
        const IndexType Dimension,
        const IndexType BlockSize)
        : BaseType(pResponseFunction),
          mAdjointSlipUtilities(Dimension, BlockSize)
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

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "SimpleSteadyAdjointScheme";
    }

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    std::vector<Matrix> mAuxMatrices;

    const FluidAdjointSlipUtilities mAdjointSlipUtilities;

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

}; /* Class SimpleSteadyAdjointScheme */

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SIMPLE_STEADY_ADJOINT_SCHEME_H_INCLUDED defined */
