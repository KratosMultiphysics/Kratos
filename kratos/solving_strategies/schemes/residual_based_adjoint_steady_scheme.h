//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_RESIDUAL_BASED_ADJOINT_STEADY_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_ADJOINT_STEADY_SCHEME_H_INCLUDED

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/schemes/residual_based_adjoint_static_scheme.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for steady adjoint equations.
/**
 * Solves the adjoint equations:
 * \f[
 *  \partial_{\mathbf{w}}\mathbf{r}^T \lambda = -\partial_{\mathbf{w}}J^{T}
 * \f]
 *
 * \f$\lambda\f$ is returned by Element::GetFirstDerivativesVector.
 * \f$\partial_{\mathbf{w}}\mathbf{r}^T\f$ is returned by Element::CalculateFirstDerivativesLHS.
 * \f$\partial_{\mathbf{w}}J^{T}\f$ is returned by ResponseFunction::CalculateFirstDerivativesGradient.
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedAdjointSteadyScheme : public ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedAdjointSteadyScheme);

    typedef ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit ResidualBasedAdjointSteadyScheme(AdjointResponseFunction::Pointer pResponseFunction)
        : ResidualBasedAdjointStaticScheme<TSparseSpace, TDenseSpace>(pResponseFunction)
    {
    }

    /// Destructor.
    ~ResidualBasedAdjointSteadyScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();
        const auto& r_const_elem_ref = rCurrentElement;
        rCurrentElement.CalculateFirstDerivativesLHS(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        r_const_elem_ref.GetFirstDerivativesVector(this->mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, this->mAdjointValues[thread_id]);

        r_const_elem_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Element::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement.CalculateFirstDerivativesLHS(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Condition::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;
        rCurrentCondition.CalculateFirstDerivativesLHS(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        this->mpResponseFunction->CalculateFirstDerivativesGradient(
            rCurrentCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        r_const_cond_ref.GetFirstDerivativesVector(this->mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, this->mAdjointValues[thread_id]);

        r_const_cond_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Condition::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition.CalculateFirstDerivativesLHS(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; /* Class ResidualBasedAdjointSteadyScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_STEADY_SCHEME_H_INCLUDED defined */
