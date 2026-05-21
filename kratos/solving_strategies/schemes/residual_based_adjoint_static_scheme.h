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

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"
#include "solving_strategies/schemes/scheme.h"
#include "response_functions/adjoint_response_function.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/// A scheme for static adjoint equations.
/**
 * Solves the adjoint equations:
 * \f[
 *  \partial_{\mathbf{u}}\mathbf{r}^T \lambda = -\partial_{\mathbf{u}}J^{T}
 * \f]
 *
 * \f$\lambda\f$ is returned by Element::GetValuesVector.
 * \f$\partial_{\mathbf{u}}\mathbf{r}^T\f$ is returned by Element::CalculateLeftHandSide.
 * \f$\partial_{\mathbf{u}}J^{T}\f$ is returned by ResponseFunction::CalculateGradient.
 */
template <class TSparseSpace, class TDenseSpace>
class ResidualBasedAdjointStaticScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedAdjointStaticScheme);

    typedef Scheme<TSparseSpace, TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType SystemMatrixType;

    typedef typename BaseType::TSystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    explicit ResidualBasedAdjointStaticScheme(AdjointResponseFunction::Pointer pResponseFunction)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        mpResponseFunction = pResponseFunction;

        int num_threads = ParallelUtilities::GetNumThreads();
        mAdjointValues.resize(num_threads);
        mLHS.resize(num_threads);
    }

    /// Destructor.
    ~ResidualBasedAdjointStaticScheme() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set the Response Function
     *
     * This sets the response function used in the sensitivity builder. This
     * is useful in cases where the LHS of the adjoint problem does not change,
     * but the RHS changes due to change in the the response function. In these
     * cases, this allows re-use of the already constructed LHS with different
     * RHSs.
     *
     * @param pResponseFunction         New Response function to be set.
     */
    void SetResponseFunction(AdjointResponseFunction::Pointer pResponseFunction)
    {
        mpResponseFunction = pResponseFunction;
    }

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        block_for_each(rModelPart.Nodes(), [](Node& rNode){
            for (auto& r_dof : rNode.GetDofs()) {
                if (r_dof->IsFree()) {
                    r_dof->GetSolutionStepValue() = 0.0;
                }
            }
        });

        BaseType::Initialize(rModelPart);

        KRATOS_CATCH("");
    }

    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Update degrees of freedom: adjoint variables associated to the
        // residual of the physical problem.
        this->mpDofUpdater->UpdateDofs(rDofSet, rDx);

        KRATOS_CATCH("");
    }

    void CalculateSystemContributions(Element& rCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHSContribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rCurrentElement;

        rCurrentElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHSContribution.size() != rLHS_Contribution.size1())
            rRHSContribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentElement, rLHS_Contribution, rRHSContribution, rCurrentProcessInfo);

        noalias(rRHSContribution) = -rRHSContribution;

        // Calculate system contributions in residual form.
        r_const_elem_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHSContribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        r_const_elem_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Element& rCurrentElement,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Element::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        int thread_id = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rCurrentElement;
        auto& lhs = mLHS[thread_id];

        rCurrentElement.CalculateLeftHandSide(lhs, rCurrentProcessInfo);

        if (rRHSContribution.size() != lhs.size1())
            rRHSContribution.resize(lhs.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentElement, lhs, rRHSContribution, rCurrentProcessInfo);

        noalias(rRHSContribution) = -rRHSContribution;

        // Calculate system contributions in residual form.
        r_const_elem_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHSContribution) -= prod(lhs, mAdjointValues[thread_id]);

        r_const_elem_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void CalculateSystemContributions(Condition& rCurrentCondition,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHSContribution,
                                      Condition::EquationIdVectorType& rEquationId,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();
        const auto& r_const_cond_ref = rCurrentCondition;
        rCurrentCondition.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHSContribution.size() != rLHS_Contribution.size1())
            rRHSContribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentCondition, rLHS_Contribution, rRHSContribution, rCurrentProcessInfo);

        noalias(rRHSContribution) = -rRHSContribution;

        // Calculate system contributions in residual form.
        r_const_cond_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHSContribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        r_const_cond_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateLHSContribution(Condition& rCurrentCondition,
                                  LocalSystemMatrixType& rLHS_Contribution,
                                  Condition::EquationIdVectorType& rEquationId,
                                  const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        rCurrentCondition.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        rCurrentCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHSContribution,
        Condition::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        int thread_id = OpenMPUtils::ThisThread();

        const auto& r_const_elem_ref = rCurrentCondition;
        auto& lhs = mLHS[thread_id];

        rCurrentCondition.CalculateLeftHandSide(lhs, rCurrentProcessInfo);

        if (rRHSContribution.size() != lhs.size1())
            rRHSContribution.resize(lhs.size1(), false);

        mpResponseFunction->CalculateGradient(
            rCurrentCondition, lhs, rRHSContribution, rCurrentProcessInfo);

        noalias(rRHSContribution) = -rRHSContribution;

        // Calculate system contributions in residual form.
        r_const_elem_ref.GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHSContribution) -= prod(lhs, mAdjointValues[thread_id]);

        r_const_elem_ref.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }

    void Clear() override
    {
        this->mpDofUpdater->Clear();
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

    AdjointResponseFunction::Pointer mpResponseFunction;
    std::vector<LocalSystemVectorType> mAdjointValues;
    std::vector<LocalSystemMatrixType> mLHS;

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

    typename TSparseSpace::DofUpdaterPointerType mpDofUpdater =
        TSparseSpace::CreateDofUpdater();

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

}; /* Class ResidualBasedAdjointStaticScheme */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/
