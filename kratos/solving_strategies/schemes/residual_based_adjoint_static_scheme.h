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

#if !defined(KRATOS_RESIDUAL_BASED_ADJOINT_STATIC_SCHEME_H_INCLUDED)
#define KRATOS_RESIDUAL_BASED_ADJOINT_STATIC_SCHEME_H_INCLUDED

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

        int num_threads = OpenMPUtils::GetNumThreads();
        mAdjointValues.resize(num_threads);
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

    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        BaseType::Initialize(rModelPart);

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

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

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            *pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        pCurrentElement->GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                    LocalSystemMatrixType& rLHS_Contribution,
                                    Element::EquationIdVectorType& rEquationId,
                                    ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
                                                LocalSystemMatrixType& rLHS_Contribution,
                                                LocalSystemVectorType& rRHS_Contribution,
                                                Condition::EquationIdVectorType& rEquationId,
                                                ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        pCurrentCondition->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        mpResponseFunction->CalculateGradient(
            *pCurrentCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // Calculate system contributions in residual form.
        pCurrentCondition->GetValuesVector(mAdjointValues[thread_id]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValues[thread_id]);

        pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
    }

    void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              Condition::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        pCurrentCondition->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("");
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

#endif /* KRATOS_RESIDUAL_BASED_ADJOINT_STATIC_SCHEME_H_INCLUDED defined */
