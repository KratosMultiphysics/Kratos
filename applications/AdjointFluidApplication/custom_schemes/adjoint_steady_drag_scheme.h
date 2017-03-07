//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ADJOINT_STEADY_DRAG_SCHEME )
#define  KRATOS_ADJOINT_STEADY_DRAG_SCHEME

// System includes
#include <cmath>
#include <iostream>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/openmp_utils.h"
#include "solving_strategies/schemes/scheme.h"

// Application includes

namespace Kratos
{

///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// A scheme for steady adjoint equations.
/**
 * The residual vector of the forward problem is:
 * \f[
 *  \mathbf{f}(\mathbf{w};\mathbf{s}) = 0
 * \f]
 *
 * The adjoint equations are:
 * \f[
 *  \partial_{\mathbf{w}}\mathbf{f}^T \lambda = -\partial_{\mathbf{w}}J^{T}
 * \f]
 *
 * with objective function \f$J=J(\mathbf{w};\mathbf{s})\f$.
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class AdjointSteadyDragScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( AdjointSteadyDragScheme );

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TSystemMatrixType SystemMatrixType;

    typedef typename BaseType::TSystemVectorType SystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef Element::IndexType IndexType;

    typedef Element::SizeType SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    AdjointSteadyDragScheme()
    : Scheme<TSparseSpace,TDenseSpace>()
      {
        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mDragFlagVector.resize(NumThreads);
        mAdjointVelocity.resize(NumThreads);
      }

    /// Destructor.
    virtual ~AdjointSteadyDragScheme() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void InitializeSolutionStep(
            ModelPart& rModelPart,
            SystemMatrixType& rA,
            SystemVectorType& rDx,
            SystemVectorType& rb
             )
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        // check for valid drag direction
        const int DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
        array_1d<double, 3>& rDragDirection = rCurrentProcessInfo[DRAG_DIRECTION];
        double magnitude = rDragDirection[0] * rDragDirection[0];
        for (int k = 1; k < DomainSize; k++)
            magnitude += rDragDirection[k] * rDragDirection[k];

        if (std::abs(magnitude - 1.0) > 1e-3)
        {
            std::cout << "Warning: DRAG_DIRECTION is not a unit vector" << std::endl;
        }

        KRATOS_CATCH("")
    }

    virtual void Update(
                ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb
                )
    {
        KRATOS_TRY

        const int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector Partition;

        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofSetBegin = rDofSet.begin() + Partition[k];
            typename DofsArrayType::iterator DofSetEnd = rDofSet.begin() + Partition[k + 1];

            for (auto itDof = DofSetBegin; itDof != DofSetEnd; itDof++)
            {
                if (itDof->IsFree())
                {
                    itDof->GetSolutionStepValue() += TSparseSpace::GetValue(rDx, itDof->EquationId());
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        int ThreadId = OpenMPUtils::ThisThread();

        this->GetDragFlagVector(mDragFlagVector[ThreadId], pCurrentElement, rCurrentProcessInfo);

        // adjoint system matrix
        pCurrentElement->Calculate(ADJOINT_MATRIX_1,rLHS_Contribution,rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size2())
            rRHS_Contribution.resize(rLHS_Contribution.size2(), false);

        noalias(rRHS_Contribution) = -prod(rLHS_Contribution,mDragFlagVector[ThreadId]);

        // residual form
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId],0);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution,mAdjointVelocity[ThreadId]);

        pCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    virtual void Calculate_LHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;

        RHS_Contribution.resize(LHS_Contribution.size1(), false);

        CalculateSystemContributions(
                pCurrentElement,
                LHS_Contribution,
                RHS_Contribution,
                EquationId,
                CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY


        KRATOS_CATCH("")
    }

    virtual void Condition_Calculate_LHS_Contribution(
            Condition::Pointer pCurrentCondition,
            LocalSystemMatrixType& LHS_Contribution,
            Condition::EquationIdVectorType& EquationId,
            ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    virtual void GetElementalDofList(
            Element::Pointer rCurrentElement,
            Element::DofsVectorType& ElementalDofList,
            ProcessInfo& CurrentProcessInfo)
        {
            rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
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

    std::vector< LocalSystemVectorType > mDragFlagVector;
    std::vector< LocalSystemVectorType > mAdjointVelocity;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void GetDragFlagVector(LocalSystemVectorType& rOutput,
            Element::Pointer pElement, ProcessInfo& rProcessInfo)
    {
        const SizeType DomainSize = static_cast<SizeType>(rProcessInfo[DOMAIN_SIZE]);
        const SizeType NumNodes = pElement->GetGeometry().PointsNumber();
        const SizeType LocalSize = (DomainSize + 1) * NumNodes;

        if (rOutput.size() != LocalSize)
        {
            rOutput.resize(LocalSize, false);
        }

        array_1d<double, 3>& rDragDirection = rProcessInfo[DRAG_DIRECTION];
        IndexType LocalIndex = 0;
        for (IndexType iNode = 0; iNode < NumNodes; ++iNode)
        {
            if (pElement->GetGeometry()[iNode].Is(STRUCTURE))
            {
                for (SizeType d = 0; d < DomainSize; d++)
                {
                    rOutput[LocalIndex++] = rDragDirection[d];
                }
            }
            else
            {
                for (SizeType d = 0; d < DomainSize; d++)
                {
                    rOutput[LocalIndex++] = 0.0;
                }
            }

            rOutput[LocalIndex++] = 0.0; // pressure dof
        }
    }

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

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

///@} // Adjoint Fluid Application group

}  /* namespace Kratos.*/

#endif /* KRATOS_ADJOINT_STEADY_DRAG_SCHEME defined */
