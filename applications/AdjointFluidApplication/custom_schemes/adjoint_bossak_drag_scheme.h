//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//   License:        BSD License
//   Kratos default license: kratos/license.txt
//
//   Project Name:        $AdjointFluidApplication        $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:          November 2016   $
//   Revision:            $Revision:                0.0   $

#if !defined(KRATOS_ADJOINT_BOSSAK_DRAG_SCHEME )
#define  KRATOS_ADJOINT_BOSSAK_DRAG_SCHEME

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

/// A scheme for unsteady adjoint equations using Bossak time discretization.
/**
 * The forward Bossak equations are:
 * \f[
 *  \mathbf{M}\dot{\mathbf{w}}^{n-\alpha} = \mathbf{f}(\mathbf{w}^{n};\mathbf{s})
 * \f]
 * \f[
 *  \dot{\mathbf{w}}^{n-\alpha} = (1 - \alpha) \dot{\mathbf{w}}^n + \alpha \dot{\mathbf{w}}^{n-1}
 * \f]
 * \f[
 *  \dot{\mathbf{w}}^n = \frac{\mathbf{w}^n - \mathbf{w}^{n-1}}{\gamma \Delta t} + \frac{\gamma - 1}{\gamma}\dot{\mathbf{w}}^{n-1}
 * \f]
 *
 * The adjoint Bossak equations are:
 * \f[
 *  \frac{1}{\gamma - 1} (\dot{\lambda}^n - \dot{\lambda}^{n+1}) + (\partial_{\mathbf{w}^n}\mathbf{f}^n - \partial_{\mathbf{w}^n}(\mathbf{M}^n\dot{\mathbf{w}}^{n-\alpha}))^T \lambda^n = -\partial_{\mathbf{w}^n}J^{nT}
 * \f]
 * \f[
 *  \frac{1}{\gamma - 1} \dot{\lambda}^n = \frac{1}{\gamma} \dot{\lambda}^{n+1} - \frac{1 - \alpha}{\gamma \Delta t}M^{nT} \lambda^n - \frac{\alpha}{\gamma \Delta t}M^{(n+1)T} \lambda^{n+1} + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{nT} + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{(n+1)T}
 * \f]
 *
 * with objective function \f$J^n=J(\mathbf{w}^n,\dot{\mathbf{w}}^n,\dot{\mathbf{w}}^{n-1};\mathbf{s})\f$.
 */
template<class TSparseSpace,
         class TDenseSpace
         >
class AdjointBossakDragScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( AdjointBossakDragScheme );

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
    AdjointBossakDragScheme(double AlphaBossak /*=-0.3*/)
    : Scheme<TSparseSpace,TDenseSpace>()
      {
        mAlphaBossak = AlphaBossak;
        mGammaNewmark = 0.5 - mAlphaBossak;
        mInvGamma = 1.0 / mGammaNewmark;
        mInvGammaMinusOne = 1.0 / (mGammaNewmark - 1.0);

        //Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mDragFlagVector.resize(NumThreads);
        mAdjointVelocity.resize(NumThreads);
        mAdjointAcceleration.resize(NumThreads);
        mMassMatrix.resize(NumThreads);
      }

    /// Destructor.
    virtual ~AdjointBossakDragScheme() {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& r_model_part)
    {
        KRATOS_TRY

        BaseType::Initialize(r_model_part);

        // this switch is used to make the discrete sensitivities exact in the
        // first time step. it is not important for most problems.
        mMass1Switch = 0.0;

        KRATOS_CATCH("")
    }

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
        double DeltaTime = -rCurrentProcessInfo[DELTA_TIME]; // DELTA_TIME < 0

        if (DeltaTime <= 0)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "detected for adjoint solution DELTA_TIME >= 0", "")
        }

        mInvDt = 1.0 / DeltaTime;

        // check for valid drag direction
        const int DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
        array_1d<double, 3>& rDragDirection = rCurrentProcessInfo[DRAG_DIRECTION];
        double magnitude = rDragDirection[0] * rDragDirection[0];
        for (int k = 1; k < DomainSize; k++)
        {
            magnitude += rDragDirection[k] * rDragDirection[k];
        }

        if (std::abs(magnitude - 1.0) > 1e-3)
        {
            std::cout << "Warning: DRAG_DIRECTION is not a unit vector" << std::endl;
        }

        const int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector Partition;
        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + Partition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + Partition[k + 1];

            for (auto it = NodesBegin; it != NodesEnd; it++)
            {
                it->GetValue(NODAL_AREA) = 0.0;
            }
        }

        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfElements(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::ElementIterator ElementsBegin = rModelPart.ElementsBegin() + Partition[k];
            ModelPart::ElementIterator ElementsEnd = rModelPart.ElementsBegin() + Partition[k + 1];

            for (auto itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                for (IndexType iNode = 0; iNode < itElem->GetGeometry().PointsNumber(); ++iNode)
                {
                    itElem->GetGeometry()[iNode].SetLock();
                    itElem->GetGeometry()[iNode].GetValue(NODAL_AREA) += 1.0;
                    itElem->GetGeometry()[iNode].UnSetLock();
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void FinalizeSolutionStep(
            ModelPart& rModelPart,
            SystemMatrixType& rA,
            SystemVectorType& rDx,
            SystemVectorType& rb)
        {
            KRATOS_TRY

            BaseType::FinalizeSolutionStep(rModelPart,rA,rDx,rb);
            mMass1Switch = 1.0;

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

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        const IndexType DomainSize = static_cast<IndexType>(rCurrentProcessInfo[DOMAIN_SIZE]);

        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + Partition[k];
            ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + Partition[k + 1];

            for (auto it = NodesBegin; it != NodesEnd; it++)
            {
                array_1d<double,3>& rCurrentAdjointAcceleration = it->FastGetSolutionStepValue(ADJOINT_ACCELERATION,0);
                const array_1d<double,3>& rOldAdjointAcceleration = it->FastGetSolutionStepValue(ADJOINT_ACCELERATION,1);
                for (IndexType d = 0; d < DomainSize; ++d)
                {
                    rCurrentAdjointAcceleration[d] = (mGammaNewmark - 1.0) * mInvGamma * rOldAdjointAcceleration[d];
                }
            }
        }

        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfElements(), NumThreads, Partition);
        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::ElementIterator ElementsBegin = rModelPart.ElementsBegin() + Partition[k];
            ModelPart::ElementIterator ElementsEnd = rModelPart.ElementsBegin() + Partition[k + 1];

            for (auto itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                this->GetDragFlagVector(mDragFlagVector[k], *(itElem.base()), rCurrentProcessInfo);

                itElem->Calculate(MASS_MATRIX_1,mMassMatrix[k],rCurrentProcessInfo);
                mMassMatrix[k] = trans(mMassMatrix[k]);

                itElem->GetFirstDerivativesVector(mAdjointVelocity[k],1);

                mAdjointAcceleration[k] = -mAlphaBossak * prod(mMassMatrix[k], mAdjointVelocity[k] + mDragFlagVector[k]);

                itElem->Calculate(MASS_MATRIX_0,mMassMatrix[k],rCurrentProcessInfo);
                mMassMatrix[k] = trans(mMassMatrix[k]);
                itElem->GetFirstDerivativesVector(mAdjointVelocity[k],0);
                noalias(mAdjointAcceleration[k]) += -(1.0 - mAlphaBossak) * prod(mMassMatrix[k], mAdjointVelocity[k] + mDragFlagVector[k]);

                noalias(mAdjointAcceleration[k]) = (mGammaNewmark - 1.0) * mInvGamma * mInvDt * mAdjointAcceleration[k];

                IndexType LocalIndex = 0;
                for (IndexType iNode = 0; iNode < itElem->GetGeometry().PointsNumber(); ++iNode)
                {
                    itElem->GetGeometry()[iNode].SetLock();
                    array_1d<double,3>& rCurrentAdjointAcceleration = itElem->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_ACCELERATION);
                    for (IndexType d = 0; d < DomainSize; ++d)
                    {
                        rCurrentAdjointAcceleration[d] += mAdjointAcceleration[k][LocalIndex++];
                    }
                    itElem->GetGeometry()[iNode].UnSetLock();
                    LocalIndex++; // pressure dof
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

        // old adjoint acceleration
        pCurrentElement->GetSecondDerivativesVector(rRHS_Contribution,1);
        IndexType LocalIndex = 0;
        const int DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
        for (IndexType iNode = 0; iNode < pCurrentElement->GetGeometry().PointsNumber(); ++iNode)
        {
            double Weight = mInvGamma * mInvGammaMinusOne / pCurrentElement->GetGeometry()[iNode].GetValue(NODAL_AREA);
            for (IndexType d = 0; d < DomainSize; ++d)
            {
                rRHS_Contribution[LocalIndex++] *= Weight;
            }
            LocalIndex++; // pressure dof
        }

        // old adjoint velocity, d (old drag) / d (primal acceleration)
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId],1);
        pCurrentElement->Calculate(MASS_MATRIX_1,mMassMatrix[ThreadId],rCurrentProcessInfo);
        mMassMatrix[ThreadId] = trans(mMassMatrix[ThreadId]);
        noalias(rRHS_Contribution) += mMass1Switch * mAlphaBossak * mInvGamma * mInvDt
                * prod(mMassMatrix[ThreadId],mAdjointVelocity[ThreadId] + mDragFlagVector[ThreadId]);

        // d (drag) / d (primal velocity)
        pCurrentElement->Calculate(ADJOINT_MATRIX_2,rLHS_Contribution,rCurrentProcessInfo);

        // d (drag) / d (primal acceleration)
        pCurrentElement->Calculate(MASS_MATRIX_0,mMassMatrix[ThreadId],rCurrentProcessInfo);
        mMassMatrix[ThreadId] = trans(mMassMatrix[ThreadId]);

        // adjoint system matrix
        noalias(rLHS_Contribution) -= mInvGamma * mInvDt * (1.0 - mAlphaBossak) * mMassMatrix[ThreadId];

        noalias(rRHS_Contribution) -= prod(rLHS_Contribution,mDragFlagVector[ThreadId]);

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

    double mAlphaBossak;
    double mGammaNewmark;
    double mInvDt;
    double mInvGamma;
    double mInvGammaMinusOne;
    double mMass1Switch;
    std::vector< LocalSystemVectorType > mDragFlagVector;
    std::vector< LocalSystemVectorType > mAdjointVelocity;
    std::vector< LocalSystemVectorType > mAdjointAcceleration;
    std::vector< LocalSystemMatrixType > mMassMatrix;

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

#endif /* KRATOS_ADJOINT_BOSSAK_DRAG_SCHEME defined */

