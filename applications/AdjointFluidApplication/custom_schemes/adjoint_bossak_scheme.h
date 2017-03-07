//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ADJOINT_BOSSAK_SCHEME)
#define KRATOS_ADJOINT_BOSSAK_SCHEME

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
#include "custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A scheme for unsteady adjoint equations using Bossak time discretization.
/**
 * The forward Bossak equations are:
 * \f[
 * \mathbf{M}\dot{\mathbf{w}}^{n-\alpha} = \mathbf{f}(\mathbf{w}^{n};\mathbf{s})
 * \f]
 * \f[
 * \dot{\mathbf{w}}^{n-\alpha}
 * = (1 - \alpha) \dot{\mathbf{w}}^n + \alpha \dot{\mathbf{w}}^{n-1}
 * \f]
 * \f[
 * \dot{\mathbf{w}}^n
 * = \frac{\mathbf{w}^n - \mathbf{w}^{n-1}}{\gamma \Delta t}
 * + \frac{\gamma - 1}{\gamma}\dot{\mathbf{w}}^{n-1}
 * \f]
 *
 * The adjoint Bossak equations are:
 * \f[
 * \frac{1}{\gamma - 1} (\dot{\lambda}^n - \dot{\lambda}^{n+1})
 * + (\partial_{\mathbf{w}^n}\mathbf{f}^n
 * -\partial_{\mathbf{w}^n}(\mathbf{M}^n\dot{\mathbf{w}}^{n-\alpha}))^T\lambda^n
 * = -\partial_{\mathbf{w}^n}J^{nT}
 * \f]
 * \f[
 * \frac{1}{\gamma - 1} \dot{\lambda}^n
 * = \frac{1}{\gamma} \dot{\lambda}^{n+1}
 * - \frac{1 - \alpha}{\gamma \Delta t}M^{nT} \lambda^n
 * - \frac{\alpha}{\gamma \Delta t}M^{(n+1)T} \lambda^{n+1}
 * + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{nT}
 * + \frac{1}{\gamma \Delta t}\partial_{\dot{\mathbf{w}}^n}J^{(n+1)T}
 * \f]
 *
 * with objective function
 *\f$J^n=J(\mathbf{w}^n,\dot{\mathbf{w}}^n,\dot{\mathbf{w}}^{n-1};\mathbf{s})\f$.
 */
template <class TSparseSpace, class TDenseSpace>
class AdjointBossakScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointBossakScheme);

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
    AdjointBossakScheme(double AlphaBossak /*=-0.3*/, ObjectiveFunction::Pointer pObjectiveFunction)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        mAlphaBossak = AlphaBossak;
        mGammaNewmark = 0.5 - mAlphaBossak;
        mInvGamma = 1.0 / mGammaNewmark;
        mInvGammaMinusOne = 1.0 / (mGammaNewmark - 1.0);
        mpObjectiveFunction = pObjectiveFunction;

        // Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mAdjointVelocity.resize(NumThreads);
        mAdjointAcceleration.resize(NumThreads);
        mObjectiveGradient.resize(NumThreads);
        mAdjointMassMatrix.resize(NumThreads);
    }

    /// Destructor.
    virtual ~AdjointBossakScheme()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY

        BaseType::Initialize(rModelPart);

        // this switch is used to make the discrete sensitivities exact in the
        // first time step. it is not important for most problems.
        mMass1Switch = 0.0;

        mpObjectiveFunction->Initialize(rModelPart);

        KRATOS_CATCH("")
    }

    virtual void InitializeSolutionStep(ModelPart& rModelPart,
                                        SystemMatrixType& rA,
                                        SystemVectorType& rDx,
                                        SystemVectorType& rb)
    {
        KRATOS_TRY

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        double DeltaTime = -rCurrentProcessInfo[DELTA_TIME]; // DELTA_TIME < 0

        if (DeltaTime <= 0)
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "detected for adjoint solution DELTA_TIME >= 0",
                               "")
        }

        mInvDt = 1.0 / DeltaTime;

#pragma omp parallel
        {
            ModelPart::NodeIterator NodesBegin;
            ModelPart::NodeIterator NodesEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Nodes(), NodesBegin, NodesEnd);

            for (auto it = NodesBegin; it != NodesEnd; ++it)
            {
                it->GetValue(NODAL_AREA) = 0.0;
            }
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElementsBegin;
            ModelPart::ElementIterator ElementsEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElementsBegin, ElementsEnd);

            for (auto it = ElementsBegin; it != ElementsEnd; ++it)
            {
                for (unsigned int iNode = 0; iNode < it->GetGeometry().PointsNumber(); ++iNode)
                {
                    itElem->GetGeometry()[iNode].SetLock();
                    itElem->GetGeometry()[iNode].GetValue(NODAL_AREA) += 1.0;
                    itElem->GetGeometry()[iNode].UnSetLock();
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void FinalizeSolutionStep(ModelPart& rModelPart,
                                      SystemMatrixType& rA,
                                      SystemVectorType& rDx,
                                      SystemVectorType& rb)
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);
        mMass1Switch = 1.0;

        KRATOS_CATCH("")
    }

    virtual void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb)
    {
        KRATOS_TRY

        const int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector Partition;

        OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, Partition);
#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            typename DofsArrayType::iterator DofSetBegin =
                rDofSet.begin() + Partition[k];
            typename DofsArrayType::iterator DofSetEnd =
                rDofSet.begin() + Partition[k + 1];

            for (auto itDof = DofSetBegin; itDof != DofSetEnd; itDof++)
            {
                if (itDof->IsFree())
                {
                    itDof->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(rDx, itDof->EquationId());
                }
            }
        }

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        const unsigned int DomainSize =
            static_cast<unsigned int>(rCurrentProcessInfo[DOMAIN_SIZE]);

        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfNodes(), NumThreads, Partition);
#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + Partition[k];
            ModelPart::NodeIterator NodesEnd =
                rModelPart.NodesBegin() + Partition[k + 1];

            for (auto it = NodesBegin; it != NodesEnd; it++)
            {
                array_1d<double, 3>& rCurrentAdjointAcceleration =
                    it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 0);
                const array_1d<double, 3>& rOldAdjointAcceleration =
                    it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 1);
                for (unsigned int d = 0; d < DomainSize; ++d)
                {
                    rCurrentAdjointAcceleration[d] =
                        (mGammaNewmark - 1.0) * mInvGamma * rOldAdjointAcceleration[d];
                }
            }
        }

        OpenMPUtils::DivideInPartitions(rModelPart.NumberOfElements(), NumThreads, Partition);
#pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();

            ModelPart::ElementIterator ElementsBegin =
                rModelPart.ElementsBegin() + Partition[k];
            ModelPart::ElementIterator ElementsEnd =
                rModelPart.ElementsBegin() + Partition[k + 1];

            for (auto it = ElementsBegin; it != ElementsEnd; ++it)
            {
                // transposed gradient of old element residual w.r.t.
                // acceleration
                it->Calculate(MASS_MATRIX_1, mAdjointMassMatrix[k], rCurrentProcessInfo);
                mAdjointMassMatrix[k] = -mAlphaBossak * trans(mAdjointMassMatrix[k]);

                // d (old objective) / d (primal acceleration)
                mpObjectiveFunction->CalculateAdjointAccelerationContribution(
                    *it, mAdjointMassMatrix[k], mObjectiveFunctionGradient[k], rCurrentProcessInfo);

                // old adjoint velocity
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 1);

                // terms depending on the old mass matrix
                mAdjointAcceleration[k] =
                    prod(mAdjointMassMatrix[k], mAdjointVelocity[k]) +
                    mObjectiveFunctionGradient[k];

                // transposed gradient of element residual w.r.t. acceleration
                it->Calculate(MASS_MATRIX_0, mAdjointMassMatrix[k], rCurrentProcessInfo);
                mAdjointMassMatrix[k] =
                    -(1.0 - mAlphaBossak) * trans(mAdjointMassMatrix[k]);

                // d (objective) / d (primal acceleration)
                mpObjectiveFunction->CalculateAdjointAccelerationContribution(
                    *it, mAdjointMassMatrix[k], mObjectiveFunctionGradient[k], rCurrentProcessInfo);

                // adjoint velocity
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 0);

                // terms depending on the mass matrix
                noalias(mAdjointAcceleration[k]) +=
                    prod(mAdjointMassMatrix[k], mAdjointVelocity[k]) +
                    mObjectiveFunctionGradient[k];

                noalias(mAdjointAcceleration[k]) = (mGammaNewmark - 1.0) * mInvGamma *
                                                   mInvDt * mAdjointAcceleration[k];

                unsigned int LocalIndex = 0;
                for (unsigned int iNode = 0; iNode < it->GetGeometry().PointsNumber(); ++iNode)
                {
                    it->GetGeometry()[iNode].SetLock();
                    array_1d<double, 3>& rCurrentAdjointAcceleration =
                        it->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_ACCELERATION);
                    for (unsigned int d = 0; d < DomainSize; ++d)
                    {
                        rCurrentAdjointAcceleration[d] +=
                            mAdjointAcceleration[k][LocalIndex++];
                    }
                    it->GetGeometry()[iNode].UnSetLock();
                    LocalIndex++; // pressure dof
                }
            }
        }

        KRATOS_CATCH("")
    }

    virtual void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int ThreadId = OpenMPUtils::ThisThread();

        // old adjoint acceleration
        pCurrentElement->GetSecondDerivativesVector(rRHS_Contribution, 1);
        unsigned int LocalIndex = 0;
        const int DomainSize = rCurrentProcessInfo[DOMAIN_SIZE];
        for (unsigned int iNode = 0; iNode < pCurrentElement->GetGeometry().PointsNumber(); ++iNode)
        {
            double Weight = mInvGamma * mInvGammaMinusOne /
                            pCurrentElement->GetGeometry()[iNode].GetValue(NODAL_AREA);
            for (int d = 0; d < DomainSize; ++d)
                rRHS_Contribution[LocalIndex++] *= Weight;
            LocalIndex++; // pressure dof
        }

        // transposed gradient of old element residual w.r.t. acceleration
        pCurrentElement->Calculate(
            MASS_MATRIX_1, mAdjointMassMatrix[ThreadId], rCurrentProcessInfo);
        mAdjointMassMatrix[ThreadId] =
            -mAlphaBossak * trans(mAdjointMassMatrix[ThreadId]);

        // d (old objective) / d (primal acceleration)
        mpObjectiveFunction->CalculateAdjointAccelerationContribution(
            *pCurrentElement,
            mAdjointMassMatrix[ThreadId],
            mObjectiveFunctionGradient[ThreadId],
            rCurrentProcessInfo);

        // old adjoint velocity
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId], 1);

        // terms depending on the old mass matrix
        noalias(rRHS_Contribution) -=
            mMass1Switch * mInvGamma * mInvDt *
            (prod(mAdjointMassMatrix[ThreadId], mAdjointVelocity[ThreadId]) +
             mObjectiveFunctionGradient[ThreadId]);

        // transposed gradient of element residual w.r.t. acceleration
        pCurrentElement->Calculate(
            MASS_MATRIX_0, mAdjointMassMatrix[ThreadId], rCurrentProcessInfo);
        mAdjointMassMatrix[ThreadId] =
            -(1.0 - mAlphaBossak) * trans(mAdjointMassMatrix[ThreadId]);

        // d (objective) / d (primal acceleration)
        mpObjectiveFunction->CalculateAdjointAccelerationContribution(
            *pCurrentElement,
            mAdjointMassMatrix[ThreadId],
            mObjectiveFunctionGradient[ThreadId],
            rCurrentProcessInfo);
        noalias(rRHS_Contribution) -=
            mInvGamma * mInvDt * mObjectiveFunctionGradient[ThreadId];

        // transposed gradient of element residual w.r.t. primal
        pCurrentElement->Calculate(ADJOINT_MATRIX_2, rLHS_Contribution, rCurrentProcessInfo);

        // d (objective) / d (primal)
        mpObjectiveFunction->CalculateAdjointVelocityContribution(
            *pCurrentElement, rLHS_Contribution, mObjectiveFunctionGradient[ThreadId], rCurrentProcessInfo);
        noalias(rRHS_Contribution) -= mObjectiveFunctionGradient[ThreadId];

        noalias(rLHS_Contribution) += mInvGamma * mInvDt * mAdjointMassMatrix[ThreadId];

        // residual form
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId], 0);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointVelocity[ThreadId]);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    virtual void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                            LocalSystemMatrixType& LHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,
                                            ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;

        RHS_Contribution.resize(LHS_Contribution.size1(), false);

        CalculateSystemContributions(
            pCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

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

    virtual void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                                      LocalSystemMatrixType& LHS_Contribution,
                                                      Condition::EquationIdVectorType& EquationId,
                                                      ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    virtual void GetElementalDofList(Element::Pointer rCurrentElement,
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
    ObjectiveFunction::Pointer mpObjectiveFunction;
    std::vector<LocalSystemVectorType> mAdjointVelocity;
    std::vector<LocalSystemVectorType> mAdjointAcceleration;
    std::vector<LocalSystemVectorType> mObjectiveGradient;
    std::vector<LocalSystemMatrixType> mAdjointMassMatrix;

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

}; /* Class Scheme */

///@}

///@name Type Definitions
///@{

///@}

///@} // Adjoint Fluid Application group

} /* namespace Kratos.*/

#endif /* KRATOS_ADJOINT_BOSSAK_SCHEME defined */
