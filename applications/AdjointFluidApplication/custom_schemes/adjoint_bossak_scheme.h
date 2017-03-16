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
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/communicator.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/kratos_parameters.h"
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
 *
 * The sensitivity is computed on the boundary model part and is defined as:
 *
 * \f[
 * d_{\mathbf{s}}\bar{J} = \frac{1}{t_{end} - t_{start}}\Sigma_{n=1}^N
 *   (\partial_{\mathbf{s}}J^n + \lambda^{nT}\partial_{\mathbf{s}}\mathbf{r}^n)
*    \Delta t
 * \f]
 *
 * with
 *
 * \f[
 *  \mathbf{r}^n = \mathbf{f}^n - \mathbf{M}\dot{\mathbf{w}}^{n-\alpha}
 * \f]
 *
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
    AdjointBossakScheme(Parameters& rParameters, ObjectiveFunction::Pointer pObjectiveFunction)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "scheme_type": "bossak",
            "boundary_model_part_name": "PLEASE_SPECIFY_MODEL_PART",
            "alpha_bossak": -0.3,
            "adjoint_start_time": 0.0,
            "adjoint_end_time": 1.0
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mBoundaryModelPartName = rParameters["boundary_model_part_name"].GetString();
        mAlphaBossak = rParameters["alpha_bossak"].GetDouble();
        mGammaNewmark = 0.5 - mAlphaBossak;
        mInvGamma = 1.0 / mGammaNewmark;
        mInvGammaMinusOne = 1.0 / (mGammaNewmark - 1.0);
        mAdjointStartTime = rParameters["adjoint_start_time"].GetDouble();
        mAdjointEndTime = rParameters["adjoint_end_time"].GetDouble();

        if (mAdjointStartTime >= mAdjointEndTime)
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters: adjoint_start_time >= adjoint_end_time",
                rParameters.PrettyPrintJsonString())

        mpObjectiveFunction = pObjectiveFunction;

        // Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mAdjointVelocity.resize(NumThreads);
        mAdjointAcceleration.resize(NumThreads);
        mObjectiveGradient.resize(NumThreads);
        mAdjointMassMatrix.resize(NumThreads);

        KRATOS_CATCH("")
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

        // check domain dimension and element
        const unsigned int WorkingSpaceDimension =
            rModelPart.Elements().begin()->WorkingSpaceDimension();

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        const unsigned int DomainSize =
            static_cast<unsigned int>(rCurrentProcessInfo[DOMAIN_SIZE]);
        if (DomainSize != 2 && DomainSize != 3)
            KRATOS_THROW_ERROR(std::runtime_error, "invalid DOMAIN_SIZE: ", DomainSize)
        if (DomainSize != WorkingSpaceDimension)
            KRATOS_THROW_ERROR(
                std::runtime_error, "DOMAIN_SIZE != WorkingSpaceDimension", "")

        if (rModelPart.HasSubModelPart(mBoundaryModelPartName) == false)
        {
            KRATOS_THROW_ERROR(
                std::runtime_error,
                "invalid parameters \"boundary_model_part_name\": ",
                mBoundaryModelPartName)
        }

        // initialize the variables to zero.
        for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            it->Set(BOUNDARY, false);
            it->FastGetSolutionStepValue(SHAPE_SENSITIVITY) = SHAPE_SENSITIVITY.Zero();
            it->FastGetSolutionStepValue(ADJOINT_VELOCITY) = ADJOINT_VELOCITY.Zero();
            it->FastGetSolutionStepValue(ADJOINT_PRESSURE) = ADJOINT_PRESSURE.Zero();
            it->FastGetSolutionStepValue(ADJOINT_ACCELERATION) =
                ADJOINT_ACCELERATION.Zero();
        }

        ModelPart& rBoundaryModelPart = rModelPart.GetSubModelPart(mBoundaryModelPartName);
        for (auto it = rBoundaryModelPart.NodesBegin(); it != rBoundaryModelPart.NodesEnd(); ++it)
            it->Set(BOUNDARY, true);

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

        for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
            it->GetValue(NODAL_AREA) = 0.0; // todo: define application variable

        for (auto it = rModelPart.ElementsBegin(); it != rModelPart.ElementsEnd(); ++it)
            for (unsigned int iNode = 0; iNode < it->GetGeometry().PointsNumber(); ++iNode)
                it->GetGeometry()[iNode].GetValue(NODAL_AREA) += 1.0;

        rModelPart.GetCommunicator().AssembleNonHistoricalData(NODAL_AREA);

        mpObjectiveFunction->InitializeSolutionStep(rModelPart);

        KRATOS_CATCH("")
    }

    virtual void FinalizeSolutionStep(ModelPart& rModelPart,
                                      SystemMatrixType& rA,
                                      SystemVectorType& rDx,
                                      SystemVectorType& rb)
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        CalculateSolutionStepSensitivityContribution(rModelPart);

        mMass1Switch = 1.0;

        mpObjectiveFunction->FinalizeSolutionStep(rModelPart);

        KRATOS_CATCH("")
    }

    /// Update adjoint and adjoint acceleration.
    virtual void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb)
    {
        KRATOS_TRY

        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();
        const unsigned int DomainSize =
            static_cast<unsigned int>(rCurrentProcessInfo[DOMAIN_SIZE]);
        Communicator& rComm = rModelPart.GetCommunicator();

        if (rComm.TotalProcesses() == 1)
        {
            for (auto it = rDofSet.begin(); it != rDofSet.end(); ++it)
                if (it->IsFree() == true)
                    it->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(rDx, it->EquationId());

            for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
            {
                array_1d<double, 3>& rCurrentAdjointAcceleration =
                    it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 0);
                const array_1d<double, 3>& rOldAdjointAcceleration =
                    it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 1);
                for (unsigned int d = 0; d < DomainSize; ++d)
                    rCurrentAdjointAcceleration[d] =
                        (mGammaNewmark - 1.0) * mInvGamma * rOldAdjointAcceleration[d];
            }
        }
        else
        {
            for (auto it = rDofSet.begin(); it != rDofSet.end(); ++it)
                if (it->GetSolutionStepValue(PARTITION_INDEX) == rComm.MyPID())
                    if (it->IsFree() == true)
                        it->GetSolutionStepValue() +=
                            TSparseSpace::GetValue(rDx, it->EquationId());

            // todo: add a function Communicator::SynchronizeDofVariables() to
            // reduce communication here.
            rComm.SynchronizeNodalSolutionStepsData();

            for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
            {
                array_1d<double, 3>& rCurrentAdjointAcceleration =
                    it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 0);

                // in the end we need to assemble so we only compute this part
                // on the process that owns the node.
                if (it->FastGetSolutionStepValue(PARTITION_INDEX) == rComm.MyPID())
                {
                    const array_1d<double, 3>& rOldAdjointAcceleration =
                        it->FastGetSolutionStepValue(ADJOINT_ACCELERATION, 1);
                    for (unsigned int d = 0; d < DomainSize; ++d)
                        rCurrentAdjointAcceleration[d] = (mGammaNewmark - 1.0) * mInvGamma *
                                                         rOldAdjointAcceleration[d];
                }
                else
                {
                    for (unsigned int d = 0; d < DomainSize; ++d)
                        rCurrentAdjointAcceleration[d] = 0.0;
                }
            }
        }

        const int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector Partition;
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
                    *it, mAdjointMassMatrix[k], mObjectiveGradient[k], rCurrentProcessInfo);

                // old adjoint velocity
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 1);

                // terms depending on the old mass matrix
                mAdjointAcceleration[k] =
                    prod(mAdjointMassMatrix[k], mAdjointVelocity[k]) +
                    mObjectiveGradient[k];

                // transposed gradient of element residual w.r.t. acceleration
                it->Calculate(MASS_MATRIX_0, mAdjointMassMatrix[k], rCurrentProcessInfo);
                mAdjointMassMatrix[k] =
                    -(1.0 - mAlphaBossak) * trans(mAdjointMassMatrix[k]);

                // d (objective) / d (primal acceleration)
                mpObjectiveFunction->CalculateAdjointAccelerationContribution(
                    *it, mAdjointMassMatrix[k], mObjectiveGradient[k], rCurrentProcessInfo);

                // adjoint velocity
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 0);

                // terms depending on the mass matrix
                noalias(mAdjointAcceleration[k]) +=
                    prod(mAdjointMassMatrix[k], mAdjointVelocity[k]) +
                    mObjectiveGradient[k];

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

        rModelPart.GetCommunicator().AssembleCurrentData(ADJOINT_ACCELERATION);

        KRATOS_CATCH("")
    }

    /// Calculate residual based element contributions to transient adjoint.
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
            *pCurrentElement, mAdjointMassMatrix[ThreadId], mObjectiveGradient[ThreadId], rCurrentProcessInfo);

        // old adjoint velocity
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId], 1);

        // terms depending on the old mass matrix
        noalias(rRHS_Contribution) -=
            mMass1Switch * mInvGamma * mInvDt *
            (prod(mAdjointMassMatrix[ThreadId], mAdjointVelocity[ThreadId]) +
             mObjectiveGradient[ThreadId]);

        // transposed gradient of element residual w.r.t. acceleration
        pCurrentElement->Calculate(
            MASS_MATRIX_0, mAdjointMassMatrix[ThreadId], rCurrentProcessInfo);
        mAdjointMassMatrix[ThreadId] =
            -(1.0 - mAlphaBossak) * trans(mAdjointMassMatrix[ThreadId]);

        // d (objective) / d (primal acceleration)
        mpObjectiveFunction->CalculateAdjointAccelerationContribution(
            *pCurrentElement, mAdjointMassMatrix[ThreadId], mObjectiveGradient[ThreadId], rCurrentProcessInfo);
        noalias(rRHS_Contribution) -= mInvGamma * mInvDt * mObjectiveGradient[ThreadId];

        // transposed gradient of element residual w.r.t. primal
        pCurrentElement->Calculate(ADJOINT_MATRIX_2, rLHS_Contribution, rCurrentProcessInfo);

        // d (objective) / d (primal)
        mpObjectiveFunction->CalculateAdjointVelocityContribution(
            *pCurrentElement, rLHS_Contribution, mObjectiveGradient[ThreadId], rCurrentProcessInfo);
        noalias(rRHS_Contribution) -= mObjectiveGradient[ThreadId];

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

    std::string mBoundaryModelPartName;
    double mAlphaBossak;
    double mGammaNewmark;
    double mInvDt;
    double mInvGamma;
    double mInvGammaMinusOne;
    double mMass1Switch;
    double mAdjointStartTime;
    double mAdjointEndTime;
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

    void CalculateSolutionStepSensitivityContribution(ModelPart& rModelPart)
    {
        KRATOS_TRY

        ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
        const int NumThreads = OpenMPUtils::GetNumThreads();
        std::vector<Vector> CoordAuxVector(NumThreads);
        std::vector<Matrix> ShapeDerivativesMatrix(NumThreads);

        double DeltaTime = -rProcessInfo[DELTA_TIME]; // DELTA_TIME < 0
        const unsigned int DomainSize =
            static_cast<unsigned int>(rProcessInfo[DOMAIN_SIZE]);

        if (DeltaTime <= 0)
        {
            KRATOS_THROW_ERROR(std::runtime_error,
                               "detected for adjoint solution DELTA_TIME >= 0",
                               "")
        }

        double Weight = DeltaTime / (mAdjointEndTime - mAdjointStartTime);

        Communicator& rComm = rModelPart.GetCommunicator();
        if (rComm.TotalProcesses() > 1)
        {
            // here we make sure we only add the old shape sensitivity once
            // when we assemble.
            for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
                if (it->FastGetSolutionStepValue(PARTITION_INDEX) != rComm.MyPID())
                    it->FastGetSolutionStepValue(SHAPE_SENSITIVITY, 0) = SHAPE_SENSITIVITY.Zero();
        }

#pragma omp parallel
        {
            ModelPart::ElementIterator ElementsBegin;
            ModelPart::ElementIterator ElementsEnd;
            OpenMPUtils::PartitionedIterators(rModelPart.Elements(), ElementsBegin, ElementsEnd);
            int k = OpenMPUtils::ThisThread();

            for (auto it = ElementsBegin; it != ElementsEnd; ++it)
            {
                Element::GeometryType& rGeom = it->GetGeometry();
                bool IsBoundary = false;
                for (unsigned int iNode = 0; iNode < rGeom.PointsNumber(); ++iNode)
                    if (rGeom[iNode].Is(BOUNDARY) == true)
                    {
                        IsBoundary = true;
                        break;
                    }

                if (IsBoundary == false) // true for most elements
                    continue;

                // transposed gradient of local element's residual w.r.t. nodal
                // coordinates
                it->Calculate(SHAPE_DERIVATIVE_MATRIX_2, ShapeDerivativesMatrix[k], rProcessInfo);

                // d (objective) / d (coordinates)
                mpObjectiveFunction->CalculateSensitivityContribution(
                    *it, ShapeDerivativesMatrix[k], mObjectiveGradient[k], rProcessInfo);

                // adjoint solution
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 0);

                if (CoordAuxVector[k].size() != ShapeDerivativesMatrix[k].size1())
                    CoordAuxVector[k].resize(ShapeDerivativesMatrix[k].size1(), false);

                noalias(CoordAuxVector[k]) =
                    prod(ShapeDerivativesMatrix[k], mAdjointVelocity[k]) +
                    mObjectiveGradient[k];

                // Carefully write results to nodal variables
                unsigned int CoordIndex = 0;
                for (unsigned int iNode = 0; iNode < rGeom.PointsNumber(); ++iNode)
                {
                    if (rGeom[iNode].Is(BOUNDARY) == true)
                    {
                        rGeom[iNode].SetLock();
                        array_1d<double, 3>& rSensitivity =
                            rGeom[iNode].FastGetSolutionStepValue(SHAPE_SENSITIVITY);
                        for (unsigned int d = 0; d < DomainSize; ++d)
                            rSensitivity[d] += Weight * CoordAuxVector[k][CoordIndex++];
                        rGeom[iNode].UnSetLock();
                    }
                    else
                    {
                        // Skip this node block.
                        CoordIndex += DomainSize;
                    }
                }
            }
        }

        rModelPart.GetCommunicator().AssembleCurrentData(SHAPE_SENSITIVITY);

        KRATOS_CATCH("")
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

} /* namespace Kratos.*/

#endif /* KRATOS_ADJOINT_BOSSAK_SCHEME defined */
