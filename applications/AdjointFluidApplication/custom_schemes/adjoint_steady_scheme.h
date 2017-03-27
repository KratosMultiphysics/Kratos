//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ADJOINT_STEADY_SCHEME)
#define KRATOS_ADJOINT_STEADY_SCHEME

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
#include "containers/variable.h"

// Application includes
#include "../../AdjointFluidApplication/custom_utilities/objective_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

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
 *
 * The sensitivity is computed on the boundary model part and is defined as:
 *
 * \f[
 * d_{\mathbf{s}}J
 * = \partial_{\mathbf{s}}J + \lambda^T\partial_{\mathbf{s}}\mathbf{f}
 * \f]
 *
 */
template <class TSparseSpace, class TDenseSpace>
class AdjointSteadyScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointSteadyScheme);

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
    AdjointSteadyScheme(Parameters& rParameters, ObjectiveFunction::Pointer pObjectiveFunction)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        KRATOS_TRY

        Parameters DefaultParams(R"(
        {
            "scheme_type": "steady",
            "boundary_model_part_name": "PLEASE_SPECIFY_MODEL_PART"
        })");

        rParameters.ValidateAndAssignDefaults(DefaultParams);

        mBoundaryModelPartName = rParameters["boundary_model_part_name"].GetString();

        mpObjectiveFunction = pObjectiveFunction;

        // Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mAdjointVelocity.resize(NumThreads);

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~AdjointSteadyScheme()
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

        // initialize the variables to zero.
        for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            it->Set(BOUNDARY, false);
            it->FastGetSolutionStepValue(SHAPE_SENSITIVITY) = SHAPE_SENSITIVITY.Zero();
            it->FastGetSolutionStepValue(ADJOINT_VELOCITY) = ADJOINT_VELOCITY.Zero();
            it->FastGetSolutionStepValue(ADJOINT_PRESSURE) = ADJOINT_PRESSURE.Zero();
        }

        ModelPart& rBoundaryModelPart = rModelPart.GetSubModelPart(mBoundaryModelPartName);
        for (auto it = rBoundaryModelPart.NodesBegin(); it != rBoundaryModelPart.NodesEnd(); ++it)
            it->Set(BOUNDARY, true);

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

        mpObjectiveFunction->FinalizeSolutionStep(rModelPart);

        KRATOS_CATCH("")
    }

    /// Update adjoint.
    virtual void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb)
    {
        KRATOS_TRY

        Communicator& rComm = rModelPart.GetCommunicator();

        if (rComm.TotalProcesses() == 1)
        {
            for (auto it = rDofSet.begin(); it != rDofSet.end(); ++it)
                if (it->IsFree() == true)
                    it->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(rDx, it->EquationId());
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
        }

        KRATOS_CATCH("")
    }

    /// Calculate residual based element contributions to steady adjoint.
    virtual void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int ThreadId = OpenMPUtils::ThisThread();

        // adjoint system matrix
        pCurrentElement->Calculate(ADJOINT_MATRIX_1, rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        // d (objective) / d (primal)
        mpObjectiveFunction->CalculateAdjointVelocityContribution(
            *pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // residual form
        pCurrentElement->GetFirstDerivativesVector(mAdjointVelocity[ThreadId]);
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
    ObjectiveFunction::Pointer mpObjectiveFunction;
    std::vector<LocalSystemVectorType> mAdjointVelocity;

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

        const unsigned int DomainSize =
            static_cast<unsigned int>(rProcessInfo[DOMAIN_SIZE]);

        Communicator& rComm = rModelPart.GetCommunicator();
        if (rComm.TotalProcesses() > 1)
        {
            // here we make sure we only add the old shape sensitivity once
            // when we assemble. this should always be the case for steady
            // adjoint since we initialize shape sensitivity to zero in
            // InitializeSolutionStep()
            for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
                if (it->FastGetSolutionStepValue(PARTITION_INDEX) != rComm.MyPID())
                    it->FastGetSolutionStepValue(SHAPE_SENSITIVITY) =
                        SHAPE_SENSITIVITY.Zero();
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
                it->Calculate(SHAPE_DERIVATIVE_MATRIX_1, ShapeDerivativesMatrix[k], rProcessInfo);

                // d (objective) / d (coordinates)
                mpObjectiveFunction->CalculateSensitivityContribution(
                    *it, ShapeDerivativesMatrix[k], CoordAuxVector[k], rProcessInfo);

                // adjoint solution
                it->GetFirstDerivativesVector(mAdjointVelocity[k], 0);

                noalias(CoordAuxVector[k]) +=
                    prod(ShapeDerivativesMatrix[k], mAdjointVelocity[k]);

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
                            rSensitivity[d] += CoordAuxVector[k][CoordIndex++];
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

#endif /* KRATOS_ADJOINT_STEADY_SCHEME defined */
