//  KratosAdjointFluidApplication
//
//  License:		 BSD License
//					 license: AdjointFluidApplication/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_ADJOINT_STEADY_VELOCITY_PRESSURE_SCHEME)
#define KRATOS_ADJOINT_STEADY_VELOCITY_PRESSURE_SCHEME

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
#include "../../AdjointFluidApplication/custom_utilities/response_function.h"

namespace Kratos
{
///@addtogroup AdjointFluidApplication
///@{

///@name Kratos Classes
///@{

/// A scheme for steady adjoint equations involving velocities and pressures.
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
 * with response function \f$J=J(\mathbf{w};\mathbf{s})\f$.
 */
template <class TSparseSpace, class TDenseSpace>
class AdjointSteadyVelocityPressureScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointSteadyVelocityPressureScheme);

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
    AdjointSteadyVelocityPressureScheme(Parameters& rParameters, ResponseFunction::Pointer pResponseFunction)
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

        mpResponseFunction = pResponseFunction;

        // Allocate auxiliary memory
        int NumThreads = OpenMPUtils::GetNumThreads();
        mAdjointValues.resize(NumThreads);

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~AdjointSteadyVelocityPressureScheme() override
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

        mpResponseFunction->Initialize();

        KRATOS_CATCH("")
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
				SystemMatrixType& rA,
				SystemVectorType& rDx,
				SystemVectorType& rb) override
    {
        KRATOS_TRY

	// Sensitivities are generally computed as a time integral. For steady
	// problems, we set the time step to -1.0 (minus because adjoint is
	// backward in time).
	  rModelPart.GetProcessInfo()[DELTA_TIME] = -1.0;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // initialize the variables to zero.
        for (auto it = rModelPart.NodesBegin(); it != rModelPart.NodesEnd(); ++it)
        {
            it->Set(BOUNDARY, false);
            it->FastGetSolutionStepValue(ADJOINT_VELOCITY) = ADJOINT_VELOCITY.Zero();
            it->FastGetSolutionStepValue(ADJOINT_PRESSURE) = ADJOINT_PRESSURE.Zero();
            it->FastGetSolutionStepValue(ACCELERATION) = ACCELERATION.Zero();
        }

        ModelPart& rBoundaryModelPart = rModelPart.GetSubModelPart(mBoundaryModelPartName);
        for (auto it = rBoundaryModelPart.NodesBegin(); it != rBoundaryModelPart.NodesEnd(); ++it)
            it->Set(BOUNDARY, true);

        mpResponseFunction->InitializeSolutionStep();

        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                                      SystemMatrixType& rA,
                                      SystemVectorType& rDx,
                                      SystemVectorType& rb) override
    {
        KRATOS_TRY

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        mpResponseFunction->FinalizeSolutionStep();

        KRATOS_CATCH("")
    }

    /// Update adjoint.
    void Update(ModelPart& rModelPart,
                        DofsArrayType& rDofSet,
                        SystemMatrixType& rA,
                        SystemVectorType& rDx,
                        SystemVectorType& rb) override
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
    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                              LocalSystemMatrixType& rLHS_Contribution,
                                              LocalSystemVectorType& rRHS_Contribution,
                                              Element::EquationIdVectorType& rEquationId,
                                              ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        int ThreadId = OpenMPUtils::ThisThread();

        // adjoint system matrix
        pCurrentElement->CalculateFirstDerivativesLHS(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        // d (response) / d (primal)
        mpResponseFunction->CalculateFirstDerivativesGradient(
            *pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        noalias(rRHS_Contribution) = -rRHS_Contribution;

        // residual form
        pCurrentElement->GetValuesVector(mAdjointValues[ThreadId]);
        noalias(rRHS_Contribution) -= prod(rLHS_Contribution, mAdjointValues[ThreadId]);

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Calculate_LHS_Contribution(Element::Pointer pCurrentElement,
                                            LocalSystemMatrixType& LHS_Contribution,
                                            Element::EquationIdVectorType& EquationId,
                                            ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        LocalSystemVectorType RHS_Contribution;

        RHS_Contribution.resize(LHS_Contribution.size1(), false);

        CalculateSystemContributions(
            pCurrentElement, LHS_Contribution, RHS_Contribution, EquationId, CurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Condition::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    void Condition_Calculate_LHS_Contribution(Condition::Pointer pCurrentCondition,
                                                      LocalSystemMatrixType& LHS_Contribution,
                                                      Condition::EquationIdVectorType& EquationId,
                                                      ProcessInfo& CurrentProcessInfo) override
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    void GetElementalDofList(Element::Pointer rCurrentElement,
                                     Element::DofsVectorType& ElementalDofList,
                                     ProcessInfo& CurrentProcessInfo) override
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
    ResponseFunction::Pointer mpResponseFunction;
    std::vector<LocalSystemVectorType> mAdjointValues;

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

#endif /* KRATOS_ADJOINT_STEADY_VELOCITY_PRESSURE_SCHEME defined */
