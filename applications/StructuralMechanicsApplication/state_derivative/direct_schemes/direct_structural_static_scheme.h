// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Kevin Braun, https://github.com/MFusseder
//

#if !defined(KRATOS_DIRECT_STRUCTURAL_STATIC_SCHEME)
#define KRATOS_DIRECT_STRUCTURAL_STATIC_SCHEME

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/openmp_utils.h"
#include "solving_strategies/schemes/scheme.h"
#include "containers/variable.h"

// Application includes
#include "state_derivative/variable_utilities/direct_sensitivity_variable.h"

namespace Kratos
{
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/// A scheme for for direct sensitivity analysis.
/**
 *
 *
 */
template <class TSparseSpace, class TDenseSpace>
class DirectStructuralStaticScheme : public Scheme<TSparseSpace, TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(DirectStructuralStaticScheme);

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
    DirectStructuralStaticScheme(Parameters rParameters, DirectSensitivityVariable::Pointer pVariable)
        : Scheme<TSparseSpace, TDenseSpace>()
    {
        KRATOS_TRY;

        Parameters default_params(R"(
        {
            "scheme_type": "direct_structural",
            "rotation_dofs": false
        })");

        rParameters.ValidateAndAssignDefaults(default_params);

        mpVariable = pVariable;

        mHasRotationDofs = rParameters["rotation_dofs"].GetBool();

        KRATOS_CATCH("");        
    }

    /// Destructor.
    ~DirectStructuralStaticScheme() override
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

        // Allocate auxiliary memory
        int num_threads = OpenMPUtils::GetNumThreads();
        mAdjointValues.resize(num_threads);

        mpVariable->Initialize();

        KRATOS_CATCH("");
    }

    void InitializeSolutionStep(ModelPart& rModelPart,
                                SystemMatrixType& rA,
                                SystemVectorType& rDx,
                                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        // Sensitivities are generally computed as a time integral. For steady
        // problems, we set the time step to -1.0 (minus because adjoint is
        // backward in time).
        rModelPart.GetProcessInfo()[DELTA_TIME] = -1.0;

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        // initialize the variables to zero.
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rModelPart.Nodes().size()); ++k)
        {
            auto it_node = rModelPart.NodesBegin() + k;
            noalias(it_node->FastGetSolutionStepValue(ADJOINT_DISPLACEMENT)) = ADJOINT_DISPLACEMENT.Zero();
        }

        if(mHasRotationDofs)
        {
            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (rModelPart.Nodes().size()); ++k)
            {
                auto it_node = rModelPart.NodesBegin() + k;
                noalias(it_node->FastGetSolutionStepValue(ADJOINT_ROTATION)) = ADJOINT_ROTATION.Zero();
            }
        }

        mpVariable->InitializeSolutionStep();

        KRATOS_CATCH("");
    }

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              SystemMatrixType& rA,
                              SystemVectorType& rDx,
                              SystemVectorType& rb) override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep(rModelPart, rA, rDx, rb);

        mpVariable->FinalizeSolutionStep();

        KRATOS_CATCH("");
    }

    /// Update adjoint.
    void Update(ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                SystemMatrixType& rA,
                SystemVectorType& rDx,
                SystemVectorType& rb) override
    {
        KRATOS_TRY;

        Communicator& r_comm = rModelPart.GetCommunicator();

        if (r_comm.TotalProcesses() == 1)
        {
            int ndofs = static_cast<int>(rDofSet.size());
            #pragma omp parallel for
            for (int i = 0; i < ndofs; ++i)
            {
                typename DofsArrayType::iterator it = rDofSet.begin() + i;
                if (it->IsFree() == true)
                    it->GetSolutionStepValue() +=
                        TSparseSpace::GetValue(rDx, it->EquationId());
            }
        }
        else
        {
            int ndofs = static_cast<int>(rDofSet.size());
            #pragma omp parallel for
            for (int i = 0; i < ndofs; ++i)
            {
                typename DofsArrayType::iterator it = rDofSet.begin() + i;
                if (it->GetSolutionStepValue(PARTITION_INDEX) == r_comm.MyPID())
                    if (it->IsFree() == true)
                        it->GetSolutionStepValue() +=
                            TSparseSpace::GetValue(rDx, it->EquationId());
            }

            // todo: add a function Communicator::SynchronizeDofVariables() to
            // reduce communication here.
            r_comm.SynchronizeNodalSolutionStepsData();
        }

        KRATOS_CATCH("");
    }

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        // check domain dimension and element
        const unsigned int working_space_dimension =
            rModelPart.Elements().begin()->GetGeometry().WorkingSpaceDimension();

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int domain_size =
            static_cast<unsigned int>(r_current_process_info[DOMAIN_SIZE]);
        KRATOS_ERROR_IF(domain_size != 2 && domain_size != 3) << "Invalid DOMAIN_SIZE: " << domain_size << std::endl;
        KRATOS_ERROR_IF(domain_size != working_space_dimension) << "DOMAIN_SIZE != WorkingSpaceDimension()" << std::endl;

        for(auto& rnode : rModelPart.Nodes())
           KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_DISPLACEMENT, rnode)
        if(mHasRotationDofs)
        {
            for(auto& rnode : rModelPart.Nodes())
                KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_ROTATION, rnode)
        }

        return BaseType::Check(rModelPart); // check elements and conditions
        KRATOS_CATCH("");
    }

    /// Calculate residual based element contributions to steady adjoint.
    void CalculateSystemContributions(Element::Pointer pCurrentElement,
                                      LocalSystemMatrixType& rLHS_Contribution,
                                      LocalSystemVectorType& rRHS_Contribution,
                                      Element::EquationIdVectorType& rEquationId,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY;

        int thread_id = OpenMPUtils::ThisThread();

        // Get element stiffness matrix
        pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);

        if (rRHS_Contribution.size() != rLHS_Contribution.size1())
            rRHS_Contribution.resize(rLHS_Contribution.size1(), false);

        // Calculate pseudo load        
        mpVariable->CalculatePseudoLoadVector(*pCurrentElement, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

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

        LocalSystemVectorType RHS_contribution;

        RHS_contribution.resize(rLHS_Contribution.size1(), false);

        CalculateSystemContributions(
            pCurrentElement, rLHS_Contribution, RHS_contribution, rEquationId, rCurrentProcessInfo);

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

        // Calculate pseudo load.
        pCurrentCondition->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        
        mpVariable->CalculatePseudoLoadVector(*pCurrentCondition, rLHS_Contribution, rRHS_Contribution, rCurrentProcessInfo);

        if( rLHS_Contribution.size1() == 0 || rLHS_Contribution.size2() == 0 )
        {
            rLHS_Contribution.resize(rRHS_Contribution.size(),rRHS_Contribution.size());
            for(IndexType i = 0; i < rLHS_Contribution.size1(); ++i)
                for(IndexType j = 0; j < rLHS_Contribution.size1(); ++j)
                    rLHS_Contribution(i,j)=0;
        } 

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
        LocalSystemVectorType RHS_contribution;

        RHS_contribution.resize(rLHS_Contribution.size1(), false);

        Condition_CalculateSystemContributions(
             pCurrentCondition, rLHS_Contribution, RHS_contribution, rEquationId, rCurrentProcessInfo);

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

    DirectSensitivityVariable::Pointer mpVariable;
    std::vector<LocalSystemVectorType> mAdjointValues;
    bool mHasRotationDofs;

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

///@} // Structural Mechanics Application group

} /* namespace Kratos.*/

#endif /* KRATOS_DIRECT_STRUCTURAL_STATIC_SCHEME defined */
