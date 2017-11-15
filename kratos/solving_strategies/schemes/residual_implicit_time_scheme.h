//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_RESIDUAL_BASED_IMPLICIT_TIME_SCHEME )
#define  KRATOS_RESIDUAL_BASED_IMPLICIT_TIME_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "utilities/openmp_utils.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"

namespace Kratos
{
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

/** @brief This is the base class for the implicit time schemes
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedImplicitTimeScheme
    : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedImplicitTimeScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                   NodesArrayType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * The implicit method method
     */
    ResidualBasedImplicitTimeScheme()
        :BaseType()
    {
        // Allocate auxiliary memory
        const unsigned int num_threads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(num_threads);
        mMatrix.D.resize(num_threads);
    }

    /** Copy Constructor.
     */
    ResidualBasedImplicitTimeScheme(ResidualBasedImplicitTimeScheme& rOther)
        :BaseType(rOther)
        ,mMatrix(rOther.mMatrix)
    {
    }

    /**
     * Clone
     */
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedImplicitTimeScheme(*this) );
    }

    /** Destructor.
     */
    ~ResidualBasedImplicitTimeScheme
    () override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * It initializes a non-linear iteration (for the element)
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        ProcessInfo& current_process_info = rModelPart.GetProcessInfo();
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i)
        {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->InitializeNonLinearIteration(current_process_info);
        }
        
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i)
        {
            auto it_elem = rModelPart.ConditionsBegin() + i;
            it_elem->InitializeNonLinearIteration(current_process_info);
        }     
        
        KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param pConditiont The condition to compute
     * @param CurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Condition::Pointer pCondition,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        pCondition->InitializeNonLinearIteration(CurrentProcessInfo);
    }
    
    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param pElement The element to compute
     * @param CurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Element::Pointer pElement,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        pElement->InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param pElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void CalculateSystemContributions(
        Element::Pointer pElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        //pElement->InitializeNonLinearIteration(CurrentProcessInfo);

        pElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        pElement->EquationIdVector(EquationId,CurrentProcessInfo);

        pElement->CalculateMassMatrix(mMatrix.M[thread],CurrentProcessInfo);

        pElement->CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        AddDynamicsToRHS(pElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Calculate_RHS_Contribution(
        Element::Pointer pElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        // pElement->InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the element considered
        pElement->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

        pElement->CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        pElement->CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        pElement->EquationIdVector(EquationId,CurrentProcessInfo);

        AddDynamicsToRHS (pElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Condition_CalculateSystemContributions(
        Condition::Pointer pCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCondition->InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the condition considered
        pCondition->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

        pCondition->EquationIdVector(EquationId,CurrentProcessInfo);

        pCondition->CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        pCondition->CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        AddDynamicsToRHS(pCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        // AssembleTimeSpaceLHS_Condition(pCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param CurrentProcessInfo The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& CurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const int thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCondition->InitializeNonLinearIteration(CurrentProcessInfo);

        // Basic operations for the condition considered
        pCondition->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

        pCondition->EquationIdVector(EquationId, CurrentProcessInfo);

        pCondition->CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

        pCondition->CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

        // Adding the dynamic contributions (static is already included)
        AddDynamicsToRHS(pCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */

    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = BaseType::Check(rModelPart);
        if(err!=0) return err;

        return 0;
        
        KRATOS_CATCH( "" );
    }


    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    struct  GeneralMatrices
    {
        std::vector<Matrix> M; // First derivative matrix  (usually mass matrix)
        std::vector<Matrix> D; // Second derivative matrix (usually damping matrix)
    };

    GeneralMatrices mMatrix; // This contains the auxiliar matrices

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * It adds the dynamic LHS contribution of the elements LHS = d(-RHS)/d(un0) = c0*c0*M + c0*D + K
     * @param LHS_Contribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToLHS" << std::endl;
    }


    /**
     * It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param rCurrentElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToRHS" << std::endl;
    }

    /**
     * It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param CurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& CurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToRHS" << std::endl;
    }

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@{

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

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
}; /* Class ResidualBasedImplicitTimeScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_IMPLICIT_TIME_SCHEME defined */
