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
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

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
     * @brief It initializes a non-linear iteration (for the element)
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
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->InitializeNonLinearIteration(current_process_info);
        }
        
        
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_elem = rModelPart.ConditionsBegin() + i;
            it_elem->InitializeNonLinearIteration(current_process_info);
        }     
        
        KRATOS_CATCH( "" );
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @param pCondition The condition to compute
     * @param rCurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Condition::Pointer pCondition,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        pCondition->InitializeNonLinearIteration(rCurrentProcessInfo);
    }
    
    /**
     * @brief It initializes a non-linear iteration (for an individual element)
     * @param pElement The element to compute
     * @param rCurrentProcessInfo The current process info instance
     */

    void InitializeNonLinearIteration(
        Element::Pointer pElement,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        pElement->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce
     * @param pElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */

    void CalculateSystemContributions(
        Element::Pointer pElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const std::size_t this_thread = OpenMPUtils::ThisThread();

        //pElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        pElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        pElement->EquationIdVector(EquationId,rCurrentProcessInfo);

        pElement->CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);

        pElement->CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(pElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param pElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */

    void Calculate_RHS_Contribution(
        Element::Pointer pElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        // pElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the element considered
        pElement->CalculateRightHandSide(RHS_Contribution,rCurrentProcessInfo);

        pElement->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pElement->CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        pElement->EquationIdVector(EquationId,rCurrentProcessInfo);

        AddDynamicsToRHS (pElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */

    void Condition_CalculateSystemContributions(
        Condition::Pointer pCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the condition considered
        pCondition->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        pCondition->EquationIdVector(EquationId,rCurrentProcessInfo);

        pCondition->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pCondition->CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(pCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        // AssembleTimeSpaceLHS_Condition(pCondition, LHS_Contribution,DampMatrix, MassMatrix,rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param pCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const std::size_t this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the condition considered
        pCondition->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);

        pCondition->EquationIdVector(EquationId, rCurrentProcessInfo);

        pCondition->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pCondition->CalculateDampingMatrix(mMatrix.D[this_thread], rCurrentProcessInfo);

        // Adding the dynamic contributions (static is already included)
        AddDynamicsToRHS(pCondition, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH( "" );
    }
    
    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        ProcessInfo current_process_info= rModelPart.GetProcessInfo();

        BaseType::InitializeSolutionStep(rModelPart, A, Dx, b);

        const double delta_time = current_process_info[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time < 1.0e-24) << "ERROR:: Detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current time step" << std::endl;
        
        KRATOS_CATCH( "" );
    }
    
    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. 
     * @details Checks can be "expensive" as the function is designed
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

    struct GeneralMatrices
    {
        std::vector< Matrix > M; /// First derivative matrix  (usually mass matrix)
        std::vector< Matrix > D; /// Second derivative matrix (usually damping matrix)
    };
    
    GeneralMatrices mMatrix;
    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief It adds the dynamic LHS contribution of the elements LHS = d(-RHS)/d(un0) = c0*c0*M + c0*D + K
     * @param LHS_Contribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToLHS(
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToLHS" << std::endl;
    }


    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param rCurrentElement The element to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToRHS(
        Element::Pointer rCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToRHS" << std::endl;
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param rCurrentCondition The condition to compute
     * @param RHS_Contribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */

    virtual void AddDynamicsToRHS(
        Condition::Pointer rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        ProcessInfo& rCurrentProcessInfo
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
