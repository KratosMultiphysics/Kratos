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

/**
 * @class ResidualBasedImplicitTimeScheme
 * @ingroup KratosCore
 * @brief This is the base class for the implicit time schemes
 * @details Other implicit schemes should derive from this one. With the use of this base scheme it is possible to reduce code duplication
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @see Scheme
 * @author Vicente Mataix Ferrandiz
 */
template<class TSparseSpace,  class TDenseSpace >
class ResidualBasedImplicitTimeScheme
    : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{
    /// Pointer definition of ResidualBasedImplicitTimeScheme
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedImplicitTimeScheme );

    /// Base class definition
    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    /// DoF array type definition
    typedef typename BaseType::DofsArrayType                 DofsArrayType;
    /// DoF vector type definition
    typedef typename Element::DofsVectorType                DofsVectorType;

    /// Data type definition
    typedef typename BaseType::TDataType                         TDataType;
    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType         TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    /// Nodes containers definition
    typedef ModelPart::NodesContainerType                   NodesArrayType;
    /// Elements containers definition
    typedef ModelPart::ElementsContainerType             ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    /// Index type definition
    typedef std::size_t                                          IndexType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     * The implicit method method
     */
    explicit ResidualBasedImplicitTimeScheme()
        :BaseType()
    {
        // Allocate auxiliary memory
        const std::size_t num_threads = OpenMPUtils::GetNumThreads();

        mMatrix.M.resize(num_threads);
        mMatrix.D.resize(num_threads);
    }

    /** Copy Constructor.
     */
    explicit ResidualBasedImplicitTimeScheme(ResidualBasedImplicitTimeScheme& rOther)
        :BaseType(rOther)
        ,mMatrix(rOther.mMatrix)
    {
    }

    /**
     * Clone
     */
    typename BaseType::Pointer Clone() override
    {
        return Kratos::make_shared<ResidualBasedImplicitTimeScheme>(*this);
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
     * @param rModelPart The model part of the problem to solve
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

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->InitializeNonLinearIteration(r_current_process_info);
        }

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->InitializeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH( "" );
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @param pCurrentCondition The condition to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Condition::Pointer pCurrentCondition,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        pCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual element)
     * @param pCurrentElement The element to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Element::Pointer pCurrentElement,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        pCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @details It "asks" the matrix needed to the element and performs the operations needed to introduce the selected time integration scheme. This function calculates at the same time the contribution to the LHS and to the RHS of the system
     * @param pCurrentElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        //pCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        pCurrentElement->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,rCurrentProcessInfo);

        pCurrentElement->EquationIdVector(EquationId,rCurrentProcessInfo);

        pCurrentElement->CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);

        pCurrentElement->CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        AddDynamicsToLHS(LHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(pCurrentElement, RHS_Contribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.CalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param pCurrentElement The element to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void Calculate_RHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current element
        // pCurrentElement->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the element considered
        pCurrentElement->CalculateRightHandSide(rRHSContribution,rCurrentProcessInfo);

        pCurrentElement->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pCurrentElement->CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        pCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

        AddDynamicsToRHS (pCurrentElement, rRHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.Calculate_RHS_Contribution");
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition The condition to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the condition considered
        pCurrentCondition->CalculateLocalSystem(rLHSContribution,rRHSContribution, rCurrentProcessInfo);

        pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);

        pCurrentCondition->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pCurrentCondition->CalculateDampingMatrix(mMatrix.D[this_thread], rCurrentProcessInfo);

        AddDynamicsToLHS(rLHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(pCurrentCondition, rRHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        // AssembleTimeSpaceLHS_Condition(pCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.Condition_CalculateSystemContributions");
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param pCurrentCondition The condition to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        // Initializing the non linear iteration for the current condition
        //pCurrentCondition->InitializeNonLinearIteration(rCurrentProcessInfo);

        // Basic operations for the condition considered
        pCurrentCondition->CalculateRightHandSide(rRHSContribution, rCurrentProcessInfo);

        pCurrentCondition->EquationIdVector(rEquationId, rCurrentProcessInfo);

        pCurrentCondition->CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        pCurrentCondition->CalculateDampingMatrix(mMatrix.D[this_thread], rCurrentProcessInfo);

        // Adding the dynamic contributions (static is already included)
        AddDynamicsToRHS(pCurrentCondition, rRHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.Condition_Calculate_RHS_Contribution");
    }

    /**
     * @brief It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart The model part of the problem to solve
     * @param rA LHS matrix
     * @param rDx Incremental update of primary variables
     * @param rb RHS Vector
     */
    void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        ) override
    {
        KRATOS_TRY;

        ProcessInfo r_current_process_info= rModelPart.GetProcessInfo();

        BaseType::InitializeSolutionStep(rModelPart, rA, rDx, rb);

        const double delta_time = r_current_process_info[DELTA_TIME];

        KRATOS_ERROR_IF(delta_time < 1.0e-24) << "ERROR:: Detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current time step" << std::endl;

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.InitializeSolutionStep");
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model part of the problem to solve
     * @return Zero means  all ok
     */
    int Check(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const int err = BaseType::Check(rModelPart);
        if(err!=0) return err;

        return 0;

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.Check");
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

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedImplicitTimeScheme";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

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
