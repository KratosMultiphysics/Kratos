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
#include "utilities/entities_utilities.h"

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
        const std::size_t num_threads = ParallelUtilities::GetNumThreads();

        mMatrix.M.resize(num_threads);
        mMatrix.D.resize(num_threads);
    }

    /**
     * @brief Constructor. The implicit method method
     * @param ThisParameters The configuration parameters
     */
    explicit ResidualBasedImplicitTimeScheme(Parameters ThisParameters)
        :ResidualBasedImplicitTimeScheme()
    {
        this->ValidateAndAssignParameters(ThisParameters);
        this->AssignSettings(ThisParameters);
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
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @details It "asks" the matrix needed to the element and performs the operations needed to introduce the selected time integration scheme. This function calculates at the same time the contribution to the LHS and to the RHS of the system
     * @param rCurrentElement The element to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Element& rCurrentElement,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateSystemContributions(rCurrentElement, rLHSContribution, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.CalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCurrentElement The element to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateRHSContribution(rCurrentElement, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.CalculateRHSContribution");
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition The condition to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateSystemContributions(rCurrentCondition, rLHSContribution, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.CalculateSystemContributions");
    }

    /**
     * @brief Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition The condition to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY;

        TCalculateRHSContribution(rCurrentCondition, rRHSContribution, rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.CalculateRHSContribution");
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

        const ProcessInfo r_current_process_info= rModelPart.GetProcessInfo();

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
    int Check(const ModelPart& rModelPart) const override
    {
        KRATOS_TRY;

        const int err = BaseType::Check(rModelPart);
        if(err!=0) return err;

        return 0;

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.Check");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "residualbased_implicit_time_scheme"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
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
     * @param rLHSContribution The dynamic contribution for the LHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void AddDynamicsToLHS(
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToLHS" << std::endl;
    }


    /**
     * @brief It adds the dynamic RHS contribution of the elements b - M*a - D*v
     * @param rCurrentElement The element to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void AddDynamicsToRHS(
        Element& rCurrentElement,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_ERROR << "YOU ARE CALLING THE BASE CLASS OF AddDynamicsToRHS" << std::endl;
    }

    /**
     * @brief It adds the dynamic RHS contribution of the condition RHS = fext - M*an0 - D*vn0 - K*dn0
     * @param rCurrentCondition The condition to compute
     * @param rRHSContribution The dynamic contribution for the RHS
     * @param D The damping matrix
     * @param M The mass matrix
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void AddDynamicsToRHS(
        Condition& rCurrentCondition,
        LocalSystemVectorType& rRHSContribution,
        LocalSystemMatrixType& D,
        LocalSystemMatrixType& M,
        const ProcessInfo& rCurrentProcessInfo
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

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @param rObject The object to compute
     * @param rLHSContribution The LHS matrix contribution
     * @param rRHSContribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    template <class TObjectType>
    void TCalculateSystemContributions(
        TObjectType& rObject,
        LocalSystemMatrixType& rLHSContribution,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rObject.CalculateLocalSystem(rLHSContribution,rRHSContribution,rCurrentProcessInfo);

        rObject.EquationIdVector(EquationId,rCurrentProcessInfo);

        rObject.CalculateMassMatrix(mMatrix.M[this_thread],rCurrentProcessInfo);

        rObject.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        AddDynamicsToLHS(rLHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        AddDynamicsToRHS(rObject, rRHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.TCalculateSystemContributions");
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rObject The object to compute
     * @param rRHSContribution The RHS vector contribution
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    template <class TObjectType>
    void TCalculateRHSContribution(
        TObjectType& rObject,
        LocalSystemVectorType& rRHSContribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY;

        const IndexType this_thread = OpenMPUtils::ThisThread();

        rObject.CalculateRightHandSide(rRHSContribution,rCurrentProcessInfo);

        rObject.CalculateMassMatrix(mMatrix.M[this_thread], rCurrentProcessInfo);

        rObject.CalculateDampingMatrix(mMatrix.D[this_thread],rCurrentProcessInfo);

        rObject.EquationIdVector(rEquationId,rCurrentProcessInfo);

        AddDynamicsToRHS(rObject, rRHSContribution, mMatrix.D[this_thread], mMatrix.M[this_thread], rCurrentProcessInfo);

        KRATOS_CATCH("ResidualBasedImplicitTimeScheme.TCalculateRHSContribution");
    }

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
