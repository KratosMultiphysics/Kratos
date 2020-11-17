//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//


#if !defined(KRATOS_SCHEME )
#define  KRATOS_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "includes/kratos_parameters.h"

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
 * @class Scheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Riccardo Rossi
 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class Scheme
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Scheme
    KRATOS_CLASS_POINTER_DEFINITION(Scheme);

    /// The definition of the current class
    typedef Scheme< TSparseSpace, TDenseSpace > ClassType;

    /// Data type definition
    typedef typename TSparseSpace::DataType TDataType;
    /// Matrix type definition
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    /// Vector type definition
    typedef typename TSparseSpace::VectorType TSystemVectorType;
    /// Local system matrix type definition
    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    /// DoF type definition
    typedef Dof<double> TDofType;
    /// DoF array type definition
    typedef ModelPart::DofsArrayType DofsArrayType;
    /// DoF iterator type definition
    typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
    /// DoF constant iterator type definition
    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;

    /// Elements containers definition
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    /// Conditions containers definition
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    explicit Scheme()
    {
        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
        mConditionsAreInitialized = false;
    }
    /**
     * @brief Constructor with Parameters
     */
    explicit Scheme(Parameters ThisParameters)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
        mConditionsAreInitialized = false;
    }

    /** Copy Constructor.
     */
    explicit Scheme(Scheme& rOther)
      :mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      ,mElementsAreInitialized(rOther.mElementsAreInitialized)
      ,mConditionsAreInitialized(rOther.mConditionsAreInitialized)
    {
    }

    /** Destructor.
     */
    virtual ~Scheme()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param ThisParameters The configuration parameters
     */
    virtual typename ClassType::Pointer Create(Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned scheme
     */
    virtual Pointer Clone()
    {
        return Kratos::make_shared<Scheme>(*this) ;
    }

    /**
     * @brief This is the place to initialize the Scheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void Initialize(ModelPart& rModelPart)
    {
        KRATOS_TRY
        mSchemeIsInitialized = true;
        KRATOS_CATCH("")
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return True if initilized, false otherwise
     */
    bool SchemeIsInitialized()
    {
        return mSchemeIsInitialized;
    }

    /**
     * @brief This method sets if the elements have been initilized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetSchemeIsInitialized(bool SchemeIsInitializedFlag = true)
    {
        mSchemeIsInitialized = SchemeIsInitializedFlag;
    }

    /**
     * @brief This method returns if the elements are initialized
     * @return True if initilized, false otherwise
     */
    bool ElementsAreInitialized()
    {
        return mElementsAreInitialized;
    }

    /**
     * @brief This method sets if the elements have been initilized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetElementsAreInitialized(bool ElementsAreInitializedFlag = true)
    {
        mElementsAreInitialized = ElementsAreInitializedFlag;
    }

    /**
     * @brief This method returns if the conditions are initialized
     * @return True if initilized, false otherwise
     */
    bool ConditionsAreInitialized()
    {
        return mConditionsAreInitialized;
    }

    /**
     * @brief This method sets if the conditions have been initilized or not (true by default)
     * @param ConditionsAreInitializedFlag If the flag must be set to true or false
     */
    void SetConditionsAreInitialized(bool ConditionsAreInitializedFlag = true)
    {
        mConditionsAreInitialized = ConditionsAreInitializedFlag;
    }

    /**
     * @brief This is the place to initialize the elements.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeElements( ModelPart& rModelPart)
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Initialize(r_current_process_info);
        }

        SetElementsAreInitialized();

        KRATOS_CATCH("")
    }

    /**
     * @brief This is the place to initialize the conditions.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeConditions(ModelPart& rModelPart)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mElementsAreInitialized) << "Before initilizing Conditions, initialize Elements FIRST" << std::endl;

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->Initialize(r_current_process_info);
        }

        SetConditionsAreInitialized();

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the beginning of each solution step.
     * @details The basic operations to be carried in there are the following:
     * - managing variables to be kept constant over the time step (for example time-Scheme constants depending on the actual time step)
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void InitializeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        // Initializes solution step for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->InitializeSolutionStep(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Initializes solution step for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->InitializeSolutionStep(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        // Initializes solution step for all of the constraints
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->InitializeSolutionStep(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        // Finalizes solution step for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->FinalizeSolutionStep(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Finalizes solution step for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->FinalizeSolutionStep(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        // Finalizes solution step for all of the constraints
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->FinalizeSolutionStep(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /************************ BEGIN FRACTIONAL STEP METHODS ****************************/
    /********************* TODO: DECIDE IF NECESSARY TO DEFINE *************************/
    /***********************************************************************************/

//     /**
//      * @brief Initializes solution step, to be used when system is not explicitely defined
//      * @details For example for fractional step strategies
//      * @warning Must be defined in derived classes
//      * @param rModelPart The model part of the problem to solve
//      */
//     virtual void InitializeSolutionStep(ModelPart& rModelPart)
//     {
//         KRATOS_TRY
//         KRATOS_CATCH("")
//     }
//
//     /**
//      * @brief Finalizes solution step, to be used when system is not explicitely defined
//      * @details For example for fractional step strategies
//      * @warning Must be defined in derived classes
//      * @param rModelPart The model part of the problem to solve
//      */
//     virtual void FinalizeSolutionStep(ModelPart& rModelPart)
//     {
//         KRATOS_TRY
//         KRATOS_CATCH("")
//     }
//
//     /**
//      * @brief Executed before each fractional step
//      * @warning Must be defined in derived classes
//      * @param rModelPart The model part of the problem to solve
//      */
//     virtual void InitializeFractionalSolutionStep(ModelPart& rModelPart)
//     {
//         KRATOS_TRY
//         KRATOS_CATCH("")
//     }
//
//     /**
//      * @brief Executed after each fractional step
//      * @warning Must be defined in derived classes
//      * @param rModelPart The model part of the problem to solve
//      */
//     virtual void FinalizeFractionalSolutionStep(ModelPart& rModelPart)
//     {
//         KRATOS_TRY
//         KRATOS_CATCH("")
//     }

    /************************ END FRACTIONAL STEP METHODS ****************************/
    /***********************************************************************************/

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param rCurrentElement The element to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void InitializeNonLinearIteration(
        Element::Pointer rCurrentElement,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param rCurrentCondition The condition to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void InitializeNonLinearIteration(
        Condition::Pointer rCurrentCondition,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        // Finalizes non-linear iteration for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->FinalizeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Finalizes non-linear iteration  for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->FinalizeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        // Finalizes non-linear iteration for all of the constraints
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->FinalizeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the update of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions to be called to prepare the data needed for the output of results.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void CalculateOutputData(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        )
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Functions that cleans the results data.
     * @warning Must be implemented in the derived classes
     */
    virtual void CleanOutputData()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is intended to be called at the end of the solution step to clean up memory storage not needed after the end of the solution step
     * @warning Must be implemented in the derived classes
     */
    virtual void Clean()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Liberate internal storage.
     * @warning Must be implemented in the derived classes
     */
    virtual void Clear()
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @details Checks can be "expensive" as the function is designed
     * @param rModelPart The model part of the problem to solve
     * @return 0 all OK, 1 otherwise
     */
    virtual int Check(const ModelPart& rModelPart) const
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Checks for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.NumberOfElements()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            const auto& r_elem = *it_elem;
            r_elem.Check(r_current_process_info);
        }

        // Checks for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.NumberOfConditions()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            const auto& r_cond = *it_cond;
            r_cond.Check(r_current_process_info);
        }

        // Checks for all of the constraints
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.NumberOfMasterSlaveConstraints()); i++) {
            auto it_constraint = rModelPart.MasterSlaveConstraintsBegin() + i;
            const auto& r_constraint = *it_constraint;
            r_constraint.Check(r_current_process_info);
        }

        return 0;
        KRATOS_CATCH("");
    }

    virtual int Check(ModelPart& rModelPart)
    {
        // calling the const version for backward compatibility
        const Scheme& r_const_this = *this;
        const ModelPart& r_const_model_part = rModelPart;
        return r_const_this.Check(r_const_model_part);
    }

    /**
     * @brief This function is designed to be called in the builder and solver to introduce the selected time integration scheme.
     * @details It "asks" the matrix needed to the element and performs the operations needed to introduce the selected time integration scheme. This function calculates at the same time the contribution to the LHS and to the RHS of the system
     * @param rElement The element to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->CalculateSystemContributions(
            Element::Pointer(&rElement),
            LHS_Contribution,
            RHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rElement.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentElement->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateSystemContributions(
        Condition& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->Condition_CalculateSystemContributions(
            Condition::Pointer(&rCondition),
            LHS_Contribution,
            RHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentCondition->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->Calculate_RHS_Contribution(
            Element::Pointer(&rElement),
            RHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rElement.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Calculate_RHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentElement->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateRHSContribution(
        Condition& rCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->Condition_Calculate_RHS_Contribution(
            Condition::Pointer(&rCondition),
            RHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentCondition->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param rElement The element to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->Calculate_LHS_Contribution(
            Element::Pointer(&rElement),
            LHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rElement.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Calculate_LHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentElement->CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void CalculateLHSContribution(
        Condition& rCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        this->Condition_Calculate_LHS_Contribution(
            Condition::Pointer(&rCondition),
            LHS_Contribution,
            rEquationIdVector,
            const_cast<ProcessInfo&>(rCurrentProcessInfo)
        ); // TODO remove this after the transition period and uncomment the following
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
        // rCondition.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Condition_Calculate_LHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentCondition->CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief This method gets the eqaution id corresponding to the current element
     * @param rElement The element to compute
     * @param rEquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationId(
        const Element& rElement,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void EquationId(
        Element::Pointer pCurrentElement,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        (pCurrentElement)->EquationIdVector(EquationId, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param rEquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationId(
        const Condition& rCondition,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        rCondition.EquationIdVector(rEquationId, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void Condition_EquationId(
        Condition::Pointer pCurrentCondition,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        (pCurrentCondition)->EquationIdVector(EquationId, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given element
     * @param pCurrentElement The element to compute
     * @param rDofList The list containing the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        const Element& rElement,
        Element::DofsVectorType& rDofList,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        rElement.GetDofList(rDofList, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void GetElementalDofList(
        Element::Pointer pCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentElement->GetDofList(ElementalDofList, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given condition
     * @param rCondition The condition to compute
     * @param rDofList The list containing the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetDofList(
        const Condition& rCondition,
        Element::DofsVectorType& rDofList,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        rCondition.GetDofList(rDofList, rCurrentProcessInfo);
    }
    KRATOS_DEPRECATED_MESSAGE("This is legacy version, please use the other overload of this function")
    virtual void GetConditionDofList(
        Condition::Pointer pCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_START_IGNORING_DEPRECATED_FUNCTION_WARNING
        pCurrentCondition->GetDofList(ConditionDofList, rCurrentProcessInfo);
        KRATOS_STOP_IGNORING_DEPRECATED_FUNCTION_WARNING
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "name" : "scheme"
        })" );
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "scheme";
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
    virtual std::string Info() const
    {
        return "Scheme";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

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

    bool mSchemeIsInitialized;      /// Flag to be used in controlling if the Scheme has been intialized or not
    bool mElementsAreInitialized;   /// Flag taking in account if the elements were initialized correctly or not
    bool mConditionsAreInitialized; /// Flag taking in account if the conditions were initialized correctly or not

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    virtual Parameters ValidateAndAssignParameters(
        Parameters ThisParameters,
        const Parameters DefaultParameters
        ) const
    {
        ThisParameters.ValidateAndAssignDefaults(DefaultParameters);
        return ThisParameters;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    virtual void AssignSettings(const Parameters ThisParameters)
    {
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

    ///@}

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class Scheme

} // namespace Kratos.

#endif /* KRATOS_SCHEME  defined */


