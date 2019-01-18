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

    /**
     * @class LocalSystemComponents
     * @brief This struct is used in the component wise calculation only is defined here and is used to declare a member variable in the component wise schemes private pointers can only be accessed by means of set and get functions
     * @details This allows to set and not copy the Element_Variables and Condition_Variables which will be asked and set by another strategy object
     */
    struct LocalSystemComponents
    {
    private:
        ///@name Member Variables
        ///@{
        // Elements
        std::vector<LocalSystemMatrixType> *mpLHS_Element_Components;
        const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Element_Variables;

        std::vector<LocalSystemVectorType> *mpRHS_Element_Components;
        const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Element_Variables;

        // Conditions
        std::vector<LocalSystemMatrixType> *mpLHS_Condition_Components;
        const std::vector< Variable< LocalSystemMatrixType > > *mpLHS_Condition_Variables;

        std::vector<LocalSystemVectorType> *mpRHS_Condition_Components;
        const std::vector< Variable< LocalSystemVectorType > > *mpRHS_Condition_Variables;
        ///@}
    public:
        ///@name Operations
        ///@{
        /**
        * @brief This method initializes the pointer of the member variables
        */
        void Initialize()
        {
            mpLHS_Element_Components = NULL;
            mpLHS_Element_Variables  = NULL;

            mpRHS_Element_Components = NULL;
            mpRHS_Element_Variables  = NULL;

            mpLHS_Condition_Components = NULL;
            mpLHS_Condition_Variables  = NULL;

            mpRHS_Condition_Components = NULL;
            mpRHS_Condition_Variables  = NULL;
        }

        /* Setting pointer variables */

        // Elements
        void SetLHS_Element_Components ( std::vector<LocalSystemMatrixType>& rLHS_Element_Components ) { mpLHS_Element_Components = &rLHS_Element_Components; };
        void SetLHS_Element_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Element_Variables ) { mpLHS_Element_Variables = &rLHS_Element_Variables; };
        void SetRHS_Element_Components ( std::vector<LocalSystemVectorType>& rRHS_Element_Components ) { mpRHS_Element_Components = &rRHS_Element_Components; };
        void SetRHS_Element_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Element_Variables ) { mpRHS_Element_Variables = &rRHS_Element_Variables; };

        bool Are_LHS_Element_Components_Set() { if( mpLHS_Element_Variables == NULL ) return false; else return true; };
        bool Are_RHS_Element_Components_Set() { if( mpRHS_Element_Variables == NULL ) return false; else return true; };

        // Conditions
        void SetLHS_Condition_Components ( std::vector<LocalSystemMatrixType>& rLHS_Condition_Components ) { mpLHS_Condition_Components = &rLHS_Condition_Components; };
        void SetLHS_Condition_Variables     ( const std::vector< Variable< LocalSystemMatrixType > >& rLHS_Condition_Variables ) { mpLHS_Condition_Variables = &rLHS_Condition_Variables; };
        void SetRHS_Condition_Components ( std::vector<LocalSystemVectorType>& rRHS_Condition_Components ) { mpRHS_Condition_Components = &rRHS_Condition_Components; };
        void SetRHS_Condition_Variables     ( const std::vector< Variable< LocalSystemVectorType > >& rRHS_Condition_Variables ) { mpRHS_Condition_Variables = &rRHS_Condition_Variables; };

        bool Are_LHS_Condition_Components_Set() { if( mpLHS_Condition_Variables == NULL ) return false; else return true; };
        bool Are_RHS_Condition_Components_Set() { if( mpRHS_Condition_Variables == NULL ) return false; else return true; };

        /* Getting pointer variables */

        // Elements
        std::vector<LocalSystemMatrixType>& GetLHS_Element_Components() { return *mpLHS_Element_Components; };
        const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Element_Variables() { return *mpLHS_Element_Variables; };
        std::vector<LocalSystemVectorType>& GetRHS_Element_Components() { return *mpRHS_Element_Components; };
        const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Element_Variables() { return *mpRHS_Element_Variables; };

        // Conditions
        std::vector<LocalSystemMatrixType>& GetLHS_Condition_Components() { return *mpLHS_Condition_Components; };
        const std::vector< Variable< LocalSystemMatrixType > >& GetLHS_Condition_Variables() { return *mpLHS_Condition_Variables; };
        std::vector<LocalSystemVectorType>& GetRHS_Condition_Components() { return *mpRHS_Condition_Components; };
        const std::vector< Variable< LocalSystemVectorType > >& GetRHS_Condition_Variables() { return *mpRHS_Condition_Variables; };

        ///@}
    };

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

    /**
     * @brief Clone method
     * @return The pointer of the cloned scheme
     */
    virtual Pointer Clone()
    {
        return Kratos::make_shared<Scheme>(*this) ;
    }

    /**
     * @brief Component wise components Get method
     * @warning Must be defined on the derived classes
     * @return The local system of components
     */
    virtual LocalSystemComponents& GetLocalSystemComponents()
    {
        KRATOS_ERROR << "Asking for Local Components to the SCHEME base class which is not component wise and not contains this member variable" << std::endl;
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

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Initialize();
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

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->Initialize();
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

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Initialize solution step for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->InitializeSolutionStep(r_current_process_info);
        }

        // Initialize solution step for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->InitializeSolutionStep(r_current_process_info);
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

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Finalizes solution step for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->FinalizeSolutionStep(r_current_process_info);
        }

        // Finalizes solution step for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->FinalizeSolutionStep(r_current_process_info);
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

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Finalizes non-linear iteration for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->FinalizeNonLinearIteration(r_current_process_info);
        }

        // Finalizes non-linear iteration  for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->FinalizeNonLinearIteration(r_current_process_info);
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
     * @brief Function to clean up "element" scratch space after each element is built.
     * @param rCurrentElement The element to compute
     */
    virtual void CleanMemory(Element::Pointer rCurrentElement)
    {
        rCurrentElement->CleanMemory();
    }

    /**
     * @brief Function to clean up "condition" scratch space after each condition is built.
     * @param rCurrentCondition The condition to compute
     */
    virtual void CleanMemory(Condition::Pointer rCurrentCondition)
    {
        rCurrentCondition->CleanMemory();
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
    virtual int Check(ModelPart& rModelPart)
    {
        KRATOS_TRY

        ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Checks for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); i++) {
            auto it_elem = rModelPart.ElementsBegin() + i;
            it_elem->Check(r_current_process_info);
        }

        // Checks for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); i++) {
            auto it_cond = rModelPart.ConditionsBegin() + i;
            it_cond->Check(r_current_process_info);
        }

        return 0;
        KRATOS_CATCH("");
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
    virtual void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentElement->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition The condition to compute
     * @param LHS_Contribution The LHS matrix contribution
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentCondition->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param pCurrentElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Calculate_RHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentElement->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentCondition->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param pCurrentElement The element to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Calculate_LHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentElement->CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition The condition to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Condition_Calculate_LHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentCondition->CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    /**
     * @brief This method gets the eqaution id corresponding to the current element
     * @param pCurrentElement The element to compute
     * @param EquationId The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void EquationId(
        Element::Pointer pCurrentElement,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        (pCurrentElement)->EquationIdVector(EquationId, rCurrentProcessInfo);
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition The condition to compute
     * @param EquationId The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void Condition_EquationId(
        Condition::Pointer pCurrentCondition,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        (pCurrentCondition)->EquationIdVector(EquationId, rCurrentProcessInfo);
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given element
     * @param pCurrentElement The element to compute
     * @param ElementalDofList The list containing the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetElementalDofList(
        Element::Pointer pCurrentElement,
        Element::DofsVectorType& ElementalDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentElement->GetDofList(ElementalDofList, rCurrentProcessInfo);
    }

    /**
     * @brief Function that returns the list of Degrees of freedom to be assembled in the system for a Given condition
     * @param pCurrentCondition The condition to compute
     * @param ConditionDofList The list containing the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    virtual void GetConditionDofList(
        Condition::Pointer pCurrentCondition,
        Element::DofsVectorType& ConditionDofList,
        ProcessInfo& rCurrentProcessInfo
        )
    {
        pCurrentCondition->GetDofList(ConditionDofList, rCurrentProcessInfo);
    }

    ///@}
    ///@name Operations
    ///@{

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


