//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "containers/csr_matrix.h"
#include "containers/system_vector.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "spaces/kratos_space.h"
#include "utilities/builtin_timer.h"
#include "utilities/dof_utilities/assembly_helper.h" //TODO: we should move this to somewhere else
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#endif

namespace Kratos::Future
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
 * @class NewScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @author Ruben Zorrilla
 */
//TODO: Think about the template parameters
template<class TMatrixType=CsrMatrix<>, class TVectorType=SystemVector<>>
class NewScheme
{
public:
    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", NewScheme, NewScheme, TMatrixType, TVectorType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", NewScheme, NewScheme, TMatrixType, TVectorType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NewScheme
    KRATOS_CLASS_POINTER_DEFINITION(NewScheme);

    /// TLS type definition
    struct ThreadLocalStorage
    {
        /// Dense space definition
        // using DenseSpaceType = UblasSpace<double, Matrix, Vector>; //TODO: We should eventually remove this and directly define the types
        using DenseSpaceType = TKratosSmpDenseSpace<double>;

        /// Local system matrix type definition
        using LocalSystemMatrixType = DenseSpaceType::MatrixType;

        /// Local system vector type definition
        using LocalSystemVectorType = DenseSpaceType::VectorType;

        // Local LHS contribution
        LocalSystemMatrixType LocalLhs;

        // Local RHS constribution
        LocalSystemVectorType LocalRhs;

        // Vector containing the localization in the system of the different terms
        Element::EquationIdVectorType LocalEqIds;
    };

    /// The definition of the current class
    using ClassType = NewScheme;

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = typename TMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TMatrixType::DataType;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Matrix pointer type definition
    using SystemMatrixPointerType = typename SystemMatrixType::Pointer;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Vector type definition
    using SystemVectorPointerType = typename SystemVectorType::Pointer;

    /// Dense space definition
    // using DenseSpaceType = UblasSpace<double, Matrix, Vector>;
    using DenseSpaceType = TKratosSmpDenseSpace<double>;

    /// Local system matrix type definition from TLS
    using LocalSystemMatrixType = typename ThreadLocalStorage::LocalSystemMatrixType;

    /// Local system vector type definition from TLS
    using LocalSystemVectorType = typename ThreadLocalStorage::LocalSystemVectorType;

    /// DoF type definition
    using DofType = Dof<double>;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// Elements containers definition
    using ElementsArrayType = ModelPart::ElementsContainerType;

    /// Conditions containers definition
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    /// Assembly helper type
    using AssemblyHelperType = AssemblyHelper<ThreadLocalStorage, TMatrixType, TVectorType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initializes the flags
     */
    explicit NewScheme() = default;

    /**
     * @brief Constructor with Parameters
     */
    explicit NewScheme(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate default parameters
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: User-definable reshaping stuff

        // Set up the assembly helper
        Parameters build_settings = ThisParameters["build_settings"];
        build_settings.AddInt("echo_level", ThisParameters["echo_level"].GetInt());
        mpAssemblyHelper = Kratos::make_unique<AssemblyHelperType>(rModelPart, build_settings);
    }

    /** Copy Constructor.
     */
    explicit NewScheme(NewScheme& rOther)
      : mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      , mElementsAreInitialized(rOther.mElementsAreInitialized)
      , mConditionsAreInitialized(rOther.mConditionsAreInitialized)
    {
        //TODO: Check this... particularly the mpAssemblyHelper pointer
    }

    /** Destructor.
     */
    virtual ~NewScheme() = default;

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
    virtual typename ClassType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned NewScheme
     */
    virtual Pointer Clone()
    {
        return Kratos::make_shared<NewScheme>(*this) ;
    }

    /**
     * @brief This is the place to initialize the NewScheme.
     * @details This is intended to be called just once when the strategy is initialized
     */
    virtual void Initialize()
    {
        KRATOS_TRY

        mSchemeIsInitialized = true;

        KRATOS_CATCH("")
    }

    /**
     * @brief This is the place to initialize the elements.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeElements()
    {
        KRATOS_TRY

        EntitiesUtilities::InitializeEntities<Element>(*mpModelPart);

        SetElementsAreInitialized(true);

        KRATOS_CATCH("")
    }

    /**
     * @brief This is the place to initialize the conditions.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    virtual void InitializeConditions()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mElementsAreInitialized) << "Before initializing Conditions, initialize Elements FIRST" << std::endl;

        EntitiesUtilities::InitializeEntities<Condition>(*mpModelPart);

        SetConditionsAreInitialized(true);

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
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Initializes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeSolutionStepAllEntities(*mpModelPart);

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
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

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
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(*mpModelPart);

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
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        // Finalizes non-linear iteration for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(*mpModelPart);

        KRATOS_CATCH("")
    }

    virtual void SetUpDofArray(DofsArrayType& rDofSet)
    {
        // Call the external utility to set up the DOFs array
        DofArrayUtilities::SetUpDofArray(*mpModelPart, rDofSet, mEchoLevel);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished DOFs array set up." << std::endl;
    }

    SizeType SetUpSystem(DofsArrayType& rDofSet)
    {
        KRATOS_ERROR_IF(rDofSet.empty()) << "DOFs set is empty. Call the 'SetUpDofArray' first." << std::endl;

        return mpAssemblyHelper->SetUpSystemIds(rDofSet);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished system set up." << std::endl;
    }

    void ResizeAndInitializeVectors(
        const DofsArrayType& rDofSet,
        SystemMatrixPointerType& rpLHS,
        SystemVectorPointerType& rpDx,
        SystemVectorPointerType& rpRHS,
        const bool CalculateReactions = false)
    {
        KRATOS_TRY

        // Call the assembly helper to allocate and initialize the required vectors
        // Note that this also allocates the required reaction vectors (e.g., elimination build)
        mpAssemblyHelper->ResizeAndInitializeVectors(rDofSet, rpLHS, rpDx, rpRHS, CalculateReactions);

        //TODO: Think on the constraints stuff!
        // ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished system initialization." << std::endl;

        KRATOS_CATCH("")
    }

    virtual void Build(
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate local LHS and RHS contributions
            ItElem->CalculateLocalSystem(rTLS.LocalLhs, rTLS.LocalRhs, rProcessInfo);
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate local LHS and RHS contributions
            ItCond->CalculateLocalSystem(rTLS.LocalLhs, rTLS.LocalRhs, rProcessInfo);
        };

        ThreadLocalStorage aux_tls;
        auto& r_assembly_helper = GetAssemblyHelper();
        r_assembly_helper.SetElementAssemblyFunction(elem_func);
        r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        r_assembly_helper.AssembleLocalSystem(rLHS, rRHS, aux_tls);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    virtual void Build(SystemVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate the RHS contributions
            ItElem->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate the RHS contributions
            ItCond->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);
        };

        ThreadLocalStorage aux_tls;
        auto& r_assembly_helper = GetAssemblyHelper();
        r_assembly_helper.SetElementAssemblyFunction(elem_func);
        r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        r_assembly_helper.AssembleRightHandSide(rRHS, aux_tls);

        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("NewScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        SystemMatrixType& rLHS,
        SystemVectorType& rRHS)
    {
        GetAssemblyHelper().ApplyDirichletConditions(rDofSet, rLHS, rRHS);
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        SystemVectorType& rRHS)
    {
        GetAssemblyHelper().ApplyDirichletConditions(rDofSet, rRHS);
    }

    virtual void CalculateReactions(
        const DofsArrayType& rDofSet,
        SystemVectorType& rRHS)
    {
        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate the RHS contributions
            ItElem->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, ThreadLocalStorage& rTLS){
            // Calculate the RHS contributions
            ItCond->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);
        };

        ThreadLocalStorage aux_tls;
        auto& r_assembly_helper = GetAssemblyHelper();
        r_assembly_helper.SetElementAssemblyFunction(elem_func);
        r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        r_assembly_helper.CalculateReactionsRightHandSide(rDofSet, rRHS, aux_tls);
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void Predict(
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
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
        DofsArrayType& rDofSet,
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
    {
        KRATOS_TRY

        //TODO: To be discussed. I think that te base scheme should do what the current IncrementalUpdateStaticScheme does that is the following
        block_for_each(rDofSet, [&Dx](DofType& rDof){
            if (rDof.IsFree()) {
                rDof.GetSolutionStepValue() += Dx[rDof.EquationId()];
            }
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Functions to be called to prepare the data needed for the output of results.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void CalculateOutputData(
        SystemMatrixType& A,
        SystemVectorType& Dx,
        SystemVectorType& b)
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

        // Reset initialization flags
        mSchemeIsInitialized = false;
        mElementsAreInitialized = false;
        mConditionsAreInitialized = false;

        // Clear the assembly helper
        GetAssemblyHelper().Clear();

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
    virtual int Check() const
    {
        KRATOS_TRY

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"({
            "name" : "new_scheme",
            "build_settings" : {
                "build_type" : "block",
                "scaling_type" : "max_diagonal"
            },
            "echo_level" : 0,
            "calculate_reactions" : false
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "new_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method sets the value of mEchoLevel
     * @param EchoLevel The value to set
     */
    void SetEchoLevel(const bool EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * @brief This method sets the flag mReshapeMatrixFlag
     * @param ReshapeMatrixFlag The flag that tells if we need to reshape the LHS matrix
     */
    void SetReshapeMatrixFlag(const bool ReshapeMatrixFlag)
    {
        mReshapeMatrixFlag = ReshapeMatrixFlag;
    }

    /**
     * @brief This method sets if the elements have been initialized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetSchemeIsInitialized(const bool SchemeIsInitializedFlag)
    {
        mSchemeIsInitialized = SchemeIsInitializedFlag;
    }

    /**
     * @brief This method sets if the elements have been initialized or not (true by default)
     * @param ElementsAreInitializedFlag If the flag must be set to true or false
     */
    void SetElementsAreInitialized(const bool ElementsAreInitializedFlag)
    {
        mElementsAreInitialized = ElementsAreInitializedFlag;
    }

    /**
     * @brief This method sets if the conditions have been initialized or not (true by default)
     * @param ConditionsAreInitializedFlag If the flag must be set to true or false
     */
    void SetConditionsAreInitialized(const bool ConditionsAreInitializedFlag)
    {
        mConditionsAreInitialized = ConditionsAreInitializedFlag;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method returns the flag mReshapeMatrixFlag
     * @return The flag that tells if we need to reshape the LHS matrix
     */
    bool GetReshapeMatrixFlag() const
    {
        return mReshapeMatrixFlag;
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return True if initialized, false otherwise
     */
    bool SchemeIsInitialized() const
    {
        return mSchemeIsInitialized;
    }

    /**
     * @brief This method returns if the elements are initialized
     * @return True if initialized, false otherwise
     */
    bool ElementsAreInitialized() const
    {
        return mElementsAreInitialized;
    }

    /**
     * @brief This method returns if the conditions are initialized
     * @return True if initialized, false otherwise
     */
    bool ConditionsAreInitialized() const
    {
        return mConditionsAreInitialized;
    }

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

    bool mSchemeIsInitialized = false;      /// Flag to be used in controlling if the Scheme has been initialized or not

    bool mElementsAreInitialized = false;   /// Flag taking in account if the elements were initialized correctly or not

    bool mConditionsAreInitialized = false; /// Flag taking in account if the conditions were initialized correctly or not

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
        const Parameters DefaultParameters) const
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
        // Set verbosity level
        mEchoLevel = ThisParameters["echo_level"].GetInt();
    }

    ///@}
    ///@name Protected  Access
    ///@{

    AssemblyHelperType& GetAssemblyHelper()
    {
        return *mpAssemblyHelper;
    }

    AssemblyHelperType& GetAssemblyHelper() const
    {
        return *mpAssemblyHelper;
    }

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

    ModelPart* mpModelPart = nullptr;

    bool mReshapeMatrixFlag = false; /// If the matrix is reshaped each step

    int mEchoLevel = 0;

    typename AssemblyHelperType::UniquePointer mpAssemblyHelper = nullptr;

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

} // namespace Kratos::Future.
