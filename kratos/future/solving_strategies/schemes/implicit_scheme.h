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
//                   Riccardo Rossi
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
#include "utilities/dof_utilities/dof_array_utilities.h"
#include "utilities/entities_utilities.h"
#include "utilities/openmp_utils.h" //TODO: SOME FILES INCLUDING scheme.h RELY ON THIS. LEAVING AS FUTURE TODO.
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/amgcl_solver.h"
#include "future/solving_strategies/schemes/assembly_helper.h"
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

/// TLS type definition

/**
 * @brief Implicit scheme TLS type definition
 * Thread Local Storage container to be used in the parallel assembly of implicit problems
 * @tparam DataType data type of the problem to be solved
 */
template<class DataType = double >
struct ImplicitThreadLocalStorage
{
    // Local LHS contribution
    DenseMatrix<DataType> LocalLhs;

    // Local RHS constribution
    DenseVector<DataType> LocalRhs;

    // Vector containing the localization in the system of the different terms
    Element::EquationIdVectorType LocalEqIds;
};

/**
 * @class ImplicitScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @author Ruben Zorrilla
 */
//TODO: Think about the template parameters

//TODO: Make base classes:
// scheme.h --> pure virtual!
// implicit_scheme.h --> the one we have in here
template<class TSparseMatrixType, class TSparseVectorType, class TSparseGraphType>
class ImplicitScheme
{
public:
    // FIXME: Does not work... ask @Charlie
    // /// Add scheme to Kratos registry
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.KratosMultiphysics", ImplicitScheme, ImplicitScheme, TSparseMatrixType, TSparseVectorType)
    // KRATOS_REGISTRY_ADD_TEMPLATE_PROTOTYPE("Schemes.All", ImplicitScheme, ImplicitScheme, TSparseMatrixType, TSparseVectorType)

    ///@name Type Definitions
    ///@{

    /// Pointer definition of ImplicitScheme
    KRATOS_CLASS_POINTER_DEFINITION(ImplicitScheme);

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType = typename TSparseMatrixType::IndexType;

    /// Data type definition
    using DataType = typename TSparseMatrixType::DataType;

    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    /// TLS type
    using TLSType = ImplicitThreadLocalStorage<DataType>;

    /// Assembly helper type
    using AssemblyHelperType = Future::AssemblyHelper<TLSType, TSparseMatrixType, TSparseVectorType, TSparseGraphType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initializes the flags
     */
    explicit ImplicitScheme() = default;

    /**
     * @brief Constructor with Parameters
     */
    explicit ImplicitScheme(
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
    explicit ImplicitScheme(ImplicitScheme& rOther)
      : mSchemeIsInitialized(rOther.mSchemeIsInitialized)
      , mSchemeSolutionStepIsInitialized(rOther.mSchemeSolutionStepIsInitialized)
    {
        //TODO: Check this... particularly the mpAssemblyHelper pointer
    }

    /** Destructor.
     */
    virtual ~ImplicitScheme() = default;

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
    virtual typename ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    /**
     * @brief Clone method
     * @return The pointer of the cloned ImplicitScheme
     */
    virtual typename ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>::Pointer Clone()
    {
        return Kratos::make_shared<ImplicitScheme<TSparseMatrixType, TSparseVectorType, TSparseGraphType>>(*this) ;
    }

    /**
     * @brief This is the place to initialize the ImplicitScheme.
     * @details This is intended to be called just once when the strategy is initialized
     */
    virtual void Initialize()
    {
        KRATOS_TRY

        // Check if the Initialize has been already performed
        if (!mSchemeIsInitialized) {
            // Initialize elements, conditions and constraints
            EntitiesUtilities::InitializeAllEntities(*mpModelPart);

            // Set the flag to avoid calling this twice
            mSchemeIsInitialized = true; //TODO: Discuss with the KTC if these should remain or not
        }

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
        DofsArrayType& rDofSet,
        typename TSparseMatrixType::Pointer& rpA,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpB,
        const bool ReformDofSet = true)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'InitializeSolutionStep\' method. Call derived class one." << std::endl;
    }

    /**
     * @brief Function called once at the end of a solution step, after convergence is reached if an iterative process is needed
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    virtual void FinalizeSolutionStep(
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        // Finalizes solution step for all of the elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(*mpModelPart);

        // Reset flags for next step
        mSchemeSolutionStepIsInitialized = false;

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
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
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
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
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

        // Set the corresponding flag
        mDofSetIsInitialized = true;

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished DOFs array set up." << std::endl;
    }

    SizeType SetUpSystemIds(DofsArrayType& rDofSet)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rDofSet.empty()) << "DOFs set is empty. Call the 'SetUpDofArray' first." << std::endl;

        const SizeType equation_system_size = mpAssemblyHelper->SetUpSystemIds(rDofSet);

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished system set up." << std::endl;

        return equation_system_size;

        KRATOS_CATCH("")
    }

    //TODO: Think on the overloads for the mass and damping matrices
    virtual void ResizeAndInitializeVectors(
        const DofsArrayType& rDofSet,
        typename TSparseMatrixType::Pointer& rpLHS,
        typename TSparseVectorType::Pointer& rpDx,
        typename TSparseVectorType::Pointer& rpRHS,
        const bool CalculateReactions = false)
    {
        KRATOS_TRY

        // Call the assembly helper to allocate and initialize the required vectors
        // Note that this also allocates the required reaction vectors (e.g., elimination build)
        (this->GetAssemblyHelper()).ResizeAndInitializeVectors(rDofSet, rpLHS, rpDx, rpRHS, CalculateReactions);

        //TODO: Think on the constraints stuff!
        // ConstructMasterSlaveConstraintsStructure(rModelPart);

        KRATOS_INFO_IF("StaticScheme", this->GetEchoLevel() >= 2) << "Finished system initialization." << std::endl;

        KRATOS_CATCH("")
    }

    virtual void Build(
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItElem->Is(ACTIVE)) {
                // Calculate local LHS and RHS contributions
                ItElem->CalculateLocalSystem(rTLS.LocalLhs, rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalLhs.resize(0,0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItCond->Is(ACTIVE)) {
                // Calculate local LHS and RHS contributions
                ItCond->CalculateLocalSystem(rTLS.LocalLhs, rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalLhs.resize(0,0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        TLSType aux_tls;
        auto& r_assembly_helper = GetAssemblyHelper();
        r_assembly_helper.SetElementAssemblyFunction(elem_func);
        r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        r_assembly_helper.AssembleLocalSystem(rLHS, rRHS, aux_tls);

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    virtual void Build(TSparseVectorType& rRHS)
    {
        Timer::Start("Build");

        const auto timer = BuiltinTimer();

        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItElem->Is(ACTIVE)) {
                // Calculate the RHS contributions
                ItElem->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItCond->Is(ACTIVE)) {
                // Calculate the RHS contributions
                ItCond->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalLhs.resize(0,0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        TLSType aux_tls;
        auto& r_assembly_helper = GetAssemblyHelper();
        r_assembly_helper.SetElementAssemblyFunction(elem_func);
        r_assembly_helper.SetConditionAssemblyFunction(cond_func);
        r_assembly_helper.AssembleRightHandSide(rRHS, aux_tls);

        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 1) << "Build time: " << timer << std::endl;
        KRATOS_INFO_IF("ImplicitScheme", mEchoLevel >= 2) << "Finished parallel building" << std::endl;

        Timer::Stop("Build");
    }

    virtual void Build(TSparseMatrixType& rLHS)
    {
        //TODO: IMPLEMENTATION
    }

    virtual void BuildMassMatrix(TSparseMatrixType& rMassMatrix)
    {
        //TODO: IMPLEMENTATION
    }

    virtual void BuildDampingMatrix(TSparseMatrixType& rDampingMatrix)
    {
        //TODO: IMPLEMENTATION
    }

    //TODO: Think about the dynamic case and the mass and damping matrices!!
    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        TSparseMatrixType& rLHS,
        TSparseVectorType& rRHS)
    {
        GetAssemblyHelper().ApplyDirichletConditions(rDofSet, rLHS, rRHS);
    }

    virtual void ApplyDirichletConditions(
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS)
    {
        GetAssemblyHelper().ApplyDirichletConditions(rDofSet, rRHS);
    }

    //TODO: Think about the dynamic case and the mass and damping matrices!!
    virtual void CalculateReactions(
        const DofsArrayType& rDofSet,
        TSparseVectorType& rRHS)
    {
        const auto elem_func = [](ModelPart::ElementConstantIterator ItElem, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItElem->Is(ACTIVE)) {
                // Calculate the RHS contributions
                ItElem->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItElem->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        const auto cond_func = [](ModelPart::ConditionConstantIterator ItCond, const ProcessInfo& rProcessInfo, TLSType& rTLS){
            if (ItCond->Is(ACTIVE)) {
                // Calculate the RHS contributions
                ItCond->CalculateRightHandSide(rTLS.LocalRhs, rProcessInfo);

                // Get the positions in the global system
                ItCond->EquationIdVector(rTLS.LocalEqIds, rProcessInfo);
            } else {
                rTLS.LocalRhs.resize(0);
                rTLS.LocalLhs.resize(0,0);
                rTLS.LocalEqIds.resize(0);
            }
        };

        TLSType aux_tls;
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
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Predict\' method. Call derived class one." << std::endl;
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
        DofsArrayType &rDofSet,
        TSparseMatrixType &A,
        TSparseVectorType &Dx,
        TSparseVectorType &b)
    {
        KRATOS_ERROR << "\'ImplicitScheme\' does not implement \'Update\' method. Call derived class one." << std::endl;
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
        TSparseMatrixType& A,
        TSparseVectorType& Dx,
        TSparseVectorType& b)
    {
        KRATOS_TRY

        //TODO: Think about creating a datavalue container in here that we can access from outside.

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
            "name" : "implicit_scheme",
            "build_settings" : {
                "build_type" : "block",
                "scaling_type" : "max_diagonal"
            },
            "echo_level" : 0,
            "move_mesh" : false
        })");

        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "implicit_scheme";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method sets the value of mMoveMesh
     * @param MoveMesh If the flag must be set to true or false
     */
    void SetMoveMesh(const bool MoveMesh)
    {
        mMoveMesh = MoveMesh;
    }

    /**
     * @brief This method sets the value of mEchoLevel
     * @param EchoLevel The value to set
     */
    void SetEchoLevel(const int EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * @brief This method sets the value of mDofSetIsInitialized
     * @param DofSetIsInitialized The value to set
     */
    void SetDofSetIsInitialized(const bool DofSetIsInitialized)
    {
        mDofSetIsInitialized = DofSetIsInitialized;
    }

    /**
     * @brief This method sets the value of mSchemeIsInitialized
     * @param SchemeIsInitialized The value to set
     */
    void SetSchemeIsInitialized(const bool SchemeIsInitialized)
    {
        mSchemeIsInitialized = SchemeIsInitialized;
    }

    /**
     * @brief This method sets the value of mSchemeIsInitialized
     * @param SchemeIsInitialized The value to set
     */
    void SetSchemeSolutionStepIsInitialized(const bool SchemeSolutionStepIsInitialized)
    {
        mSchemeSolutionStepIsInitialized = SchemeSolutionStepIsInitialized;
    }


    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    }

    ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This method returns if the mesh has to be updated
     * @return bool True if to be moved, false otherwise
     */
    int GetMoveMesh() const
    {
        return mMoveMesh;
    }

    /**
     * @brief This method returns the echo level value (verbosity level)
     * @return int Echo level value
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /**
     * @brief This method returns if the DOF set is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetDofSetIsInitialized() const
    {
        return mDofSetIsInitialized;
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetSchemeIsInitialized() const
    {
        return mSchemeIsInitialized;
    }

    /**
     * @brief This method returns if the scheme is initialized
     * @return bool True if initialized, false otherwise
     */
    bool GetSchemeSolutionStepIsInitialized() const
    {
        return mSchemeSolutionStepIsInitialized;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ImplicitScheme";
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This function is designed to move the mesh
     * @note It considers DISPLACEMENT as the variable storing the motion. Derive it to adapt to your own strategies.
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(mpModelPart->HasNodalSolutionStepVariable(DISPLACEMENT_X))
            << "It is impossible to move the mesh since the DISPLACEMENT variable is not in the ModelPart. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;

        block_for_each(mpModelPart->Nodes(), [](Node& rNode){
            noalias(rNode.Coordinates()) = rNode.GetInitialPosition().Coordinates();
            noalias(rNode.Coordinates()) += rNode.FastGetSolutionStepValue(DISPLACEMENT);
        });

        KRATOS_INFO_IF("SolvingStrategy", mEchoLevel != 0) << "Mesh moved." << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method validate and assign default parameters
     * @param rParameters Parameters to be validated
     * @param DefaultParameters The default parameters
     * @return Returns validated Parameters
     */
    Parameters ValidateAndAssignParameters(
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
        mMoveMesh = ThisParameters["move_mesh"].GetBool();
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

    int mEchoLevel = 0;

    bool mMoveMesh = false; /// Flag to activate the mesh motion from the DISPLACEMENT variable

    bool mDofSetIsInitialized = false; /// Flag to be used in controlling if the DOF set has been already set

    bool mSchemeIsInitialized = false; /// Flag to be used in controlling if the Scheme has been initialized or not

    bool mSchemeSolutionStepIsInitialized = false; /// Flag to be used in controlling if the Scheme solution step has been initialized or not

    ModelPart* mpModelPart = nullptr;

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
