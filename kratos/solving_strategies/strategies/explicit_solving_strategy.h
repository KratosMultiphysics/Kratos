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
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
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

/** @brief Explicit solving strategy base class
 * @details This is the base class from which we will derive all the explicit strategies (FE, RK4, ...)
 */
template <class TSparseSpace, class TDenseSpace>
class ExplicitSolvingStrategy
{
public:
    ///@name Type Definitions
    ///@{

    // The explicit builder and solver definition
    typedef ExplicitBuilder<TSparseSpace, TDenseSpace> ExplicitBuilderType;

    // The explicit builder and solver pointer definition
    typedef typename ExplicitBuilderType::Pointer ExplicitBuilderPointerType;

    // The DOF type from the explicit builder and solver class
    typedef typename ExplicitBuilderType::DofType DofType;

    /// The definition of the current class
    typedef ExplicitSolvingStrategy<TSparseSpace, TDenseSpace> ClassType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (empty)
     */
    explicit ExplicitSolvingStrategy()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategy(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
        mpExplicitBuilder = Kratos::make_unique<ExplicitBuilder<TSparseSpace, TDenseSpace>>();
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilder The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     * @param RebuildLevel The flag to set if the DOF set is rebuild or not
     */
    explicit ExplicitSolvingStrategy(
        ModelPart &rModelPart,
        typename ExplicitBuilderType::Pointer pExplicitBuilder,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : mpModelPart(&rModelPart),
          mpExplicitBuilder(pExplicitBuilder)
    {
        SetMoveMeshFlag(MoveMeshFlag);
        SetRebuildLevel(RebuildLevel);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     * @param RebuildLevel The flag to set if the DOF set is rebuild or not
     */
    explicit ExplicitSolvingStrategy(
        ModelPart &rModelPart,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : mpModelPart(&rModelPart)
    {
        SetMoveMeshFlag(MoveMeshFlag);
        SetRebuildLevel(RebuildLevel);
        mpExplicitBuilder = Kratos::make_unique<ExplicitBuilder<TSparseSpace, TDenseSpace>>();
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategy(const ExplicitSolvingStrategy &Other) = delete;

    /** Destructor.
     */
    virtual ~ExplicitSolvingStrategy()
    {
        mpModelPart =  nullptr;
    }

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    virtual typename ClassType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
     * values of the solution step of interest are assumed equal to the old values
     */
    virtual void Predict()
    {
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize()
    {
        // Initialize elements, conditions and constraints
        InitializeContainer(GetModelPart().Elements());
        InitializeContainer(GetModelPart().Conditions());
        InitializeContainer(GetModelPart().MasterSlaveConstraints());

        // Set the explicit DOFs rebuild level
        if (mRebuildLevel != 0) {
            mpExplicitBuilder->SetResetDofSetFlag(true);
        }

        // If the mesh is updated at each step, we require to accordingly update the lumped mass at each step
        if (mMoveMeshFlag) {
            mpExplicitBuilder->SetResetLumpedMassVectorFlag(true);
        }

        // Call the explicit builder and solver initialize (Set up DOF set and lumped mass vector)
        mpExplicitBuilder->Initialize(*mpModelPart);

        // Initialize the solution values
        InitializeDofSetValues();
    }

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear()
    {
        // This clears the DOF set and lumped mass vector
        mpExplicitBuilder->Clear();
    }

    /**
     * @brief This should be considered as a "post solution" convergence check which is useful for coupled analysis
     * @details The convergence criteria used is the one used inside the "solve" step
     */
    virtual bool IsConverged()
    {
        return true;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep()
    {
        // InitializeSolutionStep elements, conditions and constraints
        InitializeSolutionStepContainer(GetModelPart().Elements());
        InitializeSolutionStepContainer(GetModelPart().Conditions());
        InitializeSolutionStepContainer(GetModelPart().MasterSlaveConstraints());

        // Call the builder and solver initialize solution step
        mpExplicitBuilder->InitializeSolutionStep(*mpModelPart);
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void FinalizeSolutionStep()
    {
        // FinalizeSolutionStep elements, conditions and constraints
        FinalizeSolutionStepContainer(GetModelPart().Elements());
        FinalizeSolutionStepContainer(GetModelPart().Conditions());
        FinalizeSolutionStepContainer(GetModelPart().MasterSlaveConstraints());

        // Call the builder and solver finalize solution step (the reactions are computed in here)
        mpExplicitBuilder->FinalizeSolutionStep(*mpModelPart);
    }

    /**
     * @brief Solves the current step.
     * The function always return true as convergence is not checked in the explicit framework
     */
    virtual bool SolveSolutionStep()
    {
        // Call the initialize non-linear iteration
        InitializeNonLinearIterationContainer(GetModelPart().Elements());
        InitializeNonLinearIterationContainer(GetModelPart().Conditions());
        InitializeNonLinearIterationContainer(GetModelPart().MasterSlaveConstraints());

        // Apply constraints
        if(mpModelPart->MasterSlaveConstraints().size() != 0) {
            mpExplicitBuilder->ApplyConstraints(*mpModelPart);
        }

        // Solve the problem assuming that a lumped mass matrix is used
        SolveWithLumpedMassMatrix();

        // If required, update the mesh with the obtained solution
        if (mMoveMeshFlag) {
            MoveMesh();
        }

        // Call the finalize non-linear iteration
        FinalizeNonLinearIterationContainer(GetModelPart().Elements());
        FinalizeNonLinearIterationContainer(GetModelPart().Conditions());
        FinalizeNonLinearIterationContainer(GetModelPart().MasterSlaveConstraints());

        return true;
    }

    /**
     * @brief This sets the level of echo for the solving strategy
     * @param Level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing advanced information
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     */
    void SetEchoLevel(const int Level)
    {
        mEchoLevel = Level;
    }

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing advanced information
     * 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
     * }
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /**
     * This sets the build level
     * @param Level The build level
     * @details
     * {
     * 0 -> Set up the DOF set just once
     * 1 -> Set up the DOF set at the beginning of each solution step
     * }
     */
    void SetRebuildLevel(int Level)
    {
        mRebuildLevel = Level;
    }

    /**
     * @brief This returns the build level
     * @details
     * {
     * 0 -> Set up the DOF set just once
     * 1 -> Set up the DOF set at the beginning of each solution step
     * }
     * @return The build level
     */
    int GetRebuildLevel() const
    {
        return mRebuildLevel;
    }

    /**
     * @brief This function sets the flag that says if the mesh is moved
     * @param Flag True if the mesh is moved, false otherwise
     */
    void SetMoveMeshFlag(bool Flag)
    {
        mMoveMeshFlag = Flag;
    }

    /**
     * @brief This function returns the flag that says if the mesh is moved
     * @return True if the mesh is moved, false otherwise
     */
    bool MoveMeshFlag() const
    {
        return mMoveMeshFlag;
    }

    /**
     * @brief This function is designed to move the mesh
     * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        auto& r_nodes_array = GetModelPart().Nodes();
        VariableUtils().UpdateCurrentPosition(r_nodes_array);

        KRATOS_INFO_IF("ExplicitSolvingStrategy", this->GetEchoLevel() > 0) << "Mesh moved." << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    };

    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    };

    /**
     * @brief Operations to get the pointer to the explicit builder and solver
     * @return mpExplicitBuilder: The explicit builder and solver
     */
    ExplicitBuilderPointerType& pGetExplicitBuilder()
    {
        return mpExplicitBuilder;
    };

    /**
     * @brief Operations to get the residual norm
     * @return The residual norm
     */
    virtual double GetResidualNorm()
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = pGetExplicitBuilder();
        auto& r_dof_set = p_explicit_bs->GetDofSet();

        // Calculate the explicit residual
        p_explicit_bs->BuildRHS(GetModelPart());

        // Calculate the residual norm
        double res_norm = 0.0;
        res_norm = block_for_each<SumReduction<double>>(
            r_dof_set,
            [](DofType &rDof){return rDof.GetSolutionStepReactionValue();}
        );

        return res_norm;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    virtual int Check() const
    {
        KRATOS_TRY

        // Check if displacement var is needed
        if (mMoveMeshFlag) {
            VariableUtils().CheckVariableExists<>(DISPLACEMENT, GetModelPart().Nodes());
        }

        // Check elements, conditions and constraints
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        for (const auto& r_elem : GetModelPart().Elements()) {
            r_elem.Check(r_process_info);
        }
        for (const auto& r_cond : GetModelPart().Conditions()) {
            r_cond.Check(r_process_info);
        }
        for (const auto& r_cons : GetModelPart().MasterSlaveConstraints()) {
            r_cons.Check(r_process_info);
        }

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "explicit_solving_strategy" : "explicit_solving_strategy",
            "move_mesh_flag"            : false,
            "rebuild_level"             : 0,
            "echo_level"                : 1,
            "explicit_builder_settings" : {
                "name": "explicit_builder"
            }
        })");
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "explicit_solving_strategy";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ExplicitSolvingStrategy";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
        rOStream << Info();
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    // Level of echo for the solving strategy
    int mEchoLevel;

    // Settings for the rebuilding of the DOF set
    int mRebuildLevel;

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
     * @brief Calculate the explicit update
     * This method is intended to implement the explicit update calculation
     * Note that it has to be implemented according to the explicit scheme in a derived class
     */
    virtual void SolveWithLumpedMassMatrix()
    {
        KRATOS_ERROR << "Calling the base ExplicitSolvingStrategy SolveWithLumpedMassMatrix(). Implement the specific explicit scheme solution update in a derived class" << std::endl;
    }

    /**
     * @brief Get the Delta Time object
     * This method returns the DELTA_TIME from the ProcessInfo container
     * @return const double
     */
    virtual inline double GetDeltaTime()
    {
        return GetModelPart().GetProcessInfo().GetValue(DELTA_TIME);
    }

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
        const bool rebuild_level = ThisParameters["rebuild_level"].GetInt();
        const bool move_mesh_flag = ThisParameters["move_mesh_flag"].GetBool();
        SetMoveMeshFlag(move_mesh_flag);
        SetRebuildLevel(rebuild_level);
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

private:
    ///@}
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;

    bool mMoveMeshFlag;

    ExplicitBuilderPointerType mpExplicitBuilder = nullptr;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Initialize the DOF set values
     * This method initializes the current value of the unknown variables in the DOF set
     */
    void InitializeDofSetValues()
    {
        // Initialize the DOF values
        auto& r_dof_set = mpExplicitBuilder->GetDofSet();
        block_for_each(
            r_dof_set,
            [](DofType& rDof){
                auto &r_value = rDof.GetSolutionStepValue();
                r_value = 0.0;
            }
        );
    }

    /**
     * @brief Auxiliary call to the Initialize()
     * For a given container, this calls the Initialize() method
     * @tparam TContainerType Container type template (e.g. elements, conditions, ...)
     * @param rContainer Reference to the container
     */
    template <class TContainerType>
    void InitializeContainer(TContainerType &rContainer)
    {
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        block_for_each(
            rContainer,
            [&](typename TContainerType::value_type& rEntity){rEntity.Initialize(r_process_info);}
        );
    }

    /**
     * @brief Auxiliary call to the InitializeSolutionStep()
     * For a given container, this calls the InitializeSolutionStep() method
     * @tparam TContainerType Container type template (e.g. elements, conditions, ...)
     * @param rContainer Reference to the container
     */
    template <class TContainerType>
    void InitializeSolutionStepContainer(TContainerType &rContainer)
    {
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        block_for_each(
            rContainer,
            [&](typename TContainerType::value_type& rEntity){rEntity.InitializeSolutionStep(r_process_info);}
        );
    }

    /**
     * @brief Auxiliary call to the InitializeNonLinearIteration()
     * For a given container, this calls the InitializeNonLinearIteration() method
     * @tparam TContainerType Container type template (e.g. elements, conditions, ...)
     * @param rContainer Reference to the container
     */
    template <class TContainerType>
    void InitializeNonLinearIterationContainer(TContainerType &rContainer)
    {
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        block_for_each(
            rContainer,
            [&](typename TContainerType::value_type& rEntity){rEntity.InitializeNonLinearIteration(r_process_info);}
        );
    }

    /**
     * @brief Auxiliary call to the FinalizeNonLinearIteration()
     * For a given container, this calls the FinalizeNonLinearIteration() method
     * @tparam TContainerType Container type template (e.g. elements, conditions, ...)
     * @param rContainer Reference to the container
     */
    template <class TContainerType>
    void FinalizeNonLinearIterationContainer(TContainerType &rContainer)
    {
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        block_for_each(
            rContainer,
            [&](typename TContainerType::value_type& rEntity){rEntity.FinalizeNonLinearIteration(r_process_info);}
        );
    }

    /**
     * @brief Auxiliary call to the FinalizeSolutionStep()
     * For a given container, this calls the FinalizeSolutionStep() method
     * @tparam TContainerType Container type template (e.g. elements, conditions, ...)
     * @param rContainer Reference to the container
     */
    template <class TContainerType>
    void FinalizeSolutionStepContainer(TContainerType &rContainer)
    {
        const auto& r_process_info = GetModelPart().GetProcessInfo();
        block_for_each(
            rContainer,
            [&](typename TContainerType::value_type& rEntity){rEntity.FinalizeSolutionStep(r_process_info);}
        );
    }

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
}; /* Class NewExplicitSolvingStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY  defined */
