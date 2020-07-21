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
    typedef ExplicitBuilderAndSolver<TSparseSpace, TDenseSpace> ExplicitBuilderAndSolverType;

    // The explicit builder and solver pointer definition
    typedef typename ExplicitBuilderAndSolverType::Pointer ExplicitBuilderAndSolverPointerType;

    // The DOF type from the explicit builder and solver class
    typedef typename ExplicitBuilderAndSolverType::DofType DofType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitSolvingStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ExplicitSolvingStrategy(
        ModelPart &rModelPart,
        Parameters ThisParameters)
        : mrModelPart(rModelPart)
    {
        const bool rebuild_level = ThisParameters.Has("rebuild_level") ? ThisParameters["rebuild_level"].GetInt() : 0;
        const bool move_mesh_flag = ThisParameters.Has("move_mesh_flag") ? ThisParameters["move_mesh_flag"].GetBool() : false;
        SetMoveMeshFlag(move_mesh_flag);
        SetRebuildLevel(rebuild_level);
        mpExplicitBuilderAndSolver = Kratos::make_unique<ExplicitBuilderAndSolver<TSparseSpace, TDenseSpace>>();
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param pExplicitBuilderAndSolver The pointer to the explicit builder and solver
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     * @param RebuildLevel The flag to set if the DOF set is rebuild or not
     */
    explicit ExplicitSolvingStrategy(
        ModelPart &rModelPart,
        typename ExplicitBuilderAndSolverType::Pointer pExplicitBuilderAndSolver,
        bool MoveMeshFlag = false,
        int RebuildLevel = 0)
        : mrModelPart(rModelPart),
          mpExplicitBuilderAndSolver(pExplicitBuilderAndSolver)
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
        : mrModelPart(rModelPart)
    {
        SetMoveMeshFlag(MoveMeshFlag);
        SetRebuildLevel(RebuildLevel);
        mpExplicitBuilderAndSolver = Kratos::make_unique<ExplicitBuilderAndSolver<TSparseSpace, TDenseSpace>>();
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategy(const ExplicitSolvingStrategy &Other) = delete;

    /** Destructor.
     */
    virtual ~ExplicitSolvingStrategy() = default;

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
        if (!GetInitializeWasPerformedFlag()) {
            Initialize();
        }

        if (!GetInitializeSolutionStepWasPerformedFlag()) {
            InitializeSolutionStep();
        }
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize()
    {
        if (!GetInitializeWasPerformedFlag()) {
            // Initialize elements, conditions and constraints
            InitializeContainer(GetModelPart().Elements());
            InitializeContainer(GetModelPart().Conditions());
            InitializeContainer(GetModelPart().MasterSlaveConstraints());

            // Set the explicit DOFs rebuild level
            if (mRebuildLevel != 0) {
                mpExplicitBuilderAndSolver->SetResetDofSetFlag(true);
            }

            // If the mesh is updated at each step, we require to accordingly update the lumped mass at each step
            if (mMoveMeshFlag) {
                mpExplicitBuilderAndSolver->SetResetLumpedMassVectorFlag(true);
            }

            // Call the explicit builder and solver initialize (Set up DOF set and lumped mass vector)
            mpExplicitBuilderAndSolver->Initialize(mrModelPart);

            // Initialize the solution values
            InitializeDofSetValues();

            // Set the mInitializeWasPerformed flag
            mInitializeWasPerformed = true;
        }
    }

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear()
    {
        // This clears the DOF set and lumped mass vector
        mpExplicitBuilderAndSolver->Clear();

        // Initialize the explicit strategy flags
        mInitializeWasPerformed = false;
        mInitializeSolutionStepWasPerformed = false;
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
     * @brief This operations should be called before printing the results when non trivial results (e.g. stresses)
     * need to be calculated given the solution of the step
     * @details This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    virtual void CalculateOutputData()
    {
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void InitializeSolutionStep()
    {
        // Check if the Initialize() has been already performed
        if (!mInitializeWasPerformed) {
            Initialize();
        }

        // InitializeSolutionStep elements, conditions and constraints
        InitializeSolutionStepContainer(GetModelPart().Elements());
        InitializeSolutionStepContainer(GetModelPart().Conditions());
        InitializeSolutionStepContainer(GetModelPart().MasterSlaveConstraints());

        // Call the builder and solver initialize solution step
        mpExplicitBuilderAndSolver->InitializeSolutionStep(mrModelPart);

        // Set the mInitializeSolutionStepWasPerformed flag
        mInitializeSolutionStepWasPerformed = true;
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
        mpExplicitBuilderAndSolver->FinalizeSolutionStep(mrModelPart);

        // Reset the mInitializeSolutionStepWasPerformed flag
        mInitializeSolutionStepWasPerformed = false;
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    virtual bool SolveSolutionStep()
    {
        // Solve the problem assuming that a lumped mass matrix is used
        SolveWithLumpedMassMatrix();

        // If required, update the mesh with the obtained solution
        if (mMoveMeshFlag) {
            MoveMesh();
        }

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
    virtual void SetEchoLevel(const int Level)
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
    virtual int GetEchoLevel()
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
    virtual void SetRebuildLevel(int Level)
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
    virtual int GetRebuildLevel()
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
    bool MoveMeshFlag()
    {
        return mMoveMeshFlag;
    }

    /**
     * @brief This function sets the flag that says that the Initialize() has been performed
     * @param Flag True if the Initialize() has been performed, false otherwise
     */
    void SetInitializeWasPerformedFlag(bool Flag)
    {
        mInitializeWasPerformed = Flag;
    }

    /**
     * @brief This function returns the flag that says that the Initialize() has been performed
     * @return True if the Initialize() has been performed, false otherwise
     */
    bool GetInitializeWasPerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This function sets the flag that says that the InitializeSolutionStep() has been performed
     * @param Flag True if the InitializeSolutionStep() has been performed, false otherwise
     */
    void SetInitializeSolutionStepWasPerformedFlag(bool Flag)
    {
        mInitializeSolutionStepWasPerformed = Flag;
    }

    /**
     * @brief This function returns the flag that says that the InitializeSolutionStep() has been performed
     * @return True if the InitializeSolutionStep() has been performed, false otherwise
     */
    bool GetInitializeSolutionStepWasPerformedFlag()
    {
        return mInitializeSolutionStepWasPerformed;
    }

    /**
     * @brief This function is designed to move the mesh
     * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        auto& r_nodes_array = GetModelPart().Nodes();
        block_for_each(
            r_nodes_array,
            [](Node<3>& rNode){
                noalias(rNode.Coordinates()) = rNode.GetInitialPosition().Coordinates();
                noalias(rNode.Coordinates()) += rNode.FastGetSolutionStepValue(DISPLACEMENT);
            }
        );

        KRATOS_INFO_IF("ExplicitSolvingStrategy", this->GetEchoLevel() > 0) << "Mesh moved." << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Operations to get the pointer to the model
     * @return mrModelPart: The model part member variable
     */
    inline ModelPart& GetModelPart()
    {
        return mrModelPart;
    };

    /**
     * @brief Operations to get the pointer to the explicit builder and solver
     * @return mpExplicitBuilderAndSolver: The explicit builder and solver
     */
    inline ExplicitBuilderAndSolverPointerType& pGetExplicitBuilderAndSolver()
    {
        return mpExplicitBuilderAndSolver;
    };

    /**
     * @brief Operations to get the residual norm
     * @return The residual norm
     */
    virtual double GetResidualNorm()
    {
        // Get the required data from the explicit builder and solver
        const auto p_explicit_bs = pGetExplicitBuilderAndSolver();
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
    virtual int Check()
    {
        KRATOS_TRY

        // Check if displacement var is needed
        if (mMoveMeshFlag == true) {
            VariableUtils().CheckVariableExists<>(DISPLACEMENT, GetModelPart().Nodes());
        }

        // Elements check
        for (const auto& r_elem :  GetModelPart().Elements()) {
            const auto& r_process_info = GetModelPart().GetProcessInfo();
            r_elem.Check(r_process_info);
        }

        // Conditions check
        for (const auto& r_cond :  GetModelPart().Conditions()) {
            const auto& r_process_info = GetModelPart().GetProcessInfo();
            r_cond.Check(r_process_info);
        }

        return 0;

        KRATOS_CATCH("")
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

    ModelPart &mrModelPart;

    bool mMoveMeshFlag;

    bool mInitializeWasPerformed = false;

    bool mInitializeSolutionStepWasPerformed = false;

    ExplicitBuilderAndSolverPointerType mpExplicitBuilderAndSolver = nullptr;

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
        auto& r_dof_set = mpExplicitBuilderAndSolver->GetDofSet();
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
