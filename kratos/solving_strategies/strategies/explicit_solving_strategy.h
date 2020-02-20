//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

#if !defined(KRATOS_EXPLICIT_SOLVING_STRATEGY)
#define KRATOS_EXPLICIT_SOLVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/calculate_nodal_area_process.h"
#include "solving_strategies/builder_and_solvers/explicit_builder_and_solver.h"

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

    //     typedef std::set<Dof::Pointer,ComparePDof>                                    DofSetType;

    typedef typename TSparseSpace::DataType TDataType;

    typedef typename TSparseSpace::MatrixType TSystemMatrixType;

    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;

    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;

    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef ExplicitBuilderAndSolver<TSparseSpace, TDenseSpace> ExplicitBuilderAndSolverType;

    typedef typename ModelPart::DofType TDofType;

    typedef typename ModelPart::DofsArrayType DofsArrayType;

    //     typedef Dof<TDataType>                                                          TDofType;

    //     typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> >            DofsArrayType;

    //     typedef PointerVectorSet<TDofType, IndexedObject>                          DofsArrayType;

    typedef typename DofsArrayType::iterator DofIteratorType;

    typedef typename DofsArrayType::const_iterator DofConstantIteratorType;

    typedef ModelPart::NodesContainerType NodesArrayType;

    typedef ModelPart::ElementsContainerType ElementsArrayType;

    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

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
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize()
    {
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
    }

    /**
     * @brief The problem of interest is solved.
     * @details
     * {
     * This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(), SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     * }
     */
    virtual double Solve()
    {
        Initialize();
        InitializeSolutionStep();
        Predict();
        SolveSolutionStep();
        FinalizeSolutionStep();

        return 0.0;
    }

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear()
    {
        mpExplicitBuilderAndSolver->Clear();
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
        // Call the builder and solver initialize solution step
        mpExplicitBuilderAndSolver->InitializeSolutionStep(mrModelPart);

    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void FinalizeSolutionStep()
    {
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    virtual bool SolveSolutionStep()
    {
        InitializeDofSetValues();

        UpdateDofSetValues();

        SolveWithLumpedMassMatrix();

        return true;
    }

    /**
     * @brief This sets the level of echo for the solving strategy
     * @param Level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic informations
     * 2 -> Printing linear solver data
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
     * 2 -> Printing linear solver data
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
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
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
     * @brief This function is designed to move the mesh
     * @note Be careful it just consider displacements, derive this method to adapt to your own strategies (ALE, FSI, etc...)
     */
    virtual void MoveMesh()
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(GetModelPart().NodesBegin()->SolutionStepsDataHas(DISPLACEMENT_X) == false) << "It is impossible to move the mesh since the DISPLACEMENT var is not in the Model Part. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;

        auto& r_nodes_array = GetModelPart().Nodes();
        const int n_nodes = static_cast<int>(r_nodes_array.size());

#pragma omp parallel for
        for (int i_node = 0; i_node < n_nodes; ++i_node) {
            auto it_node = r_nodes_array.begin() + i_node;
            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
        }

        KRATOS_INFO_IF("ExplicitSolvingStrategy", this->GetEchoLevel() != 0 && GetModelPart().GetCommunicator().MyPID() == 0) << " MESH MOVED " << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Operations to get the pointer to the model
     * @return mrModelPart: The model part member variable
     */
    inline ModelPart &GetModelPart()
    {
        return mrModelPart;
    };

    /**
     * @brief Operations to get the residual norm
     * @return The residual norm
     */
    virtual double GetResidualNorm()
    {
        return 0.0;
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
            for (const auto& r_node : GetModelPart().Nodes()) {
                KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(DISPLACEMENT)) << "ERROR:: Problem on node with Id " << r_node.Id() <<
                "\nIt is impossible to move the mesh since the DISPLACEMENT var is not in the rModelPart. Either use (False) or add DISPLACEMENT to the list of variables" << std::endl;
            }
        }

        // Check NODAL_AREA variable. It is needed to do the solution update
            for (const auto& r_node : GetModelPart().Nodes()) {
                KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(NODAL_AREA)) << "ERROR:: Problem on node with Id " << r_node.Id() <<
                "\nIt is impossible to do the solution update since the NODAL_AREA var is not in the rModelPart. Add NODAL_AREA to the list of variables" << std::endl;
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

    virtual void UpdateDofSetValues()
    {
        KRATOS_WATCH("To be particularized in each derived strategy")

        // Get the current delta time
        const auto &r_process_info = mrModelPart.GetProcessInfo();
        const double dt = r_process_info.GetValue(DELTA_TIME);

        // Perform the solution update
        // TODO: This has to be moved to a derived RK4 class


        // Divide by the lumped mass matrix values --> THIS CAN'T BE DONE THROUGH THE DOFs... HOW DO WE DO THIS?
        // SHOULD THE RHS (REACTION) BE ALREADY SUBDIVIDED BY THE MASS?
        // WE CAN DO IT BY SAVING THE NODAL_AREA IN DE HISTORICAL DATA BASE AND it_dof->GetSolutionStepValue(NODAL_AREA)
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

    typename ExplicitBuilderAndSolverType::Pointer mpExplicitBuilderAndSolver = nullptr;

    bool mMoveMeshFlag;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void InitializeDofSetValues()
    {
        // Initialize the DOF values
        auto& r_dof_set = mpExplicitBuilderAndSolver->GetDofSet();
        const unsigned int n_dofs = r_dof_set.size();
#pragma omp parallel for
        for (int i_dof = 0; i_dof < n_dofs; ++i_dof) {
            auto it_dof = r_dof_set.begin() + i_dof;
            auto &r_value = it_dof->GetSolutionStepValue();
            r_value = 0.0;
        }
    }

    void SolveWithLumpedMassMatrix()
    {

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
