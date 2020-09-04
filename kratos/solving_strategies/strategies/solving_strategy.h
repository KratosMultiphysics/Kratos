//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_SOLVING_STRATEGY )
#define  KRATOS_SOLVING_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
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

/** @brief Solving strategy base class
 * @details This is the base class from which we will derive all the strategies (line-search, NR, etc...)
 */

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class SolvingStrategy
{
public:
    ///@name Type Definitions
    ///@{

//     typedef std::set<Dof::Pointer,ComparePDof>                                    DofSetType;

    typedef typename TSparseSpace::DataType                                        TDataType;

    typedef typename TSparseSpace::MatrixType                              TSystemMatrixType;

    typedef typename TSparseSpace::VectorType                              TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType                TSystemMatrixPointerType;

    typedef typename TSparseSpace::VectorPointerType                TSystemVectorPointerType;

    typedef typename TDenseSpace::MatrixType                           LocalSystemMatrixType;

    typedef typename TDenseSpace::VectorType                           LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace>                                    TSchemeType;

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>              ClassType;

    typedef typename ModelPart::DofType                                             TDofType;

    typedef typename ModelPart::DofsArrayType                                  DofsArrayType;

//     typedef Dof<TDataType>                                                          TDofType;

//     typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> >            DofsArrayType;

//     typedef PointerVectorSet<TDofType, IndexedObject>                          DofsArrayType;

    typedef typename DofsArrayType::iterator                                 DofIteratorType;

    typedef typename DofsArrayType::const_iterator                   DofConstantIteratorType;

    typedef ModelPart::NodesContainerType                                     NodesArrayType;

    typedef ModelPart::ElementsContainerType                               ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                           ConditionsArrayType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(SolvingStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit SolvingStrategy() { }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit SolvingStrategy(ModelPart& rModelPart, Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param MoveMeshFlag The flag to set if the mesh is moved or not
     */
    explicit SolvingStrategy(
        ModelPart& rModelPart,
        bool MoveMeshFlag = false
        ) : mpModelPart(&rModelPart)
    {
        SetMoveMeshFlag(MoveMeshFlag);
    }

    /** Destructor.
     */
    virtual ~SolvingStrategy(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    virtual typename ClassType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

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
     * 0 -> Build StiffnessMatrix just once
     * 1 -> Build StiffnessMatrix at the beginning of each solution step
     * 2 -> build StiffnessMatrix at each iteration
     * }
     */
    virtual void SetRebuildLevel(int Level)
    {
        mRebuildLevel = Level;
        mStiffnessMatrixIsBuilt = false;
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

        NodesArrayType& NodesArray = GetModelPart().Nodes();
        const int numNodes = static_cast<int>(NodesArray.size());

        #pragma omp parallel for
        for(int i = 0; i < numNodes; ++i) {
            auto it_node = NodesArray.begin() + i;

            noalias(it_node->Coordinates()) = it_node->GetInitialPosition().Coordinates();
            noalias(it_node->Coordinates()) += it_node->FastGetSolutionStepValue(DISPLACEMENT);
        }

        KRATOS_INFO_IF("SolvingStrategy", this->GetEchoLevel() != 0 && GetModelPart().GetCommunicator().MyPID() == 0) <<" MESH MOVED "<<std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Operations to get the pointer to the model
     * @return mpModelPart: The model part member variable
     */
    inline ModelPart& GetModelPart()
    {
        KRATOS_ERROR_IF_NOT(mpModelPart) << "ModelPart in the SolvingStrategy is not initialized" << std::endl;
        return *mpModelPart;
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
        if (mMoveMeshFlag == true)
        {
            for (ModelPart::NodesContainerType::iterator itNode = GetModelPart().NodesBegin();
                 itNode != GetModelPart().NodesEnd(); itNode++)
            {
                if (itNode->SolutionStepsDataHas(DISPLACEMENT) == false)
                {
                    KRATOS_ERROR << "ERROR:: Problem on node with Id " << itNode->Id() << "\nIt is impossible to move the mesh since the DISPLACEMENT var is not in the rModelPart. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables" << std::endl;
                }
            }
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
            "name"                         : "solving_strategy",
            "move_mesh_flag"               : false,
            "echo_level"                   : 1,
            "build_level"                  : 2
        })");
        return default_parameters;
    }
  
    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    virtual TSystemMatrixType& GetSystemMatrix()
    {
        KRATOS_TRY

        KRATOS_ERROR << "GetSystemMatrix not implemented in base SolvingStrategy" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    virtual TSystemVectorType& GetSystemVector()
    {
        KRATOS_TRY

        KRATOS_ERROR << "GetSystemVector not implemented in base SolvingStrategy" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    virtual TSystemVectorType& GetSolutionVector()
    {
        KRATOS_TRY

        KRATOS_ERROR << "GetSolutionVector not implemented in base SolvingStrategy" << std::endl;

        KRATOS_CATCH("");
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "solving_strategy";
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "SolvingStrategy";
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

protected:
    ///@name Protected static Member Variables
    ///@{

    // Level of echo for the solving strategy
    int mEchoLevel;

    // Settings for the rebuilding of the stiffness matrix
    int mRebuildLevel;
    bool mStiffnessMatrixIsBuilt;

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
        // By default mesh is not moved
        mMoveMeshFlag = ThisParameters["move_mesh_flag"].GetBool();

        // Be default the minimal information is shown
        mEchoLevel = ThisParameters["echo_level"].GetInt();

        // By default the matrices are rebuilt at each iteration
        mRebuildLevel = ThisParameters["build_level"].GetInt();
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

    /** Copy constructor.
     */
    SolvingStrategy(const SolvingStrategy& Other);


    ///@}

}; /* Class NewSolvingStrategy */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_SOLVING_STRATEGY  defined */
