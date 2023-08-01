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
#include "solving_strategy.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "solving_strategies/builder_and_solvers/explicit_builder.h"
#include "includes/kratos_parameters.h"
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

/** @brief Explicit solving strategy base class
 * @details This is the base class from which we will derive all the explicit strategies (FE, RK4, ...)
 */
template <class TSparseSpace, class TDenseSpace>
class ExplicitSolvingStrategy : public SolvingStrategy<TSparseSpace, TDenseSpace>
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

    /// The definition of the base class
    typedef SolvingStrategy<TSparseSpace, TDenseSpace> BaseType;

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
        : BaseType(rModelPart, ThisParameters)
    {
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
        : BaseType(rModelPart, MoveMeshFlag)
        , mpExplicitBuilder(pExplicitBuilder)
    {
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
        : BaseType(rModelPart, MoveMeshFlag)
    {
        SetRebuildLevel(RebuildLevel);
        mpExplicitBuilder = Kratos::make_unique<ExplicitBuilder<TSparseSpace, TDenseSpace>>();
    }

    /** Copy constructor.
     */
    ExplicitSolvingStrategy(const ExplicitSolvingStrategy &Other) = delete;

    /** Destructor.
     */
    virtual ~ExplicitSolvingStrategy() {}

    /**
     * @brief Create method
     * @param rModelPart The model part to be computed
     * @param ThisParameters The configuration parameters
     */
    typename BaseType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
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
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        // Initialize elements, conditions and constraints
        EntitiesUtilities::InitializeAllEntities(BaseType::GetModelPart());

        // Set the explicit DOFs rebuild level
        if (mRebuildLevel != 0) {
            mpExplicitBuilder->SetResetDofSetFlag(true);
        }

        // If the mesh is updated at each step, we require to accordingly update the lumped mass at each step
        if (BaseType::GetMoveMeshFlag()) {
            mpExplicitBuilder->SetResetLumpedMassVectorFlag(true);
        }

        // Call the explicit builder and solver initialize (Set up DOF set and lumped mass vector)
        mpExplicitBuilder->Initialize(BaseType::GetModelPart());

    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        // This clears the DOF set and lumped mass vector
        mpExplicitBuilder->Clear();
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        // InitializeSolutionStep elements, conditions and constraints
        EntitiesUtilities::InitializeSolutionStepAllEntities(BaseType::GetModelPart());

        // Call the builder and solver initialize solution step
        mpExplicitBuilder->InitializeSolutionStep(BaseType::GetModelPart());
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        // FinalizeSolutionStep elements, conditions and constraints
        EntitiesUtilities::FinalizeSolutionStepAllEntities(BaseType::GetModelPart());

        // Call the builder and solver finalize solution step (the reactions are computed in here)
        mpExplicitBuilder->FinalizeSolutionStep(BaseType::GetModelPart());
    }

    /**
     * @brief Solves the current step.
     * The function always return true as convergence is not checked in the explicit framework
     */
    bool SolveSolutionStep() override
    {
        // Call the initialize non-linear iteration
        EntitiesUtilities::InitializeNonLinearIterationAllEntities(BaseType::GetModelPart());

        // Apply constraints
        if(BaseType::GetModelPart().MasterSlaveConstraints().size() != 0) {
            mpExplicitBuilder->ApplyConstraints(BaseType::GetModelPart());
        }

        // Solve the problem assuming that a lumped mass matrix is used
        SolveWithLumpedMassMatrix();

        // If required, update the mesh with the obtained solution
        if (BaseType::GetMoveMeshFlag()) {
            BaseType::MoveMesh();
        }

        // Call the finalize non-linear iteration
        EntitiesUtilities::FinalizeNonLinearIterationAllEntities(BaseType::GetModelPart());

        return true;
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
    void SetRebuildLevel(int Level) override
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
    int GetRebuildLevel() const override
    {
        return mRebuildLevel;
    }

    /**
     * @brief Operations to get the pointer to the explicit builder and solver
     * @return mpExplicitBuilder: A pointer to the explicit builder and solver
     */
    ExplicitBuilderPointerType pGetExplicitBuilder()
    {
        return mpExplicitBuilder;
    };

    /**
     * @brief Operations to get the explicit builder and solver
     * @return The explicit builder and solver
     */
    ExplicitBuilderType& GetExplicitBuilder()
    {
        KRATOS_ERROR_IF(mpExplicitBuilder == nullptr) << "Asking for builder and solver when it is empty" << std::endl;
        return *mpExplicitBuilder;
    };

    /**
     * @brief Operations to get the explicit builder and solver
     * @return The explicit builder and solver
     */
    const ExplicitBuilderType& GetExplicitBuilder() const
    {
        KRATOS_ERROR_IF(mpExplicitBuilder == nullptr) << "Asking for builder and solver when it is empty" << std::endl;
        return *mpExplicitBuilder;
    };

    /**
     * @brief Operations to get the residual norm
     * @return The residual norm
     */
    double GetResidualNorm() override
    {
        // Get the required data from the explicit builder and solver
        auto& r_explicit_bs = GetExplicitBuilder();
        auto& r_dof_set = r_explicit_bs.GetDofSet();

        // Calculate the explicit residual
        r_explicit_bs.BuildRHS(BaseType::GetModelPart());

        // Calculate the residual norm
        double res_norm = 0.0;
        res_norm = block_for_each<SumReduction<double>>(
            r_dof_set,
            [](DofType &rDof){return rDof.GetSolutionStepReactionValue();}
        );

        return res_norm;
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "explicit_solving_strategy" : "explicit_solving_strategy",
            "rebuild_level"             : 0,
            "explicit_builder_settings" : {
                "name": "explicit_builder"
            }
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);

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
    std::string Info() const override
    {
        return "ExplicitSolvingStrategy";
    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

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
        return BaseType::GetModelPart().GetProcessInfo().GetValue(DELTA_TIME);
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Add base strategy settings
        BaseType::AssignSettings(ThisParameters);

        const bool rebuild_level = ThisParameters["rebuild_level"].GetInt();
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

    ExplicitBuilderPointerType mpExplicitBuilder = nullptr;

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
}; /* Class NewExplicitSolvingStrategy */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos.*/

#endif /* KRATOS_EXPLICIT_SOLVING_STRATEGY  defined */
