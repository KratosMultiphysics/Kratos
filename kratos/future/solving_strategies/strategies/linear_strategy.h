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
#include "includes/define.h"
#include "utilities/builtin_timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_solvers/linear_solver.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#endif

namespace Kratos::Future
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class LinearStrategy
 * @ingroup KratosCore
 * @brief This is a very simple strategy to solve linearly the problem
 * @details As a linear strategy the check on the convergence is not done and just one non linear iteration will be performed
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
//TODO: Add a base class
//FIXME: Remove the linear solver template
template<class TMatrixType, class TVectorType>
class LinearStrategy
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(LinearStrategy);

    /// The definition of the current class
    using SolvingStrategyType = LinearStrategy<TMatrixType, TVectorType>; //TODO: This will be the eventual base class

    /// The definition of the current class
    using ClassType = LinearStrategy<TMatrixType, TVectorType>;

    // Scheme type definition
    using SchemeType = ImplicitScheme<TMatrixType, TVectorType>;

    // Scheme pointer type definition
    using SchemePointerType = typename SchemeType::Pointer;

    // Linear solver pointer definition
    using LinearSolverPointerType = typename Future::LinearSolver<TMatrixType, TVectorType>::Pointer;

    /// DoF array type definition
    using DofsArrayType = typename SchemeType::DofsArrayType;

    /// Matrix type definition
    using SystemMatrixType = TMatrixType;

    /// Matrix pointer type definition
    using SystemMatrixPointerType = typename TMatrixType::Pointer;

    /// Vector type definition
    using SystemVectorType = TVectorType;

    /// Vector type definition
    using SystemVectorPointerType = typename SystemVectorType::Pointer;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit LinearStrategy() = default;

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit LinearStrategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: In here we should leverage the Registry to construct the scheme and linear solver from the json input settings
    }

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param CalculateNormDxFlag The flag sets if the norm of Dx is computed
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit LinearStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        LinearSolverPointerType pLinearSolver,
        bool ComputeReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool CalculateNormDxFlag = false,
        bool MoveMeshFlag = false)
        : mpModelPart(&rModelPart)
        , mpScheme(pScheme)
        , mpLinearSolver(pLinearSolver)
        , mReformDofsAtEachStep(ReformDofSetAtEachStep)
        , mComputeReactions(ComputeReactions)
    {
        KRATOS_TRY

        // Set EchoLevel to the default value (only time is displayed) //TODO: Get this from input settings
        this->SetEchoLevel(1);

        KRATOS_CATCH("")
    }

    /** Copy constructor.
     */
    LinearStrategy(const LinearStrategy &Other) = delete;

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    virtual ~LinearStrategy()
    {
        // Clear current class
        this->Clear();
    }

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
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize()
    {
        KRATOS_TRY

        // Initialize scheme (this has to be done once)
        pGetScheme()->Initialize();

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Boost dependencies should be replaced by std equivalent
     */
    void InitializeSolutionStep()
    {
        KRATOS_TRY

        // Call the scheme InitializeSolutionStep
        pGetScheme()->InitializeSolutionStep(mDofSet, mpA, mpdx, mpb, mReformDofsAtEachStep);

        KRATOS_CATCH("")
    }

    /**
     * @brief Operation to predict the solution
     * This provides a prediction for the current solution step solve. If not called a trivial predictor is used
     * in which the values of the solution step of interest are assumed equal to the old values.
     */
    void Predict()
    {
        KRATOS_TRY

        // Call the time scheme predict (note that this also updates the mesh if needed)
        pGetScheme()->Predict(*mpA, *mpdx, *mpb);

        //TODO: Constraints implementation
        // // Applying constraints if needed
        // const auto &r_comm = GetModelPart().GetCommunicator().GetDataCommunicator();
        // auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        // const int local_number_of_constraints = r_constraints_array.size();
        // const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        // if(global_number_of_constraints != 0) {
        //     const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.ResetSlaveDofs(r_process_info);
        //     });
        //     block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
        //         rConstraint.Apply(r_process_info);
        //     });

        //     //the following is needed since we need to eventually compute time derivatives after applying
        //     //Master slave relations
        //     TSparseSpace::SetToZero(rDx);
        //     this->pGetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        // }

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep()
    {
        // Get scheme pointer
        auto p_scheme = pGetScheme();

        // Initialize non-linear iteration (once as this is a linear strategy)
	    p_scheme->InitializeNonLinIteration(*mpA, *mpdx, *mpb);

        if (!mStiffnessMatrixIsBuilt) {
            // Initialize values
            mpA->SetValue(0.0);
            mpdx->SetValue(0.0);
            mpb->SetValue(0.0);

            // Build the local system and apply the Dirichlet conditions
            p_scheme->Build(*mpA, *mpb);
            p_scheme->ApplyDirichletConditions(mDofSet, *mpA, *mpb);
            mStiffnessMatrixIsBuilt = true;
        } else {
            // Initialize values
            mpdx->SetValue(0.0);
            mpb->SetValue(0.0);

            // Build the RHS and apply the Dirichlet conditions
            p_scheme->Build(*mpb);
            p_scheme->ApplyDirichletConditions(mDofSet, *mpb);
        }

        // Solve the system
        if(pGetLinearSolver()->AdditionalPhysicalDataIsNeeded()) {
            pGetLinearSolver()->ProvideAdditionalData(*mpA, *mpdx, *mpb, mDofSet, GetModelPart());
        }
        pGetLinearSolver()->Solve(*mpA, *mpdx, *mpb);

        // Debugging info
        EchoInfo();

        // Update results (note that this also updates the mesh if needed)
        p_scheme->Update(mDofSet, *mpA, *mpdx, *mpb);

        // Finalize current (unique) non linear iteration
        p_scheme->FinalizeNonLinIteration(*mpA, *mpdx, *mpb);

        // Calculate reactions if required
        if (mComputeReactions) {
            p_scheme->CalculateReactions(mDofSet, *mpb);
        }

        return true;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep()
    {
        KRATOS_TRY;

        // Finalisation of the solution step (operations to be done after achieving convergence)
        pGetScheme()->FinalizeSolutionStep(*mpA, *mpdx, *mpb);

        // Reset flags for next step
        mStiffnessMatrixIsBuilt = false; // By default we always rebuilt (if not implement a derived strategy)

        // Clear if needed (note that this deallocates the system arrays)
        if (mReformDofsAtEachStep) {
            Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear()
    {
        KRATOS_TRY;

        // Clear the linear solver
        if (mpLinearSolver != nullptr) {
            mpLinearSolver->Clear();
        }

        // Clearing the system of equations
        if (mpA != nullptr) {
            mpA->Clear();
        }
        if (mpdx != nullptr) {
            mpdx->Clear();
        }
        if (mpb != nullptr) {
            mpb->Clear();
        }

        // Clearing scheme
        // Note that this resets the DOF set
        if (pGetScheme() != nullptr) {
            pGetScheme()->Clear();
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check()
    {
        KRATOS_TRY

        // BaseType::Check(); //TODO: Activate when we have a base class

        pGetScheme()->Check();

        //TODO: think if we do this in debug mode only
        //TODO: this should probably go to the Scheme check (based on previous comment)
        // If reactions are to be calculated, we check if all the dofs have reactions defined
        if (mComputeReactions) {
            for (auto& r_dof : mDofSet) {
                KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
                    << "Reaction variable not set for the following: " << std::endl
                    << "- Node: "<< r_dof.Id() << std::endl
                    << "- DOF: " << r_dof << std::endl;
            }
        }

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step
     *@details This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    void CalculateOutputData()
    {
        pGetScheme()->CalculateOutputData(GetModelPart(), *mpA, *mpdx, *mpb);
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                         : "linear_strategy",
            "move_mesh"                    : false,
            "reform_dofs_at_each_step"     : false,
            "compute_reactions"            : false,
            "linear_solver_settings"       : {},
            "scheme_settings"              : {}
        })");

        // Getting base class default parameters
        // const Parameters base_default_parameters = BaseType::GetDefaultParameters(); //TODO: Once we have a base class...
        // default_parameters.RecursivelyAddMissingParameters(base_default_parameters); //TODO: Once we have a base class...
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "linear_strategy";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */

     void SetEchoLevel(const int EchoLevel) //TODO: I'd also unify the echo level with the scheme
     {
         mEchoLevel = EchoLevel;
         pGetScheme()->SetEchoLevel(EchoLevel);
     }

    /**
     * @brief This method sets the flag mComputeReactions
     * @param ComputeReactions The flag that tells if the reactions are computed
     */
    void SetComputeReactions(const bool ComputeReactions)
    {
        mComputeReactions = ComputeReactions;
    }

    /**
     * @brief This method sets the flag mReformDofsAtEachStep
     * @param ReformDofsAtEachStep The flag that indicates if the DOF set is reformed at each step
     */
    void SetReformDofsAtEachStep(const bool ReformDofsAtEachStep)
    {
        mReformDofsAtEachStep = ReformDofsAtEachStep;
    }

    //TODO: Move to the future base class
    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    };

    //TODO: Move to the future base class
    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    };

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(SchemePointerType pScheme)
    {
        mpScheme = pScheme;
    }

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    SchemePointerType pGetScheme()
    {
        return mpScheme;
    }

    /**
     * @brief Set method for the linear solver
     * @param pLinearSolver The pointer to the linear solver
     */
    void SetLinearSolver(LinearSolverPointerType pLinearSolver)
    {
        mpLinearSolver = pLinearSolver;
    }

    /**
     * @brief Get method for the linear solver
     * @return mpLinearSolver The pointer to the linear solver
     */
    LinearSolverPointerType pGetLinearSolver()
    {
        return mpLinearSolver;
    }

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    SystemMatrixType& GetSystemMatrix()
    {
        return *mpA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    SystemVectorType& GetSystemVector()
    {
        return *mpb;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    SystemVectorType& GetSolutionVector()
    {
        return *mpdx;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    DofsArrayType& GetDofSet()
    {
        return mDofSet;
    }

    /**
     * @brief It allows to get the list of Dofs from the element
     */
    const DofsArrayType& GetDofSet() const
    {
        return mDofSet;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * @brief This returns the level of echo for the solving strategy
     * @details
     * {
     * 0 -> Mute... no echo at all
     * 1 -> Printing time and basic information
     * 2 -> Printing linear solver data
     * 3 -> Print of debug information: Echo of stiffness matrix, Dx, b...
     * }
     * @return Level of echo for the solving strategy
     */
    int GetEchoLevel() const
    {
        return mEchoLevel;
    }

    /**
     * @brief This method returns the flag mComputeReactions
     * @return The flag that tells if the reactions are computed
     */
    bool GetComputeReactions() const
    {
        return mComputeReactions;
    }

    /**
     * @brief This method returns the flag mReformDofsAtEachStep
     * @return The flag that indicates if the DOF set is reformed at each step
     */
    bool GetReformDofsAtEachStep() const
    {
        return mReformDofsAtEachStep;
    }

    /**
     * @brief This method returns the residual norm
     * @return The residual norm
     */
    double GetResidualNorm()
    {
        if (mpb->size() != 0) {
            return mpb->Norm();
        } else {
            return 0.0;
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "LinearStrategy";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
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

    LinearSolverPointerType mpLinearSystemSolver = nullptr; /// Pointer to the linear solver

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters)
    {
        // BaseType::AssignSettings(ThisParameters); // TODO: Once we have a base class
        mEchoLevel = ThisParameters["echo_level"].GetInt();
        mComputeReactions = ThisParameters["compute_reactions"].GetBool();
        mReformDofsAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();

        // Saving the scheme
        if (ThisParameters["scheme_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }
    }

    ///@}
    ///@name Protected  Access
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

    SchemePointerType mpScheme = nullptr; /// The pointer to the scheme

    LinearSolverPointerType mpLinearSolver = nullptr; /// The pointer to the linear solver

    DofsArrayType mDofSet; /// The set containing the DoF of the system

    //TODO: Should these be unique_ptr?
    SystemVectorPointerType mpdx; /// The incremement in the solution //TODO: use naming convention

    //TODO: Should these be unique_ptr?
    SystemVectorPointerType mpb; /// The RHS vector of the system of equations //TODO: use naming convention (mpRHS)

    //TODO: Should these be unique_ptr?
    SystemMatrixPointerType mpA; /// The LHS matrix of the system of equations //TODO: use naming convention (mpLHS)

    int mEchoLevel;

    bool mReformDofsAtEachStep = false;

    bool mStiffnessMatrixIsBuilt = false;

    bool mComputeReactions; // Flag to activate the reactions calculation from the residual vector

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     */
    virtual void EchoInfo()
    {
        SystemMatrixType& rA  = *mpA;
        SystemVectorType& rDx = *mpdx;
        SystemVectorType& rb  = *mpb;

        if (GetEchoLevel() == 3) { //if it is needed to print the debug info
            KRATOS_INFO("LHS") << "LHS = " << rA << std::endl;
            KRATOS_INFO("RHS") << "RHS = " << rb << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
        }

        if (this->GetEchoLevel() == 4) { //print to matrix market file
            //FIXME: Impllement the MatrixMarket output

            // std::stringstream matrix_market_name;
            // matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] <<  ".mm";
            // TSparseSpace::WriteMatrixMarketMatrix((char*) (matrix_market_name.str()).c_str(), rA, false);

            // std::stringstream matrix_market_vectname;
            // matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << ".mm.rhs";
            // TSparseSpace::WriteMatrixMarketVector((char*) (matrix_market_vectname.str()).c_str(), rb);
        }
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
}; /* Class LinearStrategy */

///@}

///@name Type Definitions */
///@{


///@}

} // namespace Kratos::Future

