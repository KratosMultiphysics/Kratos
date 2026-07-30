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
#include "utilities/builtin_timer.h"

#ifdef KRATOS_USE_FUTURE
#include "future/linear_operators/sparse_matrix_linear_operator.h"
#include "future/linear_solvers/linear_solver.h"
#include "future/solving_strategies/convergence_criteria/convergence_criteria.h"
#include "future/solving_strategies/schemes/implicit_scheme.h"
#include "future/solving_strategies/strategies/implicit_strategy.h"
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
 * @class NewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is a Newton-Raphson strategy for the solution of non-linear problems.
 * @details //TODO: Add more details about the class
 * @tparam TLinearAlgebra The linear algebra type
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
template <class TLinearAlgebra>
class NewtonRaphsonStrategy : public ImplicitStrategy<TLinearAlgebra>
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of NewtonRaphsonStrategy */
    KRATOS_CLASS_POINTER_DEFINITION(NewtonRaphsonStrategy);

    /// The definition of the current class base type
    using BaseType = ImplicitStrategy<TLinearAlgebra>;

    // Scheme pointer type definition
    using SchemePointerType = typename BaseType::SchemePointerType;

    // Linear solver pointer definition
    using LinearSolverPointerType = typename BaseType::LinearSolverPointerType;

    // Convergence criteria pointer definition
    using ConvergenceCriteriaPointerType = typename BaseType::ConvergenceCriteriaPointerType;

    /// Dense vector tag type definition
    using DenseVectorTag = typename LinearSystemTags::DenseVectorTag;

    /// Sparse matrix tag type definition
    using SparseMatrixTag = typename LinearSystemTags::SparseMatrixTag;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit NewtonRaphsonStrategy() = default;

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit NewtonRaphsonStrategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : BaseType(rModelPart, ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: In here we should leverage the Registry to construct the scheme and linear solver from the json input settings
    }

    /**
     * @brief Model - Parameters constructor
     * @param rModel The model container of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit NewtonRaphsonStrategy(
        Model &rModel,
        Parameters ThisParameters)
        : BaseType(rModel, ThisParameters)
    {
        // Validate and assign defaults
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        //TODO: In here we should leverage the Registry to construct the scheme and linear solver from the json input settings
    }

    //TODO: Temporary constructor to be used until we fix the registry stuff
    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactionFlag The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit NewtonRaphsonStrategy(
        ModelPart& rModelPart,
        SchemePointerType pScheme,
        LinearSolverPointerType pLinearSolver,
        ConvergenceCriteriaPointerType pConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : BaseType(rModelPart, pScheme, pLinearSolver, pConvergenceCriteria, CalculateReactions, ReformDofSetAtEachStep, false, MoveMeshFlag)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /** Copy constructor.
     */
    NewtonRaphsonStrategy(const NewtonRaphsonStrategy &Other) = delete;

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~NewtonRaphsonStrategy() override
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

    typename Strategy::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters) const override
    {
        return Kratos::make_shared<NewtonRaphsonStrategy<TLinearAlgebra>>(rModelPart, ThisParameters);
    }

    void Initialize() override
    {
        // Call base class Initialize
        // Note that this calls the Initialize of the scheme
        BaseType::Initialize();

        // Initialize convergence criteria
        this->pGetConvergenceCriteria()->Initialize(this->GetImplicitStrategyData());
    }

    void InitializeSolutionStep() override
    {
        // Call base class InitializeSolutionStep
        // Note that this calls the InitializeSolutionStep of the scheme
        BaseType::InitializeSolutionStep();

        // Calculate residual in case the convergence criteria is residual-based
        if (this->pGetConvergenceCriteria()->RequiresResidual()) {
            // Get scheme pointer
            auto p_scheme = this->pGetScheme();

            // Initialize residual vector
            auto p_linear_system = this->GetImplicitStrategyData().pGetLinearSystem();
            auto& r_rhs = *(p_linear_system->pGetVector(DenseVectorTag::RHS));
            r_rhs.SetValue(0.0);

            // Build the residual and apply Dirichlet conditions and constraints to it
            p_scheme->Build(r_rhs);
            p_scheme->BuildLinearSystemConstraints(this->GetImplicitStrategyData());
            p_scheme->ApplyLinearSystemConstraints(this->GetImplicitStrategyData(), true); // The true flag skips the LHS application
        }

        // Initialize convergence criteria step
        this->pGetConvergenceCriteria()->InitializeSolutionStep(this->GetImplicitStrategyData());
    }

    bool SolveSolutionStep() override
    {
        // Get scheme and convergence criteria pointers
        auto p_scheme = this->pGetScheme();
        auto p_conv_crit = this->pGetConvergenceCriteria();

        // Get system data
        auto& r_dof_set = this->GetDofSet();
        auto& r_eff_dof_set = this->GetEffectiveDofSet();
        auto& r_strategy_data_container = this->GetImplicitStrategyData();

        // Get system arrays
        auto p_linear_system = r_strategy_data_container.pGetLinearSystem();
        auto& r_dx = *(p_linear_system->pGetVector(DenseVectorTag::Dx));
        auto& r_lhs = *(p_linear_system->pGetMatrix(SparseMatrixTag::LHS));
        auto& r_rhs = *(p_linear_system->pGetVector(DenseVectorTag::RHS));
        auto p_constraints_T = r_strategy_data_container.pGetConstraintsT();
        auto p_constraints_q = r_strategy_data_container.pGetConstraintsQ();

        // Build the linear system constraints
        // Note that the constraints are built once as they are not expected to change during the solution step
        p_scheme->BuildLinearSystemConstraints(r_strategy_data_container);

        // Newton-Raphson cycle
        bool is_converged = false;
        unsigned int iteration_number = 0;
        while (!is_converged && iteration_number++ <= mMaxIteration) {

            // Initialize current iteration
            this->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            p_scheme->InitializeNonLinIteration(r_strategy_data_container);
            p_conv_crit->InitializeNonLinearIteration(r_strategy_data_container);

            // Build the system
            if (!(this->GetStiffnessMatrixIsBuilt())) {
                // Initialize values
                r_lhs.SetValue(0.0);
                r_rhs.SetValue(0.0);
                r_dx.SetValue(0.0);

                // Build the local system and apply the Dirichlet conditions
                if (iteration_number == 1 && mUseOldStiffnessInFirstIteration) {
                    KRATOS_ERROR << "The option 'use_old_stiffness_in_first_iteration' is not yet implemented in the new strategy." << std::endl;
                } else {
                    p_scheme->Build(r_lhs, r_rhs);
                    // p_scheme->BuildLinearSystemConstraints(r_strategy_data_container);
                    p_scheme->ApplyLinearSystemConstraints(r_strategy_data_container);
                    this->SetStiffnessMatrixIsBuilt(true);
                }
            } else {
                //FIXME: Do the RHS-only one!!!!

                // // Initialize values
                // p_rhs->SetValue(0.0);
                // p_dx->SetValue(0.0);

                // // Build the RHS and apply the Dirichlet conditions
                // p_scheme->Build(*p_rhs);
                // p_scheme->ApplyMasterSlaveConstraints(p_rhs, p_eff_rhs, p_dx, p_eff_dx, r_T, r_b);
                // p_scheme->ApplyDirichletConditions(r_eff_dof_set, r_eff_dof_map, *p_eff_rhs);
            }

            // Get the effective arrays to solve the system
            auto p_eff_lin_sys = r_strategy_data_container.pGetEffectiveLinearSystem();
            auto p_eff_dx = p_eff_lin_sys->pGetVector(DenseVectorTag::Dx);
            auto p_eff_lhs = p_eff_lin_sys->pGetMatrix(SparseMatrixTag::LHS);
            auto p_eff_rhs = p_eff_lin_sys->pGetVector(DenseVectorTag::RHS);
            auto p_eff_lhs_lin_op = Kratos::make_shared<SparseMatrixLinearOperator<TLinearAlgebra>>(p_eff_lhs);

            // Solve the system
            const auto& rp_linear_solver = this->pGetLinearSolver();
            if (rp_linear_solver->RequiresAdditionalData()) {
                p_eff_lin_sys->SetAdditionalData(this->GetModelPart(), r_eff_dof_set);
                rp_linear_solver->PrepareAdditionalData(*p_eff_lin_sys);
            }
            rp_linear_solver->Initialize(*p_eff_lin_sys); //TODO: Check if this is needed at each iteration (only if the effective linear system changes the size, e.g., changing the Dirichlet BCs in elimination or the constraints set)
            rp_linear_solver->InitializeSolutionStep(*p_eff_lin_sys); //TODO: Check if this is needed at each iteration (only if the effective linear system changes the size, e.g., changing the Dirichlet BCs in elimination or the constraints set)
            rp_linear_solver->PerformSolutionStep(*p_eff_lin_sys);
            rp_linear_solver->FinalizeSolutionStep(*p_eff_lin_sys); //TODO: Check if this is needed at each iteration (only if the effective linear system changes the size, e.g., changing the Dirichlet BCs in elimination or the constraints set)

            // Debugging info
            this->EchoInfo();

            // Update results (note that this also updates the mesh if needed)
            p_scheme->Update(r_strategy_data_container);

            // Finalize current (unique) non linear iteration
            p_scheme->FinalizeNonLinIteration(r_strategy_data_container);
            p_conv_crit->FinalizeNonLinearIteration(r_strategy_data_container);

            // Check convergence
            // Note that the residual is computed again with current solution if the convergence criteria requires it
            if (p_conv_crit->RequiresResidual()) {

                r_rhs.SetValue(0.0);
                p_scheme->Build(r_rhs);
                // p_scheme->BuildLinearSystemConstraints(r_strategy_data_container);
                p_scheme->ApplyLinearSystemConstraints(r_strategy_data_container, true); // The true flag skips the LHS application
            }
            is_converged = p_conv_crit->IsConverged(r_strategy_data_container);

            // Set the rebuilt stifness matrix flag
            if (iteration_number == 1 && mUseOldStiffnessInFirstIteration) {
                this->SetStiffnessMatrixIsBuilt(false);
            } else if (mUpdateStiffnessAtEachIteration) {
                this->SetStiffnessMatrixIsBuilt(false);
            } else {
                this->SetStiffnessMatrixIsBuilt(true);
            }
        }

        // Plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIteration) {
            KRATOS_WARNING_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Maximum iterations (" << mMaxIteration << ") exceeded" << std::endl;
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategy", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / " << mMaxIteration << " iterations" << std::endl;
        }

        // Calculate reactions if required //TODO: Think on the constraints in here!!!
        if (this->GetComputeReactions()) {
            p_scheme->CalculateReactions(r_dof_set, r_rhs);
        }

        //FIXME: Free the effective arrays memory if p_lhs != p_eff_lhs

        return is_converged;
    }

    int Check() override
    {
        KRATOS_TRY

        int check = BaseType::Check();

        return check;

        KRATOS_CATCH("")
    }

    Parameters GetDefaultParameters() const override
    {
        // Current class default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "newton_raphson_strategy",
            "max_iteration" : 10,
            "use_old_stiffness_in_first_iteration" : false,
        })");

        // Add base class default parameters
        default_parameters.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());

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

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    std::string Info() const override
    {
        return "NewtonRaphsonStrategy";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    void PrintData(std::ostream& rOStream) const override
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
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        // Call base class AssignSettings
        BaseType::AssignSettings(ThisParameters);

        // Assign current Newton-Raphson settings
        mMaxIteration = ThisParameters["max_iteration"].GetInt();
        mUpdateStiffnessAtEachIteration = ThisParameters["update_stiffness_at_each_iteration"].GetBool();
        mUseOldStiffnessInFirstIteration = ThisParameters["use_old_stiffness_in_first_iteration"].GetBool();
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

    unsigned int mMaxIteration; // Maximum number of iterations

    bool mUpdateStiffnessAtEachIteration = true; // Indicates if the stiffness matrix is updated at each iteration

    bool mUseOldStiffnessInFirstIteration = false; // Indicates if the first iteration is linearized on the previous time step solution

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
}; /* Class NewtonRaphsonStrategy */

///@}
///@name Type Definitions */
///@{

///@}
} // namespace Kratos::Future

