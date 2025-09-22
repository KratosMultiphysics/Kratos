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
 * @class LinearStrategy
 * @ingroup KratosCore
 * @brief This is a very simple strategy to solve linearly the problem
 * @details As a linear strategy the check on the convergence is not done and just one non linear iteration will be performed
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
//TODO: Check if TSparseGraphType is really needed
template <class TMatrixType, class TVectorType, class TSparseGraphType>
class LinearStrategy : public ImplicitStrategy<TMatrixType, TVectorType, TSparseGraphType>
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(LinearStrategy);

    /// The definition of the current class
    using BaseType = ImplicitStrategy<TMatrixType, TVectorType, TSparseGraphType>;

    // Scheme pointer type definition
    using SchemePointerType = typename BaseType::SchemePointerType;

    // Linear solver pointer definition
    using LinearSolverPointerType = typename BaseType::LinearSolverPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit LinearStrategy() = default;

    //TODO: Create a Model - Parameters constructor

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit LinearStrategy(
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
        : BaseType(rModelPart, pScheme, pLinearSolver, ComputeReactions, ReformDofSetAtEachStep, CalculateNormDxFlag, MoveMeshFlag)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /** Copy constructor.
     */
    LinearStrategy(const LinearStrategy &Other) = delete;

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~LinearStrategy() override
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
        return Kratos::make_shared<LinearStrategy<TMatrixType, TVectorType, TSparseGraphType>>(rModelPart, ThisParameters);
    }

    bool SolveSolutionStep() override
    {
        // Get scheme pointer
        auto p_scheme = this->pGetScheme();

        // Get system data
        auto& r_dof_set = this->GetDofSet();
        auto& r_eff_dof_set = this->GetEffectiveDofSet();
        auto& r_linear_system_container = this->GetLinearSystemContainer();

        // Get system arrays
        auto p_dx = r_linear_system_container.pDx;
        auto p_rhs = r_linear_system_container.pRhs;
        auto p_lhs = r_linear_system_container.pLhs;
        auto p_constraints_T = r_linear_system_container.pConstraintsT;
        auto p_constraints_q = r_linear_system_container.pConstraintsQ;

        // Initialize non-linear iteration (once as this is a linear strategy)
        p_scheme->InitializeNonLinIteration(r_linear_system_container);

        if (!(this->GetStiffnessMatrixIsBuilt())) {
            // Initialize values
            p_lhs->SetValue(0.0);
            p_rhs->SetValue(0.0);
            p_dx->SetValue(0.0);

            // Build the local system and apply the Dirichlet conditions
            p_scheme->Build(*p_lhs, *p_rhs);
            p_scheme->BuildLinearSystemConstraints(r_linear_system_container);
            p_scheme->ApplyLinearSystemConstraints(r_linear_system_container);
            this->SetStiffnessMatrixIsBuilt(true);
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
        auto p_eff_dx = r_linear_system_container.pEffectiveDx;
        auto p_eff_rhs = r_linear_system_container.pEffectiveRhs;
        auto p_eff_lhs = r_linear_system_container.pEffectiveLhs;

        // Solve the system
        const auto& rp_linear_solver = this->pGetLinearSolver();
        if (rp_linear_solver->AdditionalPhysicalDataIsNeeded()) {
            rp_linear_solver->ProvideAdditionalData(*p_eff_lhs, *p_eff_dx, *p_eff_rhs, r_eff_dof_set, this->GetModelPart());
        }
        rp_linear_solver->Solve(*p_eff_lhs, *p_eff_dx, *p_eff_rhs);

        // Debugging info
        this->EchoInfo();

        // Update results (note that this also updates the mesh if needed)
        p_scheme->Update(r_linear_system_container);

        // Finalize current (unique) non linear iteration
        p_scheme->FinalizeNonLinIteration(r_linear_system_container);

        // Calculate reactions if required //TODO: Think on the constraints in here!!!
        if (this->GetComputeReactions()) {
            p_scheme->CalculateReactions(r_dof_set, *p_rhs);
        }

        //FIXME: Free the effective arrays memory if p_lhs != p_eff_lhs

        return true;
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
            "name" : "linear_strategy"
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
        return "LinearStrategy";
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
}; /* Class LinearStrategy */

///@}
///@name Type Definitions */
///@{

///@}
} // namespace Kratos::Future

