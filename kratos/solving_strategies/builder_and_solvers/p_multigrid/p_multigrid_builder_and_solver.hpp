//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#pragma once

// Project includes
#include "solving_strategies/builder_and_solvers/builder_and_solver.h" // BuilderAndSolver
#include "includes/smart_pointers.h" // KRATOS_CLASS_POINTER_DEFINITION
#include "includes/model_part.h" // ModelPart
#include "includes/code_location.h" // KRATOS_CODE_LOCATION
#include "linear_solvers/linear_solver.h" // LinearSolver

// STL Includes
#include <memory> // unique_ptr


namespace Kratos {


///@name Kratos Classes
///@{

/** @brief Class used as a standard @ref BuilderAndSolver with a built-in optional preconditioner for high-order elements. See @ref P-Multigrid.
 *
 *  @details
 *  @section Structure
 *           PMultigridBuilderAndSolver consists of 3 main components:
 *           - The top-level grid, which is responsible for assembling the main linear system
 *             of equations from @ref Element "elements" and @ref Condition "conditions".
 *           - The top-level @ref ConstraintAssembler "constraint assembler" responsible for assembling
 *             and imposing constraint equations defined by @ref MasterSlaveConstraint "multifreedom constraints".
 *           - A hierarchy of coarser grid levels, each with its own @ref LinearSolver "linear solver" and
 *             @ref ConstraintAssembler "constraint assembler".
 *
 *          @copydoc PMultigridBuilderAndSolver::GetDefaultParameters
 *
 *  @subsection Top-Level-Grid Top-Level Grid
 *              Like any other grid level, the fine grid has its own linear solver and constraint assembler, but
 *              they can be configured separately from coarse grid levels.
 *              Settings:
 *              @code
 *              {
 *                  "name"              : "p_multigrid",
 *                  "diagonal_scaling"  : "max",
 *                  "max_iterations"    : 1e2,
 *                  "tolerance"         : 1e-8,
 *                  "verbosity"         : 1,
 *                  ...
 *              }
 *              @endcode
 *              - @p "name" Name of this builder-and-solver to refer to in @ref PythonSolver solvers.
 *              - @p "diagonal_scaling" Expression to scale diagonal entries of the LHS that are
 *                                      related to @ref Dof "DoFs" constrained by Dirichlet conditions.
 *                                      It can either be a numeric literal, or a string that defines
 *                                      the scaling factor as a function of:
 *                  - @p "max" Infinity norm (maximum absolute value) of the LHS' main diagonal.
 *                  - @p "norm" 2-norm of the LHS' main diagonal.
 *              - @p "max_iterations" Max number of multigrid iterations to carry out while attempting
 *                                      to reduce the residual to the requested value.
 *              - @p "tolerance" Target relative residual norm to achieve.
 *              - @p "verbosity" Level of verbosity. Every level includes lower verbosity levels and
 *                                 adds new events to report:
 *                  - @p 0 No messages, not even in case of failure.
 *                  - @p 1 Warnings and failure messages.
 *                  - @p 2 Aggregated status reports.
 *                  - @p 3 Per-iteration status reports.
 *                  - @p 4 Output system matrices and vectors.
 *                  - @p 5 Write a VTU output with the solution and residuals at each iteration.
 *
 *  @subsection Top-Level-Constraints Top-Level Constraint Assembler
 *              The fine-level constraint assembler is responsible for building and imposing
 *              constraint equations defined by @ref MasterSlaveConstraint "multifreedom constraints" of
 *              the input @ref ModelPart "model part".
 *              Settings:
 *              @code
 *              {
 *                  ...
 *                  "constraint_imposition_settings" {
 *                      "method" "master_slave"
 *                  }
 *                  ...
 *              }
 *              @endcode
 *              - @p "constraint_imposition_settings" Settings of the top-level constraint assembler.
 *                                                    Refer to @ref ConstraintAssemblerFactory for more
 *                                                    information.
 *
 *  @subsection Coarse-Hierarchy Coarse Hierarchy
 *              Coarse grid levels have their own linear solvers as well as constraint assemblers. However,
 *              all of them share the same configuration (except the coarses grid's linear solver, which is
 *              defined by @p "linear_solver_settings" ). Coarse grids and their constraint assemblers do
 *              not actually construct their systems by assembly, but by applying their restriction operators
 *              to their parent grid's matrices/vectors. Refer to @ref PGrid for more information.
 *              Default settings:
 *              @code
 *              {
 *                  ...
 *                  "coarse_hierarchy_settings" {}
 *                  ...
 *              }
 *              @endcode
 *              - @p "coarse_hierarchy_settings" Coarse grid hierarchy configuration.
 *
 *              @copydoc PGrid
 *
 *  @subsection Solvers-Smoother Solvers and Smoothers
 *              PMultigridBuilderAndSolver uses two linear solvers in general. A smoother defined by
 *              @p "smoother_settings" for all grid levels except the coarses, and a proper linear solver
 *              on the coarsest grid, defined by @p "linear_solver_settings".
 *              Note that if @p "coarse_hierarchy_settings.max_depth" is @p 0, the coarsest grid level is
 *              the top-level one, in which case the smoother is not used. Both solvers are constructed by
 *              @ref LinearSolverFactory, so check there for more information.
 *              Settings:
 *              @code
 *              {
 *                  ...
 *                  "smoother_settings" {
 *                      "solver_type" : ""
 *                  },
 *                  "linear_solver_settings" : {
 *                      "solver_type" : "amgcl"
 *                  }
 *              }
 *              @endcode
 */
template<class TSparseSpace,
         class TDenseSpace>
class KRATOS_API(KRATOS_CORE) PMultigridBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace,TDenseSpace,LinearSolver<TSparseSpace,TDenseSpace>>
{
private:
    using LinearSolverType = LinearSolver<TSparseSpace,TDenseSpace>;

    using Interface = BuilderAndSolver<TSparseSpace,TDenseSpace,LinearSolverType>;

public:
    /// @internal
    KRATOS_CLASS_POINTER_DEFINITION(PMultigridBuilderAndSolver);

    /// @details Required by PIMPL.
    /// @internal
    ~PMultigridBuilderAndSolver() override;

    PMultigridBuilderAndSolver();

    /// @brief Construct from a linear solver and parameters.
    /// @note The provided linear solver is not used. Check the input settings
    ///       on how to configure linear solvers for this class.
    /// @see PMultigridBuilderAndSolver::GetDefaultParameters
    PMultigridBuilderAndSolver(
        const typename LinearSolverType::Pointer& pSolver,
        Parameters Settings);

    PMultigridBuilderAndSolver(PMultigridBuilderAndSolver&& rOther) noexcept = default;

    PMultigridBuilderAndSolver& operator=(PMultigridBuilderAndSolver&& rRhs) noexcept = default;

    /// @brief Construct from a linear solver and parameters.
    /// @see PMultigridBuilderAndSolver::GetDefaultParameters
    typename BuilderAndSolver<TSparseSpace,TDenseSpace,LinearSolverType>::Pointer
    Create(typename LinearSolverType::Pointer pNewLinearSystemSolver,
           Parameters ThisParameters) const override;

    /// @name Allocation and Initialization
    /// @{

    /// @copydoc BuilderAndSolver::SetUpDofSet
    void SetUpDofSet(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart) override;

    /// @brief Assign equation indices to @ref Dof "degrees of freedom".
    void SetUpSystem(ModelPart& rModelPart) override;

    /// @copydoc BuilderAndSolver::ResizeAndInitializeVectors
    void ResizeAndInitializeVectors(
        typename Interface::TSchemeType::Pointer pScheme,
        typename Interface::TSystemMatrixPointerType& rpLhs,
        typename Interface::TSystemVectorPointerType& rpSolution,
        typename Interface::TSystemVectorPointerType& rpRhs,
        ModelPart& rModelPart) override;

    /// @}
    /// @name Hooks
    /// @{

    void InitializeSolutionStep(
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Assembly
    /// @{

    /// @copydoc BuilderAndSolver::BuildLHS
    void BuildLHS(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs) override;

    /// @copydoc BuilderAndSolver::BuildRHS
    void BuildRHS(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemVectorType& rRhs) override;

    void BuildRHSAndSolve(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @copydoc BuilderAndSolver::Build
    void Build(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Constraint Imposition
    /// @{

    /// @copydoc BuilderAndSolver::ApplyDirichletConditions
    void ApplyDirichletConditions(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @copydoc BuilderAndSolver::ApplyConstraints
    void ApplyConstraints(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Solution
    /// @{

    /// @copydoc BuilderAndSolver::BuildAndSolve
    void BuildAndSolve(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& A,
        typename Interface::TSystemVectorType& Dx,
        typename Interface::TSystemVectorType& b) override;

    /// @copydoc BuilderAndSolver::SystemSolve(typename Interface::TSystemMatrixType&, typename Interface::TSystemVectorType&, typename Interface::TSystemVectorType&)
    void SystemSolve(
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Postprocessing
    /// @{

    /// @copydoc BuilderAndSolver::CalculateReactions
    void CalculateReactions(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Misc
    /// @{

    /// @copydoc BuilderAndSolver::Clear
    void Clear() override;

    /// @brief Default configuration.
    /// @details @code
    ///          {
    ///              "name"              : "p_multigrid",
    ///              "diagonal_scaling"  : "max",
    ///              "max_iterations"    : 1e2,
    ///              "tolerance"         : 1e-8,
    ///              "verbosity"         : 1,
    ///              "constraint_imposition_settings" : {
    ///                  "method" : "master_slave"
    ///              },
    ///              "coarse_hierarchy_settings" : {}
    ///              "smoother_settings" : {
    ///                  "solver_type" : ""
    ///              },
    ///              "linear_solver_settings" : {
    ///                  "solver_type" : "amgcl"
    ///              },
    ///          }
    ///          @endcode
    Parameters GetDefaultParameters() const override;

    /// @copydoc BuilderAndSolver::Info
    std::string Info() const override;

    /// @}
    /// @name Debug and Analysis
    /// @{

    /// @brief Project state variables and residuals from a coarse grid to the root grid.
    /// @param GridLevel Grid to project.
    /// @param rRootLhs LHS matrix of the root grid.
    /// @param rRootSolution Current state vector of the root grid.
    /// @param rRootRhs RHS vector of the root grid.
    /// @details State variables are stored in the corresponding @ref Dof "DoFs"' values,
    ///          while residuals are written to reactions.
    /// @warning This function changes DoFs' values and reactions, so they should be copied
    ///          before invoking this function and restored afterwards.
    void ProjectGrid(
        int GridLevel,
        const typename TSparseSpace::MatrixType& rRootLhs,
        const typename TSparseSpace::VectorType& rRootSolution,
        const typename TSparseSpace::VectorType& rRootRhs);

    /// @}
    /// @name Unsupported Interface
    /// @{

    /// @warning Not implemented.
    void ApplyRHSConstraints(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemVectorType& rRhs) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void BuildLHS_CompleteOnFreeRows(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& A) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void BuildAndSolveLinearizedOnPreviousIteration(
        typename Interface::TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        typename Interface::TSystemMatrixType& rLhs,
        typename Interface::TSystemVectorType& rSolution,
        typename Interface::TSystemVectorType& rRhs,
        const bool MoveMesh) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @}

protected:
    void AssignSettings(const Parameters Settings) override;

private:
    std::size_t GetEquationSystemSize() const noexcept;

    LinearSolverType& GetRootGridSolver() noexcept;

    PMultigridBuilderAndSolver(const PMultigridBuilderAndSolver& rOther) = delete;

    PMultigridBuilderAndSolver& operator=(const PMultigridBuilderAndSolver& rRhs) = delete;

    /// @details MSVC desperately wants to know about the destructor of @p Impl when trying to
    ///          use PIMPL with an @p std::unique_ptr, up to the point that it ICE-s if I try
    ///          dodging it with a lambda. Hopefully a raw pointer will get the message
    ///          across its thick skull that it does not need to know what the destructor
    ///          does if it only ever sees a bloody pointer.
    /// @todo Change @p Impl* to @p std::unique_ptr<Impl> when (if) MSVC sorts its shit out.
    struct Impl;
    Impl* mpImpl;
}; // class PMultigridBuilderAndSolver


} // namespace Kratos
