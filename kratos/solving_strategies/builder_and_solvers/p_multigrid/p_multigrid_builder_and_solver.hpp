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

// STL Includes
#include <memory> // unique_ptr


namespace Kratos {



/** @page P-Multigrid
 *
 *  P-Multigrid refers to a multigrid method whose coarsening strategy is based on shape functions of high order elements (quadratic elements already qualify in this case). This implementation also supports models with multifreedom constraints.
 *
 *  @section Overview
 *
 *  The system's entry point is the @ref PMultigridBuilderAndSolver. As the name suggests, the main responsibilities of this class are
 *  - allocating the left hand side matrix, as well as the right hand side vector and the solution vector
 *  - assembling the left hand side matrix, right hand side vector
 *  - solving the resulting linear system of equations.
 *
 *  Due to the current design of @ref Scheme and @ref BuilderAndSolver, there is a number of other tasks as well that either tie in the primary purpose or help in pre- and postprocessing:
 *  - allocation, assembly and imposiotion of @ref MasterSlaveConstraint "multifreedom constraints"
 *  - applying Dirichlet conditions
 *  - partial reassembly of the linear system
 *  - computing reactions.
 *
 *  The other exposed family of classes from this system is @ref MultifreedomConstraint. Its main difference compared to
 *  @ref MasterSlaveConstraint (which it derives from) is that it only models the constraint equations and does not impose
 *  a partition of slave- and master @ref Dof "DoFs".
 *  @see LinearMultifreedomConstraint.
 *  @see LinkConstraint
 *
 *  @note The rest of the classes in this system are meant for internal use and must not be exposed to the rest of Kratos.
 *
 *  Although the primary purpose of this system is to exploit the structure arising from a model using high order elements,
 *  the multigrid feature can be completely disabled to use it as a standard @ref BuilderAndSolver. One reason for
 *  doing this would be taking advantage of different constraint imposition methods, such as
 *  @ref AugmentedLagrangeConstraintAssembler "augmented Lagrange".
 *
 *  @section Coarse Hierarchy
 *
 *  @note The current implementation only supports a two-grid method, since the selection of high order elements
 *        in Kratos is limited, and does not yet justify an arbitrary depth. That said, adding support for it
 *        should be possible without major interface changes, but would involve minor changes in
 *        @ref PMultigridBuilderAndSolver and @ref PGrid, as well as major changes in @ref MakePRestrictionOperator.
 *
 *  The root grid (i.e.: the finest level) is stored in and represented by @ref PMultigridBuilderAndSolver, while
 *  coarse grids are represented by @ref PGrid in a linked list. The reason for this difference is the additional
 *  set of responsibilities of the root grid, namely the allocation and assembly of the finest level. Coarse grids
 *  do not perform assembly, but construct restriction operators that they then apply on the parent grid to compute
 *  their own system.
 *
 *  Another important distinction is that the coarse grids can have floating point types different than the root grid.
 *  This can be useful when the user has access to accelerator hardware (i.e.: GPUs). Coarse grids need not solve their
 *  own problems with high precision so they might as well use single precision floating point numbers to save VRAM,
 *  for example.
 *
 *  @section Constraints
 *
 *  Constraint assembly and imposition is extracted through @ref ConstraintAssembler "a dedicated interface" that
 *  currently supports @ref MasterSlaveConstraintAssembler "master-slave elimination",
 *  @ref AugmentedLagrangeConstraintAssembler "augmented Lagrange", and @ref NoOpConstraintAssembler "a dummy"
 *  for debugging.
 *
 *  @note If the multigrid feature is enabled, the current implementation only supports augmented Lagrange imposition.
 *
 *  @section Linear Solvers
 *
 *  Unlike other @ref BuilderAndSolver "BuilderAndSolvers", @ref PMultigridBuilderAndSolver does not use the linear
 *  solver provided to it from the python layer. The reason is that it must construct a solver for each grid level
 *  separately. These fall into two categories
 *  - smoothers for the finer grids (including the root grid)
 *  - linear solver for the coarsest grid (usually an AMG solver).
 *
 *  Instead of passing linear solver instances, the user must provide two sets of parameters for the two different
 *  solver categories, after which the grids take care of constructing their own instances.
 */


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
 *                      "method" "master_slave_elimination"
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
         class TDenseSpace,
         class TLinearSolver>
class KRATOS_API(KRATOS_CORE) PMultigridBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>
{
private:
    using Interface = BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>;

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
    PMultigridBuilderAndSolver(const typename TLinearSolver::Pointer& pSolver,
                               Parameters Settings);

    PMultigridBuilderAndSolver(PMultigridBuilderAndSolver&& rOther) noexcept = default;

    PMultigridBuilderAndSolver& operator=(PMultigridBuilderAndSolver&& rRhs) noexcept = default;

    /// @brief Construct from a linear solver and parameters.
    /// @see PMultigridBuilderAndSolver::GetDefaultParameters
    typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer
    Create(typename TLinearSolver::Pointer pNewLinearSystemSolver,
           Parameters ThisParameters) const override;

    /// @name Allocation and Initialization
    /// @{

    /// @copydoc BuilderAndSolver::SetUpDofSet
    void SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme,
                     ModelPart& rModelPart) override;

    /// @brief Assign equation indices to @ref Dof "degrees of freedom".
    void SetUpSystem(ModelPart& rModelPart) override;

    /// @copydoc BuilderAndSolver::ResizeAndInitializeVectors
    void ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme,
                                    typename Interface::TSystemMatrixPointerType& rpLhs,
                                    typename Interface::TSystemVectorPointerType& rpSolution,
                                    typename Interface::TSystemVectorPointerType& rpRhs,
                                    ModelPart& rModelPart) override;

    /// @}
    /// @name Hooks
    /// @{

    void InitializeSolutionStep(ModelPart& rModelPart,
                                typename Interface::TSystemMatrixType& rLhs,
                                typename Interface::TSystemVectorType& rSolution,
                                typename Interface::TSystemVectorType& rRhs) override;

    void FinalizeSolutionStep(ModelPart& rModelPart,
                              typename Interface::TSystemMatrixType& rLhs,
                              typename Interface::TSystemVectorType& rSolution,
                              typename Interface::TSystemVectorType& rRhs) override;

    /// @name Assembly
    /// @{

    /// @copydoc BuilderAndSolver::BuildLHS
    void BuildLHS(typename Interface::TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  typename Interface::TSystemMatrixType& rLhs) override;

    /// @copydoc BuilderAndSolver::BuildRHS
    void BuildRHS(typename Interface::TSchemeType::Pointer pScheme,
                  ModelPart& rModelPart,
                  typename Interface::TSystemVectorType& rRhs) override;

    /// @copydoc BuilderAndSolver::Build
    void Build(typename Interface::TSchemeType::Pointer pScheme,
               ModelPart& rModelPart,
               typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Constraint Imposition
    /// @{

    /// @copydoc BuilderAndSolver::ApplyDirichletConditions
    void ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme,
                                  ModelPart& rModelPart,
                                  typename Interface::TSystemMatrixType& rLhs,
                                  typename Interface::TSystemVectorType& rSolution,
                                  typename Interface::TSystemVectorType& rRhs) override;

    /// @copydoc BuilderAndSolver::ApplyConstraints
    void ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                          ModelPart& rModelPart,
                          typename Interface::TSystemMatrixType& rLhs,
                          typename Interface::TSystemVectorType& rRhs) override;

    /// @}
    /// @name Compound Assembly and Solution
    /// @{

    /// @copydoc BuilderAndSolver::BuildAndSolve
    void BuildAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                       ModelPart& rModelPart,
                       typename Interface::TSystemMatrixType& A,
                       typename Interface::TSystemVectorType& Dx,
                       typename Interface::TSystemVectorType& b) override;

    /// @}
    /// @name Postprocessing
    /// @{

    /// @copydoc BuilderAndSolver::CalculateReactions
    void CalculateReactions(typename Interface::TSchemeType::Pointer pScheme,
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
    ///                  "method" : "master_slave_elimination"
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
    void ProjectGrid(int GridLevel,
                     const typename TSparseSpace::MatrixType& rRootLhs,
                     const typename TSparseSpace::VectorType& rRootSolution,
                     const typename TSparseSpace::VectorType& rRootRhs);

    /// @}
    /// @name Unsupported Interface
    /// @{

    /// @warning Not implemented.
    void SystemSolve(typename Interface::TSystemMatrixType& rLhs,
                     typename Interface::TSystemVectorType& rSolution,
                     typename Interface::TSystemVectorType& rRhs) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void BuildRHSAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                          ModelPart& rModelPart,
                          typename Interface::TSystemMatrixType& rLhs,
                          typename Interface::TSystemVectorType& rSolution,
                          typename Interface::TSystemVectorType& rRhs) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void ApplyRHSConstraints(typename Interface::TSchemeType::Pointer pScheme,
                             ModelPart& rModelPart,
                             typename Interface::TSystemVectorType& rRhs) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void BuildLHS_CompleteOnFreeRows(typename Interface::TSchemeType::Pointer pScheme,
                                     ModelPart& rModelPart,
                                     typename Interface::TSystemMatrixType& A) override
    {KRATOS_ERROR << KRATOS_CODE_LOCATION.CleanFunctionName() << " is not implemented\n";}

    /// @warning Not implemented.
    void BuildAndSolveLinearizedOnPreviousIteration(typename Interface::TSchemeType::Pointer pScheme,
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

    TLinearSolver& GetLinearSolver() noexcept;

    PMultigridBuilderAndSolver(const PMultigridBuilderAndSolver& rOther) = delete;

    PMultigridBuilderAndSolver& operator=(const PMultigridBuilderAndSolver& rRhs) = delete;

    struct Impl;
    std::unique_ptr<Impl,std::function<void(Impl*)>> mpImpl;
}; // class PMultigridBuilderAndSolver


} // namespace Kratos
