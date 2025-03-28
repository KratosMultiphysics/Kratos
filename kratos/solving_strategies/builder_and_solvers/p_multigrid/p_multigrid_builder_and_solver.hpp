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
 *  The system's entry point is the @ref PMultigridBuilderAndSolver, which is also the only exposed class to the rest of Kratos (other classes/functions/utilities should remain inside the P-Multigrid system and are not allowed to leak out). As the name suggests, the main responsibilities of this class are
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

template<class TSparseSpace,
         class TDenseSpace,
         class TLinearSolver>
class PMultigridBuilderAndSolver
    : public BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>
{
private:
    using Interface = BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>;

public:
    KRATOS_CLASS_POINTER_DEFINITION(PMultigridBuilderAndSolver);

    /// @details Required by PIMPL.
    ~PMultigridBuilderAndSolver();

    PMultigridBuilderAndSolver();

    /// @brief Construct from a linear solver and parameters.
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

    /// @details @code
    ///          {
    ///             "name"              : "p_multigrid",
    ///             "diagonal_scaling"  : "abs_max",
    ///             "max_depth"         : -1,
    ///             "verbosity"         : 0
    ///          }
    ///          @endcode
    Parameters GetDefaultParameters() const override;

    /// @copydoc BuilderAndSolver::Info
    std::string Info() const override;

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
    std::unique_ptr<Impl> mpImpl;
}; // class PMultigridBuilderAndSolver


} // namespace Kratos
