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
#include "includes/define.h" // KRATOS_CLASS_POINTER_DEFINITION
#include "solving_strategies/builder_and_solvers/builder_and_solver.h" // BuilderAndSolver
#include "includes/model_part.h" // ModelPart
#include "includes/code_location.h" // KRATOS_CODE_LOCATION

// STL Includes
#include <memory> // unique_ptr


namespace Kratos {

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
    typename Interface::Pointer Create(typename TLinearSolver::Pointer pNewLinearSystemSolver,
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
