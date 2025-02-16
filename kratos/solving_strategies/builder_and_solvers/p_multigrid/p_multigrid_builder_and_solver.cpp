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

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.hpp" // PMultigridBuilderAndSolver
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler_factory.hpp" // ConstraintAssemblerFactory
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // MakeSparseTopology
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling, ParseDiagonalScaling, GetDiagonalScalingFactor
#include "solving_strategies/builder_and_solvers/p_multigrid/p_grid.hpp" // PGrid
#include "includes/model_part.h" // ModelPart
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "factories/linear_solver_factory.h" // LinearSolverFactory
#include "includes/kratos_components.h" // KratosComponents
#include "utilities/dof_utilities/block_build_dof_array_utility.h" // BlockBuildDofArrayUtility
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI

// System includes
#include <optional> // std::optional
#include <unordered_set> // std::unordered_set
#include <variant> // std::variant


namespace Kratos {


// --------------------------------------------------------- //
// PIMPL
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
struct PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Impl
{
    using Interface = PMultigridBuilderAndSolver<TSparse,TDense,TSolver>;

    // --------------------------------------------------------- //
    // Member Variables
    // --------------------------------------------------------- //

    Interface* mpInterface;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    typename TSparse::DataType mDiagonalScaleFactor;

    std::optional<std::variant<
        PGrid<TUblasSparseSpace<double>,TUblasDenseSpace<double>>,
        PGrid<TUblasSparseSpace<float>,TUblasDenseSpace<double>>
    >> mMaybeHierarchy;

    DiagonalScaling mDiagonalScaling;

    int mMaxIterations;

    typename TSparse::DataType mTolerance;

    int mMaxDepth;

    int mVerbosity;

    // --------------------------------------------------------- //
    // Special Member Functions
    // --------------------------------------------------------- //

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mpConstraintAssembler(),
          mDiagonalScaleFactor(1),
          mMaybeHierarchy(),
          mDiagonalScaling(DiagonalScaling::None),
          mMaxIterations(0),
          mTolerance(std::numeric_limits<typename TSparse::DataType>::max()),
          mMaxDepth(-1),
          mVerbosity(0)
    {}

    // --------------------------------------------------------- //
    // Solution
    // --------------------------------------------------------- //

    /// @brief Initialize the linear solver and solve the provided system.
    bool Solve(typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rSolution,
               typename Interface::TSystemVectorType& rRhs,
               ModelPart& rModelPart,
               Interface& rInterface)
    {
        KRATOS_TRY

        // Prepare and initialize members.
        if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
            rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                               rSolution,
                                                               rRhs,
                                                               rInterface.GetDofSet(),
                                                               rModelPart);
        }

        if (mMaybeHierarchy.has_value()) {
            std::visit([&rModelPart, &rLhs, &rSolution, &rRhs](auto& r_grid){
                            r_grid.template Initialize<TSparse>(
                                rModelPart,
                                rLhs,
                                rSolution,
                                rRhs);
                       },
                       mMaybeHierarchy.value());
        } // if mMaybeHierarchy

        KRATOS_CATCH("")

        bool linear_solver_status = false; //< Indicates whether the linear solver converged.
        typename ConstraintAssembler<TSparse,TDense>::Status constraint_status {/*converged=*/false, /*finished=*/false};

        KRATOS_TRY
        auto& r_linear_solver = rInterface.GetLinearSolver();
        std::size_t i_constraint_iteration = 0ul;

        typename TSparse::VectorType residual = rRhs,
                                     residual_update(rRhs.size()),
                                     solution_update(rSolution.size());
        TSparse::SetToZero(solution_update);

        // Compute the initial residual norm if it will be used.
        typename TSparse::DataType initial_residual_norm = 1;
        if (0 < mTolerance or 3 <= mVerbosity) {
            initial_residual_norm = TSparse::TwoNorm(residual);
            initial_residual_norm = initial_residual_norm ? initial_residual_norm : 1;
        }

        // Outer loop for constraints.
        do {
            ++i_constraint_iteration;
            mpConstraintAssembler->InitializeSolutionStep(rLhs, rSolution, rRhs, i_constraint_iteration);

            // Inner loop for multigrid.
            std::size_t i_multigrid_iteration = 0ul;

            while (++i_multigrid_iteration <= static_cast<std::size_t>(mMaxIterations)) {
                // Solve the coarse grid and apply its correction.
                if (mMaybeHierarchy.has_value()) {
                    std::visit([&solution_update, &residual](auto& r_hierarchy){
                                    return r_hierarchy.template ApplyCoarseCorrection<TSparse>(
                                        solution_update,
                                        residual);},
                               mMaybeHierarchy.value());

                    TSparse::UnaliasedAdd(rSolution, 1.0, solution_update);
                    TSparse::SetToZero(residual_update);
                    TSparse::Mult(rLhs, solution_update, residual_update);
                    TSparse::UnaliasedAdd(residual, -1.0, residual_update);
                } // if mMaybeHierarchy

                // Perform smoothing on the fine grid.
                TSparse::SetToZero(solution_update); //< do I need this?
                r_linear_solver.InitializeSolutionStep(rLhs, solution_update, residual);
                linear_solver_status = r_linear_solver.Solve(rLhs, solution_update, residual);
                r_linear_solver.FinalizeSolutionStep(rLhs, solution_update, residual);

                // Update the fine residual.
                TSparse::SetToZero(residual_update);
                TSparse::Mult(rLhs, solution_update, residual_update);
                TSparse::UnaliasedAdd(residual, -1.0, residual_update);

                // Emit status and check for convergence.
                if (0 < mTolerance or 3 <= mVerbosity) {
                    const auto relative_residual_norm = TSparse::TwoNorm(residual) / initial_residual_norm;
                    if (3 <= mVerbosity) {
                        std::stringstream residual_stream;
                        residual_stream << std::setprecision(8) << std::scientific<< relative_residual_norm;
                        std::cout << mpInterface->Info() << ": root grid: "
                                  << "multigrid iteration " << i_multigrid_iteration
                                  << " residual " << residual_stream.str()
                                  << "\n";
                    }
                    if (relative_residual_norm < mTolerance) break;
                } // 0 < mTolerance or 3 <= mVerbosity

                // Update the fine solution.
                TSparse::UnaliasedAdd(rSolution, 1.0, solution_update);
            } // while i_constraint_iteration < mMaxIterations

            constraint_status = mpConstraintAssembler->FinalizeSolutionStep(rLhs, rSolution, rRhs, i_constraint_iteration);
        } while (not constraint_status.finished);

        // Compute the final residual's norm.
        const typename TSparse::DataType final_residual_norm = TSparse::TwoNorm(residual) / initial_residual_norm;

        // Emit status.
        if (1 <= mVerbosity and (not linear_solver_status or not constraint_status.converged)) {
            std::stringstream residual_stream;
            residual_stream << std::setprecision(8) << std::scientific<< final_residual_norm;

            if (2 <= mVerbosity)
                std::cout << mpInterface->Info() << ": residual " << residual_stream.str() << "\n";
            else if (mTolerance < final_residual_norm)
                std::cerr << mpInterface->Info() << ": failed to converge "
                          << "(residual " << residual_stream.str() << ")\n";

            if (not constraint_status.converged)
                std::cerr << mpInterface->Info() << ": constraints failed to converge\n";
        }

        // Cleanup.
        mpConstraintAssembler->Finalize(rLhs,
                                        rSolution,
                                        rRhs,
                                        mpInterface->GetDofSet());

        if (mMaybeHierarchy.has_value()) {
            std::visit([&rModelPart, &rLhs, &rSolution, &rRhs](auto& r_grid){
                            r_grid.template Finalize<TSparse>(
                                rModelPart,
                                rLhs,
                                rSolution,
                                rRhs);
                       },
                       mMaybeHierarchy.value());
        } // if mMaybeHierarchy

        return final_residual_norm <= mTolerance and constraint_status.converged;
        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Topology
    // --------------------------------------------------------- //


    template <class TProxy>
    static void CollectDoFs(const TProxy& rEntities,
                            const ProcessInfo& rProcessInfo,
                            typename Interface::TSchemeType& rScheme,
                            LockObject* pLockBegin,
                            [[maybe_unused]] LockObject* pLockEnd,
                            std::unordered_set<std::size_t>* pRowSetBegin,
                            [[maybe_unused]] std::unordered_set<std::size_t>* pRowSetEnd)
    {
        KRATOS_TRY
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        using TLS = Element::EquationIdVectorType;
        block_for_each(rEntities,
                       TLS(),
                       [&rScheme, &rProcessInfo, pLockBegin, pRowSetBegin](const auto& rEntity, TLS& rTls){
            if (rEntity.GetEntity().IsActive()) {
                rScheme.EquationId(rEntity.GetEntity(), rTls, rProcessInfo);
                for (const auto equation_id : rTls) {
                    [[maybe_unused]] std::scoped_lock<LockObject> lock(pLockBegin[equation_id]);
                    auto& r_row_indices = pRowSetBegin[equation_id];
                    r_row_indices.insert(rTls.begin(), rTls.end());
                } // for equation_id in rTls
            } // if rEntity.IsActive
        });
        KRATOS_CATCH("")
    }


    void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme,
                         typename Interface::TSystemMatrixType& rLhs,
                         ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        std::vector<std::unordered_set<std::size_t> > indices(mpInterface->GetEquationSystemSize());

        {
            std::vector<LockObject> mutexes(mpInterface->GetEquationSystemSize());

            // Collect DoFs from elements.
            Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Element>(rModelPart),
                              rModelPart.GetProcessInfo(),
                              *rpScheme,
                              mutexes.data(),
                              mutexes.data() + mutexes.size(),
                              indices.data(),
                              indices.data() + indices.size());

            // Collect DoFs from conditions.
            Impl::CollectDoFs(MakeProxy<Globals::DataLocation::Condition>(rModelPart),
                              rModelPart.GetProcessInfo(),
                              *rpScheme,
                              mutexes.data(),
                              mutexes.data() + mutexes.size(),
                              indices.data(),
                              indices.data() + indices.size());
        }

        // Compute and allocate LHS topology.
        MakeSparseTopology<false,typename TSparse::DataType>(indices,
                                                             indices.size(),
                                                             rLhs);

        // Construct the coarse hierarhy's topology.
        if (mMaybeHierarchy.has_value()) {
            std::visit([&rModelPart, &rLhs, this](auto& r_grid){
                            r_grid.template MakeLhsTopology<TSparse>(rModelPart,
                                                                     rLhs,
                                                                     *mpConstraintAssembler,
                                                                     mpInterface->GetDofSet());
                       },
                       mMaybeHierarchy.value());
        } // if mMaybeHierarchy

        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Assembly
    // --------------------------------------------------------- //


    /// @brief Compute and assemble local contributions from elements and conditions into the global system.
    /// @details This function body mainly handles looping over elements/conditions and its parallelization.
    ///          The actual local system calculation, as well as assembly, is deferred to @ref MapEntityContribution.
    template <bool AssembleLHS,
              bool AssembleRHS>
    void Assemble(ModelPart& rModelPart,
                  typename Interface::TSchemeType& rScheme,
                  std::optional<typename Interface::TSystemMatrixType*> pMaybeLhs,
                  std::optional<typename Interface::TSystemVectorType*> pMaybeRhs)
    {
        KRATOS_TRY

        // Sanity checks.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pMaybeLhs.has_value() && pMaybeLhs.value() != nullptr)
                << "Requested an assembly of the left hand side, but no matrix is provided to assemble into.";
        }

        if constexpr (AssembleRHS) {
            KRATOS_ERROR_IF_NOT(pMaybeRhs.has_value() && pMaybeRhs.value() != nullptr)
                << "Requested an assembly of the right hand side, but no vector is provided to assemble into.";
        }

        // Assemble constraints.
        // Constraint assemble MUST happen before assembling the unconstrained system,
        // because constraints may have to propagate dirichlet conditions, and dirichlet
        // conditions are partly imposed during element assembly.
        mpConstraintAssembler->Assemble(
            rModelPart.MasterSlaveConstraints(),
            rModelPart.GetProcessInfo(),
            mpInterface->GetDofSet(),
            AssembleLHS,
            AssembleRHS);

        // Function-wide variables.
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        typename Interface::TSystemMatrixType* pLhs = pMaybeLhs.has_value() ? pMaybeLhs.value() : nullptr;
        typename Interface::TSystemVectorType* pRhs = pMaybeRhs.has_value() ? pMaybeRhs.value() : nullptr;
        auto p_locks = std::make_unique<std::vector<LockObject>>(AssembleLHS ? pLhs->size1() : 0ul);

        const int element_count = rModelPart.Elements().size();
        const int condition_count = rModelPart.Conditions().size();
        const auto it_element_begin = rModelPart.Elements().begin();
        const auto it_condition_begin = rModelPart.Conditions().begin();

        // Collect contributions from constitutive entities.
        #pragma omp parallel
        {
            // Thread-local variables.
            Element::EquationIdVectorType equation_indices;
            typename Interface::LocalSystemMatrixType lhs_contribution;
            typename Interface::LocalSystemVectorType rhs_contribution;

            // Collect contributions from elements.
            #pragma omp for schedule(guided, 512) nowait
            for (int i_entity=0; i_entity<element_count; ++i_entity) {
                MapEntityContribution<TSparse,TDense,AssembleLHS,AssembleRHS>(
                    *(it_element_begin + i_entity),
                    rScheme,
                    r_process_info,
                    equation_indices,
                    &lhs_contribution,
                    &rhs_contribution,
                    pLhs,
                    pRhs,
                    p_locks->data());
            } // pragma omp for

            // Collect contributions from conditions.
            #pragma omp for schedule(guided, 512)
            for (int i_entity=0; i_entity<condition_count; ++i_entity) {
                MapEntityContribution<TSparse,TDense,AssembleLHS,AssembleRHS>(
                    *(it_condition_begin + i_entity),
                    rScheme,
                    r_process_info,
                    equation_indices,
                    &lhs_contribution,
                    &rhs_contribution,
                    pLhs,
                    pRhs,
                    p_locks->data());
            } // pragma omp for
        } // pragma omp parallel

        // Assemble coarse hierarchy.
        if (this->mMaybeHierarchy.has_value()) {
            std::visit([&rModelPart, &pMaybeLhs, &pMaybeRhs, this](auto& r_grid){
                            r_grid.template Assemble<AssembleLHS,AssembleRHS,TSparse>(
                                rModelPart,
                                pMaybeLhs.has_value() ? pMaybeLhs.value() : nullptr,
                                pMaybeRhs.has_value() ? pMaybeRhs.value() : nullptr,
                                *this->mpConstraintAssembler);
                        },
                        this->mMaybeHierarchy.value());
        } // if mMaybeHierarchy

        KRATOS_CATCH("")
    }

}; // class PMultigridBuilderAndSolver::Impl


// --------------------------------------------------------- //
// Constructors
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::~PMultigridBuilderAndSolver() = default;


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver()
    : Interface(),
      mpImpl(new Impl(this))
{
}


template <class TSparse, class TDense, class TSolver>
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::PMultigridBuilderAndSolver(const typename TSolver::Pointer& pSolver,
                                                                               Parameters Settings)
    : Interface(pSolver),
      mpImpl(new Impl(this))
{
    Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());
    this->AssignSettings(Settings);
}


template <class TSparse, class TDense, class TSolver>
typename PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Interface::Pointer
PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Create(typename TSolver::Pointer pSolver,
                                                           Parameters Settings) const
{
    KRATOS_TRY
    return typename Interface::Pointer(new PMultigridBuilderAndSolver(pSolver, Settings));
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Allocation and Initialization
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme,
                                                                     ModelPart& rModelPart)
{
    KRATOS_TRY;
    BlockBuildDofArrayUtility::SetUpDofArray(rModelPart,
                                             this->GetDofSet(),
                                             this->GetEchoLevel(),
                                             this->GetCalculateReactionsFlag());
    KRATOS_CATCH("");
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::SetUpSystem(ModelPart& rModelPart)
{
    KRATOS_TRY
    this->mEquationSystemSize = this->GetDofSet().size();

    // Set equation indices of DoFs.
    IndexPartition<std::size_t>(this->GetDofSet().size()).for_each([&, this](std::size_t Index){
        (this->GetDofSet().begin() + Index)->SetEquationId(Index);
    });
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme,
                                                                                    typename Interface::TSystemMatrixPointerType& rpLhs,
                                                                                    typename Interface::TSystemVectorPointerType& rpSolution,
                                                                                    typename Interface::TSystemVectorPointerType& rpRhs,
                                                                                    ModelPart& rModelPart)
{
    KRATOS_TRY

    // Construct empty containers if necessary.
    if (!rpLhs)
        rpLhs.reset(new typename Interface::TSystemMatrixType);

    if (!rpSolution)
        rpSolution.reset(new typename Interface::TSystemVectorType);
    if (rpSolution->size() != this->mEquationSystemSize)
        rpSolution->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpSolution);

    if (!rpRhs)
        rpRhs.reset(new typename Interface::TSystemVectorType);
    if (rpRhs->size() != this->mEquationSystemSize)
        rpRhs->resize(this->mEquationSystemSize, false);
    TSparse::SetToZero(*rpRhs);

    // Construct LHS topology if necessary or requested.
    if (rpLhs->size1() == 0 || this->GetReshapeMatrixFlag() == true) {
        rpLhs->resize(this->mEquationSystemSize, this->mEquationSystemSize, false);
        mpImpl->MakeLhsTopology(pScheme, *rpLhs, rModelPart);

        // Make constraint topology.
        mpImpl->mpConstraintAssembler->Allocate(rModelPart.MasterSlaveConstraints(),
                                                rModelPart.GetProcessInfo(),
                                                *rpLhs,
                                                *rpSolution,
                                                *rpRhs,
                                                this->GetDofSet());
    } else {
        if (rpLhs->size1() != this->mEquationSystemSize || rpLhs->size2() != this->mEquationSystemSize) {
            KRATOS_ERROR <<"The equation system size has changed during the simulation. This is not permitted."<<std::endl;
        }
    }

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Hooks
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse, TDense, TSolver>::InitializeSolutionStep(ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    //
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Assembly
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Build(typename Interface::TSchemeType::Pointer pScheme,
                                                               ModelPart& rModelPart,
                                                               typename Interface::TSystemMatrixType& rLhs,
                                                               typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/true>(rModelPart,
                                                                         *pScheme,
                                                                         &rLhs,
                                                                         &rRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildLHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemMatrixType& rLhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    KRATOS_ERROR_IF(!pScheme) << "missing scheme";
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/false>(rModelPart,
                                                                          *pScheme,
                                                                          &rLhs,
                                                                          nullptr);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildRHS(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->GetDofSet(), [&rRhs](Dof<typename TDense::DataType>& rDof){
        if (rDof.IsFixed())
            rRhs[rDof.EquationId()] = 0.0;
    });

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Constraint Imposition
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme,
                                                                                  ModelPart& rModelPart,
                                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                                  [[maybe_unused]] typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    // Early exit on an empty system.
    if (this->GetDofSet().empty()) return;

    const auto it_dof_set_begin = this->GetDofSet().begin();
    const auto it_dof_set_end = this->GetDofSet().begin() + this->GetDofSet().size();

    mpImpl->mDiagonalScaleFactor = GetDiagonalScaleFactor<TSparse>(rLhs, mpImpl->mDiagonalScaling);
    Kratos::ApplyDirichletConditions<TSparse,TDense>(rLhs,
                                                     rRhs,
                                                     it_dof_set_begin,
                                                     it_dof_set_end,
                                                     mpImpl->mDiagonalScaleFactor);
    if (mpImpl->mMaybeHierarchy.has_value()) {
        std::visit([this](auto& r_grid){
                       r_grid.ApplyDirichletConditions(mpImpl->mDiagonalScaling);
                   },
                   mpImpl->mMaybeHierarchy.value());
    } // if mMaybeHierarchy
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                          ModelPart& rModelPart,
                                                                          typename Interface::TSystemMatrixType& rLhs,
                                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    mpImpl->mpConstraintAssembler->Initialize(rLhs, rRhs);
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Compound Assembly and Solution
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::BuildAndSolve(typename Interface::TSchemeType::Pointer pScheme,
                                                                       ModelPart& rModelPart,
                                                                       typename Interface::TSystemMatrixType& rLhs,
                                                                       typename Interface::TSystemVectorType& rSolution,
                                                                       typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY

    // Assemble unconstrained system.
    Build(pScheme, rModelPart, rLhs, rRhs);

    // Apply multifreedom constraints.
    ApplyConstraints(pScheme, rModelPart, rLhs, rRhs);

    // Apply Dirichlet conditions.
    ApplyDirichletConditions(pScheme, rModelPart, rLhs, rSolution, rRhs);

    // Solve constrained assembled system.
    if (not mpImpl->Solve(rLhs, rSolution, rRhs, rModelPart, *this) and 1 <= mpImpl->mVerbosity) {
        std::cerr << this->Info() << ": root grid: failed to solve the assembled system\n";
    }

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Postprocessing
// --------------------------------------------------------- //



template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::CalculateReactions(typename Interface::TSchemeType::Pointer pScheme,
                                                                            ModelPart& rModelPart,
                                                                            typename Interface::TSystemMatrixType& rLhs,
                                                                            typename Interface::TSystemVectorType& rSolution,
                                                                            typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    TSparse::SetToZero(rRhs);
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->GetDofSet(), [&rRhs](Dof<typename TDense::DataType>& rDof){
        rDof.GetSolutionStepReactionValue() = -rRhs[rDof.EquationId()];
    });
}


// --------------------------------------------------------- //
// Misc
// --------------------------------------------------------- //


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::AssignSettings(const Parameters Settings)
{
    KRATOS_TRY
    Parameters(Settings).ValidateAndAssignDefaults(this->GetDefaultParameters());
    KRATOS_CATCH("")

    KRATOS_TRY
    Interface::AssignSettings(Settings);
    KRATOS_CATCH("")

    // Set the scaling strategy for the diagonal entries of constrained DoFs.
    mpImpl->mDiagonalScaling = ParseDiagonalScaling(Settings);

    // Set the relative residual tolerance and max W cycles.
    mpImpl->mMaxIterations = Settings["max_iterations"].Get<int>();
    mpImpl->mTolerance = Settings["tolerance"].Get<double>();

    // Set multifreedom constraint imposition strategy.
    mpImpl->mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(Settings["constraint_imposition_settings"]);

    // Construct the top level smoother.
    Parameters smoother_settings = Settings["smoother_settings"];
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(smoother_settings.Has("solver_type"));
    const std::string solver_name = smoother_settings["solver_type"].Get<std::string>();
    using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;
    KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
        << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
        << "Make sure you imported the application it is defined in and that the spelling is correct.";
    const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
    this->mpLinearSystemSolver = r_factory.Create(smoother_settings);
    KRATOS_CATCH("")

    // Construct the coarse hierarchy.
    Parameters leaf_solver_settings = Settings["linear_solver_settings"];
    Parameters coarse_hierarchy_settings = Settings["coarse_hierarchy_settings"];
    std::string coarse_build_precision = coarse_hierarchy_settings.Has("precision") ?
                                         coarse_hierarchy_settings["precision"].Get<std::string>() :
                                         std::string {"double"};
    if (coarse_build_precision == "double") {
        using CoarseSparseSpace = TUblasSparseSpace<double>;
        using CoarseDenseSpace = TUblasDenseSpace<double>;
        using GridType = PGrid<CoarseSparseSpace,CoarseDenseSpace>;

        coarse_hierarchy_settings.ValidateAndAssignDefaults(GridType().GetDefaultParameters());
        if (0 < coarse_hierarchy_settings["max_depth"].Get<int>()) {
            mpImpl->mMaybeHierarchy = GridType(coarse_hierarchy_settings,
                                               smoother_settings,
                                               leaf_solver_settings);
        }
    } else if (coarse_build_precision == "single") {
        using CoarseSparseSpace = TUblasSparseSpace<float>;
        using CoarseDenseSpace = TUblasDenseSpace<double>;
        using GridType = PGrid<CoarseSparseSpace,CoarseDenseSpace>;

        coarse_hierarchy_settings.ValidateAndAssignDefaults(GridType().GetDefaultParameters());
        if (0 < coarse_hierarchy_settings["max_depth"].Get<int>()) {
            mpImpl->mMaybeHierarchy = GridType(coarse_hierarchy_settings,
                                               smoother_settings,
                                               leaf_solver_settings);
        }
    } else {
        KRATOS_ERROR << "unsupported coarse precision: \"" << coarse_build_precision << "\". "
                     << "Options are:\n"
                     << "\t\"single\""
                     << "\t\"double\"";
    }

    // Other settings.
    KRATOS_TRY
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Clear()
{
    Interface::Clear();
    mpImpl->mpConstraintAssembler->Clear();

    if (mpImpl->mMaybeHierarchy.has_value()) {
        std::visit([](auto& r_grid) {r_grid.Clear();},
                   mpImpl->mMaybeHierarchy.value());
    } // if mMaybeHierarchy
}


template <class TSparse, class TDense, class TSolver>
Parameters PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetDefaultParameters() const
{
    Parameters parameters = Parameters(R"({
"name"              : "p_multigrid",
"diagonal_scaling"  : "abs_max",
"verbosity"         : 1,
"max_iterations"    : 1e2,
"tolerance"         : 1e-8,
"smoother_settings" : {
    "solver_type" : "gauss_seidel"
},
"linear_solver_settings" : {
    "solver_type" : "amgcl"
},
"constraint_imposition_settings" : {
    "method" : "master_slave_elimination"
},
"coarse_hierarchy_settings" : {}
})");
    parameters.SetValue("coarse_hierarchy_settings", PGrid<TUblasSparseSpace<double>,TUblasDenseSpace<double>>().GetDefaultParameters());
    parameters.RecursivelyAddMissingParameters(Interface::GetDefaultParameters());
    return parameters;
}


template <class TSparse, class TDense, class TSolver>
std::string PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Info() const
{
    return "PMultigridBuilderAndSolver";
}


template <class TSparse, class TDense, class TSolver>
std::size_t PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetEquationSystemSize() const noexcept
{
    return Interface::mEquationSystemSize;
}

template <class TSparse, class TDense, class TSolver>
TSolver& PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetLinearSolver() noexcept
{
    return *Interface::mpLinearSystemSolver;
}


// --------------------------------------------------------- //
// Template Instantiations
// --------------------------------------------------------- //


template class PMultigridBuilderAndSolver<TUblasSparseSpace<double>,
                                          TUblasDenseSpace<double>,
                                          LinearSolver<TUblasSparseSpace<double>,
                                                       TUblasDenseSpace<double>>>;

template class PMultigridBuilderAndSolver<TUblasSparseSpace<float>,
                                          TUblasDenseSpace<double>,
                                          LinearSolver<TUblasSparseSpace<float>,
                                                       TUblasDenseSpace<double>>>;


} // namespace Kratos
