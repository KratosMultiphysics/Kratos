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
#include "solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp" // PMGStatusStream
#include "includes/model_part.h" // ModelPart
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "factories/linear_solver_factory.h" // LinearSolverFactory
#include "includes/kratos_components.h" // KratosComponents
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI
#include "utilities/builtin_timer.h" // BuiltinTimer

// System includes
#include <optional> // std::optional
#include <unordered_set> // std::unordered_set
#include <variant> // std::variant


namespace Kratos {


// --------------------------------------------------------- //
// PIMPL
// --------------------------------------------------------- //


template <class TSparse, class TDense>
struct PMultigridBuilderAndSolver<TSparse,TDense>::Impl
{
    using Interface = PMultigridBuilderAndSolver<TSparse,TDense>;

    // --------------------------------------------------------- //
    // Member Variables
    // --------------------------------------------------------- //

    Interface* mpInterface;

    std::optional<ModelPart*> mpMaybeModelPart;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    std::optional<std::variant<
        PGrid<TUblasSparseSpace<double>,TUblasDenseSpace<double>>,
        PGrid<TUblasSparseSpace<float>,TUblasDenseSpace<double>>
    >> mMaybeHierarchy;

    std::unique_ptr<Scaling> mpDiagonalScaling;

    int mMaxIterations;

    typename TSparse::DataType mTolerance;

    int mMaxDepth;

    int mVerbosity;

    // --------------------------------------------------------- //
    // Special Member Functions
    // --------------------------------------------------------- //

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mpMaybeModelPart(),
          mpConstraintAssembler(),
          mMaybeHierarchy(),
          mpDiagonalScaling(),
          mMaxIterations(0),
          mTolerance(std::numeric_limits<typename TSparse::DataType>::max()),
          mMaxDepth(-1),
          mVerbosity(0)
    {}

    // --------------------------------------------------------- //
    // Solution
    // --------------------------------------------------------- //


    void ExecuteMultigridLoop(typename Interface::TSystemMatrixType& rLhs,
                              typename Interface::TSystemVectorType& rSolution,
                              const typename Interface::TSystemVectorType& rRhs,
                              typename Interface::TSystemVectorType& rSolutionUpdate,
                              typename Interface::TSystemVectorType& rResidual,
                              const typename TSparse::DataType InitialResidualNorm,
                              PMGStatusStream& rStream,
                              PMGStatusStream::Report& rReport)
    {
        KRATOS_TRY
        // Inner loop for multigrid.
        rReport.multigrid_converged = false;
        rReport.multigrid_iteration = 0ul;
        rReport.multigrid_residual = std::numeric_limits<typename TSparse::DataType>::max();

        while (   rReport.multigrid_iteration < static_cast<std::size_t>(mMaxIterations)
               && !rReport.multigrid_converged) {
            // Solve the coarse grid and apply its correction.
            if (mMaybeHierarchy.has_value()) {
                std::visit([&rSolutionUpdate, &rResidual, &rStream, this](auto& r_hierarchy){
                                return r_hierarchy.template ApplyCoarseCorrection<TSparse>(
                                    rSolutionUpdate,
                                    rResidual,
                                    *mpConstraintAssembler,
                                    rStream);},
                           mMaybeHierarchy.value());

                // Update the fine solution.
                TSparse::UnaliasedAdd(rSolution, 1.0, rSolutionUpdate);

                // Update the fine residual.
                BalancedProduct<TSparse,TSparse,TSparse>(rLhs, rSolutionUpdate, rResidual, static_cast<typename TSparse::DataType>(-1));
            } // if mMaybeHierarchy

            // Perform smoothing on the fine grid.
            TSparse::SetToZero(rSolutionUpdate); //< do I need this?
            mpInterface->GetLinearSystemSolver()->InitializeSolutionStep(rLhs, rSolutionUpdate, rResidual);
            mpInterface->GetLinearSystemSolver()->Solve(rLhs, rSolutionUpdate, rResidual);
            mpInterface->GetLinearSystemSolver()->FinalizeSolutionStep(rLhs, rSolutionUpdate, rResidual);

            // Update the fine solution.
            TSparse::UnaliasedAdd(rSolution, 1.0, rSolutionUpdate);

            // Update the fine residual.
            BalancedProduct<TSparse,TSparse,TSparse>(rLhs, rSolutionUpdate, rResidual, static_cast<typename TSparse::DataType>(-1));

            // Emit status and check for convergence.
            rReport.multigrid_residual = TSparse::TwoNorm(rResidual) / InitialResidualNorm;
            rReport.multigrid_converged = rReport.multigrid_residual < mTolerance;
            if (   rReport.multigrid_iteration + 1 < static_cast<std::size_t>(mMaxIterations)
                && !rReport.multigrid_converged)
                rStream.Submit(rReport.Tag(3), mVerbosity);

            ++rReport.multigrid_iteration;
        } // while i_multigrid_iteration <= mMaxIterations

        if (rReport.multigrid_iteration) rReport.multigrid_iteration -= 1;
        KRATOS_CATCH("")
    }


    [[nodiscard]] PMGStatusStream::Report
    ExecuteConstraintLoop(typename Interface::TSystemMatrixType& rLhs,
                          typename Interface::TSystemVectorType& rSolution,
                          typename Interface::TSystemVectorType& rRhs,
                          PMGStatusStream& rStream)
    {
        KRATOS_TRY
        bool constraint_status = false;
        PMGStatusStream::Report status_report {
            /*grid_level=*/                 0ul,
            /*multigrid_converged=*/        false,
            /*multigrid_iteration=*/        0ul,
            /*multigrid_residual=*/         1.0,
            /*constraints_converged=*/      false,
            /*constraint_iteration=*/       0ul,
            /*maybe_constraint_residual=*/  {}
        };

        typename TSparse::VectorType residual(rRhs.size()),
                                     residual_update(rRhs.size()),
                                     solution_update(rSolution.size());

        // Compute the initial residual norm if it will be used.
        typename TSparse::DataType initial_residual_norm = 1;
        if (0 < mTolerance || 3 <= mVerbosity) {
            initial_residual_norm = TSparse::TwoNorm(rRhs);
            initial_residual_norm = initial_residual_norm ? initial_residual_norm : 1;
        }

        // Outer loop for constraints.
        do {
            // Initialize the constraint assembler and update residuals.
            mpConstraintAssembler->InitializeSolutionStep(rLhs, rSolution, rRhs);
            TSparse::Copy(rRhs, residual);
            BalancedProduct<TSparse,TSparse,TSparse>(rLhs, rSolution, residual, static_cast<typename TSparse::DataType>(-1));

            // Get an update on the solution with respect to the current residual.
            this->ExecuteMultigridLoop(rLhs,
                                       rSolution,
                                       rRhs,
                                       solution_update,
                                       residual,
                                       initial_residual_norm,
                                       rStream,
                                       status_report);

            // Check for constraint convergence.
            constraint_status = mpConstraintAssembler->FinalizeSolutionStep(rLhs,
                                                                            rSolution,
                                                                            rRhs,
                                                                            status_report,
                                                                            rStream);
            if (constraint_status) {
                rStream.Submit(status_report.Tag(2), mVerbosity);
                status_report.maybe_constraint_residual.reset();
                break;
            } else {
                rStream.Submit(status_report.Tag(3), mVerbosity);
                status_report.maybe_constraint_residual.reset();
            }

            ++status_report.constraint_iteration;
        } while (true);

        return status_report;
        KRATOS_CATCH("")
    }


    /// @brief Initialize the linear solver and solve the provided system.
    bool Solve(typename Interface::TSystemMatrixType& rLhs,
               typename Interface::TSystemVectorType& rSolution,
               typename Interface::TSystemVectorType& rRhs,
               ModelPart& rModelPart)
    {
        // Prepare and initialize members.
        KRATOS_TRY
        if (mpInterface->GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
            mpInterface->GetLinearSolver().ProvideAdditionalData(rLhs,
                                                                 rSolution,
                                                                 rRhs,
                                                                 mpInterface->GetDofSet(),
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

        std::optional<PMGStatusStream> status_stream = PMGStatusStream(
            /*rStream=*/            std::cout,
            /*rBuilderAndSolver=*/  *mpInterface,
            /*rModelPart=*/         rModelPart,
            /*rRootLhs=*/           rLhs,
            /*rRootSolution=*/      rSolution,
            /*rRootRhs=*/           rRhs,
            /*UseAnsiColors=*/      true);
        const auto status_report = this->ExecuteConstraintLoop(
            rLhs,
            rSolution,
            rRhs,
            status_stream.value());
        status_stream.reset();

        // Clean up.
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

        mpConstraintAssembler->Finalize(rLhs,
                                        rSolution,
                                        rRhs,
                                        mpInterface->GetDofSet());

        return status_report.multigrid_converged && status_report.constraints_converged;
        KRATOS_CATCH("")
    }


    // --------------------------------------------------------- //
    // Topology
    // --------------------------------------------------------- //


    template <class TProxy, class TIndexSet>
    static void CollectDoFs(const TProxy& rEntities,
                            const ProcessInfo& rProcessInfo,
                            typename Interface::TSchemeType& rScheme,
                            LockObject* pLockBegin,
                            [[maybe_unused]] LockObject* pLockEnd,
                            TIndexSet* pRowSetBegin,
                            [[maybe_unused]] TIndexSet* pRowSetEnd)
    {
        KRATOS_TRY
        KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
        using TLS = Element::EquationIdVectorType;
        block_for_each(rEntities,
                       TLS(),
                       [&rScheme, &rProcessInfo, pLockBegin, pRowSetBegin](const auto& r_entity, TLS& r_tls){
            if (r_entity.GetEntity().IsActive()) {
                rScheme.EquationId(r_entity.GetEntity(), r_tls, rProcessInfo);
                for (const auto equation_id : r_tls) {
                    [[maybe_unused]] std::scoped_lock<LockObject> lock(pLockBegin[equation_id]);
                    auto& r_row_indices = pRowSetBegin[equation_id];
                    r_row_indices.insert(r_tls.begin(), r_tls.end());
                } // for equation_id in r_tls
            } // if r_entity.IsActive
        });
        KRATOS_CATCH("")
    }


    void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme,
                         typename Interface::TSystemMatrixType& rLhs,
                         ModelPart& rModelPart)
    {
        KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
        KRATOS_TRY

        using IndexSet = std::unordered_set<std::size_t>;
        std::vector<IndexSet> indices(mpInterface->GetEquationSystemSize());

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
                                                             rLhs,
                                                             /*EnsureDiagonal=*/true);

        // Construct the coarse hierarchy's topology.
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
        // Constraint assembly MUST happen before assembling the unconstrained system,
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
                                *this->mpConstraintAssembler,
                                mpInterface->GetDofSet());
                        },
                        this->mMaybeHierarchy.value());
        } // if mMaybeHierarchy

        //if constexpr (AssembleLHS && AssembleRHS)
        //    NormalizeSystem<TSparse>(*pMaybeLhs.value(),
        //                             *pMaybeRhs.value(),
        //                             GetDiagonalScaleFactor<TSparse>(*pMaybeLhs.value(),
        //                                                             DiagonalScaling::AbsMax));
        KRATOS_CATCH("")
    }

}; // class PMultigridBuilderAndSolver::Impl


// --------------------------------------------------------- //
// Constructors
// --------------------------------------------------------- //


template <class TSparse, class TDense>
PMultigridBuilderAndSolver<TSparse,TDense>::~PMultigridBuilderAndSolver() = default;


template <class TSparse, class TDense>
PMultigridBuilderAndSolver<TSparse,TDense>::PMultigridBuilderAndSolver()
    : Interface(),
      mpImpl(new Impl(this))
{
}


template <class TSparse, class TDense>
PMultigridBuilderAndSolver<TSparse,TDense>::PMultigridBuilderAndSolver(const typename LinearSolverType::Pointer& pSolver,
                                                                       Parameters Settings)
    : Interface(pSolver),
      mpImpl(new Impl(this))
{
    this->AssignSettings(Settings);
}


template <class TSparse, class TDense>
typename BuilderAndSolver<TSparse,TDense,LinearSolver<TSparse,TDense>>::Pointer
PMultigridBuilderAndSolver<TSparse,TDense>::Create(typename LinearSolverType::Pointer pSolver,
                                                   Parameters Settings) const
{
    KRATOS_TRY
    return typename Interface::Pointer(new PMultigridBuilderAndSolver(pSolver, Settings));
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Allocation and Initialization
// --------------------------------------------------------- //


std::optional<std::size_t> FindNodeIndex(std::size_t NodeId,
                                         Kratos::ModelPart::NodesContainerType::const_iterator itNodeBegin,
                                         Kratos::ModelPart::NodesContainerType::const_iterator itNodeEnd) noexcept
{
    const auto it = std::lower_bound(itNodeBegin,
                                     itNodeEnd,
                                     NodeId,
                                     [](const Kratos::Node& r_node, std::size_t target_id){
                                        return r_node.Id() < target_id;
                                     });
    if (it == itNodeEnd || it->Id() != NodeId) {
        return {};
    } else {
        return std::distance(itNodeBegin, it);
    }
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::SetUpDofSet(typename Interface::TSchemeType::Pointer pScheme,
                                                             ModelPart& rModelPart)
{
    KRATOS_TRY;

    using DofsVectorType = ModelPart::DofsVectorType;
    this->GetDofSet().clear();
    const auto& r_process_info = rModelPart.GetProcessInfo();

    std::vector<std::atomic<std::uint8_t>> hanging_nodes(rModelPart.Nodes().size());
    std::fill(hanging_nodes.begin(),
              hanging_nodes.end(),
              static_cast<std::uint8_t>(1));

    // Fill the global DOF set array
    #pragma omp parallel
    {
        DofsVectorType tls_dofs, tls_constraint_dofs;

        // We create the temporal set in current thread and we reserve some space on it
        std::unordered_set<Node::DofType::Pointer, DofPointerHasher> dofs_tmp_set;
        dofs_tmp_set.reserve(20000);

        // Add the DOFs from the model part elements
        const auto& r_elements_array = rModelPart.Elements();
        const int number_of_elements = static_cast<int>(r_elements_array.size());

        #pragma omp for schedule(guided, 512) nowait
        for (int i_entity=0; i_entity<number_of_elements; ++i_entity) {
            const auto& r_entity = *(r_elements_array.begin() + i_entity);
            if (r_entity.IsActive()) {
                r_entity.GetDofList(tls_dofs, r_process_info);
                dofs_tmp_set.insert(tls_dofs.begin(), tls_dofs.end());
            } // r_entity.IsActive()
        } // for i_entity in range(number_of_elements)

        // Collect DoFs from constraints.
        const int number_of_constraints = static_cast<int>(rModelPart.MasterSlaveConstraints().size());

        #pragma omp for schedule(guided, 512) nowait
        for (int i_entity=0; i_entity<number_of_constraints; ++i_entity) {
            // Get current constraint iterator
            const auto it_entity = rModelPart.MasterSlaveConstraints().begin() + i_entity;

            // Trigger Dof index updates.
            it_entity->GetDofList(tls_dofs, tls_constraint_dofs, r_process_info);

            // Different behavior depending on whether MasterSlaveConstraint or MultifreedomConstraint.
            if (tls_dofs.empty() && !tls_constraint_dofs.empty()) {
                // Multifreedom constraint.
                dofs_tmp_set.insert(tls_constraint_dofs.begin(), tls_constraint_dofs.end());
            } else {
                // Master-slave constraint.
                dofs_tmp_set.insert(tls_dofs.begin(), tls_dofs.end());
                dofs_tmp_set.insert(tls_constraint_dofs.begin(), tls_constraint_dofs.end());
            }
        } // for i_entity in range(number_of_constraints)

        // Merge all the sets in one thread
        #pragma omp critical
        {
            this->GetDofSet().insert(dofs_tmp_set.begin(), dofs_tmp_set.end());
        }
    } // pragma omp parallel

    // Make sure that conditions act exclusively on collected DoFs.
    #ifdef KRATOS_DEBUG
    block_for_each(rModelPart.Conditions().begin(),
                   rModelPart.Conditions().end(),
                   DofsVectorType(),
                   [&r_process_info, this](const Condition& r_condition, auto& r_tls_dofs){
                        if (r_condition.IsActive()) {
                            r_condition.GetDofList(r_tls_dofs, r_process_info);
                            for (Dof<typename TDense::DataType>* p_dof : r_tls_dofs) {
                                KRATOS_ERROR_IF_NOT(this->GetDofSet().count(*p_dof))
                                    << "condition " << r_condition.Id() << " "
                                    << "acts on " << p_dof->GetVariable().Name() << " "
                                    << "of node " << p_dof->Id() << ", "
                                    << "which is not part of any element or constraint.";
                            } // for p_dof in r_tls_dofs
                        } // if r_condition.IsActive()
                   });
    #endif

    #ifdef KRATOS_DEBUG
    // If reactions are to be calculated, we check if all the dofs have reactions defined
    // This is to be done only in debug mode
    for (const auto& r_dof : this->GetDofSet()) {
        KRATOS_ERROR_IF_NOT(r_dof.HasReaction())
            << "Reaction variable not set for "
            << "DoF " << r_dof << " "
            << "of node "<< r_dof.Id();
    }
    #endif

    KRATOS_CATCH("");
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::SetUpSystem(ModelPart& rModelPart)
{
    KRATOS_TRY
    this->mEquationSystemSize = this->GetDofSet().size();

    // Set equation indices of DoFs.
    IndexPartition<std::size_t>(this->GetDofSet().size()).for_each([&, this](std::size_t Index){
        (this->GetDofSet().begin() + Index)->SetEquationId(Index);
    });
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme,
                                                                            typename Interface::TSystemMatrixPointerType& rpLhs,
                                                                            typename Interface::TSystemVectorPointerType& rpSolution,
                                                                            typename Interface::TSystemVectorPointerType& rpRhs,
                                                                            ModelPart& rModelPart)
{
    KRATOS_TRY

    // Construct empty containers if necessary.
    if (!rpLhs)
        rpLhs.reset(new typename Interface::TSystemMatrixType);
    TSparse::SetToZero(*rpLhs);

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

    this->SetDofSetIsInitializedFlag(true);

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Hooks
// --------------------------------------------------------- //


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::InitializeSolutionStep(ModelPart& rModelPart,
                                                                        typename Interface::TSystemMatrixType& rLhs,
                                                                        typename Interface::TSystemVectorType& rSolution,
                                                                        typename Interface::TSystemVectorType& rRhs)
{
    mpImpl->mpMaybeModelPart = &rModelPart;
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse, TDense>::FinalizeSolutionStep(ModelPart& rModelPart,
                                                                       typename Interface::TSystemMatrixType& rLhs,
                                                                       typename Interface::TSystemVectorType& rSolution,
                                                                       typename Interface::TSystemVectorType& rRhs)
{
    mpImpl->mpMaybeModelPart.reset();
}


// --------------------------------------------------------- //
// Assembly
// --------------------------------------------------------- //


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::Build(typename Interface::TSchemeType::Pointer pScheme,
                                                       ModelPart& rModelPart,
                                                       typename Interface::TSystemMatrixType& rLhs,
                                                       typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme" << std::endl;
    KRATOS_TRY
    mpImpl->mpMaybeModelPart = &rModelPart;
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/true>(rModelPart,
                                                                         *pScheme,
                                                                         &rLhs,
                                                                         &rRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::BuildLHS(typename Interface::TSchemeType::Pointer pScheme,
                                                          ModelPart& rModelPart,
                                                          typename Interface::TSystemMatrixType& rLhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme";
    KRATOS_TRY
    mpImpl->template Assemble</*AssembleLHS=*/true,/*AssembleRHS=*/false>(rModelPart,
                                                                          *pScheme,
                                                                          &rLhs,
                                                                          nullptr);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::BuildRHS(typename Interface::TSchemeType::Pointer pScheme,
                                                          ModelPart& rModelPart,
                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_ERROR_IF(!pScheme) << "missing scheme";
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


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme,
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

    mpImpl->mpDiagonalScaling->template Cache<TSparse>(rLhs);
    const auto diagonal_scale = mpImpl->mpDiagonalScaling->Evaluate();
    Kratos::ApplyDirichletConditions<TSparse,TDense>(rLhs,
                                                     rRhs,
                                                     it_dof_set_begin,
                                                     it_dof_set_end,
                                                     diagonal_scale);
    if (mpImpl->mMaybeHierarchy.has_value()) {
        std::visit([this](auto& r_grid){
                       r_grid.ApplyDirichletConditions(this->GetDofSet().begin(),
                                                       this->GetDofSet().end());
                   },
                   mpImpl->mMaybeHierarchy.value());
    } // if mMaybeHierarchy

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                  ModelPart& rModelPart,
                                                                  typename Interface::TSystemMatrixType& rLhs,
                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    mpImpl->mpConstraintAssembler->Initialize(rLhs,
                                              rRhs,
                                              this->GetDofSet().begin(),
                                              this->GetDofSet().end());
    if (mpImpl->mMaybeHierarchy.has_value()) {
        std::visit([](auto& r_hierarchy){r_hierarchy.ApplyConstraints();},
                   mpImpl->mMaybeHierarchy.value());
    }
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Compound Assembly and Solution
// --------------------------------------------------------- //


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::BuildAndSolve(typename Interface::TSchemeType::Pointer pScheme,
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

    this->SystemSolve(rLhs, rSolution, rRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::SystemSolve(typename Interface::TSystemMatrixType& rLhs,
                                                             typename Interface::TSystemVectorType& rSolution,
                                                             typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_TRY
    // Solve constrained assembled system.
    BuiltinTimer timer;
    if (!mpImpl->Solve(rLhs, rSolution, rRhs, *mpImpl->mpMaybeModelPart.value()) && 1 <= mpImpl->mVerbosity) {
        std::cerr << this->Info() << "Grid 0: failed to solve the assembled system\n";
    }
    KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity)
        << ": Linear solution took " << timer << "\n";
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Postprocessing
// --------------------------------------------------------- //



template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::CalculateReactions(typename Interface::TSchemeType::Pointer pScheme,
                                                                    ModelPart& rModelPart,
                                                                    typename Interface::TSystemMatrixType& rLhs,
                                                                    typename Interface::TSystemVectorType& rSolution,
                                                                    typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);

    KRATOS_TRY
    BuiltinTimer timer;

    TSparse::SetToZero(rRhs);
    mpImpl->template Assemble</*AssembleLHS=*/false,/*AssembleRHS=*/true>(rModelPart,
                                                                          *pScheme,
                                                                          nullptr,
                                                                          &rRhs);
    block_for_each(this->GetDofSet(), [&rRhs](Dof<typename TDense::DataType>& rDof){
        rDof.GetSolutionStepReactionValue() = -rRhs[rDof.EquationId()];
    });

    KRATOS_INFO_IF(this->Info(), 2 <= mpImpl->mVerbosity)
        << "Computing reactions took " << timer << "\n";
    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Misc
// --------------------------------------------------------- //


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::AssignSettings(const Parameters Settings)
{
    // Mutable parameters.
    Parameters settings = Settings;

    KRATOS_TRY
    // Parse diagonal scaling and validate other settings.
    Parameters default_parameters = this->GetDefaultParameters();
    Parameters default_diagonal_scaling = default_parameters["diagonal_scaling"].Clone();
    default_parameters.RemoveValue("diagonal_scaling");
    std::optional<Parameters> maybe_diagonal_scaling = settings.Has("diagonal_scaling") ? settings["diagonal_scaling"].Clone() : std::optional<Parameters>();
    if (maybe_diagonal_scaling.has_value()) settings.RemoveValue("diagonal_scaling");
    settings.ValidateAndAssignDefaults(default_parameters);
    settings.AddValue("diagonal_scaling", maybe_diagonal_scaling.has_value() ? maybe_diagonal_scaling.value() : default_diagonal_scaling);
    mpImpl->mpDiagonalScaling = std::make_unique<Scaling>(settings["diagonal_scaling"]);
    KRATOS_CATCH("")

    KRATOS_TRY
    Interface::AssignSettings(settings);
    KRATOS_CATCH("")

    // Set the relative residual tolerance and max W cycles.
    mpImpl->mMaxIterations = settings["max_iterations"].Get<int>();
    mpImpl->mTolerance = settings["tolerance"].Get<double>();

    // Set multifreedom constraint imposition strategy.
    mpImpl->mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(settings["constraint_imposition_settings"],
                                                                               "Grid 0 constraints");

    // Construct smoother, solver and coarse hierarchy.
    Parameters smoother_settings = settings["smoother_settings"];
    Parameters coarse_hierarchy_settings = settings["coarse_hierarchy_settings"];
    Parameters leaf_solver_settings = settings["linear_solver_settings"];

    if (coarse_hierarchy_settings["max_depth"].Get<int>()) {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(smoother_settings.Has("solver_type"));
        const std::string solver_name = smoother_settings["solver_type"].Get<std::string>();
        using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;
        KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
            << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
            << "Make sure you imported the application it is defined in and that the spelling is correct.";
        const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
        this->mpLinearSystemSolver = r_factory.Create(smoother_settings);

        // Construct the coarse hierarchy.
        const std::string coarse_build_precision = coarse_hierarchy_settings["precision"].Get<std::string>();
        if (coarse_build_precision == "double") {
            using CoarseSparseSpace = TUblasSparseSpace<double>;
            using CoarseDenseSpace = TUblasDenseSpace<double>;
            using GridType = PGrid<CoarseSparseSpace,CoarseDenseSpace>;

            coarse_hierarchy_settings.ValidateAndAssignDefaults(GridType().GetDefaultParameters());
            if (coarse_hierarchy_settings["max_depth"].Get<int>()) {
                mpImpl->mMaybeHierarchy = GridType(coarse_hierarchy_settings,
                                                   smoother_settings,
                                                   leaf_solver_settings,
                                                   settings["diagonal_scaling"]);
            }
        } /* if coarse_build_precision == "double" */ else if (coarse_build_precision == "single") {
            using CoarseSparseSpace = TUblasSparseSpace<float>;
            using CoarseDenseSpace = TUblasDenseSpace<double>;
            using GridType = PGrid<CoarseSparseSpace,CoarseDenseSpace>;

            coarse_hierarchy_settings.ValidateAndAssignDefaults(GridType().GetDefaultParameters());
            mpImpl->mMaybeHierarchy = GridType(coarse_hierarchy_settings,
                                               smoother_settings,
                                               leaf_solver_settings,
                                               settings["diagonal_scaling"]);
        } /* elif coarse_build_precision == "single" */ else {
            KRATOS_ERROR << "unsupported coarse precision: \"" << coarse_build_precision << "\". "
                         << "Options are:\n"
                         << "\t\"single\""
                         << "\t\"double\"";
        }
        KRATOS_CATCH("")
    } /* if coarse_hierarchy_settings["max_depth"].Get<int>() */ else {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(leaf_solver_settings.Has("solver_type"));
        const std::string solver_name = leaf_solver_settings["solver_type"].Get<std::string>();
        using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;
        KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
            << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
            << "Make sure you imported the application it is defined in and that the spelling is correct.";
        const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
        this->mpLinearSystemSolver = r_factory.Create(leaf_solver_settings);
        KRATOS_CATCH("")
    }

    // Other settings.
    KRATOS_TRY
    mpImpl->mVerbosity = settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::Clear()
{
    Interface::Clear();
    mpImpl->mpConstraintAssembler->Clear();

    if (mpImpl->mMaybeHierarchy.has_value()) {
        std::visit([](auto& r_grid) {r_grid.Clear();},
                   mpImpl->mMaybeHierarchy.value());
    } // if mMaybeHierarchy

    this->SetDofSetIsInitializedFlag(false);
}


template <class TSparse, class TDense>
Parameters PMultigridBuilderAndSolver<TSparse,TDense>::GetDefaultParameters() const
{
    Parameters parameters = Parameters(R"({
"name"              : "p_multigrid",
"diagonal_scaling"  : "max",
"max_iterations"    : 1e2,
"tolerance"         : 1e-8,
"verbosity"         : 1,
"smoother_settings" : {
    "solver_type" : ""
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


template <class TSparse, class TDense>
std::string PMultigridBuilderAndSolver<TSparse,TDense>::Info() const
{
    return "PMultigridBuilderAndSolver";
}


template <class TSparse, class TDense>
std::size_t PMultigridBuilderAndSolver<TSparse,TDense>::GetEquationSystemSize() const noexcept
{
    return Interface::mEquationSystemSize;
}

template <class TSparse, class TDense>
LinearSolver<TSparse,TDense>& PMultigridBuilderAndSolver<TSparse,TDense>::GetLinearSolver() noexcept
{
    return *Interface::mpLinearSystemSolver;
}


template <class TSparse, class TDense>
void PMultigridBuilderAndSolver<TSparse,TDense>::ProjectGrid(int GridLevel,
                                                             const typename TSparse::MatrixType& rRootLhs,
                                                             const typename TSparse::VectorType& rRootSolution,
                                                             const typename TSparse::VectorType& rRootRhs)
{
    KRATOS_TRY

    typename TSparse::VectorType projected;

    // Project residuals.
    if (0 < GridLevel) {
        KRATOS_ERROR_IF_NOT(mpImpl->mMaybeHierarchy.has_value());
        const auto i_grid_level = GridLevel - 1;

        // Construct a flat vector of coarse grids.
        std::variant<
            std::vector<const PGrid<TUblasSparseSpace<double>,TUblasDenseSpace<double>>*>,
            std::vector<const PGrid<TUblasSparseSpace<float>,TUblasDenseSpace<double>>*>
        > coarse_grids;

        std::visit([&coarse_grids](const auto& r_coarse_grid){
            using GridType = std::remove_cv_t<std::remove_reference_t<decltype(r_coarse_grid)>>;
            coarse_grids = std::vector<const GridType*>();
        }, mpImpl->mMaybeHierarchy.value());

        const auto vector_fill_visitor = [&coarse_grids](const auto& r_coarse_grid) {
            std::visit([&r_coarse_grid](auto& r_vector){
                using GridType = std::remove_cv_t<std::remove_reference_t<decltype(r_coarse_grid)>>;
                using ValueType = typename std::remove_cv_t<std::remove_reference_t<decltype(r_vector)>>::value_type;
                if constexpr (std::is_same_v<ValueType,const GridType*>) {
                    const GridType* p_grid = &r_coarse_grid;
                    do {
                        r_vector.push_back(p_grid);
                        const auto p_maybe_child = p_grid->GetChild();
                        p_grid = p_maybe_child.has_value() ? p_maybe_child.value() : nullptr;
                    } while (p_grid);
                }
            }, coarse_grids);
        };

        std::visit(vector_fill_visitor, mpImpl->mMaybeHierarchy.value());

        // Initialize projected solution vector.
        std::visit([&projected, i_grid_level](const auto& r_vector){
            KRATOS_ERROR_IF_NOT(i_grid_level < r_vector.size());
            const auto& r_solution = r_vector[i_grid_level]->GetSolution();
            projected.resize(r_solution.size(), false);
            std::fill(projected.begin(), projected.end(), static_cast<typename TSparse::DataType>(0));
        }, coarse_grids);

        // Project solution to the root grid.
        {
            typename TSparse::VectorType tmp;
            for (int i_grid=i_grid_level; 0<=i_grid; --i_grid) {
                std::visit([&projected, &tmp, i_grid](const auto& r_vector) {
                    const auto& r_solution = r_vector[i_grid]->GetSolution();
                    IndexPartition<std::size_t>(r_solution.size()).for_each([&r_solution, &projected](std::size_t i_component){
                        projected[i_component] += r_solution[i_component];
                    });
                    r_vector[i_grid]->template Prolong<TSparse>(projected, tmp);
                }, coarse_grids);
                tmp.swap(projected);
            }
        }

    } /*if 0 < GridLevel*/ else {
        projected.resize(rRootRhs.size(), false);
        TSparse::SetToZero(projected);
    }

    // Add diff of the root grid.
    TSparse::UnaliasedAdd(projected, static_cast<typename TSparse::DataType>(1), rRootSolution);

    // Compute residual.
    auto residual = rRootRhs;
    BalancedProduct<TSparse,TSparse,TSparse>(rRootLhs,
                                             projected,
                                             residual,
                                             static_cast<typename TSparse::DataType>(-1));

    // Update DoFs.
    IndexPartition<std::size_t>(projected.size()).for_each([this, &projected, &residual](std::size_t i_dof) {
        auto& r_dof = *(this->GetDofSet().begin() + i_dof);
        (this->GetDofSet().begin() + i_dof)->GetSolutionStepReactionValue() = residual[i_dof];
        r_dof.GetSolutionStepValue() += projected[i_dof];
    });

    KRATOS_CATCH("")
}


// --------------------------------------------------------- //
// Template Instantiations
// --------------------------------------------------------- //

template class KRATOS_API(KRATOS_CORE) PMultigridBuilderAndSolver<TUblasSparseSpace<double>,
                                                                  TUblasDenseSpace<double>>;

template class KRATOS_API(KRATOS_CORE) PMultigridBuilderAndSolver<TUblasSparseSpace<float>,
                                                                  TUblasDenseSpace<double>>;


} // namespace Kratos
