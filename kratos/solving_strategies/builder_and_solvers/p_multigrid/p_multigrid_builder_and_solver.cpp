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
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp" // MakePRestrictionOperator, MakeSparseTopology
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // MakeSparseTopology
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling, ParseDiagonalScaling, GetDiagonalScalingFactor
#include "includes/model_part.h" // ModelPart
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#include "utilities/dof_utilities/block_build_dof_array_utility.h" // BlockBuildDofArrayUtility
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/proxies.h" // MakeProxy
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE, KRATOS_PROFILE_SCOPE_MILLI

// System includesPa
#include <optional> // std::optional
#include <limits> // std::numeric_limits
#include <unordered_set> // std::unordered_set


namespace Kratos {


// --------------------------------------------------------- //
// PMG Hierarchy
// --------------------------------------------------------- //


template <class TSparse,
          class TDense,
          class TSolver>
class PGrid
{
public:
    PGrid(const unsigned Level)
        : mpRestrictionOperator(new typename TSparse::MatrixType),
          mpLhs(new typename TSparse::MatrixType),
          mpRhs(new typename TSparse::VectorType)
    {}


    PGrid()
        : PGrid(0u)
    {}


    template <class TParentSparse>
    void MakeLhsTopology(const ModelPart& rModelPart,
                         const typename TParentSparse::MatrixType& rParentLhs)
    {
        KRATOS_TRY
        // The restriction operator immediately constructs the linear equivalent,
        // because one-level coarsening strategies are not supported yet. The problem
        // is that when going from a generic polynomial level q to some other lower polynomial
        // level p!=1, new nodes must be introduced that do not exist in the original fine mesh.
        // This wouldn't be a huge problem by itself, but deciding on DoF indices on the coarse
        // level would be quite painful and require to keep track of the coarse grid's topology
        // in some manner.
        MakePRestrictionOperator<
            std::numeric_limits<unsigned>::max(),
            typename TSparse::DataType>(rModelPart,
                                        *mpRestrictionOperator,
                                        rParentLhs.size1());
        KRATOS_CATCH("")
    }


    template <bool AssembleLHS,
              bool AssembleRHS,
              class TParentSparse>
    void Assemble(const ModelPart& rModelPart,
                  const typename TParentSparse::MatrixType* pParentLhs,
                  const typename TParentSparse::VectorType* pParentRhs)
    {
        KRATOS_TRY

        // Assemble the coarse LHS matrix if requested.
        if constexpr (AssembleLHS) {
            KRATOS_ERROR_IF_NOT(pParentLhs);

            const auto check_matrix_multiplication = [](const typename TSparse::MatrixType& r_left,
                                                        const typename TSparse::MatrixType& r_right,
                                                        const typename TSparse::MatrixType& r_out) {
                KRATOS_ERROR_IF((r_left.size2() != r_right.size1()) ||
                                (r_out.size1() != r_left.size1()) ||
                                (r_out.size2() != r_right.size2()))
                    << "invalid sizes for matrix product: "
                    << "(" << r_left.size1() << "x" << r_left.size2() << ")"
                    << " x "
                    << "(" << r_right.size1() << "x" << r_right.size2() << ")"
                    << " => "
                    << "(" << r_out.size1() << "x" << r_out.size2() << ")";
            };

            typename TSparse::MatrixType prolongation_operator(mpRestrictionOperator->size2(), mpRestrictionOperator->size1()),
                                         left_multiplied_lhs(mpRestrictionOperator->size1(), pParentLhs->size2());

            check_matrix_multiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(*mpRestrictionOperator, *pParentLhs, left_multiplied_lhs);
            SparseMatrixMultiplicationUtility::TransposeMatrix(prolongation_operator, *mpRestrictionOperator, 1.0);

            mpLhs->resize(mpRestrictionOperator->size1(), mpRestrictionOperator->size1(), false);
            check_matrix_multiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
            SparseMatrixMultiplicationUtility::MatrixMultiplication(left_multiplied_lhs, prolongation_operator, *mpLhs);
        } // if AssembleLHS

        // Assemble the coarse RHS vector if requested.
        if constexpr (AssembleRHS) {
            // Sanity checks
            KRATOS_ERROR_IF_NOT(pParentRhs);
            KRATOS_ERROR_IF_NOT(mpRestrictionOperator->size2() == pParentRhs->size())
                << "expecting an RHS vector of size " << mpRestrictionOperator->size2()
                << " but got " << pParentRhs->size();

            mpRhs->resize(mpRestrictionOperator->size1(), false);
            TSparse::Mult(*mpRestrictionOperator, *pParentRhs, *mpRhs);
        } // if AssembleRHS

        KRATOS_CATCH("")
    }

private:
    std::shared_ptr<typename TSparse::MatrixType> mpRestrictionOperator;

    std::shared_ptr<typename TSparse::MatrixType> mpLhs;

    std::shared_ptr<typename TSparse::VectorType> mpRhs;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;
}; // class PGrid


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

    PGrid<TSparse,TDense,TSolver> mHierarchy;

    DiagonalScaling mDiagonalScaling;

    int mMaxDepth;

    int mVerbosity;

    // --------------------------------------------------------- //
    // Special Member Functions
    // --------------------------------------------------------- //

    Impl(Interface* pInterface)
        : mpInterface(pInterface),
          mpConstraintAssembler(),
          mDiagonalScaleFactor(1),
          mHierarchy(),
          mDiagonalScaling(DiagonalScaling::None),
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
        if (rInterface.GetLinearSolver().AdditionalPhysicalDataIsNeeded()) {
            rInterface.GetLinearSolver().ProvideAdditionalData(rLhs,
                                                               rSolution,
                                                               rRhs,
                                                               rInterface.GetDofSet(),
                                                               rModelPart);
        }
        KRATOS_CATCH("")

        std::size_t i_iteration = 0ul;
        typename ConstraintAssembler<TSparse,TDense>::Status constraint_status {/*converged=*/false, /*finished=*/false};
        bool linear_solver_status = false; //< Indicates whether the linear solver converged.
        auto& r_linear_solver = rInterface.GetLinearSolver();

        do {
            mpConstraintAssembler->InitializeSolutionStep(rModelPart.MasterSlaveConstraints(),
                                                          rModelPart.GetProcessInfo(),
                                                          rLhs,
                                                          rSolution,
                                                          rRhs,
                                                          rInterface.GetDofSet(),
                                                          i_iteration);
            r_linear_solver.InitializeSolutionStep(rLhs, rSolution, rRhs);
            linear_solver_status = r_linear_solver.Solve(rLhs, rSolution, rRhs);
            r_linear_solver.FinalizeSolutionStep(rLhs, rSolution, rRhs);
            constraint_status = mpConstraintAssembler->FinalizeSolutionStep(rModelPart.MasterSlaveConstraints(),
                                                                 rModelPart.GetProcessInfo(),
                                                                 rLhs,
                                                                 rSolution,
                                                                 rRhs,
                                                                 rInterface.GetDofSet(),
                                                                 ++i_iteration);
        } while (not constraint_status.finished);

        if (1 <= mVerbosity and not constraint_status.converged) {
            std::cerr << "PMultigridBuilderAndSolver: failed to converge constraints\n";
        }

        mpConstraintAssembler->Finalize(rModelPart.MasterSlaveConstraints(),
                                        rModelPart.GetProcessInfo(),
                                        rLhs,
                                        rSolution,
                                        rRhs,
                                        mpInterface->GetDofSet());

        return linear_solver_status && constraint_status.converged;
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
        mHierarchy.template MakeLhsTopology<TSparse>(rModelPart, rLhs);
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

        // Global variables.
        const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();
        typename Interface::TSystemMatrixType* pLhs = pMaybeLhs.has_value() ? pMaybeLhs.value() : nullptr;
        typename Interface::TSystemVectorType* pRhs = pMaybeRhs.has_value() ? pMaybeRhs.value() : nullptr;
        auto p_locks = std::make_unique<std::vector<LockObject>>(AssembleLHS ? pLhs->size1() : 0ul);

        [[maybe_unused]] const int element_count = rModelPart.Elements().size();
        [[maybe_unused]] const int condition_count = rModelPart.Conditions().size();
        [[maybe_unused]] const auto it_element_begin = rModelPart.Elements().begin();
        [[maybe_unused]] const auto it_condition_begin = rModelPart.Conditions().begin();

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
        mHierarchy.template Assemble<AssembleLHS, AssembleRHS,TSparse>(rModelPart,
                                                                       pLhs,
                                                                       pRhs);

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
    block_for_each(this->GetDofSet(), [&rRhs](Dof<double>& rDof){
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
                                                                                  typename Interface::TSystemVectorType& rSolution,
                                                                                  typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    const std::size_t system_size = rLhs.size1();
    Vector scaling_factors (system_size);

    const auto it_dof_iterator_begin = this->GetDofSet().begin();

    // NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
    IndexPartition<std::size_t>(this->GetDofSet().size()).for_each([&](std::size_t Index){
        auto it_dof_iterator = it_dof_iterator_begin + Index;
        if (it_dof_iterator->IsFixed()) {
            scaling_factors[Index] = 0.0;
        } else {
            scaling_factors[Index] = 1.0;
        }
    });

    //GetDiagonalScaleFactor<TSparse>(mpImpl->mDiagonalScaleFactor,
    //                                rLhs,
    //                                mpImpl->mDiagonalScaling,
    //                                rModelPart.GetProcessInfo());

    double* Avalues = rLhs.value_data().begin();
    std::size_t* Arow_indices = rLhs.index1_data().begin();
    std::size_t* Acol_indices = rLhs.index2_data().begin();

    IndexPartition<std::size_t>(system_size).for_each([&](std::size_t Index){
        const std::size_t col_begin = Arow_indices[Index];
        const std::size_t col_end = Arow_indices[Index+1];
        const double k_factor = scaling_factors[Index];
        if (k_factor == 0.0) {
            // Zero out the whole row, except the diagonal
            for (std::size_t j = col_begin; j < col_end; ++j)
                if (Acol_indices[j] != Index )
                    Avalues[j] = 0.0;
            // Zero out the RHS
            rRhs[Index] = 0.0;
        } else {
            // Zero out the column which is associated with the zero'ed row
            for (std::size_t j = col_begin; j < col_end; ++j)
                if(scaling_factors[ Acol_indices[j] ] == 0 )
                    Avalues[j] = 0.0;
        }
    });
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::ApplyConstraints(typename Interface::TSchemeType::Pointer pScheme,
                                                                          ModelPart& rModelPart,
                                                                          typename Interface::TSystemMatrixType& rLhs,
                                                                          typename Interface::TSystemVectorType& rRhs)
{
    KRATOS_PROFILE_SCOPE_MILLI(KRATOS_CODE_LOCATION);
    KRATOS_TRY
    mpImpl->mpConstraintAssembler->Assemble(
        rModelPart.MasterSlaveConstraints(),
        rModelPart.GetProcessInfo(),
        this->GetDofSet());
    mpImpl->mpConstraintAssembler->Initialize(
        rModelPart.MasterSlaveConstraints(),
        rModelPart.GetProcessInfo(),
        rLhs,
        rRhs,
        this->GetDofSet());
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
        std::cerr << this->Info() << ": failed to solve the assembled system\n";
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
    block_for_each(this->GetDofSet(), [&rRhs](Dof<double>& rDof){
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

    // Set multifreedom constraint imposition strategy.
    mpImpl->mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(Settings["constraint_imposition"]);

    // Other settings.
    KRATOS_TRY
    mpImpl->mMaxDepth = Settings["max_depth"].Get<int>();
    mpImpl->mVerbosity = Settings["verbosity"].Get<int>();
    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TSolver>
void PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::Clear()
{
    Interface::Clear();
    mpImpl->mpConstraintAssembler->Clear();
}


template <class TSparse, class TDense, class TSolver>
Parameters PMultigridBuilderAndSolver<TSparse,TDense,TSolver>::GetDefaultParameters() const
{
    Parameters parameters = Parameters(R"({
"name"              : "p_multigrid",
"diagonal_scaling"  : "abs_max",
"max_depth"         : -1,
"smoother_settings" : {
    "solver_type" : "gauss_seidel"
},
"solver_settings"   : {
    "solver_type" : "amgcl"
},
"constraint_imposition" : {
    "method" : "master_slave_elimination"
},
"verbosity"         : 1
})");
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


} // namespace Kratos
