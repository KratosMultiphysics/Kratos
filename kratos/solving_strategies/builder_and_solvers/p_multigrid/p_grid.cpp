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

// Boost's UBLAS pollutes macros into Kratos, breaking MSVC (of course, what else)
// when including <span>, because it manipulates iterator debugging macros (namely,
// _ITERATOR_DEBUG_LEVEL). They already realized this and have an open issue about it,
// but after 6 years it still stays unsolved https://github.com/boostorg/ublas/issues/77.
// The only viable fix is to make sure that <span> is always included BEFORE any UBLAS
// header. Of course this can't be guaranteed easily, so <span> must precede boost includes
// in all headers. Wow.
#include <span> // std::span

// Project includes
#include "solving_strategies/builder_and_solvers/p_multigrid/p_grid.hpp" // PGrid
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp" // MakePRestrictionOperator
#include "solving_strategies/builder_and_solvers/p_multigrid/constraint_assembler_factory.hpp" // ConstraintAssemblerFactory
#include "solving_strategies/builder_and_solvers/p_multigrid/sparse_utilities.hpp" // ApplyDirichletConditions
#include "solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp"
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "factories/linear_solver_factory.h" // LinearSolverFactory
#include "includes/kratos_components.h" // KratosComponents
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility
#include "utilities/builtin_timer.h" // BuiltinTimer
#include "solving_strategies/builder_and_solvers/p_multigrid/diagonal_scaling.hpp" // DiagonalScaling
#include "solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities.hpp" // detail::DofData

// System includes
#include <limits> // std::numeric_limits
#include <vector> // std::vector


namespace Kratos {


template <class TSparse, class TDense>
struct PGrid<TSparse,TDense>::Impl {
    typename TSparse::MatrixType mRestrictionOperator;

    typename TSparse::MatrixType mProlongationOperator;

    typename TSparse::MatrixType mLhs;

    typename TSparse::VectorType mSolution;

    typename TSparse::VectorType mRhs;

    VariablesList::Pointer mpVariableList;

    /// @details Array of @ref Dof "DoFs" unique to the current grid level.
    ///          DoFs need a pointer to a @ref NodalData object, which is why
    ///          pairs of @ref NodalData and @ref Dof are stored instead of just
    ///          DoFs.
    std::vector<detail::DofData> mDofSet;

    IndirectDofSet mIndirectDofSet;

    std::vector<std::size_t> mDofMap;

    std::shared_ptr<ConstraintAssembler<TSparse,TDense>> mpConstraintAssembler;

    std::unique_ptr<Scaling> mpDiagonalScaling;

    typename LinearSolverType::Pointer mpSolver;

    std::optional<std::unique_ptr<PGrid>> mMaybeChild;

    int mVerbosity;

    unsigned mDepth;
}; // struct PGrid::Impl


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid(
    Parameters Settings,
    const unsigned CurrentDepth,
    Parameters SmootherSettings,
    Parameters LeafSolverSettings,
    Parameters DiagonalScalingSettings)
        : mpImpl(new Impl {
            .mRestrictionOperator = {},
            .mProlongationOperator = {},
            .mLhs = {},
            .mSolution = {},
            .mRhs = {},
            .mDofSet = {},
            .mIndirectDofSet = {},
            .mDofMap = {},
            .mpConstraintAssembler = {},
            .mpDiagonalScaling = {},
            .mpSolver = {},
            .mMaybeChild = {},
            .mVerbosity = {},
            .mDepth = CurrentDepth}) {
                // Sanity checks.
                KRATOS_ERROR_IF_NOT(mpImpl->mDepth) << "PGrid can only have positive depths (the original system is at depth=0)";

                KRATOS_TRY
                Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

                // Check float precision.
                using value_type = typename TSparse::DataType;
                const std::string precision = Settings["precision"].Get<std::string>();
                if constexpr (std::is_same_v<value_type,double>) {
                    KRATOS_ERROR_IF_NOT(precision == "double")
                        << "attempting to construct a PGrid with inconsistent sparse space type "
                        << "(requested precision: \"" << precision << "\" "
                        << " build precision: \"double\")";
                } else if constexpr (std::is_same_v<value_type,float>) {
                    KRATOS_ERROR_IF_NOT(precision == "single")
                        << "attempting to construct a PGrid with inconsistent sparse space type "
                        << "(requested precision: \"" << precision << "\" "
                        << " build precision: \"single\")";
                } else {
                    static_assert(!std::is_same_v<value_type,value_type>, "unhandled sparse space type");
                }

                mpImpl->mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(
                    Settings["constraint_imposition_settings"],
                    "Grid " + std::to_string(mpImpl->mDepth) + " constraints");
                mpImpl->mpDiagonalScaling = std::make_unique<Scaling>(DiagonalScalingSettings);
                mpImpl->mVerbosity = Settings["verbosity"].Get<int>();

                const int max_depth = Settings["max_depth"].Get<int>();
                KRATOS_ERROR_IF_NOT(0 <= max_depth) << Settings << "\n\"max_depth\" must be non-negative";
                if (mpImpl->mDepth < static_cast<unsigned>(max_depth)) {
                    KRATOS_TRY
                    KRATOS_ERROR_IF_NOT(SmootherSettings.Has("solver_type"));
                    const std::string solver_name = SmootherSettings["solver_type"].Get<std::string>();
                    using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;
                    KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
                        << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
                        << "Make sure you imported the application it is defined in and that the spelling is correct.";
                    const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
                    mpImpl->mpSolver = r_factory.Create(SmootherSettings);
                    mpImpl->mMaybeChild = std::unique_ptr<PGrid>(
                        new PGrid(
                            Settings,
                            mpImpl->mDepth + 1u,
                            SmootherSettings,
                            LeafSolverSettings,
                            DiagonalScalingSettings));
                    KRATOS_CATCH(std::to_string(mpImpl->mDepth + 1u));
                } else {
                    KRATOS_ERROR_IF_NOT(LeafSolverSettings.Has("solver_type"));
                    const std::string solver_name = LeafSolverSettings["solver_type"].Get<std::string>();
                    using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;

                    if (!SolverFactoryRegistry::Has(solver_name)) {
                        std::stringstream message;
                        message << "PMultigridBuilderAndSolver: "
                                << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
                                << "Make sure you imported the application it is defined in and that the spelling is correct. "
                                << "Registered options are:\n";
                        for ([[maybe_unused]] const auto& [r_name, r_entry] : SolverFactoryRegistry::GetComponents()) {
                            message << "\t" << r_name << "\n";
                        }
                        KRATOS_ERROR << message.str();
                    } // if not SolverFactoryRegistry::Has(solver_name)

                    const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
                    mpImpl->mpSolver = r_factory.Create(LeafSolverSettings);
                }
                KRATOS_CATCH("")
}


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid(Parameters Settings,
                             Parameters SmootherSettings,
                             Parameters LeafSolverSettings,
                             Parameters DiagonalScalingSettings)
    : PGrid(Settings,
            1u,
            SmootherSettings,
            LeafSolverSettings,
            DiagonalScalingSettings)
{
}


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid()
    : PGrid(Parameters(),
            Parameters(R"({"solver_type" : "gauss_seidel"})"),
            Parameters(R"({"solver_type" : "amgcl"})"),
            Parameters(R"("norm")"))
{
}


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid(PGrid&&) noexcept = default;


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::~PGrid() = default;


template <class TSparse, class TDense>
PGrid<TSparse,TDense>& PGrid<TSparse,TDense>::operator=(PGrid&&) noexcept = default;


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::MakeLhsTopology(ModelPart& rModelPart,
                                            const typename TParentSparse::MatrixType& rParentLhs,
                                            [[maybe_unused]] const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                                            const IndirectDofSet& rParentDofSet)
{
    // ConstraintAssembler::Allocate is deliberately not invoked here because the relation matrix
    // and constraint gap vector are passed from the fine level in PGrid::Assemble, and mapped
    // to this level by the restriction operator.
}


template <class TSparse, class TDense>
template <bool AssembleLHS,
          bool AssembleRHS,
          class TParentSparse>
void PGrid<TSparse,TDense>::Assemble(ModelPart& rModelPart,
                                     const typename TParentSparse::MatrixType* pParentLhs,
                                     const typename TParentSparse::VectorType* pParentRhs,
                                     const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                                     IndirectDofSet& rParentDofSet)
{
    KRATOS_TRY
    using SparseUtils = SparseMatrixMultiplicationUtility;

    // Assemble LHS matrix.
    if constexpr (AssembleLHS) {
        BuiltinTimer timer;
        KRATOS_ERROR_IF_NOT(pParentLhs);

        // The restriction operator immediately constructs the linear equivalent,
        // because arbitrary-level coarsening strategies are not supported yet. The problem
        // is that when going from a generic polynomial level q to some other lower polynomial
        // level p!=1, new nodes must be introduced that do not exist in the original fine mesh.
        // This wouldn't be a huge problem by itself, but deciding on DoF indices on the coarse
        // level would be quite painful and require keeping track of the coarse grid's topology
        // in some manner.
        //
        // Note:
        // One might argue that the construction of the restriction operator could happen
        // during the allocation stage, but the problem is that some constraint imposition
        // methods (such as lagrange or its augmented version) may mutate the DoFs during their
        // assembly. As a result, the restriction operator implicitly depends on the constraint
        // imposition method of the fine grid, meaning it must be constructed AFTER constraint
        // assembly.
        MakePRestrictionOperator<std::numeric_limits<unsigned>::max(),typename TSparse::DataType>(
            rModelPart,
            pParentLhs->size1(),
            rParentDofSet,
            mpImpl->mRestrictionOperator,
            mpImpl->mpVariableList,
            mpImpl->mDofSet,
            mpImpl->mIndirectDofSet,
            mpImpl->mDofMap);

        // Compute the coarse LHS matrix.
        typename TSparse::MatrixType left_multiplied_lhs;
        SparseUtils::MatrixMultiplication(mpImpl->mRestrictionOperator, *pParentLhs, left_multiplied_lhs);
        SparseUtils::TransposeMatrix(mpImpl->mProlongationOperator, mpImpl->mRestrictionOperator, 1.0);
        SparseUtils::MatrixMultiplication(left_multiplied_lhs, mpImpl->mProlongationOperator, mpImpl->mLhs);

        KRATOS_INFO_IF("Grid " + std::to_string(mpImpl->mDepth), 2 <= this->mpImpl->mVerbosity)
            << ": grid construction took " << timer << "\n";
    } // if AssembleLHS

    // Assemble RHS vector.
    if constexpr (AssembleRHS) {
        // Sanity checks
        KRATOS_ERROR_IF_NOT(pParentRhs);
        KRATOS_ERROR_IF_NOT(mpImpl->mRestrictionOperator.size2() == pParentRhs->size())
            << "expecting an RHS vector of size " << mpImpl->mRestrictionOperator.size2()
            << " but got " << pParentRhs->size();

        // Allocate the coarse vectors.
        mpImpl->mSolution.resize(mpImpl->mRestrictionOperator.size1(), false);
        mpImpl->mRhs.resize(mpImpl->mRestrictionOperator.size1(), false);
    } // if AssembleRHS

    // Compute coarse constraints.
    // How this is actually implemented unfortunately depends on what the imposition
    // method is. The original master-slave imposition stores a different version of
    // the relation matrix, while the noop imposition stores nothing. Since there are
    // two grid levels here that may use different imposition methods, all combinations
    // must be handled separately.
    switch (rParentConstraintAssembler.GetImposition()) {
        // The parent does not impose constraints, so neither must the current level.
        case ConstraintImposition::None:
            KRATOS_ERROR_IF_NOT(mpImpl->mpConstraintAssembler->GetImposition() == ConstraintImposition::None)
                << "PMultigridBuilderAndSolver: grid " << mpImpl->mDepth
                << " imposes constraints (" << mpImpl->mpConstraintAssembler->GetValue(mpImpl->mpConstraintAssembler->GetImpositionVariable()) << ")"
                << " but its parent does not";
            break;

        // The parent imposes constraints via master-slave elimination.
        // Only linear constraints can be imposed in this manner, so there's
        // no point in computing the Hessian.
        case ConstraintImposition::MasterSlave:
            if (AssembleLHS) {
                SparseUtils::MatrixMultiplication(rParentConstraintAssembler.GetRelationMatrix(),
                                                  mpImpl->mProlongationOperator,
                                                  mpImpl->mpConstraintAssembler->GetRelationMatrix());
            }

            if (AssembleRHS) {
                mpImpl->mpConstraintAssembler->GetConstraintGapVector() = rParentConstraintAssembler.GetConstraintGapVector();
            }
            break;

        // The parent imposes constraints via augmented lagrange multipliers.
        // The only supported imposition on the current level in this case is also
        // augmented lagrange (for similar reasons why master-slave doesn't support anything).
        // There's also a theoretical reason behind the choice of dropping support for
        // master-slave elimination on this level though. Picking augmented lagrange on the
        // parent grid probably means that the user expects that constraints may become
        // linearly dependent at some point. If that's the case, master-slave elimination will
        // definitely break the linear solver while augmented lagrange might provide acceptable
        // results if configured for penalty.
        case ConstraintImposition::AugmentedLagrange:
            if (AssembleLHS) {
                SparseUtils::MatrixMultiplication(rParentConstraintAssembler.GetRelationMatrix(),
                                                  mpImpl->mProlongationOperator,
                                                  mpImpl->mpConstraintAssembler->GetRelationMatrix());

                typename TSparse::MatrixType tmp;
                SparseUtils::MatrixMultiplication(mpImpl->mRestrictionOperator,
                                                  rParentConstraintAssembler.GetHessian(),
                                                  tmp);
                SparseUtils::MatrixMultiplication(tmp,
                                                  mpImpl->mProlongationOperator,
                                                  mpImpl->mpConstraintAssembler->GetHessian());
            }

            if (AssembleRHS) {
                mpImpl->mpConstraintAssembler->GetConstraintGapVector() = rParentConstraintAssembler.GetConstraintGapVector();
            }
            break;

        // No other impositions are supported for now.
        default:
            KRATOS_ERROR << "PMultigridBuilderAndSolver: unsupported constraint imposition at depth " << mpImpl->mDepth
                         << " (parent: " << rParentConstraintAssembler.GetValue(rParentConstraintAssembler.GetImpositionVariable())
                         << " child: " << mpImpl->mpConstraintAssembler->GetValue(mpImpl->mpConstraintAssembler->GetImpositionVariable()) << ")";
    } // switch rParentConstraintAssembler.GetImposition()
    KRATOS_CATCH("")

    KRATOS_TRY
    mpImpl->mpConstraintAssembler->AllocateSystem(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs,
        mpImpl->mIndirectDofSet);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ApplyDirichletConditions(typename IndirectDofSet::const_iterator itParentDofBegin,
                                                     [[maybe_unused]] typename IndirectDofSet::const_iterator itParentDofEnd)
{

    if (mpImpl->mIndirectDofSet.empty()) return;

    // Apply dirichlet conditions on the restriction operator.
    KRATOS_TRY
    block_for_each(mpImpl->mIndirectDofSet.begin(),
                   mpImpl->mIndirectDofSet.end(),
                   [this, itParentDofBegin](const Dof<double>& r_dof){
        const std::size_t i_dof = r_dof.EquationId();
        const typename TSparse::IndexType i_entry_begin = mpImpl->mRestrictionOperator.index1_data()[i_dof];
        const typename TSparse::IndexType i_entry_end = mpImpl->mRestrictionOperator.index1_data()[i_dof + 1];

        if (r_dof.IsFixed()) {
            // Zero out the whole row, except the entry related to the dof on the fine grid.
            const auto i_fine_dof = mpImpl->mDofMap[i_dof];
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = mpImpl->mRestrictionOperator.index2_data()[i_entry];
                if (i_column == i_fine_dof) {
                    mpImpl->mRestrictionOperator.value_data()[i_entry] = static_cast<typename TSparse::DataType>(1);
                } else {
                    mpImpl->mRestrictionOperator.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)
        } /*if r_dof.IsFixed()*/ else {
            // Zero out the column which is associated with the zero'ed row.
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = mpImpl->mRestrictionOperator.index2_data()[i_entry];
                const auto it_column_dof = itParentDofBegin + i_column;
                if (it_column_dof->IsFixed()) {
                    mpImpl->mRestrictionOperator.value_data()[i_entry] = 0.0;
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)
        } /*not r_dof.IsFixed()*/
    });
    KRATOS_CATCH("")

    // Apply dirichlet conditions on the prolongation operator.
    // @todo make this more efficient.
    KRATOS_TRY
    mpImpl->mProlongationOperator = decltype(mpImpl->mProlongationOperator)();
    SparseMatrixMultiplicationUtility::TransposeMatrix(
        mpImpl->mProlongationOperator,
        mpImpl->mRestrictionOperator,
        1.0);
    KRATOS_CATCH("")

    // Apply dirichlet conditions on the LHS.
    KRATOS_TRY
    mpImpl->mpDiagonalScaling->template Cache<TSparse>(mpImpl->mLhs);
    const auto diagonal_scale = mpImpl->mpDiagonalScaling->Evaluate();
    Kratos::ApplyDirichletConditions<TSparse,TDense>(
        mpImpl->mLhs,
        mpImpl->mRhs,
        mpImpl->mIndirectDofSet.begin(),
        mpImpl->mIndirectDofSet.end(),
        diagonal_scale);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ApplyConstraints()
{
    KRATOS_TRY
    mpImpl->mpConstraintAssembler->Initialize(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs,
        mpImpl->mIndirectDofSet);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Initialize(ModelPart& rModelPart,
                                       const typename TParentSparse::MatrixType&,
                                       const typename TParentSparse::VectorType&,
                                       const typename TParentSparse::VectorType& rRhs)
{
    KRATOS_TRY
    if (4 <= mpImpl->mVerbosity) {
        const std::string file_name = "dofs_grid_" + std::to_string(this->mpImpl->mDepth) + ".csv";
        KRATOS_INFO("PMultigridBuilderAndSolver") << "writing DoF data to " << file_name << "\n";
        std::ofstream file(file_name);
        file << "Equation Id,Variable Name,constrained,Node Id,Initial Value,RHS\n";
        for (const Dof<double>& r_dof : mpImpl->mIndirectDofSet)
            file
                << r_dof.EquationId() << ','
                << r_dof.GetVariable().Name() << ','
                << (r_dof.IsFixed() ? 1 : 0) << ','
                << r_dof.Id() << ','
                << r_dof.GetSolutionStepValue() << ','
                << rRhs[r_dof.EquationId()] << '\n';
    }
    if (mpImpl->mpSolver->AdditionalPhysicalDataIsNeeded())
        mpImpl->mpSolver->ProvideAdditionalData(
            mpImpl->mLhs,
            mpImpl->mSolution,
            mpImpl->mRhs,
            mpImpl->mIndirectDofSet,
            rModelPart);
    mpImpl->mpSolver->InitializeSolutionStep(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ExecuteMultigridLoop(PMGStatusStream& rStream,
                                                 PMGStatusStream::Report& rReport)
{
    KRATOS_TRY

    rReport.multigrid_absolute_converged = false;
    rReport.multigrid_relative_converged = false;
    rReport.multigrid_iteration = 0ul;
    rReport.maybe_multigrid_absolute_residual.reset();
    rReport.maybe_multigrid_relative_residual.reset();

    // The multigrid hierarchy depth is currently capped at 1,
    // so the linear solver is used here instead of invoking
    // lower grids.
    rReport.multigrid_relative_converged = mpImpl->mpSolver->PerformSolutionStep(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ExecuteConstraintLoop(PMGStatusStream& rStream,
                                                  PMGStatusStream::Report& rReport)
{
    bool constraints_finished = false;
    auto initial_residual = TSparse::TwoNorm(mpImpl->mRhs);
    initial_residual = initial_residual ? initial_residual : 1;

    // Impose constraints and solve the coarse system.
    KRATOS_TRY
    do {
        rReport.constraints_absolute_converged = false;
        rReport.constraints_relative_converged = false;
        rReport.maybe_constraint_absolute_residual.reset();
        rReport.maybe_constraint_relative_residual.reset();

        // Initialize the constraint assembler.
        mpImpl->mpConstraintAssembler->InitializeConstraintIteration(
            mpImpl->mLhs,
            mpImpl->mSolution,
            mpImpl->mRhs,
            mpImpl->mIndirectDofSet.begin(),
            mpImpl->mIndirectDofSet.end());

        // Get an update on the solution with respect to the current right hand side.
        this->ExecuteMultigridLoop(rStream, rReport);
        constraints_finished = mpImpl->mpConstraintAssembler->FinalizeConstraintIteration(
            mpImpl->mLhs,
            mpImpl->mSolution,
            mpImpl->mRhs,
            mpImpl->mIndirectDofSet.begin(),
            mpImpl->mIndirectDofSet.end(),
            rReport,
            rStream);

        // Update state log.
        if (!constraints_finished) {
            rStream.Submit(rReport.Tag(3), mpImpl->mVerbosity);
            ++rReport.constraint_iteration;
        }
    } while (!constraints_finished);

    // Update the residual.
    BalancedProduct<TSparse,TSparse,TSparse>(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs,
        static_cast<typename TSparse::DataType>(-1));

    // Update state log.
    rReport.maybe_multigrid_absolute_residual = TSparse::TwoNorm(mpImpl->mRhs);
    rReport.maybe_multigrid_relative_residual = rReport.maybe_multigrid_absolute_residual.value() / initial_residual;
    rStream.Submit(rReport.Tag(2), mpImpl->mVerbosity);

    mpImpl->mpConstraintAssembler->Finalize(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs,
        mpImpl->mIndirectDofSet);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
bool PGrid<TSparse,TDense>::ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution,
                                                  const typename TParentSparse::VectorType& rParentRhs,
                                                  const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                                                  PMGStatusStream& rStream)
{
    #ifndef NDEBUG
    KRATOS_TRY
    CheckMatrix<typename TSparse::DataType,
                MatrixChecks::RowsAreSorted
              | MatrixChecks::ColumnsAreSorted>(mpImpl->mRestrictionOperator);
    KRATOS_CATCH("")
    KRATOS_TRY
    CheckMatrix<typename TSparse::DataType,
                MatrixChecks::RowsAreSorted
              | MatrixChecks::ColumnsAreSorted>(mpImpl->mProlongationOperator);
    KRATOS_CATCH("")
    #endif

    PMGStatusStream::Report status_report {
        /*grid_level=*/                 static_cast<std::size_t>(mpImpl->mDepth),
        /*multigrid_converged=*/        false,
        /*multigrid_iteration=*/        0ul,
        /*multigrid_residual=*/         {},
        /*constraints_converged=*/      false,
        /*constraint_iteration=*/       0ul,
        /*maybe_constraint_residual=*/  {}
    };

    TSparse::SetToZero(mpImpl->mSolution);
    this->Restrict<TParentSparse>(mpImpl->mRhs, rParentRhs, rParentConstraintAssembler);
    this->ExecuteConstraintLoop(rStream, status_report);
    this->Prolong<TParentSparse>(rParentSolution, mpImpl->mSolution, rParentConstraintAssembler);

    return (status_report.multigrid_absolute_converged || status_report.multigrid_relative_converged) && (status_report.constraints_absolute_converged || status_report.constraints_relative_converged);
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Finalize(ModelPart& rModelPart,
                                     const typename TParentSparse::MatrixType&,
                                     const typename TParentSparse::VectorType&,
                                     const typename TParentSparse::VectorType&)
{
    mpImpl->mpConstraintAssembler->Finalize(
        mpImpl->mLhs,
        mpImpl->mSolution,
        mpImpl->mRhs,
        mpImpl->mIndirectDofSet);
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Restrict(typename TSparse::VectorType& rCoarseIndependentResidual,
                                     const typename TParentSparse::VectorType& rFineIndependentResidual,
                                     const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler) const
{
    KRATOS_TRY
    // Transform fine residual from independent to dependent space.
    typename TParentSparse::VectorType fine_dependent_residual = rFineIndependentResidual;
    rParentConstraintAssembler.ComputeDependentResidual(fine_dependent_residual);

    // Use the restriction operator to transform the dependent residual to the coarse grid.
    rCoarseIndependentResidual.resize(mpImpl->mRestrictionOperator.size1(), false);
    TSparse::SetToZero(rCoarseIndependentResidual);
    BalancedProduct<TSparse,TParentSparse,TSparse>(
        mpImpl->mRestrictionOperator,
        fine_dependent_residual,
        rCoarseIndependentResidual);

    // Transform the coarse residual from dependent to independent space.
    mpImpl->mpConstraintAssembler->ComputeIndependentResidual(rCoarseIndependentResidual);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Prolong(typename TParentSparse::VectorType& rFineIndependentSolution,
                                    const typename TSparse::VectorType& rCoarseIndependentSolution,
                                    const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler) const
{
    KRATOS_TRY
    // Transform the coarse residual from independent space to dependent space.
    typename TSparse::VectorType dependent_solution = rCoarseIndependentSolution;
    mpImpl->mpConstraintAssembler->ComputeDependentSolution(dependent_solution);

    // Use the prolongation operator to transform the dependent residual to the fine grid.
    rFineIndependentSolution.resize(mpImpl->mProlongationOperator.size1(), false);
    TParentSparse::SetToZero(rFineIndependentSolution);
    BalancedProduct<TSparse,TSparse,TParentSparse>(mpImpl->mProlongationOperator,
                                                   rCoarseIndependentSolution,
                                                   rFineIndependentSolution);

    // Transform the fine residual from dependent to independent space.
    rParentConstraintAssembler.ComputeIndependentSolution(rFineIndependentSolution);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::Clear()
{
    mpImpl->mRestrictionOperator = decltype(mpImpl->mRestrictionOperator)();
    mpImpl->mProlongationOperator = decltype(mpImpl->mProlongationOperator)();
    mpImpl->mLhs = decltype(mpImpl->mLhs)();
    mpImpl->mSolution = decltype(mpImpl->mSolution)();
    mpImpl->mRhs = decltype(mpImpl->mRhs)();
    mpImpl->mIndirectDofSet = decltype(mpImpl->mIndirectDofSet)();
    mpImpl->mDofSet = decltype(mpImpl->mDofSet)();
    mpImpl->mpConstraintAssembler.reset();
    mpImpl->mpSolver.reset();

    if (mpImpl->mMaybeChild.has_value()) {
        mpImpl->mMaybeChild.value().reset();
    }
}


template <class TSparse, class TDense>
Parameters PGrid<TSparse,TDense>::GetDefaultParameters()
{
    Parameters output = Parameters(R"({
        "max_depth" : 0,
        "verbosity" : 1,
        "precision" : "",
        "constraint_imposition_settings" : {
            "method" : "augmented_lagrange",
            "max_iterations" : 1
        }
    })");

    using value_type = typename TSparse::DataType;
    if constexpr (std::is_same_v<value_type, double>) {
        output["precision"].SetString("double");
    } else if constexpr (std::is_same_v<value_type, float>) {
        output["precision"].SetString("single");
    } else {
        static_assert(!std::is_same_v<value_type,value_type>, "unhandled sparse space type");
    }

    return output;
}


template <class TSparse, class TDense>
std::optional<const PGrid<TSparse,TDense>*> PGrid<TSparse,TDense>::GetChild() const {
    return mpImpl->mMaybeChild.has_value()
        ? mpImpl->mMaybeChild.value().get()
        : std::optional<const PGrid*>();
}


template <class TSparse, class TDense>
const typename TSparse::VectorType& PGrid<TSparse,TDense>::GetSolution() const {
    return mpImpl->mSolution;
}


#define KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TParentSparse)                                                        \
    template void PGrid<TSparse,TDense>::MakeLhsTopology<TParentSparse>(ModelPart&,                                             \
                                                                        const typename TParentSparse::MatrixType&,              \
                                                                        const ConstraintAssembler<TParentSparse,TDense>&,       \
                                                                        const PointerVectorSet<Dof<double>>&);                  \
    template void PGrid<TSparse,TDense>::Assemble<false,false,TParentSparse>(ModelPart&,                                        \
                                                                             const typename TParentSparse::MatrixType*,         \
                                                                             const typename TParentSparse::VectorType*,         \
                                                                             const ConstraintAssembler<TParentSparse,TDense>&,  \
                                                                             PointerVectorSet<Dof<double>>&);                   \
    template void PGrid<TSparse,TDense>::Assemble<true,false,TParentSparse>(ModelPart&,                                         \
                                                                            const typename TParentSparse::MatrixType*,          \
                                                                            const typename TParentSparse::VectorType*,          \
                                                                            const ConstraintAssembler<TParentSparse,TDense>&,   \
                                                                            PointerVectorSet<Dof<double>>&);                    \
    template void PGrid<TSparse,TDense>::Assemble<false,true,TParentSparse>(ModelPart&,                                         \
                                                                            const typename TParentSparse::MatrixType*,          \
                                                                            const typename TParentSparse::VectorType*,          \
                                                                            const ConstraintAssembler<TParentSparse,TDense>&,   \
                                                                            PointerVectorSet<Dof<double>>&);                    \
    template void PGrid<TSparse,TDense>::Assemble<true,true,TParentSparse>(ModelPart&,                                          \
                                                                           const typename TParentSparse::MatrixType*,           \
                                                                           const typename TParentSparse::VectorType*,           \
                                                                           const ConstraintAssembler<TParentSparse,TDense>&,    \
                                                                           PointerVectorSet<Dof<double>>&);                     \
    template void PGrid<TSparse,TDense>::Initialize<TParentSparse>(ModelPart&,                                                  \
                                                                   const TParentSparse::MatrixType&,                            \
                                                                   const TParentSparse::VectorType&,                            \
                                                                   const TParentSparse::VectorType&);                           \
    template bool PGrid<TSparse,TDense>::ApplyCoarseCorrection<TParentSparse>(TParentSparse::VectorType&,                       \
                                                                              const TParentSparse::VectorType&,                 \
                                                                              const ConstraintAssembler<TParentSparse,TDense>&, \
                                                                              PMGStatusStream&);                                \
    template void PGrid<TSparse,TDense>::Finalize<TParentSparse>(ModelPart&,                                                    \
                                                                 const TParentSparse::MatrixType&,                              \
                                                                 const TParentSparse::VectorType&,                              \
                                                                 const TParentSparse::VectorType&);                             \
    template void PGrid<TSparse,TDense>::Restrict<TParentSparse>(typename TSparse::VectorType&,                                 \
                                                                 const typename TParentSparse::VectorType&,                     \
                                                                 const ConstraintAssembler<TParentSparse,TDense>&) const;       \
    template void PGrid<TSparse,TDense>::Prolong<TParentSparse>(typename TParentSparse::VectorType&,                            \
                                                                const typename TSparse::VectorType&,                            \
                                                                const ConstraintAssembler<TParentSparse,TDense>&) const

#define KRATOS_INSTANTIATE_PGRID(TSparse, TDense)                                   \
    template class PGrid<TSparse,TDense>;                                           \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<double>);   \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<float>)

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<double>, TUblasDenseSpace<double>);

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<float>, TUblasDenseSpace<double>);

#undef KRATOS_INSTANTIATE_GRID
#undef KRATOS_INSTANTIATE_GRID_MEMBERS


} // namespace Kratos
