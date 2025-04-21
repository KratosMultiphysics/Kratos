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

// System includes
#include <limits> // std::numeric_limits


namespace Kratos {


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid(Parameters Settings,
                             const unsigned CurrentDepth,
                             Parameters SmootherSettings,
                             Parameters LeafSolverSettings,
                             Parameters DiagonalScalingSettings)
    : mRestrictionOperator(),
      mProlongationOperator(),
      mLhs(),
      mSolution(),
      mRhs(),
      mDofSet(),
      mIndirectDofSet(),
      mDofMap(),
      mpConstraintAssembler(),
      mpDiagonalScaling(),
      mpSolver(),
      mMaybeChild(),
      mVerbosity(),
      mDepth(CurrentDepth)
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(mDepth) << "PGrid can only have positive depths (the original system is at depth=0)";

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

    mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(Settings["constraint_imposition_settings"],
                                                                       "Grid " + std::to_string(mDepth) + " constraints");
    mpDiagonalScaling = std::make_unique<Scaling>(DiagonalScalingSettings);
    mVerbosity = Settings["verbosity"].Get<int>();

    const int max_depth = Settings["max_depth"].Get<int>();
    KRATOS_ERROR_IF_NOT(0 <= max_depth) << Settings << "\n\"max_depth\" must be non-negative";
    if (mDepth < static_cast<unsigned>(max_depth)) {
        KRATOS_TRY
        KRATOS_ERROR_IF_NOT(SmootherSettings.Has("solver_type"));
        const std::string solver_name = SmootherSettings["solver_type"].Get<std::string>();
        using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;
        KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
            << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
            << "Make sure you imported the application it is defined in and that the spelling is correct.";
        const auto& r_factory = SolverFactoryRegistry::Get(solver_name);
        mpSolver = r_factory.Create(SmootherSettings);
        mMaybeChild = std::unique_ptr<PGrid>(new PGrid(Settings,
                                                       mDepth + 1u,
                                                       SmootherSettings,
                                                       LeafSolverSettings,
                                                       DiagonalScalingSettings));
        KRATOS_CATCH(std::to_string(mDepth + 1u));
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
        mpSolver = r_factory.Create(LeafSolverSettings);
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
            mRestrictionOperator,
            mpVariableList,
            mDofSet,
            mIndirectDofSet,
            mDofMap);

        // Compute the coarse LHS matrix.
        typename TSparse::MatrixType left_multiplied_lhs;
        SparseUtils::MatrixMultiplication(mRestrictionOperator, *pParentLhs, left_multiplied_lhs);
        SparseUtils::TransposeMatrix(mProlongationOperator, mRestrictionOperator, 1.0);
        SparseUtils::MatrixMultiplication(left_multiplied_lhs, mProlongationOperator, mLhs);

        KRATOS_INFO_IF("Grid " + std::to_string(mDepth), 2 <= this->mVerbosity)
            << ": grid construction took " << timer << "\n";
    } // if AssembleLHS

    // Assemble RHS vector.
    if constexpr (AssembleRHS) {
        // Sanity checks
        KRATOS_ERROR_IF_NOT(pParentRhs);
        KRATOS_ERROR_IF_NOT(mRestrictionOperator.size2() == pParentRhs->size())
            << "expecting an RHS vector of size " << mRestrictionOperator.size2()
            << " but got " << pParentRhs->size();

        // Allocate the coarse vectors.
        mSolution.resize(mRestrictionOperator.size1(), false);
        mRhs.resize(mRestrictionOperator.size1(), false);
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
            KRATOS_ERROR_IF_NOT(mpConstraintAssembler->GetImposition() == ConstraintImposition::None)
                << "PMultigridBuilderAndSolver: grid " << mDepth
                << " imposes constraints (" << mpConstraintAssembler->GetValue(mpConstraintAssembler->GetImpositionVariable()) << ")"
                << " but its parent does not";
            break;

        // The parent imposes constraints via master-slave elimination. The original
        // implementation stores a transformation matrix instead of the relation matrix,
        // so coarsening it is complicated and expensive => no support for it for now.
        /// @todo implement constraint imposition when the parent on a p-grid uses
        //        master-slave imposition but the child grid uses augmented lagrange.
        case ConstraintImposition::MasterSlave:
            KRATOS_ERROR << "PMultigridBuilderAndSolver: constraint imposition using master-slave elimination "
                         << "is only supported if there is no coarse hierarchy (depth=0). Consider imposition "
                         << "using augmented lagrange multipliers.";
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
                                                  mProlongationOperator,
                                                  mpConstraintAssembler->GetRelationMatrix());

                typename TSparse::MatrixType tmp;
                SparseUtils::MatrixMultiplication(mRestrictionOperator,
                                                  rParentConstraintAssembler.GetHessian(),
                                                  tmp);
                SparseUtils::MatrixMultiplication(tmp,
                                                  mProlongationOperator,
                                                  mpConstraintAssembler->GetHessian());
            }

            if (AssembleRHS) {
                mpConstraintAssembler->GetConstraintGapVector() = rParentConstraintAssembler.GetConstraintGapVector();
            }
            break;

        // No other impositions are supported for now.
        default:
            KRATOS_ERROR << "PMultigridBuilderAndSolver: unsupported constraint imposition at depth " << mDepth
                         << " (parent: " << rParentConstraintAssembler.GetValue(rParentConstraintAssembler.GetImpositionVariable())
                         << " child: " << mpConstraintAssembler->GetValue(mpConstraintAssembler->GetImpositionVariable()) << ")";
    } // switch rParentConstraintAssembler.GetImposition()

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ApplyDirichletConditions(typename IndirectDofSet::const_iterator itParentDofBegin,
                                                     [[maybe_unused]] typename IndirectDofSet::const_iterator itParentDofEnd)
{

    if (mIndirectDofSet.empty()) return;

    // Apply dirichlet conditions on the restriction operator.
    KRATOS_TRY
    block_for_each(mIndirectDofSet.begin(),
                   mIndirectDofSet.end(),
                   [this, itParentDofBegin](const Dof<double>& r_dof){
        const std::size_t i_dof = r_dof.EquationId();
        const typename TSparse::IndexType i_entry_begin = mRestrictionOperator.index1_data()[i_dof];
        const typename TSparse::IndexType i_entry_end = mRestrictionOperator.index1_data()[i_dof + 1];

        if (r_dof.IsFixed()) {
            // Zero out the whole row, except the entry related to the dof on the fine grid.
            const auto i_fine_dof = mDofMap[i_dof];
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = mRestrictionOperator.index2_data()[i_entry];
                if (i_column == i_fine_dof) {
                    mRestrictionOperator.value_data()[i_entry] = static_cast<typename TSparse::DataType>(1);
                } else {
                    mRestrictionOperator.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)
        } /*if r_dof.IsFixed()*/ else {
            // Zero out the column which is associated with the zero'ed row.
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = mRestrictionOperator.index2_data()[i_entry];
                const auto it_column_dof = itParentDofBegin + i_column;
                if (it_column_dof->IsFixed()) {
                    mRestrictionOperator.value_data()[i_entry] = 0.0;
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)
        } /*not r_dof.IsFixed()*/
    });
    KRATOS_CATCH("")

    // Apply dirichlet conditions on the prolongation operator.
    // @todo make this more efficient.
    KRATOS_TRY
    mProlongationOperator = decltype(mProlongationOperator)();
    SparseMatrixMultiplicationUtility::TransposeMatrix(mProlongationOperator, mRestrictionOperator, 1.0);
    KRATOS_CATCH("")

    // Apply dirichlet conditions on the LHS.
    KRATOS_TRY
    mpDiagonalScaling->Cache<TSparse>(mLhs);
    const auto diagonal_scale = mpDiagonalScaling->Evaluate();
    Kratos::ApplyDirichletConditions<TSparse,TDense>(
        mLhs,
        mRhs,
        mIndirectDofSet.begin(),
        mIndirectDofSet.end(),
        diagonal_scale);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ApplyConstraints()
{
    KRATOS_TRY
    mpConstraintAssembler->Initialize(mLhs,
                                      mRhs,
                                      mIndirectDofSet.begin(),
                                      mIndirectDofSet.end());
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Initialize(ModelPart& rModelPart,
                                       const typename TParentSparse::MatrixType&,
                                       const typename TParentSparse::VectorType&,
                                       const typename TParentSparse::VectorType&)
{
    KRATOS_TRY
    if (mpSolver->AdditionalPhysicalDataIsNeeded())
        mpSolver->ProvideAdditionalData(mLhs,
                                        mSolution,
                                        mRhs,
                                        mIndirectDofSet,
                                        rModelPart);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ExecuteMultigridLoop(PMGStatusStream& rStream,
                                                 PMGStatusStream::Report& rReport)
{
    KRATOS_TRY

    rReport.multigrid_converged = false;
    rReport.multigrid_iteration = 0ul;
    rReport.multigrid_residual = std::numeric_limits<typename TSparse::DataType>::max();

    // The multigrid hierarchy depth is currently capped at 1,
    // so the linear solver is used here instead of invoking
    // lower grids.
    mpSolver->InitializeSolutionStep(mLhs, mSolution, mRhs);
    rReport.multigrid_converged = mpSolver->Solve(mLhs, mSolution, mRhs);
    mpSolver->FinalizeSolutionStep(mLhs, mSolution, mRhs);

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::ExecuteConstraintLoop(PMGStatusStream& rStream,
                                                  PMGStatusStream::Report& rReport)
{
    bool constraints_finished = false;
    auto initial_residual = TSparse::TwoNorm(mRhs);
    initial_residual = initial_residual ? initial_residual : 1;

    // Impose constraints and solve the coarse system.
    KRATOS_TRY
    do {
        // Initialize the constraint assembler.
        mpConstraintAssembler->InitializeSolutionStep(mLhs, mSolution, mRhs);

        // Get an update on the solution with respect to the current right hand side.
        this->ExecuteMultigridLoop(rStream, rReport);

        constraints_finished = mpConstraintAssembler->FinalizeSolutionStep(mLhs, mSolution, mRhs, rReport, rStream);

        // Update state log.
        if (!constraints_finished)
            rStream.Submit(rReport.Tag(3), mVerbosity);
    } while (!constraints_finished);

    // Update the residual.
    BalancedProduct<TSparse,TSparse,TSparse>(mLhs, mSolution, mRhs, static_cast<typename TSparse::DataType>(-1));

    // Update state log.
    rReport.multigrid_residual = TSparse::TwoNorm(mRhs) / initial_residual;
    rStream.Submit(rReport.Tag(2), mVerbosity);

    mpConstraintAssembler->Finalize(mLhs, mSolution, mRhs, mIndirectDofSet);
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
                MatrixChecks::All>(mRestrictionOperator);
    CheckMatrix<typename TSparse::DataType,
                MatrixChecks::RowsAreSorted
              | MatrixChecks::ColumnsAreSorted>(mProlongationOperator);
    KRATOS_CATCH("")
    #endif

    PMGStatusStream::Report status_report {
        /*grid_level=*/                 static_cast<std::size_t>(mDepth),
        /*multigrid_converged=*/        false,
        /*multigrid_iteration=*/        0ul,
        /*multigrid_residual=*/         1.0,
        /*constraints_converged=*/      false,
        /*constraint_iteration=*/       0ul,
        /*maybe_constraint_residual=*/  {}
    };

    TSparse::SetToZero(mSolution);
    this->Restrict<TParentSparse>(rParentRhs, mRhs);
    this->ExecuteConstraintLoop(rStream, status_report);
    this->Prolong<TParentSparse>(mSolution, rParentSolution);

    return status_report.multigrid_converged && status_report.constraints_converged;
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Finalize(ModelPart& rModelPart,
                                     const typename TParentSparse::MatrixType&,
                                     const typename TParentSparse::VectorType&,
                                     const typename TParentSparse::VectorType&)
{
    mpConstraintAssembler->Finalize(mLhs, mSolution, mRhs, mIndirectDofSet);
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Restrict(const typename TParentSparse::VectorType& rParentVector,
                                     typename TSparse::VectorType& rOutput) const
{
    KRATOS_TRY
    rOutput.resize(mRestrictionOperator.size1(), false),
    TSparse::SetToZero(rOutput);
    BalancedProduct<TSparse,TParentSparse,TSparse>(mRestrictionOperator,
                                                   rParentVector,
                                                   rOutput);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Prolong(const typename TSparse::VectorType& rInput,
                                    typename TParentSparse::VectorType& rParentVector) const
{
    KRATOS_TRY
    rParentVector.resize(mProlongationOperator.size1(), false);
    TParentSparse::SetToZero(rParentVector);
    BalancedProduct<TSparse,TSparse,TParentSparse>(mProlongationOperator,
                                                   rInput,
                                                   rParentVector);
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void PGrid<TSparse,TDense>::Clear()
{
    mRestrictionOperator = decltype(mRestrictionOperator)();
    mProlongationOperator = decltype(mProlongationOperator)();
    mLhs = decltype(mLhs)();
    mSolution = decltype(mSolution)();
    mRhs = decltype(mRhs)();
    mIndirectDofSet = decltype(mIndirectDofSet)();
    mDofSet = decltype(mDofSet)();
    mpConstraintAssembler.reset();
    mpSolver.reset();

    if (mMaybeChild.has_value()) {
        mMaybeChild.value().reset();
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
    template void PGrid<TSparse,TDense>::Restrict<TParentSparse>(const typename TParentSparse::VectorType&,                     \
                                                                 typename TSparse::VectorType&) const;                          \
    template void PGrid<TSparse,TDense>::Prolong<TParentSparse>(const typename TSparse::VectorType&,                            \
                                                                 typename TParentSparse::VectorType&) const

#define KRATOS_INSTANTIATE_PGRID(TSparse, TDense)                                   \
    template class PGrid<TSparse,TDense>;                                           \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<double>);   \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<float>)

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<double>, TUblasDenseSpace<double>);

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<float>, TUblasDenseSpace<double>);

#undef KRATOS_INSTANTIATE_GRID
#undef KRATOS_INSTANTIATE_GRID_MEMBERS


} // namespace Kratos
