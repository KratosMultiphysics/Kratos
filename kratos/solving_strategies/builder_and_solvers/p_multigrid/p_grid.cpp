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
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "factories/linear_solver_factory.h" // LinearSolverFactory
#include "includes/kratos_components.h" // KratosComponents
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE
#include "utilities/sparse_matrix_multiplication_utility.h" // SparseMatrixMultiplicationUtility

// System includes
#include <limits> // std::numeric_limits


namespace Kratos {


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid(Parameters Settings,
                             const unsigned CurrentDepth,
                             Parameters SmootherSettings,
                             Parameters LeafSolverSettings)
    : mRestrictionOperator(),
      mProlongationOperator(),
      mLhs(),
      mSolution(),
      mRhs(),
      mDofSet(),
      mIndirectDofSet(),
      mpConstraintAssembler(),
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

    mpConstraintAssembler = ConstraintAssemblerFactory<TSparse,TDense>(Settings["constraint_imposition_settings"]);
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
                                                       LeafSolverSettings));
        KRATOS_CATCH(std::to_string(mDepth + 1u));
    } else {
        KRATOS_ERROR_IF_NOT(LeafSolverSettings.Has("solver_type"));
        const std::string solver_name = LeafSolverSettings["solver_type"].Get<std::string>();
        using SolverFactoryRegistry = KratosComponents<LinearSolverFactory<TSparse,TDense>>;

        if (not SolverFactoryRegistry::Has(solver_name)) {
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
                             Parameters LeafSolverSettings)
    : PGrid(Settings,
            1u,
            SmootherSettings,
            LeafSolverSettings)
{
}


template <class TSparse, class TDense>
PGrid<TSparse,TDense>::PGrid()
    : PGrid(Parameters(),
            Parameters(R"({"solver_type" : "gauss_seidel"})"),
            Parameters(R"({"solver_type" : "amgcl"})"))
{
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::MakeLhsTopology(ModelPart& rModelPart,
                                            const typename TParentSparse::MatrixType& rParentLhs,
                                            [[maybe_unused]] const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler,
                                            const IndirectDofSet& rParentDofSet)
{
    KRATOS_TRY

    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);
    // The restriction operator immediately constructs the linear equivalent,
    // because arbitrary-level coarsening strategies are not supported yet. The problem
    // is that when going from a generic polynomial level q to some other lower polynomial
    // level p!=1, new nodes must be introduced that do not exist in the original fine mesh.
    // This wouldn't be a huge problem by itself, but deciding on DoF indices on the coarse
    // level would be quite painful and require keeping track of the coarse grid's topology
    // in some manner.
    MakePRestrictionOperator<std::numeric_limits<unsigned>::max(),typename TSparse::DataType>(
        rModelPart,
        rParentLhs.size1(),
        rParentDofSet,
        mRestrictionOperator,
        mpVariableList,
        mDofSet,
        mIndirectDofSet);

    // ConstraintAssembler::Allocate is deliberately not invoked here because the relation matrix
    // and constraint gap vector are passed from the fine level in PGrid::Assemble, and mapped
    // to this level by the restriction operator.

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <bool AssembleLHS,
          bool AssembleRHS,
          class TParentSparse>
void PGrid<TSparse,TDense>::Assemble(const ModelPart& rModelPart,
                                     const typename TParentSparse::MatrixType* pParentLhs,
                                     const typename TParentSparse::VectorType* pParentRhs,
                                     const ConstraintAssembler<TParentSparse,TDense>& rParentConstraintAssembler)
{
    KRATOS_TRY
    using SparseUtils = SparseMatrixMultiplicationUtility;

    if constexpr (AssembleLHS) {
        KRATOS_ERROR_IF_NOT(pParentLhs);

        // Compute the coarse LHS matrix.
        typename TSparse::MatrixType left_multiplied_lhs;
        SparseUtils::MatrixMultiplication(mRestrictionOperator, *pParentLhs, left_multiplied_lhs);
        SparseUtils::TransposeMatrix(mProlongationOperator, mRestrictionOperator, 1.0);
        SparseUtils::MatrixMultiplication(left_multiplied_lhs, mProlongationOperator, mLhs);
    } // if AssembleLHS


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
                                                     typename IndirectDofSet::const_iterator itParentDofEnd,
                                                     const DiagonalScaling& rDiagonalScaling)
{
    KRATOS_TRY

    if (mIndirectDofSet.empty()) return;

    Kratos::ApplyDirichletConditions<TSparse,TDense>(
        mLhs,
        mRhs,
        mIndirectDofSet.begin(),
        mIndirectDofSet.end(),
        GetDiagonalScaleFactor<TSparse>(mLhs, rDiagonalScaling));

    Kratos::ApplyDirichletConditions<TSparse,TDense>(
        mRestrictionOperator,
        mRhs,
        mIndirectDofSet.begin(),
        mIndirectDofSet.end(),
        itParentDofBegin,
        itParentDofEnd,
        static_cast<typename TDense::DataType>(1));

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
template <class TParentSparse>
bool PGrid<TSparse,TDense>::ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution,
                                                  const typename TParentSparse::VectorType& rParentRhs)
{
    #ifndef NDEBUG
    KRATOS_TRY
    CheckMatrix<typename TSparse::DataType,MatrixChecks::All>(mRestrictionOperator);
    CheckMatrix<typename TSparse::DataType,MatrixChecks::RowsAreSorted|MatrixChecks::ColumnsAreSorted>(mProlongationOperator);
    KRATOS_CATCH("")
    #endif

    // Restrict the residual from the fine grid to the coarse one (this grid).
    KRATOS_TRY
    // This is a matrix-vector product of potentially different value types. In its current state,
    // sparse spaces do not support computing the products of arguments with different value types,
    // so I'm directly invoking the UBLAS template that does support it.
    TSparse::SetToZero(mRhs);
    TSparse::SetToZero(mSolution);
    axpy_prod(mRestrictionOperator, rParentRhs, mRhs, boost::numeric::ublas::row_major_tag());
    //for (std::size_t i_dof=0ul; i_dof<mIndirectDofSet.size(); ++i_dof) {
    //    const auto& r_dof = *(mIndirectDofSet.begin() + i_dof);
    //    if (r_dof.IsFixed()) {
    //        const auto equation_id = r_dof.EquationId();
    //        mRhs[equation_id] = 0;
    //    }
    //}
    KRATOS_CATCH("")

    // Impose constraints and solve the coarse system.
    std::size_t i_iteration = 0ul;
    typename ConstraintAssembler<TSparse,TDense>::Status constraint_status {/*converged=*/false, /*finished=*/false};
    bool linear_solver_status = false; //< Indicates whether the linear solver converged.

    KRATOS_TRY
    do {
        ++i_iteration;
        mpConstraintAssembler->InitializeSolutionStep(mLhs, mSolution, mRhs, i_iteration);
        mpSolver->InitializeSolutionStep(mLhs, mSolution, mRhs);
        linear_solver_status = mpSolver->Solve(mLhs, mSolution, mRhs);
        mpSolver->FinalizeSolutionStep(mLhs, mSolution, mRhs);
        constraint_status = mpConstraintAssembler->FinalizeSolutionStep(mLhs, mSolution, mRhs, i_iteration);
    } while (not constraint_status.finished);

    // Emit status.
    if (1 <= mVerbosity) {
        if (not linear_solver_status)
            std::cerr << "PMultigridBuilderAndSolver: grid " << mDepth << ": failed to converge\n";

        if (not constraint_status.converged)
            std::cerr << "PMultigridBuilderAndSolver: grid " << mDepth << ": failed to converge constraints\n";
    }

    mpConstraintAssembler->Finalize(mLhs, mSolution, mRhs, mIndirectDofSet);
    KRATOS_CATCH("")

    // Prolong the coarse solution to the fine grid.
    KRATOS_TRY
    // This is a matrix-vector product of potentially different value types. In its current state,
    // sparse spaces do not support computing the products of arguments with different value types,
    // so I'm directly invoking the UBLAS template that does support it.
    TParentSparse::SetToZero(rParentSolution);
    axpy_prod(mProlongationOperator, mSolution, rParentSolution, boost::numeric::ublas::row_major_tag());
    KRATOS_CATCH("")

    return linear_solver_status and constraint_status.converged;
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
    template void PGrid<TSparse,TDense>::Assemble<false,false,TParentSparse>(const ModelPart&,                                  \
                                                                             const typename TParentSparse::MatrixType*,         \
                                                                             const typename TParentSparse::VectorType*,         \
                                                                             const ConstraintAssembler<TParentSparse,TDense>&); \
    template void PGrid<TSparse,TDense>::Assemble<true,false,TParentSparse>(const ModelPart&,                                   \
                                                                            const typename TParentSparse::MatrixType*,          \
                                                                            const typename TParentSparse::VectorType*,          \
                                                                            const ConstraintAssembler<TParentSparse,TDense>&);  \
    template void PGrid<TSparse,TDense>::Assemble<false,true,TParentSparse>(const ModelPart&,                                   \
                                                                            const typename TParentSparse::MatrixType*,          \
                                                                            const typename TParentSparse::VectorType*,          \
                                                                            const ConstraintAssembler<TParentSparse,TDense>&);  \
    template void PGrid<TSparse,TDense>::Assemble<true,true,TParentSparse>(const ModelPart&,                                    \
                                                                           const typename TParentSparse::MatrixType*,           \
                                                                           const typename TParentSparse::VectorType*,           \
                                                                           const ConstraintAssembler<TParentSparse,TDense>&);   \
    template void PGrid<TSparse,TDense>::Initialize<TParentSparse>(ModelPart&,                                                  \
                                                                   const TParentSparse::MatrixType&,                            \
                                                                   const TParentSparse::VectorType&,                            \
                                                                   const TParentSparse::VectorType&);                           \
    template bool PGrid<TSparse,TDense>::ApplyCoarseCorrection<TParentSparse>(TParentSparse::VectorType&,                       \
                                                                              const TParentSparse::VectorType&);                \
    template void PGrid<TSparse,TDense>::Finalize<TParentSparse>(ModelPart&,                                                    \
                                                                 const TParentSparse::MatrixType&,                              \
                                                                 const TParentSparse::VectorType&,                              \
                                                                 const TParentSparse::VectorType&)

#define KRATOS_INSTANTIATE_PGRID(TSparse, TDense)                                   \
    template class PGrid<TSparse,TDense>;                                           \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<double>);   \
    KRATOS_INSTANTIATE_PGRID_MEMBERS(TSparse, TDense, TUblasSparseSpace<float>)

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<double>, TUblasDenseSpace<double>);

KRATOS_INSTANTIATE_PGRID(TUblasSparseSpace<float>, TUblasDenseSpace<double>);

#undef KRATOS_INSTANTIATE_GRID
#undef KRATOS_INSTANTIATE_GRID_MEMBERS


} // namespace Kratos
