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
        KRATOS_ERROR_IF_NOT(SolverFactoryRegistry::Has(solver_name))
            << "\"" << solver_name << "\" is not a valid linear solver name in the registry. "
            << "Make sure you imported the application it is defined in and that the spelling is correct.";
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
    MakePRestrictionOperator<
        std::numeric_limits<unsigned>::max(),
        typename TSparse::DataType>(rModelPart,
                                    rParentLhs.size1(),
                                    rParentDofSet,
                                    mRestrictionOperator,
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

        // Compute the coarse relation matrix.
        SparseUtils::MatrixMultiplication(rParentConstraintAssembler.GetRelationMatrix(),
                                          mProlongationOperator,
                                          mpConstraintAssembler->GetRelationMatrix());
    } // if AssembleLHS

    if constexpr (AssembleRHS) {
        // Sanity checks
        KRATOS_ERROR_IF_NOT(pParentRhs);
        KRATOS_ERROR_IF_NOT(mRestrictionOperator.size2() == pParentRhs->size())
            << "expecting an RHS vector of size " << mRestrictionOperator.size2()
            << " but got " << pParentRhs->size();

        // Compute the coarse constraint gaps.
        mpConstraintAssembler->GetConstraintGapVector() = rParentConstraintAssembler.GetConstraintGapVector();

        // Allocate the coarse vectors.
        mSolution.resize(mRestrictionOperator.size1(), false);
        mRhs.resize(mRestrictionOperator.size1(), false);
    } // if AssembleRHS

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Initialize(ModelPart& rModelPart,
                                       const typename TParentSparse::MatrixType& rParentLhs,
                                       const typename TParentSparse::VectorType& rParentSolution,
                                       const typename TParentSparse::VectorType& rParentRhs)
{
    KRATOS_TRY
    mpConstraintAssembler->Initialize(mLhs, mRhs);
    if (mpSolver->AdditionalPhysicalDataIsNeeded()) {
        mpSolver->ProvideAdditionalData(mLhs,
                                        mSolution,
                                        mRhs,
                                        mIndirectDofSet,
                                        rModelPart);
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
template <class TParentSparse>
bool PGrid<TSparse,TDense>::ApplyCoarseCorrection(typename TParentSparse::VectorType& rParentSolution,
                                                  const typename TParentSparse::VectorType& rParentRhs)
{
    // Restrict the residual from the fine grid to the coarse one (this grid).
    KRATOS_TRY
    // This is a matrix-vector product of potentially different value types. In its current state,
    // sparse spaces do not support computing the products of arguments with different value types,
    // so I'm directly invoking the UBLAS template that does support it.
    TSparse::SetToZero(mRhs);
    TSparse::SetToZero(mSolution);
    axpy_prod(mRestrictionOperator, rParentRhs, mRhs, true);
    axpy_prod(mRestrictionOperator, rParentSolution, mSolution, true);
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
    axpy_prod(mProlongationOperator, mSolution, rParentSolution, true);
    KRATOS_CATCH("")

    return linear_solver_status && constraint_status.converged;
}


template <class TSparse, class TDense>
template <class TParentSparse>
void PGrid<TSparse,TDense>::Finalize([[maybe_unused]] ModelPart& rModelPart,
                                     [[maybe_unused]] const typename TParentSparse::MatrixType& rParentLhs,
                                     [[maybe_unused]] const typename TParentSparse::VectorType& rParentSolution,
                                     [[maybe_unused]] const typename TParentSparse::VectorType& rParentRhs)
{
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
