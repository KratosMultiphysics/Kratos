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

// External includes
#include "mcgs/mcgs.hpp" // mcgs::Partition, mcgs::CSRAdaptor, mcgs::ColorSettings, mcgs::SolveSettings

// Project includes
#include "custom_solvers/gauss_seidel.h" // GaussSeidelSmoother
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace

// System includes
#include <memory> // std::unique_ptr


namespace Kratos {


template <class TSparse, class TDense>
struct GaussSeidelSmoother<TSparse,TDense>::Impl
{
    using CSRAdaptor = mcgs::CSRAdaptor<typename TSparse::IndexType,
                                        typename TSparse::DataType>;

    using Partition = mcgs::Partition<typename TSparse::IndexType>;

    Impl()
        : mpRowPartition(nullptr),
          mpReorderedRowPartition(nullptr),
          mColorSettings(),
          mSolveSettings(),
          mReformDofsAtEachStep(false)
    {
    }

    ~Impl()
    {
        if (mpReorderedRowPartition) {
            mcgs::destroyPartition(mpReorderedRowPartition);
        }

        if (mpRowPartition) {
            mcgs::destroyPartition(mpRowPartition);
        }
    }

    /// @brief Pointer to a partition of rows in the LHS matrix with respect to a coloring.
    /// @details This pointer is allocated by MCGS and must be deallocated by invoking
    ///          @p mcgs::destroyPartition.
    Partition* mpRowPartition;

    /// @brief Pointer to a reordered partition of rows in the LHS matrix with respect to a coloring.
    /// @details This pointer is allocated by MCGS and must be deallocated by invoking
    ///          @p mcgs::destroyPartition.
    Partition* mpReorderedRowPartition;

    mcgs::ColorSettings<typename TSparse::DataType> mColorSettings;

    mcgs::SolveSettings<typename TSparse::IndexType,typename TSparse::DataType> mSolveSettings;

    bool mReformDofsAtEachStep;
}; // struct GaussSeidelSmoother::Impl


template <class TSparse, class TDense>
GaussSeidelSmoother<TSparse,TDense>::GaussSeidelSmoother()
    : GaussSeidelSmoother(GetDefaultParameters())
{
}


template <class TSparse, class TDense>
GaussSeidelSmoother<TSparse,TDense>::GaussSeidelSmoother(Parameters Settings)
    : mpImpl(new Impl)
{
    KRATOS_TRY

    Settings.RecursivelyValidateAndAssignDefaults(this->GetDefaultParameters());

    // Set coloring settings.
    const Parameters color_parameters = Settings["coloring_settings"];
    mcgs::ColorSettings<typename TSparse::DataType> color_settings;
    color_settings.shrinkingFactor = color_parameters["shrinking_factor"].Get<int>();
    color_settings.maxStallCount   = color_parameters["max_stall_count"].Get<int>();
    color_settings.tolerance       = color_parameters["tolerance"].Get<double>();
    color_settings.verbosity       = color_parameters["verbosity"].Get<int>();
    mpImpl->mColorSettings = color_settings;

    // Set smoothing settings.
    mcgs::SolveSettings<typename TSparse::IndexType,typename TSparse::DataType> solve_settings;
    solve_settings.residualAbsoluteTolerance = Settings["absolute_tolerance"].Get<double>();
    solve_settings.residualRelativeTolerance = Settings["relative_tolerance"].Get<double>();
    solve_settings.maxIterations             = Settings["max_iterations"].Get<int>();
    solve_settings.relaxation                = Settings["relaxation_factor"].Get<double>();
    solve_settings.verbosity                 = Settings["verbosity"].Get<int>();

    const std::string parallelization = Settings["parallelization"].Get<std::string>();
    if (parallelization == "row-wise") {
        solve_settings.parallelization = mcgs::Parallelization::RowWise;
    } else if (parallelization == "entrywise") {
        solve_settings.parallelization = mcgs::Parallelization::EntryWise;
    } else if (parallelization == "none") {
        solve_settings.parallelization = mcgs::Parallelization::None;
    } else {
        KRATOS_ERROR
            << "invalid parallelization strategy for GassSeidelSmoother: "
            << '\"' << parallelization << "\"."
            << " Options are " << "\"row-wise\", \"entrywise\" or \"none\".";
    }
    mpImpl->mSolveSettings = solve_settings;

    // Flag indicating whether the topology changes between solution steps.
    mpImpl->mReformDofsAtEachStep = Settings["reform_dofs_at_each_step"].Get<bool>();

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
GaussSeidelSmoother<TSparse,TDense>::~GaussSeidelSmoother() = default;


template <class TSparse, class TDense>
void GaussSeidelSmoother<TSparse,TDense>::Initialize(SparseMatrixType& rLhs,
                                                     VectorType& rSolution,
                                                     VectorType& rRhs)
{
    // Compute a coloring of the LHS matrix, but only if parallelization is requested.
    KRATOS_TRY
    if (mpImpl->mSolveSettings.parallelization != mcgs::Parallelization::None) {
        // Destroy the partition if it exists already.
        /// @todo Do not destroy existing partitions if the linear system is not
        ///       rebuilt at each iteration.
        if (mpImpl->mpRowPartition) {
            mcgs::destroyPartition(mpImpl->mpRowPartition);
        }

        if (mpImpl->mpReorderedRowPartition) {
            mcgs::destroyPartition(mpImpl->mpReorderedRowPartition);
        }

        // Construct a view of the left hand side matrix.
        typename Impl::CSRAdaptor lhs_adaptor;
        lhs_adaptor.rowCount        = rLhs.size1();
        lhs_adaptor.columnCount     = rLhs.size2();
        lhs_adaptor.entryCount      = rLhs.nnz();
        lhs_adaptor.pRowExtents     = rLhs.index1_data().begin();
        lhs_adaptor.pColumnIndices  = rLhs.index2_data().begin();
        lhs_adaptor.pEntries        = rLhs.value_data().begin();

        // Compute coloring.
        std::vector<unsigned> colors(lhs_adaptor.rowCount);
        if (mcgs::color(colors.data(), lhs_adaptor, mpImpl->mColorSettings) != MCGS_SUCCESS) {
            KRATOS_ERROR << "coloring for multicolor gauss-seidel smoother failed";
        } // if coloring failed

        // Compute partition based on coloring.
        mpImpl->mpRowPartition = mcgs::makePartition(colors.data(), lhs_adaptor.rowCount);
        if (!mpImpl->mpRowPartition) {
            KRATOS_ERROR << "computing partition for multicolor gauss-seidel smoother failed";
        }

        // Reorder the system with respect to the computed coloring.
        if (!(mpImpl->mpReorderedRowPartition = mcgs::reorder(rLhs.size1(),
                                                              rLhs.size2(),
                                                              rLhs.nnz(),
                                                              &*rLhs.index1_data().begin(),
                                                              &*rLhs.index2_data().begin(),
                                                              &*rLhs.value_data().begin(),
                                                              &*rSolution.begin(),
                                                              &*rRhs.begin(),
                                                              mpImpl->mpRowPartition))) {
            KRATOS_ERROR << "reordering for multicolor gauss-seidel failed";
        }
    } // if parallelization != None
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void GaussSeidelSmoother<TSparse,TDense>::InitializeSolutionStep(SparseMatrixType& rLhs,
                                                                 VectorType& rSolution,
                                                                 VectorType& rRhs)
{
    KRATOS_TRY
    bool needs_initialization = mpImpl->mReformDofsAtEachStep;
    needs_initialization |= (!mpImpl->mpRowPartition || !mpImpl->mpReorderedRowPartition)
                            && mpImpl->mSolveSettings.parallelization != mcgs::Parallelization::None;

    if (needs_initialization) {
        this->Initialize(rLhs, rSolution, rRhs);
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
bool GaussSeidelSmoother<TSparse,TDense>::Solve(SparseMatrixType& rLhs,
                                                VectorType& rSolution,
                                                VectorType& rRhs)
{
    // Sanity checks.
    if (rLhs.size1() == 0 || rLhs.size2() == 0) {
        return true;
    }

    KRATOS_TRY
    // Check for convergence only if valid tolerances were provided.
    const bool check_convergence =    0 < mpImpl->mSolveSettings.residualAbsoluteTolerance
                                   || 0 < mpImpl->mSolveSettings.residualRelativeTolerance;
    [[maybe_unused]] typename TSparse::DataType initial_residual = 1.0;
    if (check_convergence) {
        VectorType residual(rRhs.size());
        TSparse::SetToZero(rRhs);
        TSparse::Mult(rLhs, rSolution, residual);
        residual -= rRhs;
        initial_residual = TSparse::TwoNorm(residual);
    } // if check_convergence

    if (initial_residual) {
        // Construct a view of the left hand side matrix.
        typename Impl::CSRAdaptor lhs_adaptor;
        lhs_adaptor.rowCount        = rLhs.size1();
        lhs_adaptor.columnCount     = rLhs.size2();
        lhs_adaptor.entryCount      = rLhs.nnz();
        lhs_adaptor.pRowExtents     = rLhs.index1_data().begin();
        lhs_adaptor.pColumnIndices  = rLhs.index2_data().begin();
        lhs_adaptor.pEntries        = rLhs.value_data().begin();

        // Perform smoothing
        if (mcgs::solve(&*rSolution.begin(),
                        lhs_adaptor,
                        &*rRhs.begin(),
                        mpImpl->mpRowPartition,
                        mpImpl->mSolveSettings) != MCGS_SUCCESS) {
            KRATOS_ERROR << "gauss-seidel smoothing failed";
        } // if mcgs::solve failed

        // Check for convergence only if (theoretically) attainable tolerances were provided.
        if (check_convergence) {
            VectorType residual(rRhs.size());
            TSparse::SetToZero(rRhs);
            TSparse::Mult(rLhs, rSolution, residual);
            residual -= rRhs;
            const typename TSparse::DataType final_residual = TSparse::TwoNorm(residual);

            // "initial_residual" is checked, so division by it is safe here.
            if (final_residual < mpImpl->mSolveSettings.residualAbsoluteTolerance
                || final_residual / initial_residual < mpImpl->mSolveSettings.residualRelativeTolerance) {
                return true;
            } else {
                return false;
            }
        } /*if check_convergence*/ else {
            // Tolerances are set to unattainable values, so we did definitely not converge.
            return false;
        } // if !check_convergence
    } /*if initial_residual*/ else {
        return true;
    } // if !initial_residual
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void GaussSeidelSmoother<TSparse,TDense>::FinalizeSolutionStep(SparseMatrixType& rLhs,
                                                               VectorType& rSolution,
                                                               VectorType& rRhs)
{
    KRATOS_TRY
    if (mpImpl->mSolveSettings.parallelization != mcgs::Parallelization::None) {
        if (mcgs::revertReorder(rLhs.size1(),
                                rLhs.size2(),
                                rLhs.nnz(),
                                &*rLhs.index1_data().begin(),
                                &*rLhs.index2_data().begin(),
                                &*rLhs.value_data().begin(),
                                &*rSolution.begin(),
                                &*rRhs.begin(),
                                mpImpl->mpRowPartition) != MCGS_SUCCESS) {
            KRATOS_ERROR << "multicolor gauss-seidel failed to revert its reordering of the linear system";
        }
    }
    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
Parameters GaussSeidelSmoother<TSparse,TDense>::GetDefaultParameters() const
{
    return Parameters(R"({
"solver_type" : "gauss_seidel",
"relaxation_factor"  :  1.0,
"absolute_tolerance" : -1,
"relative_tolerance" : -1,
"max_iterations"     :  1,
"reform_dofs_at_each_step" : false,
"parallelization"    : "none",
"verbosity"          : 1,
"coloring_settings" : {
    "shrinking_factor"  : -1,
    "max_stall_count"   :  1000,
    "tolerance"         :  0,
    "verbosity"         :  1
}
})");
}


template <class TSparse, class TDense>
typename LinearSolver<TSparse,TDense>::Pointer
GaussSeidelSmootherFactory<TSparse,TDense>::CreateSolver(Parameters Settings) const
{
    return typename LinearSolver<TSparse,TDense>::Pointer(new GaussSeidelSmoother<TSparse,TDense>(Settings));
}


template class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmoother<TUblasSparseSpace<double>,
                                                                         TUblasDenseSpace<double>>;
template class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmoother<TUblasSparseSpace<float>,
                                                                         TUblasDenseSpace<double>>;

template class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmootherFactory<TUblasSparseSpace<double>,
                                                                                TUblasDenseSpace<double>>;
template class KRATOS_API(LINEARSOLVERS_APPLICATION) GaussSeidelSmootherFactory<TUblasSparseSpace<float>,
                                                                                TUblasDenseSpace<double>>;

} // namespace Kratos
