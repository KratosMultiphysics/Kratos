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
#include "gauss_seidel_relaxation.h"
#include "spaces/ublas_space.h"
#include "utilities/profiler.h"

// System includes
#include <optional>


namespace Kratos {


template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mRelaxation;

    double mTolerance;

    int mVerbosity;

    std::size_t mMaxIterations;

    bool mBackward;
}; // struct GaussSeidelRelaxation::Impl


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::GaussSeidelRelaxation(Parameters parameters)
    : mpImpl(new Impl)
{
    KRATOS_TRY

    Parameters default_parameters = this->GetDefaultParameters();
    parameters.ValidateAndAssignDefaults(default_parameters);

    KRATOS_ERROR_IF_NOT(parameters["solver_type"].GetString() == "gauss_seidel")
        << "Requested a(n) '" << parameters["solver_type"].GetString() << "' solver,"
        << " but constructing a GaussSeidelRelaxation";

    mpImpl->mRelaxation = parameters["relaxation"].Get<double>();
    mpImpl->mTolerance = parameters["tolerance"].Get<double>();
    mpImpl->mVerbosity = parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = parameters["max_iterations"].Get<int>();
    mpImpl->mBackward = parameters["backward"].Get<bool>();

    KRATOS_CATCH("")
}


// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::~GaussSeidelRelaxation()
{
}


namespace {
template <class TIterator>
void GaussSeidelSweep(TIterator itRow,
                      const TIterator itRowEnd,
                      Vector& rX,
                      const Vector& rB,
                      const double relaxation)
{
    for (; itRow!=itRowEnd; ++itRow) {
        const std::size_t i_row = itRow.index1();
        double value = rB[i_row];
        double diagonal = 1.0;

        const auto it_column_end = itRow.end();
        for (auto it_column=itRow.begin(); it_column!=it_column_end; ++it_column) {
            const auto i_column = it_column.index2();
            if (i_column == i_row) {
                diagonal = *it_column;
            } else {
                value -= *it_column * rX[i_column];
            }
        }

        rX[i_row] += relaxation * (value / diagonal - rX[i_row]);
    }
}
} // unnamed namespace


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                       Vector& rX,
                                                                       Vector& rB)
{
    KRATOS_TRY

    const double relaxation = mpImpl->mRelaxation;
    Vector residual(rX.size());

    // Compute residuals after each iteration only if a
    // non-negative tolerance is requested. A negative
    // tolerance makes no sense anyway, so let the solver
    // chug on until the max iteration count without ever
    // computing residuals.
    const bool residual_criterion = 0 <= this->mpImpl->mTolerance;

    for (std::size_t i_relax=0ul; i_relax<mpImpl->mMaxIterations; ++i_relax) {
        if (this->mpImpl->mBackward) {
            GaussSeidelSweep(rA.rbegin1(), rA.rend1(), rX, rB, relaxation);
        } else {
            GaussSeidelSweep(rA.begin1(), rA.end1(), rX, rB, relaxation);
        }

        if (residual_criterion || 3 <= mpImpl->mVerbosity) {
            TSparseSpace::Mult(rA, rX, residual);
            TSparseSpace::ScaleAndAdd(1.0, rB, -1.0, residual);
            const double residual_norm = TSparseSpace::TwoNorm(residual);

            if (3 <= mpImpl->mVerbosity) {
                KRATOS_INFO("GaussSeidelRelaxation")
                    << "iteration " << i_relax
                    << " residual " << residual_norm << "\n";
            }

            if (residual_norm <= mpImpl->mTolerance) {
                break;
            }
        }
    }

    TSparseSpace::Mult(rA, rX, residual);
    TSparseSpace::ScaleAndAdd(1.0, rB, -1.0, residual);
    const double residual_norm = TSparseSpace::TwoNorm(residual);
    return residual_norm <= mpImpl->mTolerance;
    KRATOS_CATCH("")
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
Parameters
GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::GetDefaultParameters()
{
    return Parameters(R"(
{
    "solver_type" : "gauss_seidel",
    "verbosity" : 0,
    "max_iterations" : 500,
    "tolerance" : -1,
    "relaxation" : 1.0,
    "backward" : false
}
    )");
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                       DenseMatrix& rX,
                                                                       DenseMatrix& rB)
{
    return false;
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::PrintInfo(std::ostream& rStream) const
{
    rStream << "GaussSeidelRelaxation";
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
void GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::PrintData(std::ostream& rStream) const
{
    rStream
        << (mpImpl->mBackward ? "backward" : "forward") << " gauss-seidel\n"
        << "tolerance     : " << mpImpl->mTolerance << "\n"
        << "max iterations: " << mpImpl->mMaxIterations << "\n"
        << "verbosity     : " << mpImpl->mVerbosity << "\n"
        ;
}


template
class GaussSeidelRelaxation<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>,
    Reorderer<
        TUblasSparseSpace<double>,
        TUblasDenseSpace<double>
    >
>;


} // namespace Kratos
