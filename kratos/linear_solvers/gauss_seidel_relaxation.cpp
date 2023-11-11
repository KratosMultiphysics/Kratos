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


namespace Kratos {


template <class TSparseSpace,
          class TDenseSpace,
          class TReorderer>
struct GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::Impl
{
    double mTolerance;

    int mVerbosity;

    std::size_t mMaxIterations;
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

    mpImpl->mTolerance = parameters["tolerance"].Get<double>();
    mpImpl->mVerbosity = parameters["verbosity"].Get<int>();
    mpImpl->mMaxIterations = parameters["max_iterations"].Get<int>();

    KRATOS_CATCH("")
}


// Necessary for PIMPL
template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::~GaussSeidelRelaxation()
{
}


template<class TSparseSpace,
         class TDenseSpace,
         class TReorderer>
bool GaussSeidelRelaxation<TSparseSpace,TDenseSpace,TReorderer>::Solve(SparseMatrix& rA,
                                                                       Vector& rX,
                                                                       Vector& rB)
{
    KRATOS_TRY

    Vector residual(rX.size());

    for (std::size_t i_relax=0ul; i_relax<mpImpl->mMaxIterations; ++i_relax) {
        const auto it_row_end = rA.end1();
        for (auto it_row=rA.begin1(); it_row!=it_row_end; ++it_row) {
            const std::size_t i_row = it_row.index1();
            double value = rB[i_row];
            double diagonal = 1.0;

            const auto it_column_end = it_row.end();
            for (auto it_column=it_row.begin(); it_column!=it_column_end; ++it_column) {
                const auto i_column = it_column.index2();
                if (i_column == i_row) {
                    diagonal = *it_column;
                } else {
                    value -= *it_column * rX[i_column];
                }
            }

            rX[i_row] = value / diagonal;
        }

        //noalias(residual) = rB - prod(rA, rX);
        //if (norm_2(residual) < mpImpl->mTolerance) {
        //    break;
        //}
    }

    noalias(residual) = rB - prod(rA, rX);
    return norm_2(residual) < mpImpl->mTolerance;
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
    "tolerance" : 1e-6
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
