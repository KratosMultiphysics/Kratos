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
#include "mkl.h" // MKL_INT

// Project includes
#include "custom_solvers/mkl_solver_base.hpp" // MKLSolverBase
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace

// STL includes
#include <optional> // std::optional
#include <tuple> // std::tuple


namespace Kratos {



/** @details Sparse Intel solvers use @p MKL_INT as their index type, which usually translates to @p int,
 *           while sparse matrices used by Kratos use @p std::size_t. Even disregarding the sign bit,
 *           these two integer types will not be binary-compatible on modern hardware because
 *           @p int takes up 4 bytes while @p std::size_t 8, meaning their alignment does not match.
 *
 *           As a result, the LHS matrix must be copied with @p MKL_INT as its index type.
 */
template <class TSparse, class TDense>
struct MKLSolverBase<TSparse,TDense>::Impl
{
    static_assert(std::is_same_v<MKL_INT,int>);

    std::optional<std::vector<MKL_INT>> mMaybeRowExtents;

    std::optional<std::vector<MKL_INT>> mMaybeColumnIndices;
}; // struct MKLSolverBase::Impl


template <class TSparse, class TDense>
MKLSolverBase<TSparse,TDense>::MKLSolverBase()
    : mpImpl(new Impl)
{
}


template <class TSparse, class TDense>
MKLSolverBase<TSparse,TDense>::~MKLSolverBase() = default;


template <class TSparse, class TDense>
bool MKLSolverBase<TSparse,TDense>::Solve(SparseMatrix& rLhs,
                                          Vector& rSolution,
                                          Vector& rRhs)
{
    if (!rLhs.nnz()) return true;

    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rSolution.size());
    KRATOS_ERROR_IF_NOT(rRhs.size());

    KRATOS_TRY
    auto [lhs_view, solution_view, rhs_view] = this->MakeSystemView(
        rLhs,
        &*rSolution.begin(),
        (&*rSolution.begin()) + rSolution.size(),
        &*rRhs.begin(),
        (&*rRhs.begin()) + rRhs.size());
    return this->Solve(lhs_view, solution_view, rhs_view);
    KRATOS_CATCH("")
}


/// Copy and transform the index arrays to MKL_INT, and a 1-based indices. ...
template <class TSparse, class TDense>
void MKLSolverBase<TSparse,TDense>::ProvideAdditionalData(SparseMatrix& rLhs,
                                                          Vector&,
                                                          Vector&,
                                                          ModelPart::DofsArrayType&,
                                                          ModelPart&)
{
    KRATOS_TRY
    mpImpl->mMaybeRowExtents.emplace(rLhs.index1_data().size());
    mpImpl->mMaybeColumnIndices.emplace(rLhs.index2_data().size());

    for (auto [p_source_array, p_target_array] : {std::make_pair(&rLhs.index1_data(), &mpImpl->mMaybeRowExtents.value()),
                                                  std::make_pair(&rLhs.index2_data(), &mpImpl->mMaybeColumnIndices.value())}) {
        auto& r_source = *p_source_array;
        auto& r_target = *p_target_array;
        IndexPartition<typename TSparse::IndexType>(p_source_array->size()).for_each(
            [&r_source, &r_target](typename TSparse::IndexType i) mutable {
                r_target[i] = static_cast<MKL_INT>(r_source[i])
                            + static_cast<MKL_INT>(1); //< intel sometimes exclusively uses 1-based indexing
            });
    } // for p_source_array, p_target_array in array_pairs

    KRATOS_CATCH("")
}


template <class TSparse, class TDense>
void MKLSolverBase<TSparse,TDense>::FinalizeSolutionStep(SparseMatrix&,
                                                         Vector&,
                                                         Vector&)
{
}


template <class TSparse, class TDense>
void MKLSolverBase<TSparse,TDense>::Clear()
{
    Base::Clear();
    mpImpl->mMaybeRowExtents.reset();
    mpImpl->mMaybeColumnIndices.reset();
}


template <class TSparse, class TDense>
std::tuple<
    typename MKLSolverBase<TSparse,TDense>::CSRView,
    typename MKLSolverBase<TSparse,TDense>::template VectorView</*IsMutable=*/true>,
    typename MKLSolverBase<TSparse,TDense>::template VectorView</*IsMutable=*/false>>
MKLSolverBase<TSparse,TDense>::MakeSystemView(const SparseMatrix& rLhs,
                                              typename TSparse::DataType* itSolutionBegin,
                                              typename TSparse::DataType* itSolutionEnd,
                                              const typename TSparse::DataType* itRhsBegin,
                                              const typename TSparse::DataType* itRhsEnd) const
{
    KRATOS_TRY

    // Sanity checks.
    KRATOS_ERROR_IF_NOT(static_cast<typename TSparse::IndexType>(std::distance(itSolutionBegin, itSolutionEnd)) == rLhs.size1());
    KRATOS_ERROR_IF_NOT(static_cast<typename TSparse::IndexType>(std::distance(itRhsBegin, itRhsEnd)) == rLhs.size1());

    // Construct views.
    return std::make_tuple(
        typename MKLSolverBase<TSparse,TDense>::CSRView {
            /*row_count=*/      static_cast<MKL_INT>(rLhs.size1()),
            /*column_count=*/   static_cast<MKL_INT>(rLhs.size2()),
            /*entry_count=*/    static_cast<MKL_INT>(rLhs.nnz()),
            /*it_row_begin=*/   mpImpl->mMaybeRowExtents.value().data(),
            /*it_column_begin=*/mpImpl->mMaybeColumnIndices.value().data(),
            /*it_entry_begin=*/ &*rLhs.value_data().begin()},
        typename MKLSolverBase<TSparse,TDense>::template VectorView</*IsMutable=*/true> {
            /*size=*/           static_cast<MKL_INT>(std::distance(itSolutionBegin, itSolutionEnd)),
            /*it_begin=*/       itSolutionBegin},
        typename MKLSolverBase<TSparse,TDense>::template VectorView</*IsMutable=*/false>{
            /*size=*/           static_cast<MKL_INT>(std::distance(itRhsBegin, itRhsEnd)),
            /*it_begin=*/       itRhsBegin
        }
    );

    KRATOS_CATCH("")
}


template class MKLSolverBase<TUblasSparseSpace<double>,TUblasDenseSpace<double>>;


template class MKLSolverBase<TUblasSparseSpace<float>,TUblasDenseSpace<double>>;


} // namespace Kratos
