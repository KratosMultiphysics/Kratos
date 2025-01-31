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

#pragma once

// Project includes
#include "solving_strategies/schemes/scheme.h" // Scheme
#include "spaces/ublas_space.h" // TUblasSparseSpace


namespace Kratos {


template <bool SortedRows,
          class TValue,
          class TRowMapContainer>
void MakeSparseTopology(const TRowMapContainer& rRows,
                        const std::size_t ColumnCount,
                        typename TUblasSparseSpace<TValue>::MatrixType& rMatrix)
{
    KRATOS_TRY

    const std::size_t row_count = rRows.size();
    std::size_t entry_count = 0ul;

    {
        auto row_extents = rMatrix.index1_data();

        // Resize row extents.
        row_extents.resize(row_count + 1, false);

        // Fill row extents and collect the total number of entries to store.
        row_extents[0] = 0;
        for (int i = 0; i < static_cast<int>(row_count); i++) {
            row_extents[i + 1] = row_extents[i] + rRows[i].size();
            entry_count += rRows[i].size();
        } // for i in range(row_count)

        // Resize the output matrix and all its containers.
        rMatrix = typename TUblasSparseSpace<TValue>::MatrixType(rRows.size(), ColumnCount, entry_count);
        rMatrix.index1_data().swap(row_extents);
    }

    auto& r_row_extents = rMatrix.index1_data();
    auto& r_column_indices = rMatrix.index2_data();

    // Copy column indices.
    IndexPartition<std::size_t>(row_count).for_each([&r_row_extents, &r_column_indices, &rRows](const std::size_t i_row){
        const unsigned i_entry_begin = r_row_extents[i_row];
        const unsigned i_entry_end = r_row_extents[i_row + 1];
        KRATOS_DEBUG_ERROR_IF(r_column_indices.size() < i_entry_end);
        std::copy(rRows[i_row].begin(),
                  rRows[i_row].end(),
                  r_column_indices.begin() + i_entry_begin);

        if constexpr (!SortedRows) {
            std::sort(r_column_indices.begin() + i_entry_begin,
                      r_column_indices.begin() + i_entry_end);
        }
    });

    KRATOS_TRY
    rMatrix.set_filled(row_count + 1, entry_count);
    KRATOS_CATCH("")

    KRATOS_CATCH("")
}


/// @brief Map LHS contributions from a local dense matrix to a sparse global matrix.
template <class TSparse, class TDense>
void MapRowContribution(typename TSparse::MatrixType& rLhs,
                        const typename TDense::MatrixType& rLocalLhs,
                        const unsigned iRow,
                        const unsigned iLocalRow,
                        const Element::EquationIdVectorType& rEquationIds) noexcept
{
    auto& r_entries = rLhs.value_data();
    const auto& r_row_extents = rLhs.index1_data();
    const auto& r_column_indices = rLhs.index2_data();

    const auto i_row_begin = r_row_extents[iRow];
    const auto i_row_end = r_row_extents[iRow + 1];

    for (unsigned int i_local_column=0; i_local_column<rEquationIds.size(); i_local_column++) {
        const unsigned int iColumn = rEquationIds[i_local_column];
        const auto it = std::lower_bound(r_column_indices.begin() + i_row_begin,
                                            r_column_indices.begin() + i_row_end,
                                            iColumn);
        r_entries[std::distance(r_column_indices.begin(), it)] += rLocalLhs(iLocalRow, i_local_column);
    }
}


/// @brief Map RHS contributions from local to global space.
template <class TSparse, class TDense>
void MapContribution(typename TSparse::VectorType& rRhs,
                     const typename TDense::VectorType& rContribution,
                     const Element::EquationIdVectorType& rEquationIds) noexcept
{
    const unsigned local_size = rContribution.size();

    for (unsigned i_local = 0; i_local < local_size; i_local++) {
        const unsigned i_global = rEquationIds[i_local];
        AtomicAdd(rRhs[i_global], rContribution[i_local]);
    }
}


/// @brief Map LHS contributions from local to global space.
template <class TSparse, class TDense>
void MapContribution(typename TSparse::MatrixType& rLhs,
                     const typename TDense::MatrixType& rContribution,
                     const Element::EquationIdVectorType& rEquationIds,
                     LockObject* pLockBegin) noexcept
{
    const std::size_t local_size = rContribution.size1();
    for (IndexType i_local = 0; i_local < local_size; i_local++) {
        const IndexType i_global = rEquationIds[i_local];
        std::scoped_lock<LockObject> lock(pLockBegin[i_global]);
        MapRowContribution<TSparse,TDense>(rLhs,
                                           rContribution,
                                           i_global,
                                           i_local,
                                           rEquationIds);
    }
}


/// @brief Map LHS and RHS contributions from local to global space.
template <class TSparse, class TDense>
void MapContribution(typename TSparse::MatrixType& rLhs,
                     typename TSparse::VectorType& rRhs,
                     const typename TDense::MatrixType& rLhsContribution,
                     const typename TDense::VectorType& rRhsContribution,
                     const Element::EquationIdVectorType& rEquationIds,
                     LockObject* pLockBegin) noexcept
{
    const unsigned local_size = rLhsContribution.size1();
    for (unsigned i_local = 0; i_local < local_size; i_local++) {
        const unsigned i_global = rEquationIds[i_local];
        std::scoped_lock<LockObject> lock(pLockBegin[i_global]);
        rRhs[i_global] += rRhsContribution[i_local];
        MapRowContribution<TSparse,TDense>(rLhs,
                                           rLhsContribution,
                                           i_global,
                                           i_local,
                                           rEquationIds);
    }
}


/// @brief Compute the local contributions of an @ref Element or @ref Condition and assemble them into the global system.
/// @details Whether from elements or conditions, whether assembling into the LHS, RHS, or both,
///          all cases are handled in this one function to reduce code duplication and have logically
///          coherent parts of the code in one place.
template <class TSparse,
          class TDense,
          bool AssembleLHS,
          bool AssembleRHS,
          class TEntity>
void MapEntityContribution(TEntity& rEntity,
                           Scheme<TSparse,TDense>& rScheme,
                           const ProcessInfo& rProcessInfo,
                           Element::EquationIdVectorType& rEquationIndices,
                           typename TDense::MatrixType* pLhsContribution,
                           typename TDense::VectorType* pRhsContribution,
                           typename TSparse::MatrixType* pLhs,
                           typename TSparse::VectorType* pRhs,
                           LockObject* pLockBegin)
{
    if constexpr (AssembleLHS || AssembleRHS) {
        if (rEntity.IsActive()) {
            if constexpr (AssembleLHS) {
                if constexpr (AssembleRHS) {
                    // Assemble LHS and RHS.
                    rScheme.CalculateSystemContributions(rEntity,
                                                         *pLhsContribution,
                                                         *pRhsContribution,
                                                         rEquationIndices,
                                                         rProcessInfo);
                    MapContribution<TSparse,TDense>(*pLhs,
                                                    *pRhs,
                                                    *pLhsContribution,
                                                    *pRhsContribution,
                                                    rEquationIndices,
                                                    pLockBegin);
                } /*if AssembleRHS*/ else {
                    // Assemble LHS only.
                    rScheme.CalculateLHSContribution(rEntity,
                                                    *pLhsContribution,
                                                    rEquationIndices,
                                                    rProcessInfo);
                    MapContribution<TSparse,TDense>(*pLhs,
                                                    *pLhsContribution,
                                                    rEquationIndices,
                                                    pLockBegin);
                } // if !AssembleRHS
            } /*if AssembleLHS*/ else {
                // Assemble RHS only.
                rScheme.CalculateRHSContribution(rEntity,
                                                 *pRhsContribution,
                                                 rEquationIndices,
                                                 rProcessInfo);
                MapContribution<TSparse,TDense>(*pRhs,
                                                *pRhsContribution,
                                                rEquationIndices);
            } // if !AssembleLHS
        } // if rEntity.IsActive
    } // if AssembleLHS or AssembleRHS
}


} // namespace Kratos
