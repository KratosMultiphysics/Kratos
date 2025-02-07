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


/// @brief Computes \f$A += c * B\f$.
/// @tparam TValue Value type of the operand matrices.
/// @param rLeft Left hand side matrix to add the scaled right hand matrix to.
/// @param rRight Right hand matrix to scale and add to the left hand side.
/// @param Coefficient Coefficient to scale the right hand side matrix by.
template <class TValue>
static void
InPlaceMatrixAdd(typename TUblasSparseSpace<TValue>::MatrixType& rLeft,
                 const typename TUblasSparseSpace<TValue>::MatrixType& rRight,
                 const TValue Coefficient = 1.0)
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rLeft.size1() == rRight.size1() && rLeft.size2() == rRight.size2())
        << "cannot add matrices off different shapes "
        << "(" << rLeft.size1() << "x" << rLeft.size2() << ") + "
        << "(" << rRight.size1() << "x" << rRight.size2() << ")";

    // Early exit on empty update matrix.
    if (not rRight.size1() or not rRight.size2() or not rRight.nnz()) return;

    KRATOS_ERROR_IF_NOT(rRight.nnz() <= rLeft.nnz());

    #define KRATOS_MATRIX_SUM(sum_operator)                                                                                         \
        IndexPartition<std::size_t>(rLeft.size1()).for_each([&rLeft, &rRight, Coefficient](const std::size_t i_row) {               \
            const auto i_right_entry_begin = rRight.index1_data()[i_row];                                                           \
            const auto i_right_entry_end   = rRight.index1_data()[i_row + 1];                                                       \
                                                                                                                                    \
            const auto i_left_entry_begin = rLeft.index1_data()[i_row];                                                             \
            const auto i_left_entry_end   = rLeft.index1_data()[i_row + 1];                                                         \
            const auto it_left_column_begin = rLeft.index2_data().begin() + i_left_entry_begin;                                     \
            const auto it_left_column_end   = rLeft.index2_data().begin() + i_left_entry_end;                                       \
                                                                                                                                    \
            auto it_left_column = it_left_column_begin;                                                                             \
            for (auto i_right_entry=i_right_entry_begin; i_right_entry!=i_right_entry_end; ++i_right_entry) {                       \
                const auto i_column = rRight.index2_data()[i_right_entry];                                                          \
                                                                                                                                    \
                /* Find the entry in the left matrix corresponding to the entry in the right one. */                                \
                it_left_column = std::lower_bound(it_left_column, it_left_column_end, i_column);                                    \
                KRATOS_ERROR_IF(it_left_column == it_left_column_end)                                                               \
                    << "left hand matrix has no entry in row " << i_row << " at column " << i_column;                               \
                                                                                                                                    \
                const auto i_left_entry = std::distance(rLeft.index2_data().begin(), it_left_column);                               \
                rLeft.value_data()[i_left_entry] sum_operator rRight.value_data()[i_right_entry];                                   \
            } /* for i_right in range(i_right_begin, i_right_end) */                                                                \
        }) /*for i_row in range(rLeft.size1())*/

    if (Coefficient == 1.0) {
        KRATOS_TRY
        KRATOS_MATRIX_SUM(+=);
        KRATOS_CATCH("")
    } /*if Coefficient == 1.0*/ else if (Coefficient == -1.0) {
        KRATOS_TRY
        KRATOS_MATRIX_SUM(-=);
        KRATOS_CATCH("")
    } /*if Coefficient == -1.0*/ else {
        KRATOS_TRY
        KRATOS_MATRIX_SUM(+= Coefficient *);
        KRATOS_CATCH("")
    } // if Coefficient != 1.0
    #undef KRATOS_MATRIX_SUM
}


/// @brief Compute the union of two sparse matrix' sparsity patterns.
template <class TValue>
static void
MergeMatrices(typename TUblasSparseSpace<TValue>::MatrixType& rLeft,
              const typename TUblasSparseSpace<TValue>::MatrixType& rRight)
{
    using SpaceType = TUblasSparseSpace<TValue>;
    using MatrixType = typename SpaceType::MatrixType;
    using IndexType = typename MatrixType::index_array_type::value_type;
    MatrixType output;

    KRATOS_TRY

    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rLeft.size1() == rRight.size1() && rLeft.size2() == rRight.size2())
        << "cannot merge incompatible matrices "
        << "(" << rLeft.size1() << "x" << rLeft.size2() << ")"
        << " and "
        << "(" << rRight.size1() << "x" << rRight.size2() << ")";

    if (not rRight.size1() or not rRight.size2() or not rRight.nnz()) return;

    // Declare new containers for the merged matrix.
    typename MatrixType::index_array_type row_extents(rLeft.index1_data().size());
    typename MatrixType::index_array_type column_indices;
    block_for_each(row_extents, [](auto& r_item){r_item = 0;});

    // Merge rows into separate containers.
    {
        std::vector<std::vector<IndexType>> rows(rLeft.size1());
        IndexPartition<IndexType>(rLeft.size1()).for_each([&rows, &rLeft, &rRight](const IndexType i_row){
            const IndexType i_left_begin = rLeft.index1_data()[i_row];
            const IndexType i_left_end = rLeft.index1_data()[i_row + 1];
            const IndexType i_right_begin = rRight.index1_data()[i_row];
            const IndexType i_right_end = rRight.index1_data()[i_row + 1];
            rows[i_row].reserve(i_left_end - i_left_begin + i_right_end - i_right_begin);

            rows[i_row].insert(rows[i_row].end(),
                               rLeft.index2_data().begin() + i_left_begin,
                               rLeft.index2_data().begin() + i_left_end);
            rows[i_row].insert(rows[i_row].end(),
                               rRight.index2_data().begin() + i_right_begin,
                               rRight.index2_data().begin() + i_right_end);
            std::sort(rows[i_row].begin(),
                        rows[i_row].end());
            rows[i_row].erase(std::unique(rows[i_row].begin(),
                                          rows[i_row].end()),
                                rows[i_row].end());
            rows[i_row].shrink_to_fit();
        });

        // Compute new row extents.
        for (IndexType i_row=0; i_row<rLeft.size1(); ++i_row) {
            row_extents[i_row + 1] = row_extents[i_row] + rows[i_row].size();
        } // for i_row in range(rLeft.size1)

        // Fill column indices and entries.
        column_indices.resize(row_extents[rLeft.size1()], false);
        IndexPartition<IndexType>(rLeft.size1()).for_each([&rows, &row_extents, &column_indices](const IndexType i_row){
            const IndexType i_begin = row_extents[i_row];
            std::copy(rows[i_row].begin(),
                        rows[i_row].end(),
                        column_indices.begin() + i_begin);
        });
    }

    // Construct the new matrix.
    rLeft = MatrixType(rLeft.size1(), rLeft.size2(), column_indices.size());
    rLeft.index1_data().swap(row_extents);
    rLeft.index2_data().swap(column_indices);
    block_for_each(rLeft.value_data(), [](auto& r_item){r_item = 0;});
    rLeft.set_filled(rLeft.size1() + 1, column_indices.size());

    KRATOS_CATCH("")
}


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
    block_for_each(rMatrix.value_data(), [](TValue& r_entry){r_entry=0;});
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
        AtomicAdd(rRhs[i_global], static_cast<typename TSparse::DataType>(rContribution[i_local]));
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
