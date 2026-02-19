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
#include "utilities/profiler.h" // KRATOS_PROFILE_SCOPE

// System includes
#include <cstdint> // std::uint8_t
#include <climits> // CHAR_BIT



// Define the hash table type used for computing CSR patterns.
// The fallback option is std::unordered_*, but a more performant
// alternative is boost::unordered_flat_* that requires boost 1.81
// or newer.
#if BOOST_VERSION < 108100
    #include <unordered_set>
    #include <unordered_map>

    namespace Kratos {
        template <class TKey, class TValue, class ...TArgs>
        using CSRHashMap = std::unordered_map<TKey,TValue,TArgs...>;

        template <class TValue, class ...TArgs>
        using CSRHashSet = std::unordered_set<TValue,TArgs...>;
    } // namespace Kratos
#else
    #include <boost/unordered/unordered_flat_map.hpp>
    #include <boost/unordered/unordered_flat_set.hpp>

    namespace Kratos {
        template <class TKey, class TValue, class ...TArgs>
        using CSRHashMap = boost::unordered_flat_map<TKey,TValue,TArgs...>;

        template <class TValue, class ...TArgs>
        using CSRHashSet = boost::unordered_flat_set<TValue,TArgs...>;
    } // namespace Kratos
#endif


namespace Kratos {


// "Hey compiler trust me, my arrays are aligned so try vectorizing my loops over them."
#if defined(__GNUC__) || defined(__clang__)
#define KRATOS_GET_ALIGNED_INDEX_ARRAY(POINTER_TYPE, POINTER_NAME, POINTER_VALUE)                                                   \
    POINTER_TYPE __restrict POINTER_NAME = (POINTER_TYPE) __builtin_assume_aligned(POINTER_VALUE, CHAR_BIT * sizeof(POINTER_TYPE));
#else
// Someone can figure out how to encourage AVX on MSVC, I won't.
#define KRATOS_GET_ALIGNED_INDEX_ARRAY(POINTER_TYPE, POINTER_NAME, POINTER_VALUE)   \
    POINTER_TYPE POINTER_NAME = POINTER_VALUE;
#endif


/** @brief Computes \f$ A += c * B \f$.
 *  @tparam TValue Value type of the operand matrices.
 *  @param rLeft Left hand side matrix to add the scaled right hand matrix to.
 *  @param rRight Right hand matrix to scale and add to the left hand side.
 *  @param Coefficient Coefficient to scale the right hand side matrix by.
 *  @warning The sparsity pattern of @p rLeft must be a superset of @p rRight.
 */
template <class TValue>
static void
InPlaceMatrixAdd(typename TUblasSparseSpace<TValue>::MatrixType& rLeft,
                 const typename TUblasSparseSpace<TValue>::MatrixType& rRight,
                 const TValue Coefficient = 1.0)
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rLeft.size1() == rRight.size1() && rLeft.size2() == rRight.size2())
        << "cannot add matrices of different shapes "
        << "(" << rLeft.size1() << "x" << rLeft.size2() << ") + "
        << "(" << rRight.size1() << "x" << rRight.size2() << ")";

    // Early exit on empty update matrix.
    if (!rRight.size1() || !rRight.size2() || !rRight.nnz()) return;

    KRATOS_ERROR_IF_NOT(rRight.nnz() <= rLeft.nnz());

    #define KRATOS_MATRIX_SUM(sum_operator)                                                                                         \
        IndexPartition<IndexType>(rLeft.size1()).for_each([&rLeft, &rRight, Coefficient](const IndexType i_row) {                   \
            (void)Coefficient; /*<== suppress unused capture warnings.*/                                                            \
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
                KRATOS_ERROR_IF(it_left_column == it_left_column_end || *it_left_column != i_column)                                \
                    << "left hand matrix has no entry in row " << i_row << " at column " << i_column;                               \
                                                                                                                                    \
                const auto i_left_entry = std::distance(rLeft.index2_data().begin(), it_left_column);                               \
                rLeft.value_data()[i_left_entry] sum_operator rRight.value_data()[i_right_entry];                                   \
            } /* for i_right_entry in range(i_right_begin, i_right_end) */                                                          \
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

    if (!rRight.size1() || !rRight.size2() || !rRight.nnz()) return;

    // Declare new containers for the merged matrix.
    typename MatrixType::index_array_type row_extents(rLeft.index1_data().size());
    typename MatrixType::index_array_type column_indices;
    typename MatrixType::value_array_type values;
    block_for_each(
        row_extents,
        [](typename MatrixType::index_array_type::value_type& r_item){
            r_item = static_cast<TValue>(0);
        });

    // Merge rows into separate containers.
    {
        std::vector<
            std::vector<
                std::pair<typename MatrixType::index_array_type::value_type,
                          typename MatrixType::value_array_type::value_type
        >>> rows(rLeft.size1());

        IndexPartition<IndexType>(rLeft.size1()).for_each([&rows, &rLeft, &rRight](const IndexType i_row){
            const IndexType i_left_begin = rLeft.index1_data()[i_row];
            const IndexType i_left_end = rLeft.index1_data()[i_row + 1];
            const IndexType i_right_begin = rRight.index1_data()[i_row];
            const IndexType i_right_end = rRight.index1_data()[i_row + 1];
            rows[i_row].reserve((i_left_end - i_left_begin) + (i_right_end - i_right_begin));

            for (IndexType i_entry=i_left_begin; i_entry<i_left_end; ++i_entry) {
                rows[i_row].emplace_back(rLeft.index2_data()[i_entry], rLeft.value_data()[i_entry]);
            }

            for (IndexType i_entry=i_right_begin; i_entry<i_right_end; ++i_entry) {
                rows[i_row].emplace_back(rRight.index2_data()[i_entry], rRight.value_data()[i_entry]);
            }

            std::stable_sort(rows[i_row].begin(),
                             rows[i_row].end(),
                             [](const auto& r_left, const auto& r_right) {return r_left.first < r_right.first;});

            rows[i_row].erase(std::unique(rows[i_row].begin(),
                                          rows[i_row].end(),
                                          [](const auto& r_left, const auto& r_right) {return r_left.first == r_right.first;}),
                              rows[i_row].end());
            rows[i_row].shrink_to_fit();
        }); // for i_row in range(rLeft.size1())

        // Compute new row extents.
        for (IndexType i_row=0; i_row<rLeft.size1(); ++i_row) {
            row_extents[i_row + 1] = row_extents[i_row] + rows[i_row].size();
        } // for i_row in range(rLeft.size1)

        // Fill column indices and entries.
        column_indices.resize(row_extents[rLeft.size1()], false);
        values.resize(row_extents[rLeft.size1()], false);
        IndexPartition<IndexType>(rLeft.size1()).for_each([&rows, &row_extents, &column_indices, &values](const IndexType i_row){
            const IndexType i_begin = row_extents[i_row];
            for (IndexType i_pair=0ul; i_pair<static_cast<IndexType>(rows[i_row].size()); ++i_pair) {
                const auto i_entry = i_begin + i_pair;
                column_indices[i_entry] = rows[i_row][i_pair].first;
                values[i_entry] = rows[i_row][i_pair].second;
            }
        }); // for i_row in range(rLeft.size1())
    }

    // Construct the new matrix.
    rLeft = MatrixType(rLeft.size1(), rLeft.size2(), column_indices.size());
    rLeft.index1_data().swap(row_extents);
    rLeft.index2_data().swap(column_indices);
    rLeft.value_data().swap(values);
    rLeft.set_filled(rLeft.size1() + 1, column_indices.size());

    KRATOS_CATCH("")
}


template <bool SortedRows,
          class TValue,
          class TRowMapContainer>
void MakeSparseTopology(TRowMapContainer& rRows,
                        const IndexType ColumnCount,
                        typename TUblasSparseSpace<TValue>::MatrixType& rMatrix,
                        bool EnsureDiagonal)
{
    KRATOS_TRY

    if (EnsureDiagonal) {
        IndexPartition<IndexType>(rRows.size()).for_each([&rRows](IndexType i_row){
            rRows[i_row].emplace(i_row);
        });
    }

    const IndexType row_count = rRows.size();
    IndexType entry_count = 0ul;

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
    IndexPartition<IndexType>(row_count).for_each([&r_row_extents, &r_column_indices, &rRows](const IndexType i_row){
        const IndexType i_entry_begin = r_row_extents[i_row];
        const IndexType i_entry_end = r_row_extents[i_row + 1];
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
                        const IndexType iRow,
                        const IndexType iLocalRow,
                        const Element::EquationIdVectorType& rEquationIds) noexcept
{
    auto& r_entries = rLhs.value_data();
    const auto& r_row_extents = rLhs.index1_data();
    const auto& r_column_indices = rLhs.index2_data();

    const auto i_row_begin = r_row_extents[iRow];
    const auto i_row_end = r_row_extents[iRow + 1];
    const typename TDense::IndexType row_size = rEquationIds.size();

    for (typename TDense::IndexType i_local_column=0; i_local_column<row_size; ++i_local_column) {
        const typename TSparse::IndexType i_column = rEquationIds[i_local_column];
        const auto it = std::lower_bound(r_column_indices.begin() + i_row_begin,
                                         r_column_indices.begin() + i_row_end,
                                         i_column);
        r_entries[std::distance(r_column_indices.begin(), it)] += rLocalLhs(iLocalRow, i_local_column);
    }
}


/// @brief Map RHS contributions from local to global space.
template <class TSparse, class TDense>
void MapContribution(typename TSparse::VectorType& rRhs,
                     const typename TDense::VectorType& rContribution,
                     const Element::EquationIdVectorType& rEquationIds) noexcept
{
    const typename TDense::IndexType local_size = rContribution.size();

    for (typename TDense::IndexType i_local = 0; i_local < local_size; i_local++) {
        const typename TSparse::IndexType i_global = rEquationIds[i_local];
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
    const typename TDense::IndexType local_size = rContribution.size1();
    for (typename TDense::IndexType i_local = 0; i_local < local_size; i_local++) {
        const IndexType i_global = rEquationIds[i_local];
        std::scoped_lock<LockObject> lock(pLockBegin[i_global]);
        MapRowContribution<TSparse,TDense>(
            rLhs,
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
    const typename TDense::IndexType local_size = rLhsContribution.size1();
    for (typename TDense::IndexType i_local = 0; i_local < local_size; i_local++) {
        const typename TSparse::IndexType i_global = rEquationIds[i_local];
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
                    rScheme.CalculateSystemContributions(
                        rEntity,
                        *pLhsContribution,
                        *pRhsContribution,
                        rEquationIndices,
                        rProcessInfo);
                    MapContribution<TSparse,TDense>(
                        *pLhs,
                        *pRhs,
                        *pLhsContribution,
                        *pRhsContribution,
                        rEquationIndices,
                        pLockBegin);
                } /*if AssembleRHS*/ else {
                    // Assemble LHS only.
                    rScheme.CalculateLHSContribution(
                        rEntity,
                        *pLhsContribution,
                        rEquationIndices,
                        rProcessInfo);
                    MapContribution<TSparse,TDense>(
                        *pLhs,
                        *pLhsContribution,
                        rEquationIndices,
                        pLockBegin);
                } // if !AssembleRHS
            } /*if AssembleLHS*/ else {
                // Assemble RHS only.
                rScheme.CalculateRHSContribution(
                    rEntity,
                    *pRhsContribution,
                    rEquationIndices,
                    rProcessInfo);
                MapContribution<TSparse,TDense>(
                    *pRhs,
                    *pRhsContribution,
                    rEquationIndices);
            } // if !AssembleLHS
        } // if rEntity.IsActive
    } // if AssembleLHS or AssembleRHS
}


template <class TSparse,
          class TDense,
          class TItRowDof,
          class TItColumnDof>
void ApplyDirichletConditions(typename TSparse::MatrixType& rLhs,
                              typename TSparse::VectorType& rRhs,
                              const TItRowDof itRowDofBegin,
                              const TItRowDof itRowDofEnd,
                              const TItColumnDof itColumnDofBegin,
                              [[maybe_unused]] const TItColumnDof itColumnDofEnd,
                              const typename TSparse::DataType DiagonalScaleFactor)
{
    // Type checks.
    static_assert(std::is_same_v<
        typename std::iterator_traits<TItRowDof>::value_type,
        Dof<typename TDense::DataType>>);

    KRATOS_TRY
    KRATOS_PROFILE_SCOPE(KRATOS_CODE_LOCATION);

    block_for_each(itRowDofBegin,
                   itRowDofEnd,
                   [&rLhs, &rRhs, itColumnDofBegin, DiagonalScaleFactor](const Dof<double>& r_dof){
        const IndexType i_dof = r_dof.EquationId();
        const typename TSparse::IndexType i_entry_begin = rLhs.index1_data()[i_dof];
        const typename TSparse::IndexType i_entry_end = rLhs.index1_data()[i_dof + 1];
        bool found_diagonal = false;

        if (r_dof.IsFixed()) {
            // Zero out the whole row, except the diagonal.
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = rLhs.index2_data()[i_entry];
                if (i_column == i_dof) {
                    found_diagonal = true;
                    rLhs.value_data()[i_entry] = DiagonalScaleFactor;
                } else {
                    rLhs.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)

            rRhs[i_dof] = static_cast<typename TSparse::DataType>(0);
        } /*if r_dof.IsFixed()*/ else {
            // Zero out the column which is associated with the zero'ed row.
            for (typename TSparse::IndexType i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = rLhs.index2_data()[i_entry];
                const auto it_column_dof = itColumnDofBegin + i_column;

                if (i_column == i_dof) {
                    found_diagonal = true;

                    // Check the entry on the main diagonal.
                    if (!rLhs.value_data()[i_entry]) {
                        // Explicit zero on the main diagonal.
                        // - If the associated entry on the RHS is also zero,
                        //   force the DoF to also vanish.
                        // - Otherwise it's impossible to solve the system,
                        //   so throw an exception.
                        KRATOS_ERROR_IF(rRhs[i_dof])
                            << "explicit zero on the main diagonal of row " << i_dof << " "
                            << "related to dof " << r_dof.GetVariable().Name() << " "
                            << "of node " << r_dof.Id()
                            << ", but the related entry on the right hand side does not vanish ("
                            << rLhs.value_data()[i_entry] << ")";
                        rLhs.value_data()[i_entry] = DiagonalScaleFactor;
                    }
                } /*if i_column == i_dof*/ else if (it_column_dof->IsFixed()) {
                    rLhs.value_data()[i_entry] = static_cast<typename TSparse::DataType>(0);
                }

            } // for i_entry in range(i_entry_begin, i_entry_end)
        } /*not r_dof.IsFixed()*/

        KRATOS_ERROR_IF_NOT(found_diagonal)
        << "implicit zero on the main diagonal of row " << i_dof << " "
        << "related to dof " << r_dof.GetVariable().Name() << " "
        << "of node " << r_dof.Id();
    });

    KRATOS_CATCH("")
}


template <class TSparse, class TDense, class TItDof>
void ApplyDirichletConditions(typename TSparse::MatrixType& rLhs,
                              typename TSparse::VectorType& rRhs,
                              const TItDof itDofBegin,
                              const TItDof itDofEnd,
                              const typename TSparse::DataType DiagonalScaleFactor)
{
    return ApplyDirichletConditions<TSparse,TDense>(
        rLhs,
        rRhs,
        itDofBegin,
        itDofEnd,
        itDofBegin,
        itDofEnd,
        DiagonalScaleFactor);
}


/// @internal
struct MatrixChecks
{
    static constexpr std::uint8_t None                    = 0;
    static constexpr std::uint8_t DiagonalExists          = 1;
    static constexpr std::uint8_t DiagonalIsNonNegative   = 1 << 1;
    static constexpr std::uint8_t DiagonalIsPositive      = 1 << 2;
    static constexpr std::uint8_t IsDiagonallyDominant    = 1 << 3;
    static constexpr std::uint8_t RowsAreSorted           = 1 << 4;
    static constexpr std::uint8_t ColumnsAreSorted        = 1 << 5;
    static constexpr std::uint8_t All                     = 1
                                                          + (1 << 1)
                                                          + (1 << 2)
                                                          + (1 << 3)
                                                          + (1 << 4)
                                                          + (1 << 5);
}; // struct MatrixCheck


/// @internal
template <class TValue, std::uint8_t Checks>
void CheckMatrix(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix)
{
    if constexpr (Checks == MatrixChecks::None) return;

    KRATOS_ERROR_IF_NOT(rMatrix.size1() + 1 == rMatrix.index1_data().size())
        << "input matrix has inconsistent row extents";
    KRATOS_ERROR_IF_NOT(rMatrix.index1_data()[rMatrix.size1()] == rMatrix.nnz())
        << "row extents of the input matrix do not consistently cover its contents";

    if constexpr (Checks & MatrixChecks::ColumnsAreSorted) {
        KRATOS_ERROR_IF_NOT(std::is_sorted(rMatrix.index1_data().begin(), rMatrix.index1_data().end()))
            << "the input matrix' columns are not sorted";
    }

    for (IndexType i_row=0ul; i_row<rMatrix.size1(); ++i_row) {
        const auto i_entry_begin = rMatrix.index1_data()[i_row];
        const auto i_entry_end = rMatrix.index1_data()[i_row + 1];

        if constexpr (Checks & MatrixChecks::RowsAreSorted)
            KRATOS_ERROR_IF_NOT(std::is_sorted(
                rMatrix.index2_data().begin() + i_entry_begin,
                rMatrix.index2_data().begin() + i_entry_end))
                    << "row " << i_row << " of the input matrix is not sorted";

        if constexpr (Checks & (  MatrixChecks::DiagonalExists
                                | MatrixChecks::DiagonalIsNonNegative
                                | MatrixChecks::DiagonalIsPositive
                                | MatrixChecks::IsDiagonallyDominant)) {
            std::optional<IndexType> maybe_diagonal;

            for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                const auto i_column = rMatrix.index2_data()[i_entry];
                if (i_column == i_row) {
                    maybe_diagonal = i_entry;
                    break;
                }
            } // for i_entry in range(i_entry_begin, i_entry_end)

            if (maybe_diagonal.has_value()) {
                const auto diagonal = rMatrix.value_data()[maybe_diagonal.value()];

                if constexpr (Checks & MatrixChecks::DiagonalIsNonNegative)
                    KRATOS_ERROR_IF_NOT(static_cast<TValue>(0) <= diagonal)
                        << "diagonal in row " << i_row << " is negative (" << diagonal << ")";

                if constexpr (Checks & MatrixChecks::DiagonalIsPositive)
                    KRATOS_ERROR_IF_NOT(static_cast<TValue>(0) < diagonal)
                        << "diagonal in row " << i_row << " is not positive (" << diagonal << ")";

                if constexpr (Checks & MatrixChecks::IsDiagonallyDominant) {
                    const auto abs_diagonal = std::abs(diagonal);
                    for (auto i_entry=i_entry_begin; i_entry<i_entry_end; ++i_entry) {
                        const auto i_column = rMatrix.index2_data()[i_entry];
                        const auto entry = rMatrix.value_data()[i_entry];
                        KRATOS_ERROR_IF(abs_diagonal <= std::abs(entry) && i_column != i_row)
                            << "row " << i_row << " of the input matrix is not diagonally dominant "
                            << "(entry in column " << i_column << " {"
                            << entry << "} has a larger magnitude than the diagonal {"
                            << diagonal << "})";
                    } // for i_entry in range(i_entry_begin, i_entry_end)
                } // if MatrixChecks::IsDiagonallyDominant
            } /*if diagonal found*/ else {
                if constexpr (Checks & (MatrixChecks::DiagonalExists | MatrixChecks::DiagonalIsPositive))
                    KRATOS_ERROR << "diagonal on row " << i_row << " does not exist";
            }
        } // if any of the diagonal checks
    } // for i_row in range(rMatrix.size1())
}


/** @brief Compute a scaled matrix-vector product and add it to the provided output vector.
 *  @details Computes @f[
 *              r += a * A @ b
 *           @f]
 *           where @p A is the input matrix, @p b the input vector, @p r the output vector and @p a the scaling coefficient.
 *           This function parallelizes on equal sized nonzero chunks instead of rows, which may lead to better
 *           performance when the matrix' rows have very different number of nonzeros.
 *  @tparam TLHSSparse Sparse space of the input matrix.
 *  @tparam TRHSSparse Sparse space of the input vector.
 *  @tparam TOutputSparse Sparse space of the output vector.
 *  @param rLhs Input matrix.
 *  @param rRhs Input vector.
 *  @param rOutput Output vector.
 *  @param Coefficient Coefficient to scale each output component by.
 */
template <class TLHSSparse, class TRHSSparse, class TOutputSparse>
void BalancedProduct(const typename TLHSSparse::MatrixType& rLhs,
                     const typename TRHSSparse::VectorType& rRhs,
                     typename TOutputSparse::VectorType& rOutput,
                     const typename TOutputSparse::DataType Coefficient = static_cast<typename TOutputSparse::DataType>(1))
{
    // Sanity checks.
    KRATOS_ERROR_IF_NOT(rLhs.size2() == rRhs.size() && rLhs.size1() == rOutput.size())
        << "incompatible matrix-vector product: "
        << "(" << rLhs.size1() << "x" << rLhs.size2() << ") "
        << "@ "
        << "(" << rRhs.size() << ") "
        << "=> (" << rOutput.size() << ")";

    KRATOS_TRY

    // Create partition for entries in the matrix.
    const auto thread_count = ParallelUtilities::GetNumThreads();
    std::vector<typename TLHSSparse::IndexType> partition(thread_count + 1);

    partition.front() = 0ul;
    const auto chunk_size = rLhs.nnz() / partition.size();
    for (typename TLHSSparse::IndexType i_end=1ul; i_end<partition.size(); ++i_end) {
        partition[i_end] = partition[i_end - 1] + chunk_size;
    } // for i_end in range(1, partition.size())
    partition.back() = rLhs.nnz();

    // Compute matrix-vector product.
    #define KRATOS_BALANCED_MATRIX_VECTOR_PRODUCT(OPERATOR)                                                                                 \
        IndexPartition<typename TLHSSparse::IndexType>(thread_count).for_each(                                                              \
        [&rLhs, &rRhs, &rOutput, &partition, Coefficient](const typename TLHSSparse::IndexType i_chunk){                                    \
            (void)Coefficient; /*<== suppress unused capture warnings*/                                                                     \
            /* Define the entry range to compute on. */                                                                                     \
            const auto i_chunk_begin = partition[i_chunk];                                                                                  \
            const auto i_chunk_end = partition[i_chunk + 1];                                                                                \
                                                                                                                                            \
            /* Find the initial row's index. */                                                                                             \
            const auto it_initial_row = std::lower_bound(rLhs.index1_data().begin(),                                                        \
                                                         rLhs.index1_data().end(),                                                          \
                                                         static_cast<typename TLHSSparse::IndexType>(i_chunk_begin));                       \
            typename TLHSSparse::IndexType i_row = std::distance(rLhs.index1_data().begin(), it_initial_row);                               \
            if (rLhs.index1_data()[i_row] != i_chunk_begin) --i_row;                                                                        \
                                                                                                                                            \
            do {                                                                                                                            \
                const auto i_row_begin = rLhs.index1_data()[i_row];                                                                         \
                const auto i_row_end = rLhs.index1_data()[i_row + 1];                                                                       \
                                                                                                                                            \
                const auto i_begin = std::max(i_row_begin, i_chunk_begin);                                                                  \
                const auto i_end = std::min(i_row_end, i_chunk_end);                                                                        \
                const typename TLHSSparse::IndexType chunk_size = i_end - i_begin;                                                          \
                                                                                                                                            \
                auto contribution = static_cast<typename TOutputSparse::DataType>(0);                                                       \
                                                                                                                                            \
                KRATOS_GET_ALIGNED_INDEX_ARRAY(const typename TLHSSparse::IndexType*, it_column, &*(rLhs.index2_data().begin() + i_begin)); \
                KRATOS_GET_ALIGNED_INDEX_ARRAY(const typename TLHSSparse::DataType*, it_entry, &*(rLhs.value_data().begin() + i_begin));    \
                KRATOS_GET_ALIGNED_INDEX_ARRAY(const typename TRHSSparse::DataType*, it_rhs, &*rRhs.begin());                               \
                                                                                                                                            \
                for (typename TLHSSparse::IndexType i=0; i<chunk_size; ++i) {                                                               \
                    const auto i_column = it_column[i];                                                                                     \
                    const auto entry = it_entry[i];                                                                                         \
                    contribution += entry * it_rhs[i_column];                                                                               \
                } /* for i_entry in range(i_begin, i_end) */                                                                                \
                                                                                                                                            \
                if constexpr (std::is_same_v<OPERATOR,std::plus<typename TOutputSparse::DataType>>) {                                       \
                    AtomicAdd(rOutput[i_row], contribution);                                                                                \
                } else if constexpr (std::is_same_v<OPERATOR,std::minus<typename TOutputSparse::DataType>>) {                               \
                    AtomicSub(rOutput[i_row], contribution);                                                                                \
                } else if constexpr (std::is_same_v<OPERATOR,std::multiplies<typename TOutputSparse::DataType>>) {                          \
                    AtomicAdd(rOutput[i_row], Coefficient * contribution);                                                                  \
                } else {                                                                                                                    \
                    static_assert(std::is_same_v<OPERATOR,void>, "unsupported operator");                                                   \
                }                                                                                                                           \
                                                                                                                                            \
                ++i_row;                                                                                                                    \
                if (i_end == i_chunk_end)                                                                                                   \
                    break;                                                                                                                  \
            } while (true);                                                                                                                 \
        }) /* for i_chunk in range(thread_count) */

    if (Coefficient == static_cast<typename TOutputSparse::DataType>(1)) {
        KRATOS_BALANCED_MATRIX_VECTOR_PRODUCT(std::plus<typename TOutputSparse::DataType>);
    } else if (Coefficient == static_cast<typename TOutputSparse::DataType>(-1)) {
        KRATOS_BALANCED_MATRIX_VECTOR_PRODUCT(std::minus<typename TOutputSparse::DataType>);
    } else {
        KRATOS_BALANCED_MATRIX_VECTOR_PRODUCT(std::multiplies<typename TOutputSparse::DataType>);
    }

    #undef KRATOS_BALANCED_MATRIX_VECTOR_PRODUCT
    KRATOS_CATCH("")
}


#undef KRATOS_GET_ALIGNED_INDEX_ARRAY


} // namespace Kratos
