#define MCGS_INTERNAL

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::reorder, mcgs::CSRAdaptor
#include "partition.hpp" // mcgs::Partition
#include "defineMacros.hpp" // MCGS_EXPORT_SYMBOL

// --- STL Includes ---
#include <vector> // std::vector
#include <algorithm> // std::copy, std::swap, std::transform
#include <numeric> // std::iota


namespace mcgs {


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
Partition<TIndex>* reorder(const unsigned long rowCount, const unsigned long columnCount, const unsigned long nonzeroCount,
                           TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
                           TValue* pSolution, TValue* pRHS,
                           const Partition<TIndex>* pPartition)
{
    // Allocate arrays for the new partition
    std::vector<TIndex> newPartitionExtents(pPartition->size() + 1);
    newPartitionExtents[0] = static_cast<TIndex>(0);

    // Reorder matrix and vectors
    {
        // Allocate temporary matrix
        const auto partitionCount = pPartition->size();
        std::vector<TIndex> newRowExtents(rowCount + 1);
        std::vector<TIndex> newColumnIndices(nonzeroCount);
        std::vector<TValue> newNonzeros(nonzeroCount);

        // Allocate temporary vectors
        std::vector<TValue> solution, rhs;
        if (pSolution != nullptr) solution.resize(columnCount);
        if (pRHS != nullptr) rhs.resize(rowCount);
        const bool reorderSolution = !solution.empty();
        const bool reorderRHS = !rhs.empty();

        newRowExtents.front() = 0;

        // Compute new extents
        {
            TIndex iNewRow = 0;

            for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
                const auto itPartitionBegin = pPartition->begin(iPartition);
                const auto partitionSize = pPartition->size(iPartition);

                for (std::remove_const_t<decltype(partitionSize)> iLocal=0; iLocal<partitionSize; ++iLocal) {
                    const TIndex iOldRow = itPartitionBegin[iLocal];
                    const auto rowSize = pRowExtents[iOldRow + 1] - pRowExtents[iOldRow];
                    newRowExtents[iNewRow + 1] = newRowExtents[iNewRow] + rowSize;
                    ++iNewRow;
                } // for iLocal in range(parititionSize)

                newPartitionExtents[iPartition + 1] = newPartitionExtents[iPartition] + partitionSize;
            } // for iPartition in range(partitionCount)
        }

        std::vector<TIndex> columnMap(columnCount);

        #ifdef MCGS_OPENMP
        #pragma omp parallel
        #endif
        {
            for (std::size_t iPartition=0; iPartition<partitionCount; ++iPartition) {
                auto itPartitionBegin = pPartition->begin(iPartition);
                const auto partitionSize = pPartition->size(iPartition);

                #ifdef MCGS_OPENMP
                #pragma omp for
                #endif
                for (int iLocal=0; iLocal<static_cast<int>(partitionSize); ++iLocal) {
                    const TIndex iOldRow = itPartitionBegin[iLocal];
                    const TIndex iNewRow = newPartitionExtents[iPartition] + iLocal;
                    columnMap[iOldRow] = iNewRow;

                    std::copy(pColumnIndices + pRowExtents[iOldRow],
                            pColumnIndices + pRowExtents[iOldRow + 1],
                            newColumnIndices.data() + newRowExtents[iNewRow]);

                    std::copy(pNonzeros + pRowExtents[iOldRow],
                            pNonzeros + pRowExtents[iOldRow + 1],
                            newNonzeros.data() + newRowExtents[iNewRow]);

                    if (reorderSolution) solution[iNewRow] = pSolution[iOldRow];
                    if (reorderRHS) rhs[iNewRow] = pRHS[iOldRow];
                } // for iLocal in range(parititionSize)
            } // for iPartition in range(partitionCount)
        } // omp parallel

        std::copy(newRowExtents.begin(), newRowExtents.end(), pRowExtents);
        std::copy(newNonzeros.begin(), newNonzeros.end(), pNonzeros);
        std::copy(solution.begin(), solution.end(), pSolution);
        std::copy(rhs.begin(), rhs.end(), pRHS);
        std::transform(newColumnIndices.begin(),
                       newColumnIndices.end(),
                       pColumnIndices,
                       [&columnMap](const TIndex iOldColumn) -> TIndex {
                               return columnMap[iOldColumn];
                       });
    }

    // Allocate arrays for the new partition
    // and assign the new row indices
    std::vector<TIndex> newRowIndices(rowCount + 1);
    std::iota(newRowIndices.begin(), newRowIndices.end(), TIndex(0));

    return new Partition<TIndex>(std::move(newPartitionExtents), std::move(newRowIndices));
}


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int revertReorder(TValue* pRHS,
                  const unsigned long columnCount,
                  const Partition<TIndex>* pPartition)
{
    std::vector<TValue> swap(columnCount);

    TIndex iNewRow = 0;
    for (typename Partition<TIndex>::size_type iPartition=0; iPartition<pPartition->size(); ++iPartition) {
        for (auto itPartition=pPartition->begin(iPartition); itPartition!=pPartition->end(iPartition); ++itPartition) {
            swap[*itPartition] = pRHS[iNewRow++];
        }
    }

    if (iNewRow != static_cast<TIndex>(columnCount)) {
        return MCGS_FAILURE;
    }

    std::copy(swap.begin(), swap.end(), pRHS);
    return MCGS_SUCCESS;
}


template <class TIndex, class TValue>
MCGS_EXPORT_SYMBOL
int revertReorder(const unsigned long rowCount, const unsigned long columnCount, const unsigned long nonzeroCount,
                  TIndex* pRowExtents, TIndex* pColumnIndices, TValue* pNonzeros,
                  TValue* pSolution, TValue* pRHS,
                  const Partition<TIndex>* pPartition)
{
    // Construct inverse partition
    std::vector<TIndex> partitionExtents {static_cast<TIndex>(0), static_cast<TIndex>(rowCount)}, rowIndices(rowCount);
    const TIndex* itPartitionBegin = pPartition->begin(0);

    #ifdef MCGS_OPENMP
    #pragma omp parallel for
    #endif
    for (int iRow=0; iRow<static_cast<int>(rowCount); ++iRow) {
        rowIndices[itPartitionBegin[iRow]] = iRow;
    } // for iRow in range(rowCount)

    Partition<TIndex> inversePartition(std::move(partitionExtents), std::move(rowIndices));

    // Reorder
    Partition<TIndex>* pDummy = reorder(rowCount, columnCount, nonzeroCount,
                                        pRowExtents, pColumnIndices, pNonzeros,
                                        pSolution, pRHS,
                                        &inversePartition);

    if (pDummy == nullptr) {
        return MCGS_FAILURE;
    } else {
        delete pDummy;
    }

    return MCGS_SUCCESS;
}


#define MCGS_INSTANTIATE_REORDER(TIndex, TValue)                                                                                            \
    template MCGS_EXPORT_SYMBOL Partition<TIndex>* reorder<TIndex,TValue>(const unsigned long, const unsigned long, const unsigned long,    \
                                                       TIndex*, TIndex*, TValue*,                                                           \
                                                       TValue*, TValue*,                                                                    \
                                                       const Partition<TIndex>*);                                                           \
    template MCGS_EXPORT_SYMBOL int revertReorder<TIndex,TValue>(TValue*,                                                                   \
                                                                 const unsigned long,                                                       \
                                                                 const Partition<TIndex>*);                                                 \
    template MCGS_EXPORT_SYMBOL int revertReorder<TIndex,TValue>(const unsigned long, const unsigned long, const unsigned long,             \
                                                                 TIndex*, TIndex*, TValue*,                                                 \
                                                                 TValue*, TValue*,                                                          \
                                                                 const Partition<TIndex>*)

MCGS_INSTANTIATE_REORDER(int, double);
MCGS_INSTANTIATE_REORDER(long, double);
MCGS_INSTANTIATE_REORDER(unsigned, double);
MCGS_INSTANTIATE_REORDER(unsigned long, double);
MCGS_INSTANTIATE_REORDER(int, float);
MCGS_INSTANTIATE_REORDER(long, float);
MCGS_INSTANTIATE_REORDER(unsigned, float);
MCGS_INSTANTIATE_REORDER(unsigned long, float);

#undef MCGS_INSTANTIATE_REORDER
#undef MCGS_INTERNAL
#include "undefineMacros.hpp"


} // namespace mcgs
