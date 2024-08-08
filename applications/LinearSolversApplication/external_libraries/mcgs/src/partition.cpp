#define MCGS_INTERNAL

// --- Internal Includes ---
#include "mcgs/mcgs.hpp" // mcgs::partition
#include "partition.hpp" // mcgs::Partition
#include "defineMacros.hpp" // MCGS_EXPORT_SYMBOL

// --- STL Includes ---
#include <cstddef> // std::size_t
#include <vector> // std::vector
#include <unordered_map> // std::unordered_map
#include <algorithm> // std::copy, std::is_sorted


namespace mcgs {


template <class TIndex>
template <class TColor>
Partition<TIndex>::Partition(const TColor* pColors, const TIndex rowCount)
    : _isContiguous(false),
      _partitionExtents(),
      _rowIndices()
{
    std::unordered_map<
        TColor,
        std::vector<TIndex>
    > partitionMap;

    for (TIndex iColumn=0; iColumn<rowCount; ++iColumn) {
        const TColor color = pColors[iColumn];
        partitionMap.emplace(color, std::vector<TIndex> {}) // <== make sure an entry is mapped to color
            .first                                          // <== iterator pointing to the entry
            ->second                                        // <== reference to the mapped vector
            .push_back(iColumn);                            // <== insert the column index into the mapped vector
    } // for iColumn in range(rowCount)

    _rowIndices.resize(rowCount + 1);
    _partitionExtents.resize(partitionMap.size() + 1);

    TIndex iColor = 0;
    TIndex rowCounter = 0;

    for ([[maybe_unused]] const auto& [color, rColumns] : partitionMap) {
        _partitionExtents[iColor] = rowCounter;
        std::copy(rColumns.begin(), rColumns.end(), _rowIndices.begin() + rowCounter);

        ++iColor;
        rowCounter += rColumns.size();
    }

    _partitionExtents.back() = rowCounter;
    _rowIndices.back() = rowCount;
    _isContiguous = std::is_sorted(_rowIndices.begin(), _rowIndices.end());
}


template <class TIndex>
Partition<TIndex>::Partition(std::vector<TIndex>&& rPartitionExtents,
                             std::vector<TIndex>&& rRowIndices) noexcept
    : _isContiguous(false),
      _partitionExtents(std::move(rPartitionExtents)),
      _rowIndices(std::move(rRowIndices))
{
    _isContiguous = std::is_sorted(_rowIndices.begin(), _rowIndices.end());
}


template <class TIndex, class TColor>
[[nodiscard]] MCGS_EXPORT_SYMBOL
Partition<TIndex>* makePartition(const TColor* pColors,
                                 const TIndex rowCount)
{
    return new Partition<TIndex>(pColors, rowCount);
}


template <class TIndex>
MCGS_EXPORT_SYMBOL
void destroyPartition(Partition<TIndex>* pPartition)
{
    delete pPartition;
}


#define MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, TColor)                                  \
    template MCGS_EXPORT_SYMBOL Partition<TIndex>* makePartition<TIndex,TColor>(const TColor*, const TIndex);

#define MCGS_INSTANTIATE_PARTITION(TIndex)                                          \
    template class Partition<TIndex>;                                               \
    MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, unsigned)                            \
    MCGS_INSTANTIATE_PARTITION_FACTORY(TIndex, unsigned long)                       \
    template MCGS_EXPORT_SYMBOL void destroyPartition<TIndex>(Partition<TIndex>*)

MCGS_INSTANTIATE_PARTITION(int);
MCGS_INSTANTIATE_PARTITION(long);
MCGS_INSTANTIATE_PARTITION(unsigned);
MCGS_INSTANTIATE_PARTITION(unsigned long);

#undef MCGS_INSTANTIATE_PARTITION_FACTORY
#undef MCGS_INSTANTIATE_PARTITION
#undef MCGS_INTERNAL
#include "undefineMacros.hpp"


} // namespace mcgs
