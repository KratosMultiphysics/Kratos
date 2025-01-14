#pragma once

// --- Internal Includes ---
#include "defineMacros.hpp"

// --- STL Includes ---
#include <vector> // std::vector
#include <cstddef> // std::size_t, std::ptrdiff_t
#include <iterator> // std::distance


namespace mcgs {


template <class TIndex>
class MCGS_EXPORT_SYMBOL Partition
{
private:
    bool _isContiguous;

    std::vector<TIndex> _partitionExtents;

    std::vector<TIndex> _rowIndices;

    Partition() = delete;
    Partition(Partition&&) = delete;
    Partition(const Partition&) = delete;
    Partition& operator=(Partition&&) = delete;
    Partition& operator=(const Partition&) = delete;

public:
    using value_type = TIndex;
    using reference = TIndex&;
    using const_reference = const TIndex&;
    using iterator = TIndex*;
    using const_iterator = const TIndex*;
    using difference_type = std::ptrdiff_t;
    using size_type = std::size_t;

    template <class TColor>
    Partition(const TColor* pColors, const TIndex rowCount);

    Partition(std::vector<TIndex>&& rPartitionExtents,
              std::vector<TIndex>&& rRowIndices) noexcept;

    size_type size() const noexcept
    {return _partitionExtents.size() - 1;}

    size_type size(const size_type iPartition) const noexcept
    {return std::distance(this->begin(iPartition), this->end(iPartition));}

    const_iterator begin(const size_type iPartition) const noexcept
    {return &_rowIndices[_partitionExtents[iPartition]];}

    const_iterator end(const size_type iPartition) const noexcept
    {return &_rowIndices[_partitionExtents[iPartition + 1]];}

    bool isContiguous() const noexcept
    {return _isContiguous;}
};


} // namespace mcgs
