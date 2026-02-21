//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "distance_matrix.h"

namespace Kratos {
///@name Kratos Classes
///@{

DistanceMatrix::DistanceMatrix()
    : mN(0)
{
}

void DistanceMatrix::Update(const TensorAdaptor<double>& rDistancesTensorAdaptor)
{
    KRATOS_TRY

    auto shape = rDistancesTensorAdaptor.Shape();

    this->mN = shape[0];
    this->mDistances.resize(this->GetEntriesSize());

    const auto data_view = rDistancesTensorAdaptor.ViewData();
    const auto dimensionality = data_view.size() / shape[0];

    IndexPartition<IndexType>(this->mDistances.size()).for_each([&](const auto Index) {
        const auto& index_pair = GetIndexPair(Index);

        const auto i_index = std::get<0>(index_pair);
        const auto i_data_begin = i_index * dimensionality;
        const auto j_index = std::get<1>(index_pair);
        const auto j_data_begin = j_index * dimensionality;

        double distance = 0.0;

        for (IndexType i_comp = 0; i_comp < dimensionality; ++i_comp) {
            const double i_value = data_view[i_data_begin + i_comp];
            const double j_value = data_view[j_data_begin + i_comp];
            distance += std::pow(i_value - j_value, 2.0);
        }

        mDistances[Index] = std::sqrt(distance);
    });

    KRATOS_CATCH("");
}

double DistanceMatrix::GetDistance(
    const IndexType iIndex,
    const IndexType jIndex) const
{
    auto ordered_i = iIndex;
    auto ordered_j = jIndex;
    if (ordered_i != ordered_j) {
        if (ordered_j < ordered_i) {
            std::swap(ordered_i, ordered_j);
        }
        return mDistances[GetEntryIndex(ordered_i, ordered_j)];
    } else {
        return 0.0;
    }
}

double DistanceMatrix::GetDistance(const IndexType EntryIndex) const
{
    return mDistances[EntryIndex];
}

IndexType DistanceMatrix::GetEntriesSize() const
{
    return this->mN * (this->mN - 1) / 2;
}

IndexType DistanceMatrix::GetNumberOfItems() const
{
    return this->mN;
}

IndexType DistanceMatrix::GetEntryIndex(
    const IndexType iIndex,
    const IndexType jIndex) const
{
    return this->mN * iIndex + jIndex - ((iIndex + 2) * (iIndex + 1)) / 2;
}

std::tuple<IndexType, IndexType> DistanceMatrix::GetIndexPair(const IndexType EntryIndex) const
{

    const IndexType i = static_cast<IndexType>(std::floor(this->mN - 0.5 - std::sqrt(std::pow(this->mN - 0.5, 2) - 2 * EntryIndex)));
    const IndexType j = static_cast<IndexType>(EntryIndex - this->mN * i + (i + 1) * (i + 2) / 2);

    return std::make_tuple(i, j);
}

} /* namespace Kratos.*/