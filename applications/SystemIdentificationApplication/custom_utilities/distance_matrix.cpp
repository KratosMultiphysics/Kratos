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

void DistanceMatrix::Update(
    std::variant<
        ContainerExpression<ModelPart::NodesContainerType>::Pointer,
        ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
        ContainerExpression<ModelPart::ElementsContainerType>::Pointer> pDistancesExpression)
{
    KRATOS_TRY

    std::visit([&](const auto& pExp) {
        this->mN = pExp->GetContainer().size();
        this->mDistances.resize(this->GetEntriesSize());

        const auto& r_expression = pExp->GetExpression();
        const auto dimensionality = r_expression.GetItemComponentCount();

        IndexPartition<IndexType>(this->mDistances.size()).for_each([&](const auto Index) {
            const auto& index_pair = GetIndexPair(Index, this->mN);

            const auto i_index = std::get<0>(index_pair);
            const auto i_data_begin = i_index * dimensionality;
            const auto j_index = std::get<1>(index_pair);
            const auto j_data_begin = j_index * dimensionality;

            double distance = 0.0;

            for (IndexType i_comp = 0; i_comp < dimensionality; ++i_comp) {
                const double i_value = r_expression.Evaluate(i_index, i_data_begin, i_comp);
                const double j_value = r_expression.Evaluate(j_index, j_data_begin, i_comp);
                distance += std::pow(i_value - j_value, 2.0);
            }

            mDistances[Index] = std::sqrt(distance);
        });
    }, pDistancesExpression);

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
        return mDistances[GetEntryIndex(ordered_i, ordered_j, this->mN)];
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
    return EntriesSize(this->mN);
}

IndexType DistanceMatrix::GetNumberOfItems() const
{
    return this->mN;
}

IndexType DistanceMatrix::GetEntryIndex(
    const IndexType iIndex,
    const IndexType jIndex,
    const IndexType N)
{
    return N * iIndex + jIndex - ((iIndex + 2) * (iIndex + 1)) / 2;
}

std::tuple<IndexType, IndexType> DistanceMatrix::GetIndexPair(
    const IndexType EntryIndex,
    const IndexType N)
{

    const IndexType i = static_cast<IndexType>(std::floor(N - 0.5 - std::sqrt(std::pow(N - 0.5, 2) - 2 * EntryIndex)));
    const IndexType j = static_cast<IndexType>(EntryIndex - N * i + (i + 1) * (i + 2) / 2);

    // #pragma omp critical
    // {
    //     KRATOS_WATCH(i)
    //     KRATOS_WATCH(j)
    //     KRATOS_WATCH(EntryIndex)
    //     KRATOS_WATCH(N);
    // }

    return std::make_tuple(i, j);
}

IndexType DistanceMatrix::EntriesSize(const IndexType N)
{
    return N * (N - 1) / 2;
}

} /* namespace Kratos.*/