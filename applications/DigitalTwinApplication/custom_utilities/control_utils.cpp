//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes
#include <cmath>
#include <algorithm>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes

// Include base h
#include "control_utils.h"

namespace Kratos {

IndexType ControlUtils::GetDistVectorSize(const IndexType N)
{
    return N * (N - 1) / 2;
}

IndexType ControlUtils::GetDistIndexFromPairIndices(
    const IndexType N,
    const IndexType I,
    const IndexType J)
{
    return N * I + J - ((I + 2) * (I + 1)) / 2;
}

std::tuple<IndexType, IndexType> ControlUtils::GetPairIndicesFromDistIndex(
    const IndexType N,
    const IndexType DistIndex)
{
    const int b = 1 - 2 * N;
    const int i = std::floor((-b - std::sqrt(b * b - 8 * DistIndex)) / 2);
    const int j = (DistIndex + i * (b + i + 2) / 2 + 1);
    return std::make_tuple(i, j);
}

template<class TContainerType>
void ControlUtils::AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer)
{
    KRATOS_TRY

    const IndexType number_of_entities = rSourceContainer.size();

    KRATOS_ERROR_IF_NOT(number_of_entities == rDestinationContainer.size())
        << "Number of entities mismatch [ rSourceContainer.size() = " << number_of_entities
        << ", rDestinationContainer.size() = " << rDestinationContainer.size() << " ].\n";

    // TODO: Temporarily fix for the pointer vector set not being ordered.
    rSourceContainer.Unique();
    rDestinationContainer.Unique();

    IndexPartition<IndexType>(number_of_entities).for_each([&rSourceContainer, &rDestinationContainer](const auto Index) {
        auto& r_destination_entity = *(rDestinationContainer.begin() + Index);
        auto p_itr = rSourceContainer.find(r_destination_entity.Id());

        KRATOS_ERROR_IF(p_itr == rSourceContainer.end()) <<
            "The entity with id = " << r_destination_entity.Id() << " not found in the source container.";

        r_destination_entity.SetProperties(p_itr->pGetProperties());
    });

    KRATOS_CATCH("");
}

// template instantiations
template void ControlUtils::AssignEquivalentProperties(ModelPart::ConditionsContainerType&, ModelPart::ConditionsContainerType&);
template void ControlUtils::AssignEquivalentProperties(ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&);

} /* namespace Kratos.*/