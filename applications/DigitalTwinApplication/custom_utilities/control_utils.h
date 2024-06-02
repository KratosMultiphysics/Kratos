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

#pragma once

// System includes
#include <tuple>

// External includes

// Project includes
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class ControlUtils
{
public:
    ///@name Public static operations
    ///@{

    static IndexType GetDistVectorSize(const IndexType N);

    static IndexType GetDistIndexFromPairIndices(
        const IndexType N,
        const IndexType I,
        const IndexType J);

    static std::tuple<IndexType, IndexType> GetPairIndicesFromDistIndex(
        const IndexType N,
        const IndexType DistIndex);

    template<class TContainerType>
    static void AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/