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

#pragma once

// System includes
#include <tuple>
#include <variant>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) DistanceMatrix
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(DistanceMatrix);

    ///@}
    ///@name Life cycle
    ///@{

    DistanceMatrix();

    ///@}
    ///@name Public operations
    ///@{

    double GetDistance(
        const IndexType iIndex,
        const IndexType jIndex) const;

    double GetDistance(const IndexType EntryIndex) const;

    IndexType GetEntriesSize() const;

    IndexType GetEntryIndex(
        const IndexType iIndex,
        const IndexType jIndex) const;

    std::tuple<IndexType, IndexType> GetIndexPair(const IndexType iEntryIndex) const;

    IndexType GetNumberOfItems() const;

    void Update(const TensorAdaptor<double>& rDistancesTensorAdaptor);

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::vector<double> mDistances;

    IndexType mN;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/