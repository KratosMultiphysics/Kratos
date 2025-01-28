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
#include "expression/container_expression.h"

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

    void Update(
        std::variant<
            ContainerExpression<ModelPart::NodesContainerType>::Pointer,
            ContainerExpression<ModelPart::ConditionsContainerType>::Pointer,
            ContainerExpression<ModelPart::ElementsContainerType>::Pointer> pDistancesExpression);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::vector<double> mDistances;

    IndexType mN;

    ///@}
};

///@} // Kratos Classes

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DistanceMatrix& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} /* namespace Kratos.*/