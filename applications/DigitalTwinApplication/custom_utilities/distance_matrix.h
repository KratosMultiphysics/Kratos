//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"


// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class DistanceMatrix
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(DistanceMatrix);

    ///@}
    ///@name Life cycle
    ///@{

    DistanceMatrix(const std::vector<double>& rCompressedDistancesMatrix);

    ///@}
    ///@name Public operations
    ///@{

    const std::vector<double>& GetCompressedDistanceMatrix() const { return mCompressedDistanceMatrix; }

    double GetDistance(
        const IndexType Index1,
        const IndexType Index2) const;

    ///@}
private:
    ///@name Private member variables
    ///@{

    IndexType mNumberOfItems;

    const std::vector<double> mCompressedDistanceMatrix;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/