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
#include <variant>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/kratos_parameters.h"
#include "geometries/geometry.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorUtils
{
public:
    ///@name Public static operations
    ///@{

    static bool IsPointInGeometry(
        const Point& rPoint,
        const Geometry<Node>& rGeometry);

    static void ReadVariableData(
        DataValueContainer& rDataValueContainer,
        Parameters VariableDataParameters);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/