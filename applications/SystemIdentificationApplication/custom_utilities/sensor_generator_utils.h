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

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorGeneratorUtils
{
public:
    ///@name Public static operations
    ///@{

    KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) static bool IsPointInGeometry(
        const Point& rPoint,
        const Geometry<Node>& rGeometry);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/