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

// External includes

// Project includes

// Application includes

// Include base h
#include "sensor_generator_utils.h"

namespace Kratos {

bool SensorGeneratorUtils::IsPointInGeometry(
    const Point& rPoint,
    const Geometry<Node>& rGeometry)
{
    Point result;
    return rGeometry.IsInside(rPoint, result);
}

} /* namespace Kratos.*/