//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"

namespace Kratos
{

    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, double, PERTURBATION_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, bool, ADAPT_PERTURBATION_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, bool, HAS_ROTATION_DOFS)

    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DIGITAL_TWIN_APPLICATION, SENSOR_DIRECTION)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, double, SENSOR_WEIGHT)

    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, std::string, SENSOR_NAME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, int, SENSOR_ID)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(DIGITAL_TWIN_APPLICATION, SENSOR_LOCATION)

    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, int, SENSOR_NODE_ID)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, int, SENSOR_ELEMENT_ID)
    KRATOS_DEFINE_APPLICATION_VARIABLE(DIGITAL_TWIN_APPLICATION, int, SENSOR_ELEMENT_NODE_INDEX)

} // namespace Kratos
