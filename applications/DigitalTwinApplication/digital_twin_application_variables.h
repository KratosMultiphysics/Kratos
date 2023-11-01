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

// External includes

// Project includes
#include "containers/variable.h"
#include "includes/define.h"

namespace Kratos
{

    KRATOS_DEFINE_VARIABLE(double, PERTURBATION_SIZE)
    KRATOS_DEFINE_VARIABLE(bool, ADAPT_PERTURBATION_SIZE)
    KRATOS_DEFINE_VARIABLE(bool, HAS_ROTATION_DOFS)
    KRATOS_DEFINE_VARIABLE(int, SENSOR_NODE_ID)
    KRATOS_DEFINE_VARIABLE(double, SENSOR_WEIGHT)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SENSOR_DIRECTION)

} // namespace Kratos
