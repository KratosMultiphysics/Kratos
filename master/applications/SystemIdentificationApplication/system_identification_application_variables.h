//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
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
#include "expression/container_expression.h"

namespace Kratos
{
    // Adjoint variables
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(SYSTEM_IDENTIFICATION_APPLICATION, ADJOINT_DISPLACEMENT)
    KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS(SYSTEM_IDENTIFICATION_APPLICATION, ADJOINT_ROTATION)

    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, double, PERTURBATION_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, bool, ADAPT_PERTURBATION_SIZE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, bool, HAS_ROTATION_DOFS)

    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, std::string, SENSOR_NAME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, std::string, TEST_ANALYSIS_NAME)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, double, SENSOR_MEASURED_VALUE)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, double, SENSOR_ERROR)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SYSTEM_IDENTIFICATION_APPLICATION, int, SENSOR_ELEMENT_ID)

} // namespace Kratos
