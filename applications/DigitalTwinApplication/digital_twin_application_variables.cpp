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

// Include base h
#include "digital_twin_application_variables.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(double, PERTURBATION_SIZE)
    KRATOS_CREATE_VARIABLE(bool, ADAPT_PERTURBATION_SIZE)
    KRATOS_CREATE_VARIABLE(bool, HAS_ROTATION_DOFS)

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SENSOR_DIRECTION)
    KRATOS_CREATE_VARIABLE(double, SENSOR_WEIGHT)

    KRATOS_CREATE_VARIABLE(std::string, SENSOR_NAME)
    KRATOS_CREATE_VARIABLE(std::string, TEST_ANALYSIS_NAME)
    KRATOS_CREATE_VARIABLE(int, SENSOR_ID)
    KRATOS_CREATE_VARIABLE(int, SENSOR_CLUSTER_ID)
    KRATOS_CREATE_VARIABLE(double, SENSOR_MEASURED_VALUE)
    KRATOS_CREATE_VARIABLE(double, SENSOR_ERROR)
    KRATOS_CREATE_VARIABLE(double, SENSOR_SENSITIVITY_NORM_INF)
    KRATOS_CREATE_VARIABLE(double, SENSOR_SENSITIVITY_NORM_L2)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SENSOR_LOCATION)

    KRATOS_CREATE_VARIABLE(int, SENSOR_NODE_ID)
    KRATOS_CREATE_VARIABLE(int, SENSOR_ELEMENT_ID)
    KRATOS_CREATE_VARIABLE(int, SENSOR_ELEMENT_NODE_INDEX)
    KRATOS_CREATE_VARIABLE(Vector, SENSOR_ENTITY_IDS)
} // namespace Kratos
