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

// System includes

// External includes

// Project includes
#include "digital_twin_application.h"
#include "digital_twin_application_variables.h"

namespace Kratos {

KratosDigitalTwinApplication::KratosDigitalTwinApplication()
    : KratosApplication("DigitalTwinApplication")
{
}

void KratosDigitalTwinApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosDigitalTwinApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(PERTURBATION_SIZE)
    KRATOS_REGISTER_VARIABLE(ADAPT_PERTURBATION_SIZE)
    KRATOS_REGISTER_VARIABLE(HAS_ROTATION_DOFS)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SENSOR_DIRECTION)
    KRATOS_REGISTER_VARIABLE(SENSOR_WEIGHT)

    KRATOS_REGISTER_VARIABLE(SENSOR_NAME)
    KRATOS_REGISTER_VARIABLE(TEST_ANALYSIS_NAME)
    KRATOS_REGISTER_VARIABLE(SENSOR_ID)
    KRATOS_REGISTER_VARIABLE(SENSOR_CLUSTER_ID)
    KRATOS_REGISTER_VARIABLE(SENSOR_MEASURED_VALUE)
    KRATOS_REGISTER_VARIABLE(SENSOR_ERROR)
    KRATOS_REGISTER_VARIABLE(SENSOR_SENSITIVITY_NORM_INF)
    KRATOS_REGISTER_VARIABLE(SENSOR_SENSITIVITY_NORM_L2)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SENSOR_LOCATION)

    KRATOS_REGISTER_VARIABLE(SENSOR_NODE_ID)
    KRATOS_REGISTER_VARIABLE(SENSOR_ELEMENT_ID)
    KRATOS_REGISTER_VARIABLE(SENSOR_ELEMENT_NODE_INDEX)
    KRATOS_REGISTER_VARIABLE(SENSOR_ENTITY_IDS)
}

} // namespace Kratos.
