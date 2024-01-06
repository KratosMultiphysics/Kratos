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
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_sensors_to_python.h"
#include "digital_twin_application.h"
#include "digital_twin_application_variables.h"

namespace Kratos::Python {

PYBIND11_MODULE(KratosDigitalTwinApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosDigitalTwinApplication, KratosDigitalTwinApplication::Pointer, KratosApplication>(
        m, "KratosDigitalTwinApplication")
        .def(py::init<>());

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomSensorsToPython(m);

    // registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERTURBATION_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADAPT_PERTURBATION_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HAS_ROTATION_DOFS)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SENSOR_DIRECTION)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_WEIGHT)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_NAME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEST_ANALYSIS_NAME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_GROUP_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_NEIGHBOUR_IDS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_LOCATION_ROBUSTNESS)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_CLUSTER_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_MEASURED_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_COVERAGE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_SENSITIVITY_NORM_INF)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_SENSITIVITY_NORM_L2)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SENSOR_LOCATION)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_NODE_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ELEMENT_ID)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ELEMENT_NODE_INDEX)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ENTITY_IDS)

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);
}

} // namespace Kratos::Python
