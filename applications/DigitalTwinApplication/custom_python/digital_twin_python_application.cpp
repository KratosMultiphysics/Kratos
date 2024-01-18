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

    AddCustomUtilitiesToPython(m);
    AddCustomSensorsToPython(m);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADJOINT_DISPLACEMENT)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, ADJOINT_ROTATION)

    // registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERTURBATION_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ADAPT_PERTURBATION_SIZE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HAS_ROTATION_DOFS)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_NAME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TEST_ANALYSIS_NAME)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_MEASURED_VALUE)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ERROR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SENSOR_ELEMENT_ID)
}

} // namespace Kratos::Python
